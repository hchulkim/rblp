#' @title Integration Configuration
#' @description Configuration for building integration nodes and weights for
#'   approximating integrals over the distribution of random coefficients.
#' @export
BLPIntegration <- R6::R6Class("BLPIntegration",
  public = list(
    #' @description Create integration configuration
    #' @param specification Character: method name
    #' @param size Integer: number of draws or quadrature level
    #' @param seed Optional RNG seed
    initialize = function(specification, size, seed = NULL) {
      valid <- c("monte_carlo", "halton", "lhs", "mlhs", "product", "grid")
      if (!specification %in% valid) {
        stop(sprintf("Invalid specification '%s'. Choose from: %s",
                     specification, paste(valid, collapse = ", ")))
      }
      private$specification_ <- specification
      private$size_ <- as.integer(size)
      private$seed_ <- seed
    },

    #' @description Build nodes and weights
    #' @param dimensions Number of integration dimensions (K2)
    #' @return List with \code{nodes} (N x K2) and \code{weights} (N x 1)
    build = function(dimensions) {
      switch(private$specification_,
        "monte_carlo" = private$build_monte_carlo(dimensions),
        "halton"      = private$build_halton(dimensions),
        "lhs"         = private$build_lhs(dimensions),
        "mlhs"        = private$build_mlhs(dimensions),
        "product"     = private$build_product(dimensions),
        "grid"        = private$build_grid(dimensions),
        stop("Unknown specification")
      )
    },

    #' @description Print integration configuration
    print = function(...) {
      cat(sprintf("BLPIntegration: %s (size=%d)\n",
                  private$specification_, private$size_))
      invisible(self)
    }
  ),
  private = list(
    specification_ = NULL,
    size_ = NULL,
    seed_ = NULL,

    # Monte Carlo integration: draw i.i.d. N(0,1) random variates for each
    # dimension of the random coefficients. Equal weights 1/N implement the law
    # of large numbers approximation E[f(v)] ~ (1/N) sum f(v_i). Simple but
    # converges slowly at rate O(N^{-1/2}).
    build_monte_carlo = function(dims) {
      if (!is.null(private$seed_)) set.seed(private$seed_)
      n <- private$size_
      nodes <- matrix(stats::rnorm(n * dims), nrow = n, ncol = dims)
      weights <- rep(1 / n, n)
      list(nodes = nodes, weights = weights)
    },

    # Halton quasi-random sequences: low-discrepancy points that fill the
    # unit hypercube more evenly than pseudo-random draws, yielding faster
    # convergence for numerical integration (roughly O(N^{-1} log(N)^d)).
    # Each dimension uses a different prime base (2, 3, 5, ...) to avoid
    # correlation across dimensions. The first 1000 points are discarded
    # because early Halton draws in high bases can exhibit strong correlation
    # patterns. After generating uniform [0,1] draws, the inverse-normal
    # transform (qnorm) maps them to N(0,1) nodes for the random coefficients.
    build_halton = function(dims) {
      if (!is.null(private$seed_)) set.seed(private$seed_)
      n <- private$size_
      discard <- 1000L

      if (requireNamespace("randtoolbox", quietly = TRUE)) {
        raw <- randtoolbox::halton(n + discard, dim = dims, normal = FALSE)
        raw <- raw[(discard + 1):(discard + n), , drop = FALSE]
      } else {
        # Fallback: van der Corput sequence in each prime base. Each index i
        # is expanded in base-b digits and "reflected" about the decimal point,
        # producing a deterministic sequence that is equidistributed in [0,1].
        primes <- c(2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47)
        raw <- matrix(0, n + discard, dims)
        for (d in seq_len(dims)) {
          base <- primes[d]
          for (i in seq_len(n + discard)) {
            f <- 1
            r <- 0
            ii <- i
            while (ii > 0) {
              f <- f / base
              r <- r + f * (ii %% base)
              ii <- ii %/% base
            }
            raw[i, d] <- r
          }
        }
        raw <- raw[(discard + 1):(discard + n), , drop = FALSE]
      }

      # Map uniform [0,1] draws to standard normal via the quantile function,
      # then clamp any non-finite values (from qnorm(0) or qnorm(1)) to zero.
      nodes <- stats::qnorm(raw)
      nodes[!is.finite(nodes)] <- 0
      weights <- rep(1 / n, n)
      list(nodes = nodes, weights = weights)
    },

    # Latin Hypercube Sampling (LHS): a stratified sampling scheme that
    # partitions each marginal dimension into N equal-probability strata and
    # draws exactly one point per stratum. This guarantees good marginal
    # coverage even at moderate sample sizes, reducing variance relative to
    # pure Monte Carlo. Each dimension is independently permuted so that the
    # joint distribution remains approximately uniform over the hypercube.
    # The uniform draws are then inverse-normal-transformed to N(0,1) nodes.
    build_lhs = function(dims) {
      if (!is.null(private$seed_)) set.seed(private$seed_)
      n <- private$size_
      if (requireNamespace("lhs", quietly = TRUE)) {
        raw <- lhs::randomLHS(n, dims)
      } else {
        # Manual LHS: for each dimension, create a random permutation of strata
        # {1,...,N}, then place a uniform draw within each stratum: (perm - U)/N.
        raw <- matrix(0, n, dims)
        for (d in seq_len(dims)) {
          perm <- sample.int(n)
          raw[, d] <- (perm - stats::runif(n)) / n
        }
      }
      nodes <- stats::qnorm(raw)
      nodes[!is.finite(nodes)] <- 0
      weights <- rep(1 / n, n)
      list(nodes = nodes, weights = weights)
    },

    # Modified Latin Hypercube Sampling (Hess, Train, Polak 2004): generates
    # multiple independent LHS draws and averages the resulting normal nodes.
    # By the central limit theorem, this averaging smooths out the randomness
    # in any single LHS draw while preserving the stratification benefit.
    # The result is a set of integration nodes with lower simulation variance
    # than standard LHS, approaching the regularity of quasi-random sequences.
    # Using 10 replications balances variance reduction against computational cost.
    build_mlhs = function(dims) {
      if (!is.null(private$seed_)) set.seed(private$seed_)
      n <- private$size_
      n_reps <- 10L
      nodes_sum <- matrix(0, n, dims)
      for (r in seq_len(n_reps)) {
        raw <- matrix(0, n, dims)
        for (d in seq_len(dims)) {
          perm <- sample.int(n)
          raw[, d] <- (perm - stats::runif(n)) / n
        }
        # Accumulate the inverse-normal-transformed draws across replications
        nodes_sum <- nodes_sum + stats::qnorm(raw)
      }
      # Average over replications to get the final MLHS nodes
      nodes <- nodes_sum / n_reps
      nodes[!is.finite(nodes)] <- 0
      weights <- rep(1 / n, n)
      list(nodes = nodes, weights = weights)
    },

    # Gauss-Hermite product rule: deterministic quadrature that is exact for
    # polynomial integrands up to degree 2*level - 1. The raw Gauss-Hermite
    # rule integrates f(x)*exp(-x^2); to evaluate E[f(v)] for v ~ N(0,1), we
    # apply the change of variables v = x*sqrt(2), yielding the scaling
    # nodes * sqrt(2) and weights / sqrt(pi). In multiple dimensions, the
    # product (tensor) rule takes all combinations of 1-D nodes, so N = level^dims
    # points total. This gives high accuracy for smooth integrands but suffers
    # from the "curse of dimensionality" -- the number of nodes grows
    # exponentially in the number of random coefficient dimensions.
    build_product = function(dims) {
      level <- private$size_
      gh <- gauss_hermite(level)
      # Transform from Hermite measure exp(-x^2) to standard normal (2*pi)^{-1/2}*exp(-x^2/2)
      gh_nodes <- gh$nodes * sqrt(2)
      gh_weights <- gh$weights / sqrt(pi)

      if (dims == 1) {
        return(list(nodes = matrix(gh_nodes, ncol = 1), weights = gh_weights))
      }

      # Build the full tensor product grid: each multidimensional node is a
      # combination of 1-D nodes, and its weight is the product of the
      # corresponding 1-D weights (because the normal distribution factors).
      grids <- rep(list(seq_len(level)), dims)
      idx <- as.matrix(expand.grid(grids))
      n <- nrow(idx)
      nodes <- matrix(0, n, dims)
      weights <- rep(1, n)
      for (d in seq_len(dims)) {
        nodes[, d] <- gh_nodes[idx[, d]]
        weights <- weights * gh_weights[idx[, d]]
      }
      list(nodes = nodes, weights = weights)
    },

    build_grid = function(dims) {
      # Sparse grid (simplified Smolyak)
      # For now, fall back to product rule with warning for high dims
      if (dims > 5) {
        warning("Sparse grid approximated by product rule for dims > 5")
      }
      private$build_product(dims)
    }
  )
)

#' Gauss-Hermite quadrature points and weights
#' @param n Number of points
#' @return List with nodes and weights for integrating f(x) * exp(-x^2)
#' @keywords internal
gauss_hermite <- function(n) {
  if (n < 1 || n > 30) stop("Gauss-Hermite level must be between 1 and 30")
  # Golub-Welsch eigenvalue method: the n-point Gauss-Hermite nodes are the
  # eigenvalues of the symmetric tridiagonal Jacobi matrix associated with the
  # Hermite polynomials. The sub-diagonal entries are b_k = sqrt(k/2), derived
  # from the three-term recurrence relation of the (physicists') Hermite
  # polynomials. This approach is numerically stable and avoids root-finding.
  i <- seq_len(n - 1)
  b <- sqrt(i / 2)
  cm <- diag(0, n)
  for (k in seq_along(b)) {
    cm[k, k + 1] <- b[k]
    cm[k + 1, k] <- b[k]
  }
  eig <- eigen(cm, symmetric = TRUE)
  nodes <- eig$values
  # The quadrature weights equal the squared first component of each
  # eigenvector times sqrt(pi). This follows from the connection between
  # eigenvectors of the Jacobi matrix and the Christoffel numbers.
  weights <- eig$vectors[1, ]^2 * sqrt(pi)
  ord <- order(nodes)
  list(nodes = nodes[ord], weights = weights[ord])
}

#' Create integration configuration
#'
#' @param specification Character: "monte_carlo", "halton", "lhs", "mlhs", "product", or "grid"
#' @param size Integer: number of draws or quadrature level
#' @param seed Optional RNG seed for reproducibility
#' @return A BLPIntegration object
#' @export
#' @examples
#' int_mc <- blp_integration("monte_carlo", 50, seed = 42)
#' int_prod <- blp_integration("product", 7)
blp_integration <- function(specification, size, seed = NULL) {
  BLPIntegration$new(specification, size, seed)
}

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

    build_monte_carlo = function(dims) {
      if (!is.null(private$seed_)) set.seed(private$seed_)
      n <- private$size_
      nodes <- matrix(stats::rnorm(n * dims), nrow = n, ncol = dims)
      weights <- rep(1 / n, n)
      list(nodes = nodes, weights = weights)
    },

    build_halton = function(dims) {
      if (!is.null(private$seed_)) set.seed(private$seed_)
      n <- private$size_
      discard <- 1000L

      if (requireNamespace("randtoolbox", quietly = TRUE)) {
        raw <- randtoolbox::halton(n + discard, dim = dims, normal = FALSE)
        raw <- raw[(discard + 1):(discard + n), , drop = FALSE]
      } else {
        # Simple Halton sequence fallback
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

      nodes <- stats::qnorm(raw)
      nodes[!is.finite(nodes)] <- 0
      weights <- rep(1 / n, n)
      list(nodes = nodes, weights = weights)
    },

    build_lhs = function(dims) {
      if (!is.null(private$seed_)) set.seed(private$seed_)
      n <- private$size_
      if (requireNamespace("lhs", quietly = TRUE)) {
        raw <- lhs::randomLHS(n, dims)
      } else {
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

    build_mlhs = function(dims) {
      # Modified LHS (Hess, Train, Polak 2004): average of multiple LHS draws
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
        nodes_sum <- nodes_sum + stats::qnorm(raw)
      }
      nodes <- nodes_sum / n_reps
      nodes[!is.finite(nodes)] <- 0
      weights <- rep(1 / n, n)
      list(nodes = nodes, weights = weights)
    },

    build_product = function(dims) {
      # Gauss-Hermite product rule
      level <- private$size_
      gh <- gauss_hermite(level)
      # Scale for standard normal: nodes * sqrt(2), weights / sqrt(pi)
      gh_nodes <- gh$nodes * sqrt(2)
      gh_weights <- gh$weights / sqrt(pi)

      if (dims == 1) {
        return(list(nodes = matrix(gh_nodes, ncol = 1), weights = gh_weights))
      }

      # Product rule: expand.grid of all combinations
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
  # Use eigenvalue method for Gauss-Hermite
  i <- seq_len(n - 1)
  b <- sqrt(i / 2)
  cm <- diag(0, n)
  for (k in seq_along(b)) {
    cm[k, k + 1] <- b[k]
    cm[k + 1, k] <- b[k]
  }
  eig <- eigen(cm, symmetric = TRUE)
  nodes <- eig$values
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

#' @title BLP Parameter Manager
#' @description Manages BLP model parameter compression, expansion, bounds, and labels.
#' @keywords internal
#
# Central bookkeeper for the nonlinear parameters of the BLP model (sigma,
# pi, rho). The GMM optimizer works on a flat "theta" vector, but the model
# needs structured matrices. This class handles the mapping between the two
# representations, tracks which elements are free vs. fixed, and stacks
# bound constraints for box-constrained optimizers like L-BFGS-B.
#
# Convention: in the initial sigma/pi/rho supplied by the user,
#   0   = parameter is FIXED at zero (not estimated),
#   any non-zero value = parameter is FREE with that value as starting point.
# This follows pyblp's convention and makes specification concise.
BLPParameters <- R6::R6Class("BLPParameters",
  public = list(
    K2 = 0L,
    D = 0L,
    H = 0L,
    K1 = 0L,
    K3 = 0L,

    initialize = function(sigma = NULL, pi = NULL, rho = NULL,
                          beta = NULL, gamma = NULL,
                          sigma_bounds = NULL, pi_bounds = NULL,
                          rho_bounds = NULL, beta_bounds = NULL,
                          gamma_bounds = NULL, rc_types = NULL) {
      # Sigma: K2 x K2 lower-triangular Cholesky factor of the random
      # coefficient covariance matrix Omega = sigma %*% t(sigma).
      # Only the lower triangle (including diagonal) can be free; the upper
      # triangle is always zero by the Cholesky structure. This
      # parameterization guarantees that Omega is positive semi-definite
      # by construction, avoiding constrained optimization on the PSD cone.
      # Diagonal elements are bounded below at 0 (standard deviations are
      # non-negative); off-diagonal elements are unconstrained (they
      # capture taste correlations and can be negative).
      if (!is.null(sigma)) {
        sigma <- as.matrix(sigma)
        self$K2 <- nrow(sigma)
        private$sigma_free_ <- (sigma != 0) & lower.tri(sigma, diag = TRUE)
        private$sigma_values_ <- sigma
        if (is.null(sigma_bounds)) {
          private$sigma_lb_ <- matrix(-Inf, self$K2, self$K2)
          private$sigma_ub_ <- matrix(Inf, self$K2, self$K2)
          diag(private$sigma_lb_) <- 0  # diagonal elements non-negative
        } else {
          private$sigma_lb_ <- sigma_bounds[[1]]
          private$sigma_ub_ <- sigma_bounds[[2]]
        }
      }

      # Pi: K2 x D matrix of demographic interaction coefficients.
      # Entry pi[k,d] governs how demographic variable d shifts the random
      # coefficient on characteristic k. Together with sigma, pi determines
      # the full individual-level taste deviation:
      #   mu_ij = X2_j' (sigma * nu_i + pi * d_i)
      # where nu_i ~ N(0,I) are unobserved taste shocks and d_i are observed
      # demographics (e.g., income, age). Pi is unconstrained by default.
      if (!is.null(pi)) {
        pi <- as.matrix(pi)
        self$D <- ncol(pi)
        private$pi_free_ <- (pi != 0)
        private$pi_values_ <- pi
        if (is.null(pi_bounds)) {
          private$pi_lb_ <- matrix(-Inf, self$K2, self$D)
          private$pi_ub_ <- matrix(Inf, self$K2, self$D)
        } else {
          private$pi_lb_ <- pi_bounds[[1]]
          private$pi_ub_ <- pi_bounds[[2]]
        }
      }

      # Rho: nesting parameters for the nested logit component (if used).
      # rho_h in [0, 1) measures within-nest correlation in the GEV error
      # structure. rho = 0 reduces to standard logit; rho -> 1 means products
      # in the same nest are nearly perfect substitutes. Bounded strictly
      # below 1 (default upper = 0.99) to ensure well-defined choice probs.
      if (!is.null(rho)) {
        rho <- as.numeric(rho)
        self$H <- length(rho)
        private$rho_free_ <- (rho != 0)
        private$rho_values_ <- rho
        if (is.null(rho_bounds)) {
          private$rho_lb_ <- rep(0, self$H)
          private$rho_ub_ <- rep(0.99, self$H)
        } else {
          private$rho_lb_ <- rho_bounds[[1]]
          private$rho_ub_ <- rho_bounds[[2]]
        }
      }

      # Beta: linear demand parameters. By default these are "concentrated
      # out" of the GMM objective (NA = concentrated), meaning they are
      # recovered analytically by IV/2SLS for each trial value of theta.
      # This dramatically reduces the dimensionality of the nonlinear
      # optimization problem from (K1 + K2*(K2+1)/2 + K2*D) to just the
      # nonlinear parameters, which is the standard BLP approach.
      if (!is.null(beta)) {
        beta <- as.numeric(beta)
        self$K1 <- length(beta)
        private$beta_free_ <- !is.na(beta) & (beta != 0)
        private$beta_concentrated_ <- is.na(beta)
        private$beta_values_ <- beta
      }

      # Gamma: supply-side parameters, also concentrated out by default
      # via IV regression of recovered marginal costs on cost shifters X3.
      if (!is.null(gamma)) {
        gamma <- as.numeric(gamma)
        self$K3 <- length(gamma)
        private$gamma_free_ <- !is.na(gamma) & (gamma != 0)
        private$gamma_concentrated_ <- is.na(gamma)
        private$gamma_values_ <- gamma
      }

      private$rc_types_ <- rc_types %||% rep("linear", self$K2)
      private$build_labels_()
    },

    # Compress: pack the free elements of sigma, pi, rho into a single flat
    # vector "theta" for the optimizer. The order is always:
    #   [sigma free entries | pi free entries | rho free entries]
    # Only elements marked as free (non-zero in the initial specification)
    # are included; fixed elements stay at zero and are never touched by
    # the optimizer. This is the vector that L-BFGS-B or other optimizers
    # see as their decision variable.
    compress = function(sigma = NULL, pi = NULL, rho = NULL,
                        beta = NULL, gamma = NULL) {
      theta <- numeric(0)
      if (!is.null(sigma %||% private$sigma_values_)) {
        s <- sigma %||% private$sigma_values_
        theta <- c(theta, s[private$sigma_free_])
      }
      if (!is.null(pi %||% private$pi_values_)) {
        p <- pi %||% private$pi_values_
        theta <- c(theta, p[private$pi_free_])
      }
      if (!is.null(rho %||% private$rho_values_)) {
        r <- rho %||% private$rho_values_
        theta <- c(theta, r[private$rho_free_])
      }
      theta
    },

    # Expand: unpack the flat theta vector back into structured matrices.
    # This is the inverse of compress(): it takes the optimizer's current
    # candidate theta, distributes values into the correct positions of
    # sigma (lower-triangular), pi (full matrix), and rho (vector), while
    # keeping fixed elements at zero. Called at every GMM objective
    # evaluation to reconstruct the model parameters from the optimizer's
    # state. The upper triangle of sigma is explicitly zeroed to maintain
    # the Cholesky structure.
    expand = function(theta) {
      idx <- 1L

      sigma <- NULL
      if (self$K2 > 0 && !is.null(private$sigma_values_)) {
        sigma <- private$sigma_values_
        sigma[sigma != 0] <- 0  # reset free elements
        n_sigma <- sum(private$sigma_free_)
        if (n_sigma > 0) {
          sigma[private$sigma_free_] <- theta[idx:(idx + n_sigma - 1)]
          idx <- idx + n_sigma
        }
        # Ensure lower triangular
        sigma[upper.tri(sigma)] <- 0
      }

      pi_mat <- NULL
      if (self$D > 0 && !is.null(private$pi_values_)) {
        pi_mat <- private$pi_values_
        pi_mat[pi_mat != 0] <- 0
        n_pi <- sum(private$pi_free_)
        if (n_pi > 0) {
          pi_mat[private$pi_free_] <- theta[idx:(idx + n_pi - 1)]
          idx <- idx + n_pi
        }
      }

      rho <- NULL
      if (self$H > 0 && !is.null(private$rho_values_)) {
        rho <- private$rho_values_
        rho[rho != 0] <- 0
        n_rho <- sum(private$rho_free_)
        if (n_rho > 0) {
          rho[private$rho_free_] <- theta[idx:(idx + n_rho - 1)]
          idx <- idx + n_rho
        }
      }

      list(sigma = sigma, pi = pi_mat, rho = rho)
    },

    # Stack the lower and upper bounds for all free parameters into vectors
    # aligned with the theta vector produced by compress(). The optimizer
    # (e.g., L-BFGS-B in optim()) uses these directly as box constraints.
    # The order matches compress(): sigma bounds first, then pi, then rho.
    # Typical constraints: sigma diagonal >= 0, rho in [0, 0.99], pi
    # unconstrained (-Inf, Inf).
    get_bounds = function() {
      lower <- upper <- numeric(0)
      if (!is.null(private$sigma_lb_)) {
        lower <- c(lower, private$sigma_lb_[private$sigma_free_])
        upper <- c(upper, private$sigma_ub_[private$sigma_free_])
      }
      if (!is.null(private$pi_lb_)) {
        lower <- c(lower, private$pi_lb_[private$pi_free_])
        upper <- c(upper, private$pi_ub_[private$pi_free_])
      }
      if (!is.null(private$rho_lb_)) {
        lower <- c(lower, private$rho_lb_[private$rho_free_])
        upper <- c(upper, private$rho_ub_[private$rho_free_])
      }
      list(lower = lower, upper = upper)
    },

    get_labels = function() private$labels_,
    get_sigma_free = function() private$sigma_free_,
    get_pi_free = function() private$pi_free_,

    n_free = function() {
      n <- 0L
      if (!is.null(private$sigma_free_)) n <- n + sum(private$sigma_free_)
      if (!is.null(private$pi_free_)) n <- n + sum(private$pi_free_)
      if (!is.null(private$rho_free_)) n <- n + sum(private$rho_free_)
      n
    },

    print = function(...) {
      cat(sprintf("BLPParameters: %d free nonlinear parameters\n", self$n_free()))
      invisible(self)
    }
  ),
  private = list(
    sigma_free_ = NULL, sigma_values_ = NULL,
    sigma_lb_ = NULL, sigma_ub_ = NULL,
    pi_free_ = NULL, pi_values_ = NULL,
    pi_lb_ = NULL, pi_ub_ = NULL,
    rho_free_ = NULL, rho_values_ = NULL,
    rho_lb_ = NULL, rho_ub_ = NULL,
    beta_free_ = NULL, beta_values_ = NULL, beta_concentrated_ = NULL,
    gamma_free_ = NULL, gamma_values_ = NULL, gamma_concentrated_ = NULL,
    rc_types_ = NULL,
    labels_ = NULL,

    # Build human-readable labels for each free parameter, matching the
    # order used by compress(). These labels appear in summary_table()
    # output and diagnostic messages, making it easy to identify which
    # element of theta corresponds to which structural parameter. Sigma
    # labels traverse column-major within the lower triangle (matching
    # R's matrix indexing convention).
    build_labels_ = function() {
      labs <- character(0)
      if (!is.null(private$sigma_free_)) {
        for (j in seq_len(self$K2)) {
          for (i in j:self$K2) {
            if (private$sigma_free_[i, j]) {
              labs <- c(labs, sprintf("sigma[%d,%d]", i, j))
            }
          }
        }
      }
      if (!is.null(private$pi_free_)) {
        for (j in seq_len(self$D)) {
          for (i in seq_len(self$K2)) {
            if (private$pi_free_[i, j]) {
              labs <- c(labs, sprintf("pi[%d,%d]", i, j))
            }
          }
        }
      }
      if (!is.null(private$rho_free_)) {
        for (i in seq_len(self$H)) {
          if (private$rho_free_[i]) {
            labs <- c(labs, sprintf("rho[%d]", i))
          }
        }
      }
      private$labels_ <- labs
    }
  )
)

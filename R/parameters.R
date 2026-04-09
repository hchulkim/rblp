#' @title BLP Parameter Manager
#' @description Manages BLP model parameter compression, expansion, bounds, and labels.
#' @keywords internal
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
      # Sigma: K2 x K2 lower triangular. 0 = fixed, non-zero = free starting value
      if (!is.null(sigma)) {
        sigma <- as.matrix(sigma)
        self$K2 <- nrow(sigma)
        private$sigma_free_ <- (sigma != 0) & lower.tri(sigma, diag = TRUE)
        private$sigma_values_ <- sigma
        # Default bounds: diagonal >= 0
        if (is.null(sigma_bounds)) {
          private$sigma_lb_ <- matrix(-Inf, self$K2, self$K2)
          private$sigma_ub_ <- matrix(Inf, self$K2, self$K2)
          diag(private$sigma_lb_) <- 0  # diagonal elements non-negative
        } else {
          private$sigma_lb_ <- sigma_bounds[[1]]
          private$sigma_ub_ <- sigma_bounds[[2]]
        }
      }

      # Pi: K2 x D. 0 = fixed, non-zero = free
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

      # Rho: nesting parameters. 0 = fixed, non-zero = free
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

      # Beta: concentrated out by default (NA = concentrated)
      if (!is.null(beta)) {
        beta <- as.numeric(beta)
        self$K1 <- length(beta)
        private$beta_free_ <- !is.na(beta) & (beta != 0)
        private$beta_concentrated_ <- is.na(beta)
        private$beta_values_ <- beta
      }

      # Gamma: concentrated out by default
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

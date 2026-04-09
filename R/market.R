#' @title BLP Market Computations
#' @description Market-level computations for BLP demand estimation.
#'   All probability, share, delta, markup, and Jacobian computations happen here.
#' @keywords internal
BLPMarket <- R6::R6Class("BLPMarket",
  public = list(
    J = 0L,
    I = 0L,

    initialize = function(products, agents, sigma = NULL, pi = NULL,
                          rho = NULL, rc_types = NULL, epsilon_scale = 1.0,
                          costs_type = "linear",
                          sigma_free = NULL, pi_free = NULL) {
      private$products_ <- products
      private$agents_ <- agents
      private$sigma_ <- sigma
      private$pi_ <- pi
      private$rho_ <- if (!is.null(rho) && any(rho != 0)) rho else NULL
      private$rc_types_ <- rc_types %||% rep("linear",
                              if (!is.null(sigma)) nrow(sigma) else 0L)
      private$epsilon_scale_ <- epsilon_scale
      private$costs_type_ <- costs_type
      private$sigma_free_ <- sigma_free
      private$pi_free_ <- pi_free

      self$J <- nrow(products$X2 %||% products$X1)
      self$I <- length(agents$weights)

      # Build nesting group indices if rho is used
      if (!is.null(private$rho_) && !is.null(products$nesting_ids)) {
        nids <- products$nesting_ids
        private$nesting_ids_ <- nids
        private$nesting_groups_ <- split(seq_len(self$J), nids)
      }
    },

    compute_random_coefficients = function() {
      if (is.null(private$sigma_) || self$I == 0) return(NULL)
      K2 <- nrow(private$sigma_)
      nodes <- private$agents_$nodes  # I x K2
      # coeffs = sigma %*% t(nodes) -> K2 x I
      coeffs <- private$sigma_ %*% t(nodes)
      if (!is.null(private$pi_) && !is.null(private$agents_$demographics)) {
        demos <- private$agents_$demographics  # I x D
        coeffs <- coeffs + private$pi_ %*% t(demos)
      }
      # Apply transforms
      for (k in seq_len(K2)) {
        if (private$rc_types_[k] == "log") {
          coeffs[k, ] <- exp(coeffs[k, ])
        } else if (private$rc_types_[k] == "logit") {
          coeffs[k, ] <- 1 / (1 + exp(-coeffs[k, ]))
        }
      }
      coeffs  # K2 x I
    },

    compute_mu = function(coefficients = NULL) {
      if (is.null(private$products_$X2)) return(matrix(0, self$J, self$I))
      if (is.null(coefficients)) coefficients <- self$compute_random_coefficients()
      if (is.null(coefficients)) return(matrix(0, self$J, self$I))
      private$products_$X2 %*% coefficients  # J x I
    },

    compute_probabilities = function(delta, mu = NULL) {
      if (is.null(mu)) mu <- matrix(0, self$J, self$I)
      eps <- private$epsilon_scale_

      # V = (delta + mu) / epsilon_scale
      V <- sweep(mu, 1, delta, "+") / eps  # J x I

      if (is.null(private$rho_) || is.null(private$nesting_ids_)) {
        return(private$compute_logit_probs(V))
      }
      private$compute_nested_logit_probs(V)
    },

    compute_shares = function(probabilities, weights = NULL) {
      if (is.null(weights)) weights <- private$agents_$weights
      as.numeric(probabilities %*% weights)
    },

    compute_delta = function(initial_delta, iteration, fp_type = "safe_linear") {
      # No agents means logit: delta is already analytic, skip contraction
      if (self$I == 0) {
        return(list(delta = initial_delta, converged = TRUE,
                    iterations = 0L, evaluations = 0L))
      }

      s_obs <- private$products_$shares
      rho_val <- if (!is.null(private$rho_)) private$rho_[1] else 0

      # Precompute mu (does NOT depend on delta, only on sigma/pi/nodes)
      mu <- self$compute_mu()

      if (fp_type %in% c("safe_linear", "linear")) {
        contraction <- function(d) {
          prob <- self$compute_probabilities(d, mu)
          if (is.list(prob)) prob <- prob$probabilities
          s_pred <- self$compute_shares(prob)
          s_pred <- pmax(s_pred, 1e-300)
          d_new <- d + (1 - rho_val) * (log(s_obs) - log(s_pred))
          d_new
        }
      } else {
        # Nonlinear contraction on exp scale
        contraction <- function(d) {
          prob <- self$compute_probabilities(d, mu)
          if (is.list(prob)) prob <- prob$probabilities
          s_pred <- self$compute_shares(prob)
          s_pred <- pmax(s_pred, 1e-300)
          log(exp(d) * (s_obs / s_pred)^(1 - rho_val))
        }
      }

      result <- iteration$iterate(initial_delta, contraction)
      list(delta = result$values, converged = result$converged,
           iterations = result$iterations, evaluations = result$evaluations)
    },

    compute_shares_by_xi_jacobian = function(probabilities, conditionals = NULL) {
      # ds/d_delta: J x J matrix
      w <- private$agents_$weights
      P <- if (is.list(probabilities)) probabilities$probabilities else probabilities
      s <- self$compute_shares(P, w)

      if (is.null(private$rho_)) {
        # Standard logit: ds_j/dd_k = sum_i w_i P_ji (1{j=k} - P_ki)
        # = diag(s) - P %*% diag(w) %*% t(P)
        return(diag(s, nrow = self$J) - P %*% (w * t(P)))
      }

      # Nested logit correction
      rho <- private$rho_[1]
      C <- if (!is.null(conditionals)) conditionals else {
        if (is.list(probabilities)) probabilities$conditionals else NULL
      }

      Jacob <- matrix(0, self$J, self$J)
      for (i in seq_len(self$I)) {
        p_i <- P[, i]
        c_i <- if (!is.null(C)) C[, i] else p_i
        for (j in seq_len(self$J)) {
          for (k in seq_len(self$J)) {
            same_group <- !is.null(private$nesting_ids_) &&
              private$nesting_ids_[j] == private$nesting_ids_[k]
            term1 <- p_i[j] * ((j == k) - p_i[k]) / (1 - rho)
            if (same_group && j != k) {
              term1 <- term1 - rho / (1 - rho) * p_i[j] * c_i[k]
            }
            if (j == k) {
              term1 <- term1 + rho / (1 - rho) * p_i[j] * (1 - c_i[j])
            }
            Jacob[j, k] <- Jacob[j, k] + w[i] * term1
          }
        }
      }
      Jacob
    },

    compute_shares_by_theta_jacobian = function(probabilities, conditionals = NULL,
                                                 sigma, pi_mat = NULL, rho = NULL) {
      P <- if (is.list(probabilities)) probabilities$probabilities else probabilities
      w <- private$agents_$weights
      X2 <- private$products_$X2
      nodes <- private$agents_$nodes
      demos <- private$agents_$demographics
      K2 <- ncol(X2)

      # Use free masks if available, otherwise infer from values
      sigma_free <- private$sigma_free_
      if (is.null(sigma_free)) sigma_free <- sigma != 0 & lower.tri(sigma, diag = TRUE)
      pi_free <- private$pi_free_
      if (is.null(pi_free) && !is.null(pi_mat)) pi_free <- pi_mat != 0

      n_sigma <- sum(sigma_free)
      n_pi <- if (!is.null(pi_free)) sum(pi_free) else 0L
      n_rho <- if (!is.null(rho)) sum(rho != 0) else 0L
      n_theta <- n_sigma + n_pi + n_rho
      if (n_theta == 0) return(matrix(0, self$J, 0))

      ds_dtheta <- matrix(0, self$J, n_theta)

      # Vectorized: eliminate loop over agents
      # For sigma[row, col]:
      #   ds_j = sum_i w_i * P_ji * nodes_i_col * (X2_j_row - sum_l P_li * X2_l_row)
      #        = X2_j_row * C_j - D_j
      # where C = P %*% (w * nodes[,col]), D = P %*% (w * nodes[,col] * B)
      # and B = t(X2[,row]) %*% P (I-vector of weighted X2 means)

      col_idx <- 0L

      # Sigma derivatives (lower triangular elements)
      for (col in seq_len(K2)) {
        node_col <- nodes[, col]  # I-vector
        wn <- w * node_col        # I-vector: w_i * node_i_col

        for (row in col:K2) {
          if (!sigma_free[row, col]) next
          col_idx <- col_idx + 1L

          x2_row <- X2[, row]  # J-vector
          # B_i = sum_j P_ji * X2_j_row  (I-vector)
          B <- as.numeric(crossprod(x2_row, P))
          # C_j = sum_i w_i * node_i_col * P_ji  (J-vector)
          C <- as.numeric(P %*% wn)
          # D_j = sum_i w_i * node_i_col * B_i * P_ji  (J-vector)
          D <- as.numeric(P %*% (wn * B))

          ds_dtheta[, col_idx] <- x2_row * C - D
        }
      }

      # Pi derivatives
      if (!is.null(pi_free) && !is.null(demos)) {
        D_dim <- ncol(demos)
        for (d in seq_len(D_dim)) {
          demo_col <- demos[, d]  # I-vector
          wd <- w * demo_col      # I-vector

          for (row in seq_len(K2)) {
            if (!pi_free[row, d]) next
            col_idx <- col_idx + 1L

            x2_row <- X2[, row]
            B <- as.numeric(crossprod(x2_row, P))
            C <- as.numeric(P %*% wd)
            D_vec <- as.numeric(P %*% (wd * B))

            ds_dtheta[, col_idx] <- x2_row * C - D_vec
          }
        }
      }

      ds_dtheta
    },

    compute_xi_by_theta_jacobian = function(probabilities, conditionals = NULL) {
      sigma <- private$sigma_
      pi_mat <- private$pi_
      rho <- private$rho_
      if (is.null(sigma) && is.null(rho)) return(NULL)

      ds_ddelta <- self$compute_shares_by_xi_jacobian(probabilities, conditionals)
      ds_dtheta <- self$compute_shares_by_theta_jacobian(
        probabilities, conditionals, sigma, pi_mat, rho)

      if (ncol(ds_dtheta) == 0) return(NULL)
      # IFT: d_xi/d_theta = -(ds/d_delta)^{-1} %*% ds/d_theta
      -approximately_solve(ds_ddelta, ds_dtheta)
    },

    compute_eta = function(probabilities, ownership = NULL) {
      P <- if (is.list(probabilities)) probabilities$probabilities else probabilities
      w <- private$agents_$weights
      s <- self$compute_shares(P, w)

      if (is.null(ownership)) ownership <- private$products_$ownership

      alpha <- private$get_alpha()

      # Vectorized: Delta = diag(A) - B (same structure as elasticities)
      v <- w * alpha
      Pv <- sweep(P, 2, v, "*")
      A <- rowSums(Pv)
      B <- tcrossprod(Pv, P)
      Delta <- diag(A, nrow = self$J) - B

      # Omega = ownership * Delta (element-wise)
      Omega <- ownership * Delta

      # eta = -Omega^{-1} %*% s
      -approximately_solve(Omega, s)
    },

    compute_costs = function(probabilities, prices, ownership = NULL) {
      eta <- self$compute_eta(probabilities, ownership)
      prices - eta
    },

    compute_equilibrium_prices = function(costs, iteration, ownership = NULL,
                                           prices_init = NULL) {
      if (is.null(prices_init)) prices_init <- costs * 1.5

      contraction <- function(p) {
        # Update mu based on new prices
        delta_adj <- private$products_$delta_base
        if (!is.null(private$products_$price_col_x1)) {
          delta_adj <- delta_adj + private$products_$beta_price * (p - private$products_$prices)
        }

        prob <- self$compute_probabilities(delta_adj)
        P <- if (is.list(prob)) prob$probabilities else prob
        eta <- self$compute_eta(P, ownership)
        costs + eta
      }

      result <- iteration$iterate(prices_init, contraction)
      result$values
    },

    compute_elasticities = function(probabilities, prices = NULL, name = "prices") {
      P <- if (is.list(probabilities)) probabilities$probabilities else probabilities
      w <- private$agents_$weights
      s <- self$compute_shares(P, w)
      if (is.null(prices)) prices <- private$products_$prices

      alpha <- private$get_alpha()

      # Vectorized: ds/dp = diag(A) - B
      # where A_j = sum_i v_i * P_ji, B_jk = sum_i v_i * P_ji * P_ki, v = w * alpha
      v <- w * alpha
      Pv <- sweep(P, 2, v, "*")  # J x I: P_ji * v_i
      A <- rowSums(Pv)           # J-vector
      B <- tcrossprod(Pv, P)     # J x J: sum_i v_i * P_ji * P_ki

      dsdp <- diag(A, nrow = self$J) - B

      # E_jk = dsdp_jk * prices[k] / s[j]
      E <- sweep(sweep(dsdp, 2, prices, "*"), 1, s, "/")
      E
    },

    compute_diversion_ratios = function(probabilities, name = "prices") {
      P <- if (is.list(probabilities)) probabilities$probabilities else probabilities
      w <- private$agents_$weights
      alpha <- private$get_alpha()

      # Vectorized ds/dp
      v <- w * alpha
      Pv <- sweep(P, 2, v, "*")
      A <- rowSums(Pv)
      B <- tcrossprod(Pv, P)
      dsdp <- diag(A, nrow = self$J) - B

      # D_jk = -dsdp[k,j] / dsdp[j,j]
      D <- -t(dsdp) / diag(dsdp)
      diag(D) <- 0
      D
    },

    compute_consumer_surplus = function(probabilities, delta, mu) {
      # CS = sum_i w_i * log(1 + sum_j exp(V_ji)) / (-alpha_i)
      w <- private$agents_$weights
      alpha <- private$get_alpha()
      V <- sweep(mu, 1, delta, "+")

      cs <- 0
      for (i in seq_len(self$I)) {
        V_i <- V[, i]
        V_max <- max(V_i)
        log_denom <- V_max + log(exp(-V_max) + sum(exp(V_i - V_max)))
        a_i <- alpha[i]
        if (abs(a_i) > 1e-300) {
          cs <- cs + w[i] * log_denom / (-a_i)
        }
      }
      cs
    },

    compute_hhi = function(shares = NULL, firm_ids = NULL) {
      if (is.null(shares)) shares <- private$products_$shares
      if (is.null(firm_ids)) firm_ids <- private$products_$firm_ids
      firm_shares <- tapply(shares, firm_ids, sum)
      sum(firm_shares^2) * 10000
    },

    compute_markups = function(prices, costs) {
      (prices - costs) / prices
    },

    compute_profits = function(prices, shares, costs) {
      (prices - costs) * shares
    }
  ),
  private = list(
    products_ = NULL,
    agents_ = NULL,
    sigma_ = NULL,
    pi_ = NULL,
    rho_ = NULL,
    rc_types_ = NULL,
    epsilon_scale_ = 1.0,
    costs_type_ = "linear",
    nesting_ids_ = NULL,
    nesting_groups_ = NULL,
    sigma_free_ = NULL,
    pi_free_ = NULL,

    compute_logit_probs = function(V) {
      # Standard logit probabilities with numerical stability
      # Include outside good (utility=0) in max for log-sum-exp trick
      V_max <- pmax(apply(V, 2, max), 0)
      exp_V <- exp(sweep(V, 2, V_max))       # exp(V_j - V_max)
      exp_outside <- exp(-V_max)              # exp(0 - V_max) = outside good
      denom <- exp_outside + colSums(exp_V)   # I vector
      probs <- sweep(exp_V, 2, denom, "/")
      list(probabilities = probs, conditionals = NULL)
    },

    compute_nested_logit_probs = function(V) {
      rho <- private$rho_[1]
      groups <- private$nesting_groups_
      J <- self$J
      I <- self$I

      probs <- matrix(0, J, I)
      conditionals <- matrix(0, J, I)

      for (i in seq_len(I)) {
        V_i <- V[, i]
        # Compute within-group conditionals and inclusive values
        n_groups <- length(groups)
        IV <- numeric(n_groups)  # inclusive values

        for (g in seq_along(groups)) {
          idx <- groups[[g]]
          V_g <- V_i[idx] / (1 - rho)
          V_g_max <- max(V_g)
          exp_V_g <- exp(V_g - V_g_max)
          sum_exp <- sum(exp_V_g)
          conditionals[idx, i] <- exp_V_g / sum_exp
          IV[g] <- V_g_max + log(sum_exp)
        }

        # Across-group probabilities
        IV_scaled <- IV * (1 - rho)
        IV_max <- max(IV_scaled)
        exp_IV <- exp(IV_scaled - IV_max)
        denom <- exp(-IV_max) + sum(exp_IV)

        for (g in seq_along(groups)) {
          idx <- groups[[g]]
          group_prob <- exp_IV[g] / denom
          probs[idx, i] <- conditionals[idx, i] * group_prob
        }
      }

      list(probabilities = probs, conditionals = conditionals)
    },

    get_alpha = function() {
      # Get individual-specific price coefficient for each agent
      # alpha_i = beta_price + sigma_price * node_i + pi_price * demo_i
      if (is.null(private$products_$X2)) {
        # No random coefficients: alpha is constant
        return(rep(private$products_$beta_price %||% -1, self$I))
      }

      K2 <- nrow(private$sigma_)
      price_col <- private$products_$price_col_x2
      if (is.null(price_col) || price_col == 0) {
        return(rep(private$products_$beta_price %||% -1, self$I))
      }

      nodes <- private$agents_$nodes
      alpha_base <- private$products_$beta_price %||% 0

      alpha <- rep(alpha_base, self$I)
      if (!is.null(private$sigma_)) {
        alpha <- alpha + as.numeric(private$sigma_[price_col, ] %*% t(nodes))
      }
      if (!is.null(private$pi_) && !is.null(private$agents_$demographics)) {
        alpha <- alpha + as.numeric(private$pi_[price_col, ] %*%
                                     t(private$agents_$demographics))
      }
      alpha
    }
  )
)

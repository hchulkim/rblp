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
      # Compute individual-specific deviations from mean taste: beta_ik - beta_k.
      # In the BLP model, consumer i's coefficient on characteristic k is:
      #   beta_ik = beta_k + sigma_k * nu_ik + pi_k' * D_i
      # where nu_ik are integration nodes (draws from the mixing distribution),
      # sigma is the Cholesky factor allowing correlated tastes, and pi*D
      # captures observed demographic heterogeneity. The matrix product
      # sigma %*% t(nodes) gives the K2 x I matrix of taste deviations.
      coeffs <- private$sigma_ %*% t(nodes)
      if (!is.null(private$pi_) && !is.null(private$agents_$demographics)) {
        demos <- private$agents_$demographics  # I x D
        coeffs <- coeffs + private$pi_ %*% t(demos)
      }
      # Optional nonlinear transforms on random coefficients: "log" ensures
      # positivity (e.g., for a price coefficient that must be negative when
      # negated), "logit" bounds coefficients to (0,1).
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
      # Compute mu_ij = X2_j' * (sigma*nu_i + pi*D_i), the individual-specific
      # portion of utility that varies across consumers. This is the
      # heterogeneous part of indirect utility: U_ij = delta_j + mu_ij + eps_ij.
      # mu does NOT depend on delta, only on the nonlinear parameters and
      # the integration nodes/demographics, so it can be precomputed once
      # per contraction mapping.
      if (is.null(private$products_$X2)) return(matrix(0, self$J, self$I))
      if (is.null(coefficients)) coefficients <- self$compute_random_coefficients()
      if (is.null(coefficients)) return(matrix(0, self$J, self$I))
      private$products_$X2 %*% coefficients  # J x I
    },

    compute_probabilities = function(delta, mu = NULL) {
      # Compute choice probabilities P_ij = Pr(consumer i chooses product j).
      # The total utility is V_ij = (delta_j + mu_ij) / epsilon_scale.
      # For the standard logit kernel, P_ij = exp(V_ij) / (1 + sum_k exp(V_ik)),
      # where the 1 in the denominator represents the outside good (utility = 0).
      # Epsilon_scale allows rescaling the variance of the Type-I EV error.
      if (is.null(mu)) mu <- matrix(0, self$J, self$I)
      eps <- private$epsilon_scale_

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
      # For a pure logit model (no random coefficients / no agents), delta
      # has the closed-form inverse: delta = log(s) - log(s0). No contraction needed.
      if (self$I == 0) {
        return(list(delta = initial_delta, converged = TRUE,
                    iterations = 0L, evaluations = 0L))
      }

      s_obs <- private$products_$shares
      rho_val <- if (!is.null(private$rho_)) private$rho_[1] else 0

      # mu = X2*(sigma*nu + pi*D) is precomputed because it only depends on
      # the current nonlinear parameters and integration nodes, not on delta.
      # This avoids redundant matrix multiplications inside each contraction step.
      mu <- self$compute_mu()

      # BLP contraction mapping (Berry 1994):
      #   delta^{h+1} = delta^h + (1-rho) * [log(s_obs) - log(s_pred(delta^h))]
      # This exploits the fact that s(delta) is a smooth, monotone function of
      # delta: if predicted shares are too low (s_pred < s_obs), delta is
      # increased to make products more attractive, and vice versa. The (1-rho)
      # factor accounts for the nesting parameter in the nested logit case.
      # Berry, Levinsohn & Pakes (1995) proved this is a contraction mapping
      # under standard conditions, guaranteeing convergence to the unique delta
      # that rationalizes the observed shares for any given theta.
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
        # Equivalent nonlinear form on the exp(delta) scale:
        # exp(delta^{h+1}) = exp(delta^h) * (s_obs/s_pred)^{(1-rho)}.
        # This can be more numerically stable for extreme share values.
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
      # Compute the J x J Jacobian ds/d_delta, which measures how predicted
      # market shares respond to changes in mean utilities. This matrix is
      # central to both the implicit function theorem (for d_xi/d_theta) and
      # the markup equation. For the mixed logit, the (j,k) entry is:
      #   ds_j/dd_k = sum_i w_i * P_ij * (1{j=k} - P_ik)
      # The diagonal is positive (own-effect) and off-diagonals are negative
      # (substitution). In matrix form: ds/d_delta = diag(s) - P*diag(w)*P'.
      w <- private$agents_$weights
      P <- if (is.list(probabilities)) probabilities$probabilities else probabilities
      s <- self$compute_shares(P, w)

      if (is.null(private$rho_)) {
        return(diag(s, nrow = self$J) - P %*% (w * t(P)))
      }

      # Nested logit share Jacobian ds/d_delta. The nesting structure introduces
      # additional terms relative to the standard logit formula. For products
      # j,k in the same nest, substitution is amplified by the within-nest
      # correlation rho: a change in delta_k shifts conditional probabilities
      # within the nest. The (j,k) entry has three components:
      # (1) standard logit term scaled by 1/(1-rho),
      # (2) within-nest cross-substitution correction (-rho/(1-rho)*P_j*C_k),
      # (3) own-product correction (+rho/(1-rho)*P_j*(1-C_j)).
      # C_ij = P(j|g,i) is the within-nest conditional probability.
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
      # Compute the J x n_theta Jacobian ds/d_theta, which measures how
      # predicted shares respond to changes in the nonlinear parameters (sigma, pi, rho).
      # For sigma[row,col], the derivative of share j w.r.t. sigma_{row,col} is:
      #   ds_j/d_sigma_{rc} = sum_i w_i * P_ij * nu_{ic} * (x2_{jr} - sum_k P_ik * x2_{kr})
      # This has the covariance structure: the share of j increases with
      # sigma_{rc} when product j's characteristic r is above the probability-
      # weighted average, scaled by the integration node for dimension c.
      # The analogous formula holds for pi with demographics replacing nodes.
      P <- if (is.list(probabilities)) probabilities$probabilities else probabilities
      w <- private$agents_$weights
      X2 <- private$products_$X2
      nodes <- private$agents_$nodes
      demos <- private$agents_$demographics
      K2 <- ncol(X2)

      # Only compute derivatives for free (estimated) parameters; fixed-at-zero
      # parameters are excluded from the compressed theta vector.
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

      # Vectorized computation avoids explicit loops over I agents.
      # Key identity: ds_j/d_sigma_{rc} = x2_{jr} * C_j - D_j, where
      #   C_j = sum_i w_i * nu_{ic} * P_ij  (weighted node-probability product)
      #   D_j = sum_i w_i * nu_{ic} * B_i * P_ij  (weighted by mean characteristic)
      #   B_i = sum_k P_ik * x2_{kr}  (probability-weighted average of x2 for agent i)

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
      # Apply the implicit function theorem (IFT) to obtain d_xi/d_theta.
      # The BLP share inversion defines delta(theta) implicitly via
      #   s(delta(theta), theta) = s_obs.
      # Total differentiation gives:
      #   (ds/d_delta)(d_delta/d_theta) + (ds/d_theta) = 0
      # Since xi = delta - X1*beta (with beta concentrated out separately),
      # d_xi/d_theta = d_delta/d_theta = -(ds/d_delta)^{-1} (ds/d_theta).
      # This avoids differentiating through the contraction mapping iterations,
      # which would be both expensive and numerically unstable.
      sigma <- private$sigma_
      pi_mat <- private$pi_
      rho <- private$rho_
      if (is.null(sigma) && is.null(rho)) return(NULL)

      ds_ddelta <- self$compute_shares_by_xi_jacobian(probabilities, conditionals)
      ds_dtheta <- self$compute_shares_by_theta_jacobian(
        probabilities, conditionals, sigma, pi_mat, rho)

      if (ncol(ds_dtheta) == 0) return(NULL)
      -approximately_solve(ds_ddelta, ds_dtheta)
    },

    compute_eta = function(probabilities, ownership = NULL) {
      # Compute Nash-Bertrand equilibrium markups eta = p - mc.
      # Under multi-product Bertrand-Nash pricing, each firm f maximizes
      # sum_{j in f} (p_j - mc_j)*s_j(p), leading to the FOC system:
      #   s + (O * Delta) * (p - mc) = 0,
      # where O is the ownership matrix (O_{jk} = 1 if j,k owned by same firm)
      # and Delta_{jk} = ds_k/dp_j is the price-derivative matrix.
      # The markup vector is: eta = p - mc = -Omega^{-1} * s,
      # where Omega = O * Delta (Hadamard product). This inverts the
      # system of first-order conditions to recover implied marginal costs.
      P <- if (is.list(probabilities)) probabilities$probabilities else probabilities
      w <- private$agents_$weights
      s <- self$compute_shares(P, w)

      if (is.null(ownership)) ownership <- private$products_$ownership

      alpha <- private$get_alpha()

      # Delta_{jk} = ds_j/dp_k = sum_i w_i * alpha_i * P_ij * (1{j=k} - P_ik).
      # This is the same structure as ds/d_delta but scaled by alpha_i (the
      # individual price coefficient), since dp enters utility as alpha*dp.
      v <- w * alpha
      Pv <- sweep(P, 2, v, "*")
      A <- rowSums(Pv)
      B <- tcrossprod(Pv, P)
      Delta <- diag(A, nrow = self$J) - B

      # Omega = O * Delta: zero out cross-firm entries so that only
      # own-firm substitution patterns enter the markup equation.
      Omega <- ownership * Delta

      # Solve the FOC: Omega * eta = -s  =>  eta = -Omega^{-1} s.
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
      # Compute the J x J matrix of price elasticities of demand.
      # Own-price elasticity (diagonal): E_{jj} = (ds_j/dp_j)*(p_j/s_j)
      #   = alpha * p_j * (1 - s_j) in the simple logit; in the mixed logit,
      #   it integrates over heterogeneous alpha_i, yielding richer patterns.
      # Cross-price elasticity (off-diagonal): E_{jk} = (ds_j/dp_k)*(p_k/s_j).
      #   In the simple logit, cross-elasticities depend only on shares (IIA),
      #   but mixed logit breaks IIA through taste heterogeneity.
      P <- if (is.list(probabilities)) probabilities$probabilities else probabilities
      w <- private$agents_$weights
      s <- self$compute_shares(P, w)
      if (is.null(prices)) prices <- private$products_$prices

      alpha <- private$get_alpha()

      # ds_j/dp_k = sum_i w_i * alpha_i * P_ij * (1{j=k} - P_ik).
      # Own-price derivatives are negative (alpha < 0, diagonal positive in
      # the P*(1-P) term), cross-price derivatives are positive (products
      # are substitutes). The mixed logit integrates over heterogeneous
      # alpha_i, generating more realistic substitution patterns than logit.
      v <- w * alpha
      Pv <- sweep(P, 2, v, "*")  # J x I: P_ji * v_i
      A <- rowSums(Pv)           # J-vector
      B <- tcrossprod(Pv, P)     # J x J: sum_i v_i * P_ji * P_ki

      dsdp <- diag(A, nrow = self$J) - B

      # Convert derivatives to elasticities: E_{jk} = (ds_j/dp_k) * (p_k / s_j).
      E <- sweep(sweep(dsdp, 2, prices, "*"), 1, s, "/")
      E
    },

    compute_diversion_ratios = function(probabilities, name = "prices") {
      # Diversion ratios D_{jk}: the fraction of product j's lost sales that
      # are captured by product k when j's price increases. Formally,
      # D_{jk} = -(ds_k/dp_j) / (ds_j/dp_j). This is a key antitrust
      # statistic: high diversion between two merging products indicates
      # strong competitive interaction and larger merger price effects.
      P <- if (is.list(probabilities)) probabilities$probabilities else probabilities
      w <- private$agents_$weights
      alpha <- private$get_alpha()

      v <- w * alpha
      Pv <- sweep(P, 2, v, "*")
      A <- rowSums(Pv)
      B <- tcrossprod(Pv, P)
      dsdp <- diag(A, nrow = self$J) - B

      D <- -t(dsdp) / diag(dsdp)
      diag(D) <- 0
      D
    },

    compute_consumer_surplus = function(probabilities, delta, mu) {
      # Small-Rosen expected consumer surplus (Small & Rosen, 1981).
      # For consumer i with Type-I extreme value errors, the expected maximum
      # utility (inclusive value) is the log-sum-exp:
      #   E[max_j U_ij] = log(1 + sum_j exp(V_ij))
      # Converting to dollar terms by dividing by the marginal utility of
      # income (-alpha_i) gives the compensating variation:
      #   CS_i = log(1 + sum_j exp(V_ij)) / (-alpha_i).
      # Aggregate surplus is the population-weighted sum: CS = sum_i w_i * CS_i.
      # The log-sum-exp trick (subtracting V_max) prevents numerical overflow.
      w <- private$agents_$weights
      alpha <- private$get_alpha()
      V <- sweep(mu, 1, delta, "+")

      cs <- 0
      for (i in seq_len(self$I)) {
        V_i <- V[, i]
        V_max <- max(V_i)
        # log(exp(-V_max) + sum(exp(V_i - V_max))) + V_max = log(1 + sum(exp(V_i)))
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
      # Standard multinomial logit choice probabilities:
      #   P_ij = exp(V_ij) / (1 + sum_k exp(V_ik))
      # The "1" in the denominator is the outside good with utility normalized
      # to zero. Direct exponentiation of V can overflow for large utilities,
      # so we use the log-sum-exp numerical trick: subtract V_max (the
      # column-wise maximum, also compared against 0 for the outside good)
      # from all utilities before exponentiating. This is mathematically
      # equivalent but prevents exp() from returning Inf.
      V_max <- pmax(apply(V, 2, max), 0)
      exp_V <- exp(sweep(V, 2, V_max))       # exp(V_j - V_max)
      exp_outside <- exp(-V_max)              # exp(0 - V_max) = outside good
      denom <- exp_outside + colSums(exp_V)   # I vector
      probs <- sweep(exp_V, 2, denom, "/")
      list(probabilities = probs, conditionals = NULL)
    },

    compute_nested_logit_probs = function(V) {
      # Nested logit (Cardell 1997 / McFadden 1978) choice probabilities.
      # Products are partitioned into nests (groups). Consumer i's probability
      # of choosing product j in nest g is:
      #   P_ij = P(j|g,i) * P(g|i),
      # where the within-nest conditional is:
      #   P(j|g,i) = exp(V_ij/(1-rho)) / sum_{k in g} exp(V_ik/(1-rho))
      # and the nest selection probability uses the inclusive value (log-sum):
      #   IV_g = (1-rho) * log(sum_{k in g} exp(V_ik/(1-rho)))
      #   P(g|i) = exp(IV_g) / (1 + sum_h exp(IV_h)).
      # rho in [0,1) is the within-nest correlation of the GEV error.
      # rho=0 collapses to standard logit; rho->1 means perfect within-nest
      # correlation (products in the same nest are near-perfect substitutes).
      rho <- private$rho_[1]
      groups <- private$nesting_groups_
      J <- self$J
      I <- self$I

      probs <- matrix(0, J, I)
      conditionals <- matrix(0, J, I)

      for (i in seq_len(I)) {
        V_i <- V[, i]
        n_groups <- length(groups)
        IV <- numeric(n_groups)

        # Step 1: Within each nest, compute conditional probabilities and
        # the inclusive value (the log-sum-exp of scaled utilities).
        for (g in seq_along(groups)) {
          idx <- groups[[g]]
          V_g <- V_i[idx] / (1 - rho)
          V_g_max <- max(V_g)
          exp_V_g <- exp(V_g - V_g_max)
          sum_exp <- sum(exp_V_g)
          conditionals[idx, i] <- exp_V_g / sum_exp
          # IV_g = log(sum exp(V/(1-rho))), stored in log-sum-exp stable form.
          IV[g] <- V_g_max + log(sum_exp)
        }

        # Step 2: Across-nest probabilities using scaled inclusive values.
        # The nest probability P(g|i) is a logit over inclusive values,
        # with the outside good (utility 0) included in the denominator.
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
      # Recover the individual-specific price coefficient alpha_i for each
      # agent. In the random coefficients model:
      #   alpha_i = beta_price + sigma_price * nu_i + pi_price' * D_i
      # alpha_i is the marginal utility of income for consumer i and is
      # needed for: (1) price elasticities, (2) markup computation (Bertrand
      # FOCs), and (3) converting the inclusive value to dollar-metric
      # consumer surplus (Small-Rosen formula).
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

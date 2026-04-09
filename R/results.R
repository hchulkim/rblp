#' @title BLP Estimation Results
#' @description Stores estimation results and provides post-estimation methods.
#' @export
BLPResults <- R6::R6Class("BLPResults",
  public = list(
    #' @field problem The originating BLPProblem
    problem = NULL,
    #' @field sigma Estimated Cholesky root of random coefficient covariance
    sigma = NULL,
    #' @field pi Estimated demographics interaction matrix
    pi = NULL,
    #' @field rho Estimated nesting parameters
    rho = NULL,
    #' @field beta Estimated demand linear parameters
    beta = NULL,
    #' @field gamma Estimated supply linear parameters
    gamma = NULL,
    #' @field delta Estimated mean utilities
    delta = NULL,
    #' @field xi Demand-side structural error
    xi = NULL,
    #' @field omega Supply-side structural error
    omega = NULL,
    #' @field objective GMM objective value
    objective = NULL,
    #' @field gradient GMM gradient at solution
    gradient = NULL,
    #' @field hessian GMM Hessian at solution
    hessian = NULL,
    #' @field se Standard errors for nonlinear parameters
    se = NULL,
    #' @field parameter_covariances Covariance matrix for nonlinear parameters
    parameter_covariances = NULL,
    #' @field W Final weighting matrix
    W = NULL,
    #' @field optimization_converged Whether optimization converged
    optimization_converged = NULL,
    #' @field optimization_iterations Number of optimization iterations
    optimization_iterations = NULL,
    #' @field optimization_evaluations Number of function evaluations
    optimization_evaluations = NULL,
    #' @field fp_converged Whether all fixed points converged
    fp_converged = NULL,
    #' @field fp_iterations Total fixed-point iterations
    fp_iterations = NULL,
    #' @field method GMM method ("1s" or "2s")
    method = NULL,
    #' @field se_type Standard error type
    se_type = NULL,

    #' @description Create results object
    initialize = function(problem, params, sigma, pi, rho, beta, gamma,
                          delta, xi, omega, objective, gradient, hessian,
                          se, parameter_covariances, W,
                          step_results, optimization_converged,
                          optimization_iterations, optimization_evaluations,
                          fp_converged, fp_iterations, method, se_type,
                          beta_se = NULL, gamma_se = NULL) {
      self$problem <- problem
      private$params_ <- params
      self$sigma <- sigma
      self$pi <- pi
      self$rho <- rho
      self$beta <- beta
      self$gamma <- gamma
      self$delta <- delta
      self$xi <- xi
      self$omega <- omega
      self$objective <- objective
      self$gradient <- gradient
      self$hessian <- hessian
      self$se <- se
      self$parameter_covariances <- parameter_covariances
      self$W <- W
      private$step_results_ <- step_results
      self$optimization_converged <- optimization_converged
      self$optimization_iterations <- optimization_iterations
      self$optimization_evaluations <- optimization_evaluations
      self$fp_converged <- fp_converged
      self$fp_iterations <- fp_iterations
      self$method <- method
      self$se_type <- se_type
      private$beta_se_ <- beta_se
      private$gamma_se_ <- gamma_se
    },

    #' @description Compute own-price elasticities for a specific market
    #' @param market_id Market identifier
    #' @return J x J elasticity matrix
    compute_elasticities = function(market_id = NULL) {
      if (is.null(market_id)) {
        # Compute for all markets
        results <- list()
        for (t in self$problem$unique_market_ids) {
          results[[as.character(t)]] <- private$compute_market_elasticities(t)
        }
        return(results)
      }
      private$compute_market_elasticities(market_id)
    },

    #' @description Compute diversion ratios for a specific market
    #' @param market_id Market identifier
    #' @return J x J diversion ratio matrix
    compute_diversion_ratios = function(market_id = NULL) {
      if (is.null(market_id)) {
        results <- list()
        for (t in self$problem$unique_market_ids) {
          results[[as.character(t)]] <- private$compute_market_diversion(t)
        }
        return(results)
      }
      private$compute_market_diversion(market_id)
    },

    #' @description Extract marginal costs
    #' @return Named list with costs per market, or full vector
    compute_costs = function() {
      costs_vec <- numeric(self$problem$N)
      for (t in self$problem$unique_market_ids) {
        mkt_data <- private$build_market_(t)
        P <- if (is.list(mkt_data$prob)) mkt_data$prob$probabilities else mkt_data$prob
        costs_t <- mkt_data$market$compute_costs(
          P, mkt_data$products$prices, mkt_data$products$ownership
        )
        costs_vec[mkt_data$indices] <- costs_t
      }
      costs_vec
    },

    #' @description Compute markups (p - c) / p
    #' @return Markup vector
    compute_markups = function() {
      costs <- self$compute_costs()
      prices <- self$problem$products$prices
      (prices - costs) / prices
    },

    #' @description Compute consumer surplus per market
    #' @return Named numeric vector of CS per market
    compute_consumer_surplus = function() {
      cs <- numeric(self$problem$T)
      names(cs) <- as.character(self$problem$unique_market_ids)
      for (i in seq_along(self$problem$unique_market_ids)) {
        t <- self$problem$unique_market_ids[i]
        mkt_data <- private$build_market_(t)
        delta_t <- self$delta[mkt_data$indices]
        mu_t <- mkt_data$market$compute_mu()
        cs[i] <- mkt_data$market$compute_consumer_surplus(
          mkt_data$prob, delta_t, mu_t
        )
      }
      cs
    },

    #' @description Compute HHI per market
    #' @return Named numeric vector of HHI per market
    compute_hhi = function() {
      hhi <- numeric(self$problem$T)
      names(hhi) <- as.character(self$problem$unique_market_ids)
      for (i in seq_along(self$problem$unique_market_ids)) {
        t <- self$problem$unique_market_ids[i]
        md <- self$problem$get_market_data(t)
        hhi[i] <- BLPMarket$new(
          products = md$products, agents = md$agents
        )$compute_hhi(md$products$shares, md$products$firm_ids)
      }
      hhi
    },

    #' @description Simulate merger and compute new equilibrium prices
    #' @param new_firm_ids Updated firm ownership (length N)
    #' @param iteration BLPIteration for price equilibrium
    #' @param costs Optional pre-computed costs; computed if NULL
    #' @return List with new_prices, new_shares, new_costs, delta_cs
    compute_merger = function(new_firm_ids, iteration = NULL, costs = NULL) {
      if (is.null(iteration)) iteration <- blp_iteration("simple", list(atol = 1e-12))
      if (is.null(costs)) costs <- self$compute_costs()

      new_prices <- numeric(self$problem$N)
      new_shares <- numeric(self$problem$N)
      cs_pre <- numeric(self$problem$T)
      cs_post <- numeric(self$problem$T)

      for (i in seq_along(self$problem$unique_market_ids)) {
        t <- self$problem$unique_market_ids[i]
        md <- self$problem$get_market_data(t)
        idx <- md$indices

        # Build new ownership matrix
        new_own <- build_ownership_matrix(
          new_firm_ids[idx],
          rep(t, length(idx))
        )

        # Attach needed info
        md$products$delta_base <- self$delta[idx]
        md$products$beta_price <- if (!is.null(self$beta) &&
          !is.null(md$products$price_col_x1)) {
          self$beta[md$products$price_col_x1]
        } else -1

        mkt <- BLPMarket$new(
          products = md$products,
          agents = md$agents,
          sigma = self$sigma,
          pi = self$pi,
          rho = self$rho,
          rc_types = self$problem$rc_types,
          epsilon_scale = self$problem$epsilon_scale,
          costs_type = self$problem$costs_type
        )

        # Pre-merger CS
        mu_t <- mkt$compute_mu()
        prob_pre <- mkt$compute_probabilities(self$delta[idx], mu_t)
        cs_pre[i] <- mkt$compute_consumer_surplus(prob_pre, self$delta[idx], mu_t)

        # Find post-merger equilibrium
        new_p <- mkt$compute_equilibrium_prices(
          costs[idx], iteration, new_own, md$products$prices
        )
        new_prices[idx] <- new_p

        # Compute new delta (adjust for price change)
        delta_post <- self$delta[idx]
        if (!is.null(md$products$price_col_x1) && !is.null(self$beta)) {
          bp <- self$beta[md$products$price_col_x1]
          delta_post <- delta_post + bp * (new_p - md$products$prices)
        }

        prob_post <- mkt$compute_probabilities(delta_post, mu_t)
        P_post <- if (is.list(prob_post)) prob_post$probabilities else prob_post
        new_shares[idx] <- mkt$compute_shares(P_post)
        cs_post[i] <- mkt$compute_consumer_surplus(prob_post, delta_post, mu_t)
      }

      list(
        new_prices = new_prices,
        new_shares = new_shares,
        costs = costs,
        cs_pre = cs_pre,
        cs_post = cs_post,
        delta_cs = cs_post - cs_pre,
        price_change = new_prices - self$problem$products$prices,
        price_change_pct = (new_prices - self$problem$products$prices) /
          self$problem$products$prices * 100
      )
    },

    #' @description Run Hansen's J-test for overidentifying restrictions
    #' @return List with statistic, df, p_value
    run_hansen_test = function() {
      stat <- self$objective
      if (self$method == "2s") {
        # J = N * min objective (at optimal W)
        stat <- self$problem$N * self$objective / self$problem$N  # already scaled
      }
      n_moments <- nrow(self$W)
      n_params <- (private$params_$n_free()) + self$problem$K1
      if (!is.null(self$gamma)) n_params <- n_params + self$problem$K3
      df <- n_moments - n_params
      p_value <- if (df > 0) 1 - stats::pchisq(stat, df) else NA_real_
      list(statistic = stat, df = df, p_value = p_value)
    },

    #' @description Run Wald test for parameter restrictions
    #' @param R Restriction matrix (q x p)
    #' @param r Restriction values (q x 1), default 0
    #' @return List with statistic, df, p_value
    run_wald_test = function(R, r = NULL) {
      if (is.null(r)) r <- rep(0, nrow(R))
      theta <- private$params_$compress(self$sigma, self$pi, self$rho)
      diff <- R %*% theta - r
      V <- self$parameter_covariances
      if (is.null(V)) stop("No parameter covariance matrix available")
      RVR <- R %*% V %*% t(R)
      stat <- as.numeric(crossprod(diff, approximately_solve(RVR, diff)))
      df <- nrow(R)
      p_value <- 1 - stats::pchisq(stat, df)
      list(statistic = stat, df = df, p_value = p_value)
    },

    #' @description Extract sigma squared (Sigma %*% Sigma')
    #' @return K2 x K2 covariance matrix
    sigma_squared = function() {
      if (is.null(self$sigma)) return(NULL)
      self$sigma %*% t(self$sigma)
    },

    #' @description Get a summary table of all estimated parameters
    #' @return Data frame with estimates and standard errors
    summary_table = function() {
      rows <- list()

      # Beta
      if (!is.null(self$beta)) {
        x1_names <- colnames(self$problem$products$X1)
        if (is.null(x1_names)) x1_names <- paste0("beta_", seq_along(self$beta))
        for (k in seq_along(self$beta)) {
          rows[[length(rows) + 1]] <- data.frame(
            parameter = x1_names[k],
            type = "linear (beta)",
            estimate = self$beta[k],
            se = if (!is.null(private$beta_se_)) private$beta_se_[k] else NA_real_,
            stringsAsFactors = FALSE
          )
        }
      }

      # Sigma
      if (!is.null(self$sigma)) {
        labels <- private$params_$get_labels()
        theta_se <- self$se
        idx <- 0L
        K2 <- nrow(self$sigma)
        for (j in seq_len(K2)) {
          for (i in j:K2) {
            if (self$sigma[i, j] != 0 || (i == j)) {
              idx <- idx + 1L
              rows[[length(rows) + 1]] <- data.frame(
                parameter = if (idx <= length(labels)) labels[idx] else
                  sprintf("sigma[%d,%d]", i, j),
                type = "nonlinear (sigma)",
                estimate = self$sigma[i, j],
                se = if (!is.null(theta_se) && idx <= length(theta_se))
                  theta_se[idx] else NA_real_,
                stringsAsFactors = FALSE
              )
            }
          }
        }
      }

      # Pi
      if (!is.null(self$pi)) {
        K2 <- nrow(self$pi)
        D <- ncol(self$pi)
        for (d in seq_len(D)) {
          for (k in seq_len(K2)) {
            if (self$pi[k, d] != 0) {
              rows[[length(rows) + 1]] <- data.frame(
                parameter = sprintf("pi[%d,%d]", k, d),
                type = "demographics (pi)",
                estimate = self$pi[k, d],
                se = NA_real_,  # Pi SE requires more complex extraction
                stringsAsFactors = FALSE
              )
            }
          }
        }
      }

      # Gamma
      if (!is.null(self$gamma)) {
        x3_names <- colnames(self$problem$products$X3)
        if (is.null(x3_names)) x3_names <- paste0("gamma_", seq_along(self$gamma))
        for (k in seq_along(self$gamma)) {
          rows[[length(rows) + 1]] <- data.frame(
            parameter = x3_names[k],
            type = "supply (gamma)",
            estimate = self$gamma[k],
            se = if (!is.null(private$gamma_se_)) private$gamma_se_[k] else NA_real_,
            stringsAsFactors = FALSE
          )
        }
      }

      # Rho
      if (!is.null(self$rho)) {
        for (h in seq_along(self$rho)) {
          rows[[length(rows) + 1]] <- data.frame(
            parameter = sprintf("rho[%d]", h),
            type = "nesting (rho)",
            estimate = self$rho[h],
            se = NA_real_,
            stringsAsFactors = FALSE
          )
        }
      }

      if (length(rows) == 0) {
        return(data.frame(parameter = character(), type = character(),
                          estimate = numeric(), se = numeric()))
      }
      result <- do.call(rbind, rows)
      result$t_stat <- result$estimate / result$se
      result$p_value <- 2 * stats::pnorm(-abs(result$t_stat))
      result
    },

    #' @description Print estimation results
    print = function(...) {
      cat("BLP Estimation Results\n")
      cat(sprintf("  Method: %s GMM\n", toupper(self$method)))
      cat(sprintf("  Objective: %.6e\n", self$objective))
      cat(sprintf("  Optimization converged: %s\n", self$optimization_converged))
      cat(sprintf("  FP converged: %s (%d total iterations)\n",
                  self$fp_converged, self$fp_iterations))
      cat("\n")

      tbl <- self$summary_table()
      if (nrow(tbl) > 0) {
        cat("Parameter Estimates:\n")
        tbl_print <- tbl[, c("parameter", "estimate", "se", "t_stat")]
        tbl_print$estimate <- sprintf("%.6f", tbl_print$estimate)
        tbl_print$se <- ifelse(is.na(tbl$se), "  ---", sprintf("%.6f", tbl$se))
        tbl_print$t_stat <- ifelse(is.na(tbl$t_stat), " ---",
                                    sprintf("%.3f", tbl$t_stat))
        print(tbl_print, row.names = FALSE, right = FALSE)
      }
      invisible(self)
    }
  ),
  private = list(
    params_ = NULL,
    step_results_ = NULL,
    beta_se_ = NULL,
    gamma_se_ = NULL,

    compute_all_se_ = function() {
      # Extract beta and gamma SEs from step_results_
      # These were computed in problem.R's compute_standard_errors
      sr <- private$step_results_
      if (length(sr) > 0) {
        last <- sr[[length(sr)]]
        # beta_se and gamma_se are stored by the problem's SE computation
        # We need to recompute here
      }
    },

    build_market_ = function(market_id) {
      md <- self$problem$get_market_data(market_id)
      md$products$beta_price <- if (!is.null(self$beta) &&
        !is.null(md$products$price_col_x1)) {
        self$beta[md$products$price_col_x1]
      } else NULL

      mkt <- BLPMarket$new(
        products = md$products,
        agents = md$agents,
        sigma = self$sigma,
        pi = self$pi,
        rho = self$rho,
        rc_types = self$problem$rc_types,
        epsilon_scale = self$problem$epsilon_scale,
        costs_type = self$problem$costs_type
      )
      mu <- mkt$compute_mu()
      prob <- mkt$compute_probabilities(self$delta[md$indices], mu)

      list(market = mkt, products = md$products, agents = md$agents,
           indices = md$indices, prob = prob, mu = mu)
    },

    compute_market_elasticities = function(market_id) {
      mkt_data <- private$build_market_(market_id)
      mkt_data$market$compute_elasticities(
        mkt_data$prob, mkt_data$products$prices
      )
    },

    compute_market_diversion = function(market_id) {
      mkt_data <- private$build_market_(market_id)
      mkt_data$market$compute_diversion_ratios(mkt_data$prob)
    }
  )
)

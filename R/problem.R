#' @title BLP Problem
#' @description The main user-facing class for BLP demand estimation.
#'   Inherits from BLPEconomy and adds the solve() interface.
#' @export
BLPProblem <- R6::R6Class("BLPProblem",
  inherit = BLPEconomy,
  public = list(
    #' @description Create a BLP Problem
    #' @param product_formulations List of BLPFormulation objects (1-3).
    #'   First: linear demand (X1). Second: nonlinear demand (X2). Third: supply (X3).
    #' @param product_data Data frame with columns: market_ids, shares, prices,
    #'   and optionally firm_ids, nesting_ids, demand_instruments*, supply_instruments*
    #' @param agent_formulation Optional BLPFormulation for demographics
    #' @param agent_data Optional data frame with columns: market_ids, weights, nodes*
    #' @param integration Optional BLPIntegration object
    #' @param rc_types Character vector of random coefficient types ("linear", "log", "logit")
    #' @param epsilon_scale Scaling for epsilon (default 1)
    #' @param costs_type "linear" or "log"
    #' @param add_exogenous Whether to add exogenous X1/X3 columns to instruments
    initialize = function(product_formulations, product_data,
                          agent_formulation = NULL, agent_data = NULL,
                          integration = NULL, rc_types = NULL,
                          epsilon_scale = 1.0, costs_type = "linear",
                          add_exogenous = TRUE) {
      super$initialize(product_formulations, product_data,
                       agent_formulation, agent_data,
                       integration, rc_types,
                       epsilon_scale, costs_type, add_exogenous)
      rblp_message(sprintf(
        "BLPProblem: %d markets, %d products, K1=%d, K2=%d, K3=%d, MD=%d, MS=%d",
        self$T, self$N, self$K1, self$K2, self$K3, self$MD, self$MS))
    },

    #' @description Solve the BLP estimation problem
    #' @param sigma Initial K2 x K2 Cholesky root (0 = fixed, non-zero = free)
    #' @param pi Initial K2 x D demographics interaction (0 = fixed)
    #' @param rho Initial nesting parameters (0 = fixed)
    #' @param beta Demand linear coefficients (NA = concentrated out)
    #' @param gamma Supply linear coefficients (NA = concentrated out)
    #' @param sigma_bounds List of lower, upper bound matrices for sigma
    #' @param pi_bounds List of lower, upper bound matrices for pi
    #' @param rho_bounds List of lower, upper bound vectors for rho
    #' @param method '1s' or '2s' for one-step or two-step GMM
    #' @param optimization BLPOptimization object
    #' @param iteration BLPIteration object (for contraction mapping)
    #' @param fp_type Fixed point type: "safe_linear", "linear", "nonlinear"
    #' @param W_type Weighting matrix type: "robust", "clustered", "unadjusted"
    #' @param se_type Standard error type: "robust", "clustered", "unadjusted"
    #' @param initial_W Optional initial weighting matrix
    #' @param scale_objective Whether to scale objective by N
    #' @param center_moments Whether to center moment conditions
    #' @param delta_behavior How to initialize delta: "first", "logit", "last"
    #' @param micro_moments Optional list of MicroMoment objects
    #' @param error_behavior "revert", "punish", or "raise"
    #' @param error_punishment Punishment scale for objective on error
    #' @param processes Number of parallel processes (1 = no parallel)
    #' @return A BLPResults object
    solve = function(sigma = NULL, pi = NULL, rho = NULL,
                     beta = NULL, gamma = NULL,
                     sigma_bounds = NULL, pi_bounds = NULL,
                     rho_bounds = NULL,
                     method = "2s",
                     optimization = NULL,
                     iteration = NULL,
                     fp_type = "safe_linear",
                     W_type = "robust", se_type = "robust",
                     initial_W = NULL,
                     initial_update = FALSE,
                     scale_objective = TRUE,
                     center_moments = TRUE,
                     delta_behavior = "first",
                     micro_moments = NULL,
                     error_behavior = "revert",
                     error_punishment = 1,
                     processes = 1L) {
      # Defaults
      if (is.null(optimization)) optimization <- blp_optimization("l-bfgs-b")
      if (is.null(iteration)) iteration <- blp_iteration("squarem")

      # Build parameter manager
      params <- BLPParameters$new(
        sigma = sigma, pi = pi, rho = rho,
        beta = beta, gamma = gamma,
        sigma_bounds = sigma_bounds, pi_bounds = pi_bounds,
        rho_bounds = rho_bounds,
        rc_types = self$rc_types
      )

      # Initialize mean utilities (delta) at the closed-form logit solution
      # delta_jt = log(s_jt) - log(s_0t), which is the analytical inverse of
      # the logit share equation. This provides a good starting point for the
      # BLP contraction mapping that will iteratively recover the full
      # random-coefficients delta from observed shares.
      rho_init <- if (!is.null(rho)) rho else NULL
      delta <- self$compute_logit_delta(rho_init)

      # Build the initial GMM weighting matrix W for the first step.
      # The default is W = (Z'Z/N)^{-1}, which corresponds to the 2SLS weighting.
      # Under homoskedasticity E[u^2 Z'Z] = sigma^2 * Z'Z, so (Z'Z)^{-1} is
      # proportional to the efficient weight. This choice makes the first GMM
      # step equivalent to standard 2SLS/IV estimation. The demand and supply
      # instrument blocks are stacked via block-diagonal so that each set of
      # moments is weighted independently.
      W <- initial_W
      if (is.null(W)) {
        ZD <- self$products$ZD
        ZS <- self$products$ZS
        Z_blocks <- list()
        if (!is.null(ZD)) Z_blocks[[length(Z_blocks) + 1]] <- ZD
        if (!is.null(ZS)) Z_blocks[[length(Z_blocks) + 1]] <- ZS
        if (length(Z_blocks) > 0) {
          Z_full <- Matrix::bdiag(Z_blocks)
          Z_full <- as.matrix(Z_full)
          S_init <- crossprod(Z_full) / self$N
          W <- compute_gmm_weights(S_init)
        } else {
          W <- diag(1)
        }
      }

      # Handle micro moments
      n_micro <- 0L
      if (!is.null(micro_moments)) {
        n_micro <- length(micro_moments)
        n_total <- nrow(W) + n_micro
        W_new <- matrix(0, n_total, n_total)
        W_new[seq_len(nrow(W)), seq_len(ncol(W))] <- W
        # Micro moment weight = N (as in pyblp)
        for (mm in seq_len(n_micro)) {
          W_new[nrow(W) + mm, ncol(W) + mm] <- self$N
        }
        W <- W_new
      }

      # Two-step GMM: Step 1 uses the 2SLS weight W = (Z'Z/N)^{-1} to get
      # consistent (but inefficient) estimates and residuals. Those residuals
      # are then used to estimate the optimal weighting matrix
      # W* = S^{-1} = [Var(Z'u)]^{-1}, which is used in Step 2 to produce
      # efficient estimates. One-step GMM ("1s") skips the weight update.
      n_steps <- if (method == "2s") 2L else 1L

      # Compress the structural parameter matrices (sigma, pi, rho) into a
      # single vector theta that the optimizer searches over. Only free (nonzero)
      # elements are included: e.g., the lower triangle of sigma (Cholesky root
      # of the random coefficient covariance) and nonzero pi entries. Fixed-at-zero
      # elements are excluded, reducing the dimensionality of the search.
      theta <- params$compress()
      bounds <- params$get_bounds()

      rblp_message(sprintf("Solving BLP with %d nonlinear parameters, %s GMM",
                           params$n_free(), method))

      # Iteration state
      last_delta <- delta
      last_objective <- Inf
      step_results <- list()

      # Initial update: evaluate at starting values to get delta and residuals,
      # then update the weighting matrix before the first GMM step. This is
      # important when using micro moments or when the initial W = (Z'Z/N)^{-1}
      # is a poor approximation (pyblp: initial_update=True).
      if (initial_update) {
        rblp_message("\n--- Initial Update ---")
        init_progress <- private$compute_progress(
          theta, params, last_delta, W,
          iteration, fp_type, scale_objective, center_moments,
          micro_moments, error_behavior, error_punishment, processes
        )
        if (!is.null(init_progress$delta)) last_delta <- init_progress$delta

        # Update W from initial residuals
        S_init <- compute_gmm_moment_covariances(
          init_progress$u_list, init_progress$Z_list,
          type = W_type,
          clustering_ids = self$products$clustering_ids
        )
        if (!is.null(micro_moments) && n_micro > 0) {
          total_agg <- nrow(S_init)
          n_total <- total_agg + n_micro
          S_new <- matrix(0, n_total, n_total)
          S_new[seq_len(total_agg), seq_len(total_agg)] <- S_init
          for (mm in seq_len(n_micro)) {
            S_new[total_agg + mm, total_agg + mm] <- 1 / self$N
          }
          S_init <- S_new
        }
        W <- compute_gmm_weights(S_init)
        rblp_message(sprintf("Initial update: objective = %.8e", init_progress$objective))
      }

      for (step in seq_len(n_steps)) {
        rblp_message(sprintf("\n--- GMM Step %d/%d ---", step, n_steps))

        # Build the GMM objective function that the nonlinear optimizer calls.
        # For each candidate theta, this: (1) runs the BLP contraction mapping
        # to recover delta(theta), (2) concentrates out beta via IV regression
        # to get xi(theta) = delta - X1*beta, (3) forms moments g = Z'xi/N,
        # and (4) returns the quadratic form g'Wg and its gradient.
        # The warm-start delta from the previous evaluation (last_delta) is
        # carried forward to reduce contraction iterations.
        objective_fn <- function(theta_vec) {
          result <- private$compute_progress(
            theta_vec, params, last_delta, W,
            iteration, fp_type, scale_objective, center_moments,
            micro_moments, error_behavior, error_punishment, processes
          )

          if (!is.null(result$delta)) last_delta <<- result$delta

          list(
            objective = result$objective,
            gradient = result$gradient
          )
        }

        # Optimize
        if (params$n_free() > 0) {
          opt_result <- optimization$optimize(theta, bounds, objective_fn)
          theta <- opt_result$values
          opt_converged <- opt_result$converged
          opt_iterations <- opt_result$iterations
          opt_evaluations <- opt_result$evaluations
        } else {
          opt_converged <- TRUE
          opt_iterations <- 0L
          opt_evaluations <- 0L
        }

        # Final evaluation at optimum
        progress <- private$compute_progress(
          theta, params, last_delta, W,
          iteration, fp_type, scale_objective, center_moments,
          micro_moments, error_behavior, error_punishment, processes
        )
        last_delta <- progress$delta

        step_results[[step]] <- list(
          theta = theta,
          progress = progress,
          converged = opt_converged,
          iterations = opt_iterations,
          evaluations = opt_evaluations
        )

        rblp_message(sprintf("Step %d: objective = %.8e, converged = %s",
                             step, progress$objective, opt_converged))

        # Between GMM steps, update the weighting matrix to the efficient
        # (optimal) weight W* = S^{-1}. The moment covariance S estimates
        # Var(g) = E[Z_i' u_i u_i' Z_i] using the Step-1 residuals.
        # Under "robust" type, S is heteroskedasticity-consistent (like HC0).
        # Under "clustered", moments are summed within clusters before forming
        # the outer product, yielding cluster-robust inference.
        # This update makes Step 2 the efficient two-step GMM estimator.
        if (step < n_steps) {
          S <- compute_gmm_moment_covariances(
            progress$u_list, progress$Z_list,
            type = W_type,
            clustering_ids = self$products$clustering_ids
          )
          if (!is.null(micro_moments) && n_micro > 0) {
            total_agg <- nrow(S)
            n_total <- total_agg + n_micro
            S_new <- matrix(0, n_total, n_total)
            S_new[seq_len(total_agg), seq_len(total_agg)] <- S
            # Micro covariance block
            for (mm in seq_len(n_micro)) {
              S_new[total_agg + mm, total_agg + mm] <- 1 / self$N
            }
            S <- S_new
          }
          W <- compute_gmm_weights(S)
        }
      }

      # Final step results
      final <- step_results[[n_steps]]
      progress <- final$progress

      # Compute standard errors
      se_result <- private$compute_standard_errors(
        progress, W, params, se_type, scale_objective
      )

      # Build and return results
      expanded <- params$expand(theta)

      BLPResults$new(
        problem = self,
        params = params,
        sigma = expanded$sigma,
        pi = expanded$pi,
        rho = expanded$rho,
        beta = progress$beta,
        gamma = progress$gamma,
        delta = last_delta,
        xi = progress$xi,
        omega = progress$omega,
        objective = progress$objective,
        gradient = progress$gradient,
        hessian = se_result$hessian,
        se = se_result$se,
        parameter_covariances = se_result$param_cov,
        W = W,
        step_results = step_results,
        optimization_converged = final$converged,
        optimization_iterations = final$iterations,
        optimization_evaluations = final$evaluations,
        fp_converged = progress$fp_converged,
        fp_iterations = progress$fp_iterations,
        method = method,
        se_type = se_type,
        beta_se = se_result$beta_se,
        gamma_se = se_result$gamma_se
      )
    },

    #' @description Print the problem
    print = function(...) {
      cat(sprintf("BLPProblem: %d markets, %d products\n", self$T, self$N))
      cat(sprintf("  Linear demand (K1):    %d\n", self$K1))
      cat(sprintf("  Nonlinear demand (K2): %d\n", self$K2))
      cat(sprintf("  Supply (K3):           %d\n", self$K3))
      cat(sprintf("  Demand instruments:    %d\n", self$MD))
      cat(sprintf("  Supply instruments:    %d\n", self$MS))
      cat(sprintf("  Agents (I):            %d\n", self$I))
      invisible(self)
    }
  ),
  private = list(
    compute_progress = function(theta, params, delta, W,
                                 iteration, fp_type, scale_objective,
                                 center_moments, micro_moments,
                                 error_behavior, error_punishment,
                                 processes) {
      # Expand the compressed theta vector back into the full parameter matrices.
      # sigma is the K2 x K2 lower-triangular Cholesky factor of the random
      # coefficient covariance (Sigma = sigma * sigma'). pi is the K2 x D
      # matrix of interactions between product characteristics and consumer
      # demographics. rho contains nesting parameters for nested logit.
      expanded <- params$expand(theta)
      sigma <- expanded$sigma
      pi_mat <- expanded$pi
      rho <- expanded$rho

      # Free masks track which elements of sigma/pi are being estimated (nonzero
      # initial values). This ensures the Jacobian columns line up with the
      # compressed theta vector: only free parameters get Jacobian columns.
      sigma_free <- params$get_sigma_free()
      pi_free <- params$get_pi_free()

      N <- self$N
      fp_iterations_total <- 0L
      fp_converged_all <- TRUE
      xi_jacobian_list <- list()
      delta_new <- delta

      # The BLP inner loop is separable across markets: each market t has its
      # own contraction mapping to invert observed shares s_t into mean
      # utilities delta_t, conditional on the current nonlinear parameters
      # theta. This per-market structure is the key computational unit of BLP.
      market_fn <- function(market_id) {
        md <- self$get_market_data(market_id)

        # Create a market-level object that holds this market's products,
        # agents (integration nodes/weights), and the current sigma/pi/rho.
        mkt <- BLPMarket$new(
          products = md$products,
          agents = md$agents,
          sigma = sigma,
          pi = pi_mat,
          rho = rho,
          rc_types = self$rc_types,
          epsilon_scale = self$epsilon_scale,
          costs_type = self$costs_type,
          sigma_free = sigma_free,
          pi_free = pi_free
        )

        # Run the BLP contraction mapping: delta^{h+1} = delta^h + log(s_obs) - log(s_pred).
        # This is the "inner loop" of BLP. For any given theta, the contraction
        # finds the unique delta(theta) such that the model-predicted shares
        # exactly match observed shares. Berry (1994) showed this map is a
        # contraction under standard regularity conditions.
        fp_result <- mkt$compute_delta(
          initial_delta = delta[md$indices],
          iteration = iteration,
          fp_type = fp_type
        )

        delta_t <- fp_result$delta

        # After convergence, recompute choice probabilities P_ijt at the
        # converged delta. These are needed for the Jacobian computation.
        mu <- mkt$compute_mu()
        prob_result <- mkt$compute_probabilities(delta_t, mu)

        # Compute the Jacobian d_xi/d_theta via the implicit function theorem.
        # Since delta(theta) is defined implicitly by s(delta, theta) = s_obs,
        # differentiating gives: d_xi/d_theta = -(ds/d_delta)^{-1} * (ds/d_theta).
        # This Jacobian is needed for the gradient of the GMM objective.
        xi_jac <- NULL
        if (params$n_free() > 0) {
          xi_jac <- mkt$compute_xi_by_theta_jacobian(prob_result)
        }

        # Supply side: compute costs and omega Jacobian
        costs <- NULL
        omega_jac <- NULL
        if (self$K3 > 0 && !is.null(md$products$ownership)) {
          P <- if (is.list(prob_result)) prob_result$probabilities else prob_result
          costs <- mkt$compute_costs(P, md$products$prices, md$products$ownership)
          if (self$costs_type == "log") costs <- log(costs)
        }

        list(
          delta = delta_t,
          indices = md$indices,
          xi_jacobian = xi_jac,
          costs = costs,
          fp_converged = fp_result$converged,
          fp_iterations = fp_result$iterations
        )
      }

      # Run markets
      if (processes > 1L && self$T > 1L) {
        market_results <- parallel::mclapply(
          self$unique_market_ids, market_fn, mc.cores = processes
        )
      } else {
        market_results <- lapply(self$unique_market_ids, market_fn)
      }

      # Aggregate results
      full_xi_jac <- if (params$n_free() > 0) {
        matrix(0, N, params$n_free())
      } else NULL
      costs_vec <- if (self$K3 > 0) numeric(N) else NULL

      for (res in market_results) {
        delta_new[res$indices] <- res$delta
        if (!is.null(res$xi_jacobian) && !is.null(full_xi_jac)) {
          full_xi_jac[res$indices, ] <- res$xi_jacobian
        }
        if (!is.null(res$costs)) {
          costs_vec[res$indices] <- res$costs
        }
        fp_iterations_total <- fp_iterations_total + res$fp_iterations
        if (!res$fp_converged) fp_converged_all <- FALSE
      }

      # Concentrate out beta via IV/GMM regression: beta = (X1'Z W Z'X1)^{-1} X1'Z W Z'delta.
      # This is the "outer loop" concentration step. Because beta enters linearly
      # in delta = X1*beta + xi, we can solve for beta analytically given delta(theta),
      # reducing the nonlinear search to theta only. The structural error is
      # xi = delta - X1*beta, and the moment condition is E[Z'xi] = 0.
      X1 <- self$products$X1
      ZD <- self$products$ZD

      # When fixed effects are absorbed via the Frisch-Waugh-Lovell (FWL) theorem,
      # X1 and ZD were already group-demeaned during economy construction. For
      # consistency, delta and the xi-Jacobian must also be demeaned by the same
      # groups before the IV regression. This is equivalent to including the FE
      # dummies in the regression but is computationally cheaper.
      delta_iv <- delta_new
      if (!is.null(private$absorb_groups_)) {
        grp <- private$absorb_groups_
        # Use pre-computed group indices for fast demeaning
        grp_idx <- split(seq_len(N), grp)
        grp_sizes <- lengths(grp_idx)

        # Demean delta: subtract group mean
        for (g in names(grp_idx)) {
          idx_g <- grp_idx[[g]]
          gm <- mean(delta_iv[idx_g])
          delta_iv[idx_g] <- delta_iv[idx_g] - gm
        }

        # The xi-Jacobian (d_xi/d_theta) must also be FWL-demeaned. Without this,
        # the gradient would be inconsistent with the demeaned moment conditions,
        # because the Jacobian feeds into G = Z' (d_xi/d_theta) / N which must
        # be computed in the same demeaned space as the moments g = Z' xi / N.
        if (!is.null(full_xi_jac)) {
          for (g in names(grp_idx)) {
            idx_g <- grp_idx[[g]]
            gm_j <- colMeans(full_xi_jac[idx_g, , drop = FALSE])
            full_xi_jac[idx_g, ] <- sweep(full_xi_jac[idx_g, , drop = FALSE], 2, gm_j)
          }
        }
      }

      # Extract the demand block of the GMM weighting matrix. The full W is
      # block-diagonal with demand (MD x MD) and supply (MS x MS) blocks,
      # reflecting that demand and supply moment conditions are independent.
      MD <- self$MD
      MS <- self$MS
      W_demand <- W[seq_len(MD), seq_len(MD), drop = FALSE]

      # IV estimation concentrates out beta: given delta(theta), solve the
      # GMM normal equations for beta, yielding xi = delta - X1*beta.
      # The residual_jacobian is d_xi/d_theta after accounting for the
      # concentration: d_xi/d_theta = d_delta/d_theta - X1 * d_beta/d_theta.
      # This "concentrated Jacobian" is what enters the gradient formula.
      demand_iv <- iv_estimate(X1, ZD, W_demand, delta_iv, full_xi_jac)
      beta <- demand_iv$parameters
      xi <- demand_iv$residuals
      xi_jac_concentrated <- demand_iv$residual_jacobian

      # Supply side IV
      gamma <- NULL
      omega <- NULL
      omega_jac_concentrated <- NULL
      if (self$K3 > 0 && !is.null(costs_vec)) {
        X3 <- self$products$X3
        ZS <- self$products$ZS
        W_supply <- W[(MD + 1):(MD + MS), (MD + 1):(MD + MS), drop = FALSE]
        supply_iv <- iv_estimate(X3, ZS, W_supply, costs_vec, NULL)
        gamma <- supply_iv$parameters
        omega <- supply_iv$residuals
      }

      # Form the sample moment vector g = (1/N) * Z' u, which is the empirical
      # analog of the population moment condition E[Z'u] = 0. For demand,
      # u = xi (unobserved product quality); for supply, u = omega (cost shock).
      # The identifying assumption is that instruments Z are uncorrelated with
      # the structural errors: E[Z_d' xi] = 0 and E[Z_s' omega] = 0.
      # Stacking demand and supply moments gives the full moment vector.
      u_list <- list(xi)
      Z_list <- list(ZD)
      if (!is.null(omega)) {
        u_list[[2]] <- omega
        Z_list[[2]] <- ZS
      }

      g_parts <- lapply(seq_along(u_list), function(a) {
        crossprod(Z_list[[a]], u_list[[a]]) / N
      })
      g <- do.call(rbind, g_parts)

      # Add micro moments
      if (!is.null(micro_moments) && length(micro_moments) > 0) {
        micro_g <- sapply(micro_moments, function(mm) {
          mm$value - mm$compute_simulated_value(self, delta_new, sigma, pi_mat, rho)
        })
        g <- rbind(g, matrix(micro_g, ncol = 1))
      }

      # Center moments
      if (center_moments && nrow(g) > 0) {
        # moments are already averaged, centering not needed at aggregate level
      }

      # The GMM objective is the quadratic form Q(theta) = g(theta)' W g(theta).
      # The optimizer minimizes this over theta. At the true parameter values,
      # g -> 0 in probability, so Q -> 0. Scaling by N gives Q*N which has
      # better numerical conditioning for optimization.
      objective <- as.numeric(crossprod(g, W %*% g))
      if (scale_objective) objective <- objective * N

      # Analytic gradient of the GMM objective: dQ/d_theta = 2 * G' W g,
      # where G = (1/N) * Z' (d_xi/d_theta) is the Jacobian of the moment
      # vector with respect to theta. The d_xi/d_theta used here is the
      # "concentrated" Jacobian that accounts for beta being re-optimized
      # at each theta. This analytic gradient dramatically speeds convergence
      # compared to numerical differencing, especially with many parameters.
      gradient <- NULL
      if (params$n_free() > 0 && !is.null(xi_jac_concentrated)) {
        G_parts <- list(crossprod(ZD, xi_jac_concentrated) / N)
        if (!is.null(omega_jac_concentrated)) {
          G_parts[[2]] <- crossprod(ZS, omega_jac_concentrated) / N
        } else if (MS > 0) {
          # Supply moments present but no supply Jacobian: add zero block
          # so G dimensions match the full moment vector g = [g_demand; g_supply]
          G_parts[[2]] <- matrix(0, MS, params$n_free())
        }
        G <- do.call(rbind, G_parts)

        if (!is.null(micro_moments) && length(micro_moments) > 0) {
          # Micro moment Jacobian: d(g_micro)/d(theta) via finite differences.
          # Each micro moment g_m = f_m(theta) - target_m depends on theta
          # through the choice probabilities. The analytic derivative is complex
          # (requiring derivatives of conditional expectations), so we use
          # central finite differences as in pyblp's default implementation.
          eps_fd <- getOption("rblp.finite_differences_epsilon",
                              sqrt(.Machine$double.eps))
          n_micro <- length(micro_moments)
          n_theta <- params$n_free()
          G_micro <- matrix(0, n_micro, n_theta)

          for (k in seq_len(n_theta)) {
            theta_plus <- theta
            theta_minus <- theta
            theta_plus[k] <- theta_plus[k] + eps_fd
            theta_minus[k] <- theta_minus[k] - eps_fd

            exp_plus <- params$expand(theta_plus)
            exp_minus <- params$expand(theta_minus)

            for (m in seq_len(n_micro)) {
              val_plus <- micro_moments[[m]]$compute_simulated_value(
                self, delta_new, exp_plus$sigma, exp_plus$pi, exp_plus$rho
              )
              val_minus <- micro_moments[[m]]$compute_simulated_value(
                self, delta_new, exp_minus$sigma, exp_minus$pi, exp_minus$rho
              )
              G_micro[m, k] <- (val_plus - val_minus) / (2 * eps_fd)
            }
          }
          G <- rbind(G, G_micro)
        }

        gradient <- as.numeric(2 * crossprod(G, W %*% g))
        if (scale_objective) gradient <- gradient * N
      }

      list(
        objective = objective,
        gradient = gradient,
        delta = delta_new,
        beta = beta,
        gamma = gamma,
        xi = xi,
        omega = omega,
        g = g,
        G = if (exists("G")) G else NULL,
        u_list = u_list,
        Z_list = Z_list,
        fp_converged = fp_converged_all,
        fp_iterations = fp_iterations_total
      )
    },

    compute_standard_errors = function(progress, W, params, se_type,
                                        scale_objective) {
      N <- self$N
      n_theta <- params$n_free()

      # Estimate the moment covariance matrix S = Var(sqrt(N)*g).
      # Under "robust" (heteroskedasticity-consistent), S = (1/N) sum_i (Z_i'u_i)(Z_i'u_i)'.
      # Under "clustered", observations within each cluster are summed before
      # forming the outer product, allowing for within-cluster correlation.
      S <- compute_gmm_moment_covariances(
        progress$u_list, progress$Z_list,
        type = se_type,
        clustering_ids = self$products$clustering_ids
      )

      # If micro moments are present, extend S to include micro blocks.
      # The micro moment covariance block is diagonal with entries 1/n_obs.
      n_S <- nrow(S)
      n_W <- nrow(W)
      if (n_W > n_S) {
        n_micro <- n_W - n_S
        S_ext <- matrix(0, n_W, n_W)
        S_ext[seq_len(n_S), seq_len(n_S)] <- S
        for (mm in seq_len(n_micro)) {
          S_ext[n_S + mm, n_S + mm] <- 1 / N
        }
        S <- S_ext
      }

      # GMM sandwich covariance for nonlinear parameters theta:
      # V(theta) = (1/N) * (G'WG)^{-1} G'W S W'G (G'WG)^{-1}.
      # Under efficient two-step GMM where W = S^{-1}, this simplifies to
      # (G'S^{-1}G)^{-1} / N. The "sandwich" form is robust even if the
      # weighting matrix is not exactly optimal.
      param_cov <- NULL
      se_vec <- NULL
      hessian <- NULL

      if (n_theta > 0 && !is.null(progress$G)) {
        G <- progress$G
        param_cov <- compute_gmm_parameter_covariances(W, S, G, se_type)
        param_cov <- param_cov / N

        se_vec <- sqrt(pmax(diag(param_cov), 0))
      }

      # Standard errors for the concentrated-out linear parameters beta and gamma.
      # Since beta was solved analytically (not searched over by the optimizer),
      # its covariance requires a separate sandwich formula. This is the standard
      # GMM/IV sandwich: V(beta) = (X'Z W Z'X)^{-1} (X'Z W S W Z'X) (X'Z W Z'X)^{-1}.
      # The "bread" is (X'Z W Z'X)^{-1} and the "meat" is X'Z W S W Z'X.
      # Under efficient GMM (W = S^{-1}), this simplifies to (X'Z S^{-1} Z'X)^{-1}.
      # Note: this treats theta as fixed at its estimate; the joint covariance
      # of (theta, beta) would require the full influence function.
      beta_se <- NULL
      gamma_se <- NULL
      if (!is.null(progress$beta)) {
        X1 <- self$products$X1
        ZD <- self$products$ZD
        MD <- self$MD
        W_d <- W[seq_len(MD), seq_len(MD), drop = FALSE]
        bread <- crossprod(X1, ZD) %*% W_d %*% crossprod(ZD, X1)
        bread_inv <- approximately_invert(bread)$inverse
        XZW <- crossprod(X1, ZD) %*% W_d
        S_d <- S[seq_len(MD), seq_len(MD), drop = FALSE]
        beta_cov <- N * bread_inv %*% XZW %*% S_d %*% t(XZW) %*% bread_inv
        beta_se <- sqrt(pmax(diag(beta_cov), 0))
      }

      if (!is.null(progress$gamma) && self$K3 > 0) {
        X3 <- self$products$X3
        ZS <- self$products$ZS
        MS <- self$MS
        MD <- self$MD
        W_s <- W[(MD + 1):(MD + MS), (MD + 1):(MD + MS), drop = FALSE]
        bread_s <- crossprod(X3, ZS) %*% W_s %*% crossprod(ZS, X3)
        bread_s_inv <- approximately_invert(bread_s)$inverse
        XZW_s <- crossprod(X3, ZS) %*% W_s
        S_s <- S[(MD + 1):(MD + MS), (MD + 1):(MD + MS), drop = FALSE]
        gamma_cov <- N * bread_s_inv %*% XZW_s %*% S_s %*% t(XZW_s) %*% bread_s_inv
        gamma_se <- sqrt(pmax(diag(gamma_cov), 0))
      }

      list(
        param_cov = param_cov,
        se = se_vec,
        hessian = hessian,
        beta_se = beta_se,
        gamma_se = gamma_se
      )
    }
  )
)

#' Create a BLP Problem
#'
#' Main entry point for BLP demand estimation.
#'
#' @param product_formulations List of BLPFormulation objects (1-3)
#' @param product_data Data frame with market_ids, shares, prices, etc.
#' @param agent_formulation Optional demographics formulation
#' @param agent_data Optional agent data frame
#' @param integration Optional BLPIntegration object
#' @param rc_types Character vector of random coefficient types
#' @param epsilon_scale Epsilon scaling (default 1)
#' @param costs_type "linear" or "log"
#' @param add_exogenous Whether to add exogenous regressors to instruments
#' @return A BLPProblem object
#' @export
#' @examples
#' \dontrun{
#' # Logit model
#' f1 <- blp_formulation(~ prices + sugar + mushy)
#' problem <- blp_problem(list(f1), product_data)
#' results <- problem$solve()
#'
#' # Random coefficients
#' f2 <- blp_formulation(~ prices + sugar + mushy)
#' problem <- blp_problem(list(f1, f2), product_data,
#'                         integration = blp_integration("product", 5))
#' results <- problem$solve(sigma = diag(3))
#' }
blp_problem <- function(product_formulations, product_data,
                          agent_formulation = NULL, agent_data = NULL,
                          integration = NULL, rc_types = NULL,
                          epsilon_scale = 1.0, costs_type = "linear",
                          add_exogenous = TRUE) {
  BLPProblem$new(product_formulations, product_data,
                  agent_formulation, agent_data,
                  integration, rc_types,
                  epsilon_scale, costs_type, add_exogenous)
}

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

      # Initialize delta
      rho_init <- if (!is.null(rho)) rho else NULL
      delta <- self$compute_logit_delta(rho_init)

      # Build initial weighting matrix
      W <- initial_W
      if (is.null(W)) {
        # Default: 2SLS weighting (Z'Z/N)^{-1}
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

      # Determine number of GMM steps
      n_steps <- if (method == "2s") 2L else 1L

      # Initial compressed theta
      theta <- params$compress()
      bounds <- params$get_bounds()

      rblp_message(sprintf("Solving BLP with %d nonlinear parameters, %s GMM",
                           params$n_free(), method))

      # Iteration state
      last_delta <- delta
      last_objective <- Inf
      step_results <- list()

      for (step in seq_len(n_steps)) {
        rblp_message(sprintf("\n--- GMM Step %d/%d ---", step, n_steps))

        # Build objective function
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

        # Update W for second step
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
      expanded <- params$expand(theta)
      sigma <- expanded$sigma
      pi_mat <- expanded$pi
      rho <- expanded$rho

      # Get free masks for consistent Jacobian dimensions
      sigma_free <- params$get_sigma_free()
      pi_free <- params$get_pi_free()

      N <- self$N
      fp_iterations_total <- 0L
      fp_converged_all <- TRUE
      xi_jacobian_list <- list()
      delta_new <- delta

      # Market-by-market computation
      market_fn <- function(market_id) {
        md <- self$get_market_data(market_id)

        # Create market object
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

        # Compute delta via contraction
        fp_result <- mkt$compute_delta(
          initial_delta = delta[md$indices],
          iteration = iteration,
          fp_type = fp_type
        )

        delta_t <- fp_result$delta

        # Compute probabilities at converged delta
        mu <- mkt$compute_mu()
        prob_result <- mkt$compute_probabilities(delta_t, mu)

        # Compute Jacobian: d_xi/d_theta
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

      # IV regression for demand side
      X1 <- self$products$X1
      ZD <- self$products$ZD

      # Apply absorb demeaning to delta (X1 and ZD already demeaned)
      delta_iv <- delta_new
      if (!is.null(private$absorb_groups_)) {
        grp <- private$absorb_groups_
        gm <- tapply(delta_iv, grp, mean)
        delta_iv <- as.numeric(delta_iv - gm[match(grp, names(gm))])
      }

      # Extract demand block of W
      MD <- self$MD
      MS <- self$MS
      W_demand <- W[seq_len(MD), seq_len(MD), drop = FALSE]

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

      # Compute moments g = Z' u / N
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

      # Objective: g' W g
      objective <- as.numeric(crossprod(g, W %*% g))
      if (scale_objective) objective <- objective * N

      # Gradient
      gradient <- NULL
      if (params$n_free() > 0 && !is.null(xi_jac_concentrated)) {
        G_parts <- list(crossprod(ZD, xi_jac_concentrated) / N)
        if (!is.null(omega_jac_concentrated)) {
          G_parts[[2]] <- crossprod(ZS, omega_jac_concentrated) / N
        }
        G <- do.call(rbind, G_parts)

        if (!is.null(micro_moments) && length(micro_moments) > 0) {
          # Micro moment gradient placeholder (zero for now)
          G_micro <- matrix(0, length(micro_moments), params$n_free())
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

      # Moment covariances
      S <- compute_gmm_moment_covariances(
        progress$u_list, progress$Z_list,
        type = se_type,
        clustering_ids = self$products$clustering_ids
      )

      # Parameter covariances
      param_cov <- NULL
      se_vec <- NULL
      hessian <- NULL

      if (n_theta > 0 && !is.null(progress$G)) {
        G <- progress$G
        param_cov <- compute_gmm_parameter_covariances(W, S, G, se_type)
        param_cov <- param_cov / N

        se_vec <- sqrt(pmax(diag(param_cov), 0))
      }

      # Extract SEs for beta and gamma
      beta_se <- NULL
      gamma_se <- NULL
      if (!is.null(progress$beta)) {
        # Beta concentrated out: need sandwich SE
        X1 <- self$products$X1
        ZD <- self$products$ZD
        MD <- self$MD
        W_d <- W[seq_len(MD), seq_len(MD), drop = FALSE]
        bread <- crossprod(X1, ZD) %*% W_d %*% crossprod(ZD, X1)
        bread_inv <- approximately_invert(bread)$inverse
        # Sandwich for beta: V = (X'ZWZ'X)^{-1} X'ZW S WZ'X (X'ZWZ'X)^{-1}
        # S_d = sum(g_i g_i')/N, so we multiply by N to get the raw sum
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

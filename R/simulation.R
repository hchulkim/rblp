#' @title BLP Simulation
#' @description Simulate BLP-type equilibrium data for testing and Monte Carlo studies.
#' @export
BLPSimulation <- R6::R6Class("BLPSimulation",
  inherit = BLPEconomy,
  public = list(
    #' @field beta True demand linear parameters
    beta = NULL,
    #' @field sigma True Cholesky root of RC covariance
    sigma = NULL,
    #' @field pi True demographics interaction
    pi = NULL,
    #' @field gamma True supply parameters
    gamma = NULL,
    #' @field rho True nesting parameters
    rho = NULL,
    #' @field xi Demand structural errors
    xi = NULL,
    #' @field omega Supply structural errors
    omega = NULL,

    #' @description Create a BLP simulation
    #' @param product_formulations List of BLPFormulation objects (1-3)
    #' @param product_data Data frame with market_ids, firm_ids, and characteristics
    #' @param beta Demand linear coefficients (required)
    #' @param sigma Optional K2 x K2 Cholesky root
    #' @param pi Optional K2 x D demographics interaction
    #' @param gamma Optional supply coefficients
    #' @param rho Optional nesting parameters
    #' @param agent_formulation Optional demographics formulation
    #' @param agent_data Optional agent data
    #' @param integration Optional BLPIntegration
    #' @param xi Optional demand errors (drawn if NULL)
    #' @param omega Optional supply errors (drawn if NULL)
    #' @param xi_variance Variance of xi (default 1)
    #' @param omega_variance Variance of omega (default 1)
    #' @param correlation Correlation between xi and omega (default 0.9)
    #' @param rc_types Character vector of RC types
    #' @param costs_type "linear" or "log"
    #' @param seed Random seed
    #
    # Sets up the data-generating process (DGP) for Monte Carlo studies.
    # The user supplies the TRUE parameter values (beta, sigma, pi, gamma)
    # that govern the structural model. The simulation then generates
    # equilibrium data consistent with these parameters, which can be
    # passed to the estimator to verify that the solver recovers the truth.
    # This is the standard "estimation roundtrip" test for BLP code.
    initialize = function(product_formulations, product_data,
                          beta, sigma = NULL, pi = NULL,
                          gamma = NULL, rho = NULL,
                          agent_formulation = NULL, agent_data = NULL,
                          integration = NULL,
                          xi = NULL, omega = NULL,
                          xi_variance = 1, omega_variance = 1,
                          correlation = 0.9,
                          rc_types = NULL, costs_type = "linear",
                          seed = NULL) {
      # Placeholder shares and prices are needed so the parent class
      # BLPEconomy can build the design matrices X1, X2, X3. The actual
      # equilibrium values will be computed later in replace_endogenous().
      if (is.null(product_data$shares)) product_data$shares <- rep(0.01, nrow(product_data))
      if (is.null(product_data$prices)) product_data$prices <- rep(1, nrow(product_data))

      super$initialize(product_formulations, product_data,
                       agent_formulation, agent_data,
                       integration, rc_types, 1.0, costs_type,
                       add_exogenous = FALSE)

      self$beta <- beta
      self$sigma <- sigma
      self$pi <- pi
      self$gamma <- gamma
      self$rho <- rho

      if (!is.null(seed)) set.seed(seed)
      N <- self$N

      # Draw the demand-side structural error xi ~ N(0, xi_variance).
      # xi captures unobserved product quality (e.g., brand reputation,
      # advertising) that enters the utility function but is not in X1.
      # This is the endogeneity source: firms observe xi when setting prices,
      # creating the E[p * xi] != 0 correlation that necessitates IV.
      if (is.null(xi)) {
        xi <- stats::rnorm(N, 0, sqrt(xi_variance))
      }
      self$xi <- xi

      if (!is.null(gamma) && self$K3 > 0) {
        if (is.null(omega)) {
          # Draw supply-side error omega correlated with xi. The correlation
          # (default 0.9) reflects the empirical regularity that products
          # with high unobserved quality (high xi) also tend to have high
          # unobserved costs (high omega) -- e.g., premium ingredients.
          # Constructed via: omega = rho*xi*(sd_omega/sd_xi) + sqrt(1-rho^2)*eps
          # which yields Corr(xi, omega) = correlation by construction.
          omega <- correlation * xi * sqrt(omega_variance / xi_variance) +
            sqrt(1 - correlation^2) * stats::rnorm(N, 0, sqrt(omega_variance))
        }
        self$omega <- omega
      }

      # Compute the initial mean utility: delta_j = X1_j' beta + xi_j.
      # This is the utility component common to all consumers, before adding
      # the individual-specific random taste shocks mu_ij. At this stage
      # prices are placeholders; they will be replaced by the equilibrium
      # prices computed in replace_endogenous().
      delta <- as.numeric(self$products$X1 %*% beta) + xi
      private$delta_ <- delta

      rblp_message(sprintf("BLPSimulation: %d markets, %d products", self$T, self$N))
    },

    #' @description Solve for equilibrium prices and shares
    #' @param iteration BLPIteration for fixed-point iteration
    #' @param constant_costs Whether costs are independent of shares
    #' @return A BLPSimulationResults object
    #
    # Computes the Nash-Bertrand equilibrium: given marginal costs (from
    # gamma, X3, omega), find prices such that each multi-product firm's
    # first-order condition is satisfied simultaneously. The equilibrium
    # is found by contraction mapping: p_{n+1} = c + markup(p_n). Once
    # prices converge, equilibrium shares follow from the demand model.
    # The resulting (prices, shares) pair is the simulated "observed data"
    # that would be taken to the BLP estimator.
    replace_endogenous = function(iteration = NULL, constant_costs = TRUE) {
      if (is.null(iteration)) iteration <- blp_iteration("simple", list(atol = 1e-12))

      N <- self$N
      prices <- numeric(N)
      shares <- numeric(N)
      costs <- numeric(N)
      delta_eq <- private$delta_

      for (t in self$unique_market_ids) {
        md <- self$get_market_data(t)
        idx <- md$indices

        # Compute marginal costs from the supply-side structural equation:
        # mc_j = X3_j' gamma + omega_j (linear costs), or
        # mc_j = exp(X3_j' gamma + omega_j) (log costs, ensuring mc > 0).
        # These are the TRUE costs used by firms in their pricing decisions.
        if (!is.null(self$gamma) && self$K3 > 0) {
          costs_t <- as.numeric(md$products$X3 %*% self$gamma)
          if (!is.null(self$omega)) costs_t <- costs_t + self$omega[idx]
          if (self$costs_type == "log") costs_t <- exp(costs_t)
        } else {
          costs_t <- rep(1, length(idx))
        }
        costs[idx] <- costs_t

        # Prepare the market for the equilibrium price search. delta_base
        # holds the non-price part of mean utility; beta_price is the
        # linear price coefficient (alpha). The price finder will update
        # delta as delta = delta_base + alpha * p at each iteration.
        md$products$delta_base <- delta_eq[idx]
        md$products$beta_price <- if (!is.null(md$products$price_col_x1)) {
          self$beta[md$products$price_col_x1]
        } else 0

        mkt <- BLPMarket$new(
          products = md$products,
          agents = md$agents,
          sigma = self$sigma,
          pi = self$pi,
          rho = self$rho,
          rc_types = self$rc_types,
          epsilon_scale = self$epsilon_scale,
          costs_type = self$costs_type
        )

        # Solve for Nash-Bertrand equilibrium prices via contraction mapping:
        # p_{n+1} = mc - Omega(p_n)^{-1} s(p_n), starting from 1.5 * costs
        # as an initial guess. Convergence gives the fixed point where each
        # firm's price equals cost plus the optimal markup given rivals' prices.
        p_eq <- mkt$compute_equilibrium_prices(
          costs_t, iteration, md$products$ownership, costs_t * 1.5
        )
        prices[idx] <- p_eq

        # Update mean utility to reflect the equilibrium price: the price
        # component of delta shifts by alpha * (p_eq - p_placeholder).
        delta_t <- delta_eq[idx]
        if (!is.null(md$products$price_col_x1)) {
          pc <- md$products$price_col_x1
          delta_t <- delta_t + self$beta[pc] * (p_eq - md$products$prices)
        }

        # Given equilibrium prices, compute choice probabilities via the
        # mixed logit: P_ij = exp(delta_j + mu_ij) / (1 + sum_k exp(...)).
        # Aggregate across the simulated consumer population (weighting by
        # agent weights) to get equilibrium market shares s_j = sum_i w_i P_ij.
        mu <- mkt$compute_mu()
        prob <- mkt$compute_probabilities(delta_t, mu)
        P <- if (is.list(prob)) prob$probabilities else prob
        shares[idx] <- mkt$compute_shares(P)
        delta_eq[idx] <- delta_t
      }

      # Replace the placeholder prices and shares with their equilibrium
      # values. The resulting product_data is a complete simulated dataset
      # that looks like real data: it contains equilibrium prices (endogenous,
      # correlated with xi) and market shares, ready for BLP estimation.
      product_data <- self$products$original_data
      product_data$prices <- prices
      product_data$shares <- shares

      BLPSimulationResults$new(
        simulation = self,
        product_data = product_data,
        delta = delta_eq,
        costs = costs,
        prices = prices,
        shares = shares
      )
    },

    #' @description Print the simulation
    print = function(...) {
      cat(sprintf("BLPSimulation: %d markets, %d products\n", self$T, self$N))
      cat(sprintf("  K1=%d, K2=%d, K3=%d\n", self$K1, self$K2, self$K3))
      invisible(self)
    }
  ),
  private = list(
    delta_ = NULL
  )
)

#' @title BLP Simulation Results
#' @description Results from a BLP simulation equilibrium computation.
#' @export
BLPSimulationResults <- R6::R6Class("BLPSimulationResults",
  public = list(
    #' @field simulation The originating BLPSimulation
    simulation = NULL,
    #' @field product_data Updated product data with equilibrium prices/shares
    product_data = NULL,
    #' @field delta Mean utilities at equilibrium
    delta = NULL,
    #' @field costs Marginal costs
    costs = NULL,
    #' @field prices Equilibrium prices
    prices = NULL,
    #' @field shares Equilibrium shares
    shares = NULL,

    #' @description Create simulation results
    initialize = function(simulation, product_data, delta, costs, prices, shares) {
      self$simulation <- simulation
      self$product_data <- product_data
      self$delta <- delta
      self$costs <- costs
      self$prices <- prices
      self$shares <- shares
    },

    #' @description Convert to a BLPProblem for estimation
    #' @param product_formulations Optional override formulations
    #' @param add_instruments Whether to auto-add BLP instruments
    #' @return A BLPProblem object
    #
    # Converts the simulated equilibrium data into a BLPProblem that can be
    # passed to solve(). This is the key step in the "estimation roundtrip":
    # simulate data with known parameters, estimate the model, and check
    # that the estimator recovers the truth. If add_instruments = TRUE,
    # BLP (1995) instruments are automatically constructed -- these are
    # sums of rival/own-firm exogenous characteristics, which provide
    # the exclusion restriction identifying the price coefficient.
    to_problem = function(product_formulations = NULL, add_instruments = TRUE) {
      if (is.null(product_formulations)) {
        product_formulations <- self$simulation$product_formulations
      }

      pd <- self$product_data

      # Auto-generate BLP instruments from exogenous characteristics if
      # the user has not supplied them. The identification logic: prices
      # are endogenous (correlated with xi), so we need instruments that
      # (a) shift markups (and thus prices) but (b) are uncorrelated with xi.
      # BLP instruments satisfy this because rival product characteristics
      # affect the competitive environment (and hence equilibrium markups)
      # but are exogenous to product j's unobserved quality xi_j.
      if (add_instruments) {
        has_demand_iv <- any(grepl("^demand_instruments", names(pd)))
        has_supply_iv <- any(grepl("^supply_instruments", names(pd)))

        # For demand: exclude price/shares (endogenous) from the instrument
        # basis, then compute own-firm and rival-firm characteristic sums.
        if (!has_demand_iv) {
          X1 <- product_formulations[[1]]$build_matrix(pd)
          x1_names <- colnames(X1)
          exog_cols <- which(!grepl("prices|shares", x1_names, ignore.case = TRUE) &
                              x1_names != "(Intercept)")
          if (length(exog_cols) > 0) {
            X_exog <- X1[, exog_cols, drop = FALSE]
            demand_iv <- build_blp_instruments(
              X_exog, pd$market_ids, pd$firm_ids
            )
            for (k in seq_len(ncol(demand_iv))) {
              pd[[paste0("demand_instruments", k - 1)]] <- demand_iv[, k]
            }
          }
        }

        # For supply: analogous instruments from cost-side characteristics
        # (X3), used to form moment conditions E[Z_s' omega] = 0.
        if (!has_supply_iv && length(product_formulations) >= 3) {
          X3 <- product_formulations[[3]]$build_matrix(pd)
          x3_names <- colnames(X3)
          exog_cols3 <- which(!grepl("prices|shares", x3_names, ignore.case = TRUE))
          if (length(exog_cols3) > 0) {
            X_exog3 <- X3[, exog_cols3, drop = FALSE]
            supply_iv <- build_blp_instruments(
              X_exog3, pd$market_ids, pd$firm_ids
            )
            for (k in seq_len(ncol(supply_iv))) {
              pd[[paste0("supply_instruments", k - 1)]] <- supply_iv[, k]
            }
          }
        }
      }

      # Pass integration from the simulation so the problem has nodes/weights
      integration <- self$simulation$.__enclos_env__$private$integration_

      blp_problem(
        product_formulations = product_formulations,
        product_data = pd,
        agent_formulation = self$simulation$agent_formulation,
        agent_data = if (!is.null(self$simulation$agents$original_data))
          self$simulation$agents$original_data else NULL,
        integration = integration,
        rc_types = self$simulation$rc_types,
        costs_type = self$simulation$costs_type,
        add_exogenous = TRUE
      )
    },

    #' @description Print simulation results
    print = function(...) {
      cat("BLPSimulationResults\n")
      cat(sprintf("  Markets: %d, Products: %d\n",
                  self$simulation$T, nrow(self$product_data)))
      cat(sprintf("  Price range: [%.4f, %.4f]\n",
                  min(self$prices), max(self$prices)))
      cat(sprintf("  Share range: [%.6f, %.6f]\n",
                  min(self$shares), max(self$shares)))
      invisible(self)
    }
  )
)

#' Create a BLP Simulation
#'
#' @param product_formulations List of BLPFormulation objects
#' @param product_data Data frame with market_ids, firm_ids, and characteristics
#' @param beta Demand linear coefficients
#' @param sigma Optional Cholesky root of RC covariance
#' @param pi Optional demographics interaction
#' @param gamma Optional supply coefficients
#' @param rho Optional nesting parameters
#' @param agent_formulation Optional demographics formulation
#' @param agent_data Optional agent data
#' @param integration Optional BLPIntegration
#' @param xi Optional demand errors
#' @param omega Optional supply errors
#' @param xi_variance Variance of xi
#' @param omega_variance Variance of omega
#' @param correlation Correlation between xi and omega
#' @param rc_types Character vector of RC types
#' @param costs_type "linear" or "log"
#' @param seed Random seed
#' @return A BLPSimulation object
#' @export
blp_simulation <- function(product_formulations, product_data,
                            beta, sigma = NULL, pi = NULL,
                            gamma = NULL, rho = NULL,
                            agent_formulation = NULL, agent_data = NULL,
                            integration = NULL,
                            xi = NULL, omega = NULL,
                            xi_variance = 1, omega_variance = 1,
                            correlation = 0.9,
                            rc_types = NULL, costs_type = "linear",
                            seed = NULL) {
  BLPSimulation$new(product_formulations, product_data, beta,
                     sigma, pi, gamma, rho,
                     agent_formulation, agent_data, integration,
                     xi, omega, xi_variance, omega_variance,
                     correlation, rc_types, costs_type, seed)
}

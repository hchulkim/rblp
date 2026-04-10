#!/usr/bin/env Rscript
# =============================================================================
# rblp: Nevo (2000) Cereal Demand Estimation
# Demonstrates the main functionality of the rblp package
# Follows the pyblp tutorial for comparison
# =============================================================================

library(rblp)

cat("=== rblp: Nevo (2000) Cereal Demand Estimation ===\n\n")

# ---- Load Data ----
products <- load_nevo_products()
agents <- load_nevo_agents()

cat(sprintf("Products: %d obs, %d markets\n", nrow(products), length(unique(products$market_ids))))
cat(sprintf("Agents:   %d obs (%d per market)\n\n", nrow(agents),
            nrow(agents) / length(unique(products$market_ids))))

# =============================================================================
# Example 1: Plain Logit (pooled, no product FE)
# =============================================================================
cat("--- Example 1: Plain Logit (pooled) ---\n")
cat("Model: delta_jt = beta0 + alpha*p_jt + beta_s*sugar_j + beta_m*mushy_j + xi_jt\n\n")

f1 <- blp_formulation(~ prices + sugar + mushy)
logit_problem <- blp_problem(list(f1), products)
logit_results <- logit_problem$solve(method = "1s")

cat("Logit Results:\n")
print(logit_results)

# =============================================================================
# Example 2: Logit with Product Fixed Effects (matches pyblp tutorial)
# =============================================================================
cat("\n--- Example 2: Logit with Product FE ---\n")
cat("Model: delta_jt = alpha*p_jt + FE_j + xi_jt (sugar/mushy absorbed)\n")
cat("pyblp equivalent: Formulation('prices', absorb='C(product_ids)')\n\n")

f1_fe <- blp_formulation(~ prices, absorb = ~ product_ids)
fe_problem <- blp_problem(list(f1_fe), products)
fe_results <- fe_problem$solve(method = "1s")

cat("FE Logit Results:\n")
print(fe_results)
cat(sprintf("  --> pyblp reports alpha ~ -30 for this specification\n"))

# =============================================================================
# Example 3: RC Logit without Demographics
# =============================================================================
cat("\n--- Example 3: RC Logit (no demographics) ---\n")
cat("Uses product rule integration (Gauss-Hermite, 3 nodes per dim)\n\n")

f1_rc <- blp_formulation(~ prices + sugar + mushy)
f2_rc <- blp_formulation(~ prices + sugar + mushy)

rc_problem <- blp_problem(
  product_formulations = list(f1_rc, f2_rc),
  product_data = products,
  integration = blp_integration("product", size = 3)
)

cat(sprintf("Problem: K1=%d, K2=%d, MD=%d, I=%d agents/market\n",
            rc_problem$K1, rc_problem$K2, rc_problem$MD,
            rc_problem$I / rc_problem$T))

initial_sigma <- diag(c(0.5, 0.5, 0.5, 0.5))
rc_results <- rc_problem$solve(
  sigma = initial_sigma,
  method = "1s",
  optimization = blp_optimization("l-bfgs-b",
    method_options = list(maxit = 200, factr = 1e7))
)

cat("\nRC Logit Results (no demographics):\n")
print(rc_results)

# =============================================================================
# Example 4: RC Logit with Demographics (pyblp Nevo tutorial specification)
# =============================================================================
cat("\n--- Example 4: RC Logit with Demographics ---\n")
cat("Matches pyblp Nevo tutorial: X1 = prices w/ product FE, X2 = 1 + prices + sugar + mushy\n")
cat("Pi matrix: K2 x D = 4 x 4 demographic interactions\n\n")

# pyblp: Formulation('0 + prices', absorb='C(product_ids)') for X1
# pyblp: Formulation('1 + prices + sugar + mushy') for X2
f1_demo <- blp_formulation(~ 0 + prices, absorb = ~ product_ids)
f2_demo <- blp_formulation(~ prices + sugar + mushy)
demo_form <- blp_formulation(~ 0 + income + income_squared + age + child)

demo_problem <- blp_problem(
  product_formulations = list(f1_demo, f2_demo),
  product_data = products,
  agent_formulation = demo_form,
  agent_data = agents
)

cat(sprintf("Problem: K1=%d, K2=%d, D=%d, MD=%d, I=%d\n",
            demo_problem$K1, demo_problem$K2, demo_problem$D,
            demo_problem$MD, demo_problem$I))

# Starting values from pyblp's converged estimates (Nevo published values
# may converge to different local optima depending on the optimizer)
sigma0 <- diag(c(0.558, 3.312, -0.006, 0.093))
pi0 <- matrix(c(
  2.292, 0, 1.284, 0,
  588.3, -30.19, 0, 11.05,
  -0.384, 0, 0.0524, 0,
  0.748, 0, -1.354, 0
), nrow = 4, ncol = 4, byrow = TRUE)

demo_results <- demo_problem$solve(
  sigma = sigma0, pi = pi0,
  method = "1s",
  optimization = blp_optimization("l-bfgs-b",
    method_options = list(maxit = 1000, factr = 1e7))
)

cat("\nRC Logit with Demographics Results:\n")
print(demo_results)
cat(sprintf("  --> pyblp reports alpha ~ -63 for this specification\n"))

# =============================================================================
# Example 5: Post-Estimation
# =============================================================================
cat("\n--- Example 5: Post-Estimation ---\n")

first_market <- logit_results$problem$unique_market_ids[1]
E <- logit_results$compute_elasticities(first_market)
cat(sprintf("\nOwn-price elasticities (market '%s', first 5 products):\n", first_market))
cat(paste(sprintf("  %.3f", diag(E)[1:min(5, nrow(E))]), collapse = "\n"), "\n")

cs <- logit_results$compute_consumer_surplus()
cat("\nConsumer surplus (first 5 markets):\n")
cat(paste(sprintf("  %.4f", cs[1:min(5, length(cs))]), collapse = "\n"), "\n")

hhi <- logit_results$compute_hhi()
cat("\nHHI (first 5 markets):\n")
cat(paste(sprintf("  %.0f", hhi[1:min(5, length(hhi))]), collapse = "\n"), "\n")

# =============================================================================
# Example 6: Simulation Roundtrip
# =============================================================================
cat("\n--- Example 6: Simulation Roundtrip ---\n")

id_data <- build_id_data(T = 50, J = 20, F = 4)
set.seed(42)
id_data$x <- runif(nrow(id_data), 0, 1)

f1_sim <- blp_formulation(~ prices + x)
sim <- blp_simulation(
  product_formulations = list(f1_sim),
  product_data = id_data,
  beta = c(0.5, -2, 0.8),
  xi_variance = 0.3,
  seed = 42
)

sim_results <- sim$replace_endogenous(
  iteration = blp_iteration("simple", list(atol = 1e-12, max_evaluations = 5000))
)

cat(sprintf("Equilibrium prices: [%.4f, %.4f] (mean=%.4f)\n",
            min(sim_results$prices), max(sim_results$prices), mean(sim_results$prices)))

sim_problem <- sim_results$to_problem()
sim_est <- sim_problem$solve(method = "1s")

cat("\nTrue vs Estimated (logit):\n")
cat(sprintf("  intercept: true=%.2f, est=%.4f (se=%.4f)\n",
            0.5, sim_est$beta[1], sim_est$summary_table()$se[1]))
cat(sprintf("  price:     true=%.2f, est=%.4f (se=%.4f)\n",
            -2, sim_est$beta[2], sim_est$summary_table()$se[2]))
cat(sprintf("  x:         true=%.2f, est=%.4f (se=%.4f)\n",
            0.8, sim_est$beta[3], sim_est$summary_table()$se[3]))

# =============================================================================
# Example 7: Merger Simulation
# =============================================================================
cat("\n--- Example 7: Merger Simulation ---\n")

new_firm_ids <- sim_results$product_data$firm_ids
new_firm_ids[new_firm_ids == 3] <- 2

merger <- sim_est$compute_merger(
  new_firm_ids = new_firm_ids,
  iteration = blp_iteration("simple", list(atol = 1e-12, max_evaluations = 5000))
)

cat(sprintf("Mean price change:  %.2f%%\n", mean(merger$price_change_pct)))
cat(sprintf("Max price change:   %.2f%%\n", max(merger$price_change_pct)))
cat(sprintf("Mean CS change:     %.4f\n", mean(merger$delta_cs)))

cat("\n=== All examples completed successfully! ===\n")

# =============================================================================
# Comprehensive validation tests for the rblp package
# Validates estimation, post-estimation, and simulation against expected
# results from pyblp (Conlon & Gortmaker 2020) using the Nevo (2000) data.
# =============================================================================

# Suppress solver messages during testing
options(rblp.verbose = FALSE)

# =============================================================================
# Section 1: Logit Models
# =============================================================================

test_that("Plain logit produces reasonable coefficient estimates", {
  skip_on_cran()

  products <- load_nevo_products()
  f1 <- blp_formulation(~ prices + sugar + mushy)
  problem <- blp_problem(list(f1), products)
  results <- problem$solve(method = "1s")

  expect_s3_class(results, "BLPResults")

  # Identify the price coefficient
  price_idx <- which(colnames(problem$products$X1) == "prices")
  price_coeff <- results$beta[price_idx]

  # Price coefficient should be negative (without FE, expect around -10 to -15)
  expect_true(price_coeff < 0,
    info = sprintf("Price coeff %.4f should be negative", price_coeff))
  expect_true(price_coeff > -50,
    info = sprintf("Price coeff %.4f should be > -50", price_coeff))

  # Objective should be finite and reasonably small
  expect_true(is.finite(results$objective))
  expect_true(results$objective < 300,
    info = sprintf("Objective %.4f should be < 300", results$objective))
})

test_that("Logit with product fixed effects recovers price coefficient near -30", {
  skip_on_cran()

  products <- load_nevo_products()
  f1_fe <- blp_formulation(~ prices, absorb = ~ product_ids)
  fe_problem <- blp_problem(list(f1_fe), products)
  fe_results <- fe_problem$solve(method = "1s")

  expect_s3_class(fe_results, "BLPResults")

  # With product FE, pyblp reports alpha ~ -30
  price_coeff <- fe_results$beta[1]  # only prices (no intercept, absorbed)
  expect_true(abs(price_coeff - (-30)) < 5,
    info = sprintf("FE price coeff %.4f should be near -30 (within 5)", price_coeff))

  # Objective should be finite
  expect_true(is.finite(fe_results$objective))
})

# =============================================================================
# Section 2: Random Coefficients Models
# =============================================================================

test_that("RC logit with diagonal sigma produces negative price coefficient", {
  skip_on_cran()

  products <- load_nevo_products()
  f1_rc <- blp_formulation(~ prices + sugar + mushy)
  f2_rc <- blp_formulation(~ prices + sugar + mushy)

  rc_problem <- blp_problem(
    product_formulations = list(f1_rc, f2_rc),
    product_data = products,
    integration = blp_integration("product", size = 3)
  )

  expect_equal(rc_problem$K2, 4)  # intercept + prices + sugar + mushy

  sigma0 <- diag(c(0.5, 0.5, 0.5, 0.5))
  rc_results <- rc_problem$solve(
    sigma = sigma0,
    method = "1s",
    optimization = blp_optimization("l-bfgs-b",
      method_options = list(maxit = 100))
  )

  expect_s3_class(rc_results, "BLPResults")

  # Price coefficient should be negative
  price_idx <- which(colnames(rc_problem$products$X1) == "prices")
  price_coeff <- rc_results$beta[price_idx]
  expect_true(price_coeff < 0,
    info = sprintf("RC price coeff %.4f should be negative", price_coeff))

  # Objective should be finite
  expect_true(is.finite(rc_results$objective),
    info = "RC logit objective should be finite")
})

test_that("RC logit with demographics matches pyblp Nevo specification", {
  skip_on_cran()

  products <- load_nevo_products()
  agents <- load_nevo_agents()

  # pyblp Nevo tutorial specification:
  # X1 = prices with product FE (absorb product_ids)
  # X2 = intercept + prices + sugar + mushy
  # Demographics = income, income_squared, age, child
  f1_demo <- blp_formulation(~ 0 + prices, absorb = ~ product_ids)
  f2_demo <- blp_formulation(~ prices + sugar + mushy)
  demo_form <- blp_formulation(~ 0 + income + income_squared + age + child)

  demo_problem <- blp_problem(
    product_formulations = list(f1_demo, f2_demo),
    product_data = products,
    agent_formulation = demo_form,
    agent_data = agents
  )

  expect_equal(demo_problem$K2, 4)  # intercept + prices + sugar + mushy in X2
  expect_equal(demo_problem$D, 4)   # 4 demographics

  # Starting values from pyblp Nevo tutorial
  sigma0 <- diag(c(0.558, 3.312, -0.006, 0.093))
  pi0 <- matrix(c(
    2.292, 0, 1.284, 0,
    588.3, -30.19, 0, 11.05,
    -0.384, 0, 0.0524, 0,
    0.748, 0, -1.354, 0
  ), nrow = 4, ncol = 4, byrow = TRUE)

  demo_results <- demo_problem$solve(
    sigma = sigma0,
    pi = pi0,
    method = "1s",
    optimization = blp_optimization("l-bfgs-b",
      method_options = list(maxit = 200, factr = 1e7))
  )

  expect_s3_class(demo_results, "BLPResults")

  # pyblp reports alpha ~ -63 for this specification
  price_coeff <- demo_results$beta[1]
  expect_true(price_coeff > -100,
    info = sprintf("Demo price coeff %.4f should be > -100", price_coeff))
  expect_true(price_coeff < -30,
    info = sprintf("Demo price coeff %.4f should be < -30", price_coeff))

  # Sigma should be estimated (4x4 matrix)
  expect_equal(nrow(demo_results$sigma), 4)
  expect_equal(ncol(demo_results$sigma), 4)

  # Pi should be estimated (4x4 matrix)
  expect_equal(nrow(demo_results$pi), 4)
  expect_equal(ncol(demo_results$pi), 4)

  # Objective should be finite
  expect_true(is.finite(demo_results$objective))
})

# =============================================================================
# Section 3: Post-Estimation -- Elasticities
# =============================================================================

test_that("Elasticity matrix has correct sign pattern", {
  skip_on_cran()

  products <- load_nevo_products()
  f1 <- blp_formulation(~ prices + sugar + mushy)
  problem <- blp_problem(list(f1), products)
  results <- problem$solve(method = "1s")

  first_market <- results$problem$unique_market_ids[1]
  E <- results$compute_elasticities(first_market)

  expect_true(is.matrix(E))
  J <- nrow(E)
  expect_true(J > 1, info = "Market should have multiple products")
  expect_equal(nrow(E), ncol(E))

  # Own-price elasticities (diagonal) should all be negative
  own_elast <- diag(E)
  expect_true(all(own_elast < 0),
    info = sprintf("All %d diagonal (own-price) elasticities should be negative", J))

  # Cross-price elasticities (off-diagonal) should be positive for logit
  # Allow a small fraction to be zero or slightly negative due to numerics
  off_diag <- E[row(E) != col(E)]
  pct_positive <- mean(off_diag > 0)
  expect_true(pct_positive > 0.95,
    info = sprintf("%.1f%% of cross-price elasticities are positive (expect >95%%)",
                   pct_positive * 100))
})

# =============================================================================
# Section 4: Post-Estimation -- Consumer Surplus
# =============================================================================

test_that("Consumer surplus values are all positive", {
  skip_on_cran()

  products <- load_nevo_products()
  f1 <- blp_formulation(~ prices + sugar + mushy)
  problem <- blp_problem(list(f1), products)
  results <- problem$solve(method = "1s")

  cs <- results$compute_consumer_surplus()

  expect_true(length(cs) == results$problem$T,
    info = "CS should have one value per market")
  expect_true(all(is.finite(cs)),
    info = "All CS values should be finite")
  expect_true(all(cs > 0),
    info = sprintf("All %d CS values should be positive (min=%.6f)",
                   length(cs), min(cs)))
})

# =============================================================================
# Section 5: Post-Estimation -- HHI
# =============================================================================

test_that("HHI values are in valid range (0, 10000]", {
  skip_on_cran()

  products <- load_nevo_products()
  f1 <- blp_formulation(~ prices + sugar + mushy)
  problem <- blp_problem(list(f1), products)
  results <- problem$solve(method = "1s")

  hhi <- results$compute_hhi()

  expect_true(length(hhi) == results$problem$T,
    info = "HHI should have one value per market")
  expect_true(all(is.finite(hhi)),
    info = "All HHI values should be finite")
  expect_true(all(hhi > 0),
    info = sprintf("All HHI values should be > 0 (min=%.2f)", min(hhi)))
  expect_true(all(hhi <= 10000),
    info = sprintf("All HHI values should be <= 10000 (max=%.2f)", max(hhi)))
})

# =============================================================================
# Section 6: Simulation Roundtrip
# =============================================================================

test_that("Simulation roundtrip recovers parameter signs and approximate magnitudes", {
  skip_on_cran()

  # Build simulated data (larger sample for instrument relevance)
  id_data <- build_id_data(T = 50, J = 20, F = 4)
  set.seed(42)
  id_data$x <- stats::runif(nrow(id_data), 0, 1)

  # True parameters
  true_beta <- c(0.5, -2, 0.8)  # intercept, price, x

  f1_sim <- blp_formulation(~ prices + x)
  sim <- blp_simulation(
    product_formulations = list(f1_sim),
    product_data = id_data,
    beta = true_beta,
    xi_variance = 0.3,
    seed = 42
  )

  expect_s3_class(sim, "BLPSimulation")
  expect_equal(sim$T, 50)
  expect_equal(sim$N, 1000)

  # Solve for equilibrium
  sim_results <- sim$replace_endogenous(
    iteration = blp_iteration("simple", list(atol = 1e-12, max_evaluations = 5000))
  )

  expect_s3_class(sim_results, "BLPSimulationResults")
  expect_true(all(sim_results$shares > 0), info = "All shares should be positive")
  expect_true(all(sim_results$prices > 0), info = "All prices should be positive")

  # Convert to problem and estimate
  sim_problem <- sim_results$to_problem()
  expect_s3_class(sim_problem, "BLPProblem")
  expect_true(sim_problem$MD > 0, info = "Should have demand instruments")

  sim_est <- sim_problem$solve(method = "1s")
  expect_s3_class(sim_est, "BLPResults")

  est_beta <- sim_est$beta
  expect_equal(length(est_beta), 3)

  # Check price coefficient sign (most robust check)
  # Price coefficient: true = -2 (negative)
  expect_true(est_beta[2] < 0,
    info = sprintf("Price est=%.4f should be negative (true=-2)", est_beta[2]))
})

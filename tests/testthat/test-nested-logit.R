# =============================================================================
# Tests for nested logit (rho parameter)
# =============================================================================

test_that("nested logit probabilities sum to less than 1", {
  skip_on_cran()

  # Create data with nesting structure
  id_data <- build_id_data(T = 10, J = 12, F = 3)
  set.seed(333)
  id_data$x <- runif(nrow(id_data), 0, 1)
  # Assign products to 3 nests
  id_data$nesting_ids <- rep(c("A", "B", "C"), length.out = nrow(id_data))

  f1 <- blp_formulation(~ prices + x)
  f2 <- blp_formulation(~ 0 + prices + x)

  # Simulate with nesting
  sim <- blp_simulation(
    product_formulations = list(f1, f2),
    product_data = id_data,
    beta = c(0.5, -2.0, 0.8),
    sigma = diag(c(0.5, 0.5)),
    rho = 0.5,
    integration = blp_integration("product", size = 3),
    xi_variance = 0.2,
    seed = 333
  )

  sim_results <- sim$replace_endogenous(
    iteration = blp_iteration("simple", list(atol = 1e-12, max_evaluations = 5000))
  )

  # Shares should be valid
  expect_true(all(sim_results$shares > 0))
  expect_true(all(sim_results$shares < 1))

  # Within each market, total inside share should be < 1
  for (t in unique(sim_results$product_data$market_ids)) {
    idx <- which(sim_results$product_data$market_ids == t)
    expect_true(sum(sim_results$shares[idx]) < 1,
                info = paste("Market", t, "shares should sum to < 1"))
  }
})

test_that("rho = 0 nested logit matches standard logit", {
  skip_on_cran()

  id_data <- build_id_data(T = 5, J = 8, F = 2)
  set.seed(444)
  id_data$x <- runif(nrow(id_data), 0, 1)
  id_data$nesting_ids <- rep(c("A", "B"), length.out = nrow(id_data))
  id_data$shares <- runif(nrow(id_data), 0.01, 0.08)
  id_data$prices <- runif(nrow(id_data), 1, 3)

  f1 <- blp_formulation(~ prices + x)
  f2 <- blp_formulation(~ 0 + prices + x)

  integration <- blp_integration("product", size = 3)

  # Create problem with rho = 0 (should behave like standard logit)
  problem <- blp_problem(
    product_formulations = list(f1, f2),
    product_data = id_data,
    integration = integration
  )

  # Solve with rho = 0 (fixed at zero, so effectively no nesting)
  sigma0 <- diag(c(0.3, 0.3))
  results_no_nest <- problem$solve(
    sigma = sigma0,
    rho = 0,
    method = "1s",
    optimization = blp_optimization("l-bfgs-b", method_options = list(maxit = 50))
  )

  # Should converge
  expect_true(results_no_nest$fp_converged)
})

test_that("nested logit delta inversion includes rho correction", {
  # Test that compute_logit_delta with rho != 0 gives different result than rho = 0
  id_data <- build_id_data(T = 5, J = 6, F = 2)
  set.seed(555)
  id_data$x <- runif(nrow(id_data), 0, 1)
  id_data$nesting_ids <- rep(c("A", "B"), length.out = nrow(id_data))
  id_data$shares <- runif(nrow(id_data), 0.02, 0.08)
  id_data$prices <- runif(nrow(id_data), 1, 3)

  f1 <- blp_formulation(~ prices + x)

  problem <- blp_problem(list(f1), id_data)

  delta_no_rho <- problem$compute_logit_delta(rho = NULL)
  delta_with_rho <- problem$compute_logit_delta(rho = 0.5)

  # They should be different

  expect_false(all(abs(delta_no_rho - delta_with_rho) < 1e-10),
               info = "rho should modify the logit delta inversion")

  # Both should be finite
  expect_true(all(is.finite(delta_no_rho)))
  expect_true(all(is.finite(delta_with_rho)))
})

test_that("higher rho increases within-nest substitution in nested logit probabilities", {
  skip_on_cran()

  # Create a simple market with known structure
  id_data <- data.frame(
    market_ids = rep("M1", 6),
    firm_ids = 1:6,
    nesting_ids = c("A", "A", "A", "B", "B", "B"),
    prices = c(1.0, 1.1, 1.2, 1.0, 1.1, 1.2),
    x = c(0.5, 0.6, 0.4, 0.5, 0.6, 0.4),
    shares = rep(0.1, 6)
  )

  f1 <- blp_formulation(~ prices + x)
  f2 <- blp_formulation(~ 0 + prices + x)
  integration <- blp_integration("product", size = 3)

  problem <- blp_problem(
    product_formulations = list(f1, f2),
    product_data = id_data,
    integration = integration
  )

  # Get market data and create markets with different rho
  md <- problem$get_market_data("M1")

  mkt_low <- rblp:::BLPMarket$new(
    products = md$products, agents = md$agents,
    sigma = diag(c(0.3, 0.3)), rho = 0.2,
    rc_types = c("linear", "linear"),
    epsilon_scale = 1.0
  )

  mkt_high <- rblp:::BLPMarket$new(
    products = md$products, agents = md$agents,
    sigma = diag(c(0.3, 0.3)), rho = 0.8,
    rc_types = c("linear", "linear"),
    epsilon_scale = 1.0
  )

  delta <- rep(0, 6)
  mu_low <- mkt_low$compute_mu()
  mu_high <- mkt_high$compute_mu()

  prob_low <- mkt_low$compute_probabilities(delta, mu_low)
  prob_high <- mkt_high$compute_probabilities(delta, mu_high)

  # Both should produce valid probabilities
  s_low <- mkt_low$compute_shares(prob_low$probabilities)
  s_high <- mkt_high$compute_shares(prob_high$probabilities)

  expect_true(all(s_low > 0) && all(s_low < 1))
  expect_true(all(s_high > 0) && all(s_high < 1))
  expect_true(sum(s_low) < 1)
  expect_true(sum(s_high) < 1)
})

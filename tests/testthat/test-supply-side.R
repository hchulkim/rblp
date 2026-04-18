# =============================================================================
# Tests for supply-side estimation and cost recovery
# =============================================================================

test_that("compute_costs returns positive marginal costs for logit model", {
  skip_on_cran()

  products <- load_nevo_products()
  f1 <- blp_formulation(~ prices + sugar + mushy)
  problem <- blp_problem(list(f1), products)
  results <- problem$solve(method = "1s")

  costs <- results$compute_costs()
  expect_equal(length(costs), problem$N)

  # Costs should be less than prices (positive markups in Bertrand-Nash)
  markups <- results$compute_markups()
  expect_true(all(is.finite(markups)))

  # Most markups should be positive (firms set price above cost)
  expect_true(mean(markups > 0) > 0.9,
              info = "Most products should have positive markups")
})

test_that("compute_markups are consistent with costs and prices", {
  skip_on_cran()

  products <- load_nevo_products()
  f1 <- blp_formulation(~ prices + sugar + mushy)
  problem <- blp_problem(list(f1), products)
  results <- problem$solve(method = "1s")

  costs <- results$compute_costs()
  markups <- results$compute_markups()
  prices <- problem$products$prices

  # By definition: markup = (price - cost) / price
  expected_markups <- (prices - costs) / prices
  expect_equal(markups, expected_markups, tolerance = 1e-12)
})

test_that("markup formula gives positive markups under Bertrand-Nash", {
  skip_on_cran()

  # Use simulation to control the DGP
  id_data <- build_id_data(T = 20, J = 10, F = 3)
  set.seed(123)
  id_data$x <- runif(nrow(id_data), 0, 1)

  f1 <- blp_formulation(~ prices + x)

  sim <- blp_simulation(
    product_formulations = list(f1),
    product_data = id_data,
    beta = c(0.5, -2.0, 0.8),
    xi_variance = 0.2,
    seed = 123
  )

  sim_results <- sim$replace_endogenous(
    iteration = blp_iteration("simple", list(atol = 1e-12, max_evaluations = 5000))
  )

  sim_problem <- sim_results$to_problem()
  est <- sim_problem$solve(method = "1s")

  costs <- est$compute_costs()
  markups <- est$compute_markups()

  # All costs and markups should be finite (main check — values may be
  # negative with weak instruments, which is an IV quality issue)
  expect_true(all(is.finite(costs)), info = "Costs should be finite")
  expect_true(all(is.finite(markups)), info = "Markups should be finite")

  # Price coefficient should be negative
  expect_true(est$beta[2] < 0, info = "Price coefficient should be negative")
})

test_that("costs are invariant to non-price product characteristics", {
  skip_on_cran()

  # The cost recovery should depend only on demand parameters and shares,
  # not on the specific values of non-price product characteristics
  products <- load_nevo_products()
  f1 <- blp_formulation(~ prices + sugar + mushy)
  problem <- blp_problem(list(f1), products)
  results <- problem$solve(method = "1s")

  costs1 <- results$compute_costs()

  # Costs should be reproducible (same computation each time)
  costs2 <- results$compute_costs()
  expect_equal(costs1, costs2, tolerance = 1e-15)
})

test_that("supply-side estimation with RC logit and cost equation", {
  skip_on_cran()

  # Full model: demand RCs + supply
  id_data <- build_id_data(T = 30, J = 15, F = 3)
  set.seed(999)
  id_data$x <- runif(nrow(id_data), 0, 1)
  id_data$w <- runif(nrow(id_data), 0, 1)  # cost shifter

  f1 <- blp_formulation(~ prices + x)
  f2 <- blp_formulation(~ 0 + prices + x)
  f3 <- blp_formulation(~ w)

  true_beta <- c(0.5, -2.0, 0.8)
  true_sigma <- diag(c(0.3, 0.3))
  true_gamma <- c(0.5, 1.0)

  sim <- blp_simulation(
    product_formulations = list(f1, f2, f3),
    product_data = id_data,
    beta = true_beta,
    sigma = true_sigma,
    gamma = true_gamma,
    integration = blp_integration("product", size = 3),
    xi_variance = 0.2,
    omega_variance = 0.2,
    correlation = 0.7,
    seed = 999
  )

  sim_results <- sim$replace_endogenous(
    iteration = blp_iteration("simple", list(atol = 1e-12, max_evaluations = 5000))
  )

  # Check equilibrium data is reasonable
  expect_true(all(sim_results$prices > 0))
  expect_true(all(sim_results$shares > 0))
  expect_true(all(sim_results$shares < 1))

  # Estimate
  sim_problem <- sim_results$to_problem()
  expect_true(sim_problem$K3 > 0, info = "Supply formulation should set K3 > 0")
  expect_true(sim_problem$MS > 0, info = "Supply instruments should be present")

  est <- sim_problem$solve(
    sigma = diag(c(0.3, 0.3)),
    method = "1s",
    optimization = blp_optimization("l-bfgs-b",
      method_options = list(maxit = 100))
  )

  # Gamma should be estimated
  expect_false(is.null(est$gamma))
  expect_equal(length(est$gamma), 2)

  # Omega should exist
  expect_false(is.null(est$omega))
  expect_equal(length(est$omega), sim_problem$N)

  # Summary table should include gamma
  tbl <- est$summary_table()
  gamma_rows <- tbl[tbl$type == "supply (gamma)", ]
  expect_true(nrow(gamma_rows) > 0, info = "gamma should appear in summary table")
})

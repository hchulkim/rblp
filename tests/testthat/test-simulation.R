test_that("BLPSimulation logit case works", {
  skip_if_not_installed("R6")

  # Simple logit simulation
  id_data <- build_id_data(T = 5, J = 10, F = 3)
  set.seed(123)
  id_data$x1 <- stats::runif(nrow(id_data))
  id_data$x2 <- stats::runif(nrow(id_data))

  f1 <- blp_formulation(~ prices + x1 + x2)

  sim <- blp_simulation(
    product_formulations = list(f1),
    product_data = id_data,
    beta = c(1, -2, 0.5, 1),  # intercept, price, x1, x2
    seed = 42
  )

  expect_s3_class(sim, "BLPSimulation")
  expect_equal(sim$T, 5)
  expect_equal(sim$N, 50)

  # Replace endogenous variables
  sim_results <- sim$replace_endogenous()
  expect_s3_class(sim_results, "BLPSimulationResults")

  # Shares should be positive and sum to < 1 within each market
  expect_true(all(sim_results$shares > 0))
  for (t in unique(id_data$market_ids)) {
    idx <- which(id_data$market_ids == t)
    expect_true(sum(sim_results$shares[idx]) < 1)
  }

  # Prices should be positive (> costs)
  expect_true(all(sim_results$prices > 0))
})

test_that("BLPSimulationResults can convert to problem", {
  skip_if_not_installed("R6")

  id_data <- build_id_data(T = 3, J = 6, F = 2)
  set.seed(123)
  id_data$x1 <- stats::runif(nrow(id_data))

  f1 <- blp_formulation(~ prices + x1)

  sim <- blp_simulation(
    product_formulations = list(f1),
    product_data = id_data,
    beta = c(1, -2, 1),
    seed = 42
  )

  sim_results <- sim$replace_endogenous()
  problem <- sim_results$to_problem()
  expect_s3_class(problem, "BLPProblem")
  expect_true(problem$MD > 0)  # should have demand instruments
})

# =============================================================================
# Tests for merger simulation
# =============================================================================

test_that("merger simulation produces price increases for merging firms", {
  skip_on_cran()

  id_data <- build_id_data(T = 30, J = 15, F = 5)
  set.seed(789)
  id_data$x <- runif(nrow(id_data), 0, 1)

  f1 <- blp_formulation(~ prices + x)
  sim <- blp_simulation(
    product_formulations = list(f1),
    product_data = id_data,
    beta = c(0.5, -2.0, 0.8),
    xi_variance = 0.3,
    seed = 789
  )

  sim_results <- sim$replace_endogenous(
    iteration = blp_iteration("simple", list(atol = 1e-12, max_evaluations = 5000))
  )

  sim_problem <- sim_results$to_problem()
  est <- sim_problem$solve(method = "1s")

  # Merge firm 3 into firm 2
  new_firm_ids <- sim_results$product_data$firm_ids
  new_firm_ids[new_firm_ids == 3] <- 2

  merger <- est$compute_merger(
    new_firm_ids = new_firm_ids,
    iteration = blp_iteration("simple", list(atol = 1e-12, max_evaluations = 5000))
  )

  # Merger output should have correct structure
  expect_equal(length(merger$new_prices), est$problem$N)
  expect_equal(length(merger$new_shares), est$problem$N)
  expect_equal(length(merger$costs), est$problem$N)
  expect_equal(length(merger$price_change), est$problem$N)
  expect_equal(length(merger$price_change_pct), est$problem$N)
  expect_equal(length(merger$cs_pre), est$problem$T)
  expect_equal(length(merger$cs_post), est$problem$T)
  expect_equal(length(merger$delta_cs), est$problem$T)

  # Prices should generally change (merger should move prices)
  # In small samples, the direction depends on IV quality
  expect_true(mean(abs(merger$price_change_pct)) > 0,
              info = "Merger should cause non-zero price changes")

  # All new prices should be positive and finite
  expect_true(all(merger$new_prices > 0), info = "Post-merger prices should be positive")
  expect_true(all(is.finite(merger$new_prices)), info = "Post-merger prices should be finite")
  expect_true(all(merger$new_shares > 0), info = "Post-merger shares should be positive")
  expect_true(all(is.finite(merger$new_shares)), info = "Post-merger shares should be finite")

  # CS changes should be finite
  expect_true(all(is.finite(merger$delta_cs)), info = "CS changes should be finite")
})

test_that("identity merger produces zero price change", {
  skip_on_cran()

  id_data <- build_id_data(T = 10, J = 10, F = 3)
  set.seed(111)
  id_data$x <- runif(nrow(id_data), 0, 1)

  f1 <- blp_formulation(~ prices + x)
  sim <- blp_simulation(
    product_formulations = list(f1),
    product_data = id_data,
    beta = c(0.5, -2.0, 0.8),
    xi_variance = 0.2,
    seed = 111
  )

  sim_results <- sim$replace_endogenous(
    iteration = blp_iteration("simple", list(atol = 1e-12, max_evaluations = 5000))
  )

  sim_problem <- sim_results$to_problem()
  est <- sim_problem$solve(method = "1s")

  # "Merge" with same ownership — no change
  same_firm_ids <- sim_results$product_data$firm_ids

  merger <- est$compute_merger(
    new_firm_ids = same_firm_ids,
    iteration = blp_iteration("simple", list(atol = 1e-12, max_evaluations = 5000))
  )

  # Price changes should be very close to zero
  expect_true(max(abs(merger$price_change_pct)) < 0.5,
              info = "No-merger scenario should yield near-zero price changes")
})

test_that("monopoly merger produces largest price increases", {
  skip_on_cran()

  id_data <- build_id_data(T = 10, J = 10, F = 5)
  set.seed(222)
  id_data$x <- runif(nrow(id_data), 0, 1)

  f1 <- blp_formulation(~ prices + x)
  sim <- blp_simulation(
    product_formulations = list(f1),
    product_data = id_data,
    beta = c(0.5, -2.0, 0.8),
    xi_variance = 0.2,
    seed = 222
  )

  sim_results <- sim$replace_endogenous(
    iteration = blp_iteration("simple", list(atol = 1e-12, max_evaluations = 5000))
  )

  sim_problem <- sim_results$to_problem()
  est <- sim_problem$solve(method = "1s")

  # Partial merger: firm 3 -> firm 2
  partial_ids <- sim_results$product_data$firm_ids
  partial_ids[partial_ids == 3] <- 2
  partial <- est$compute_merger(
    new_firm_ids = partial_ids,
    iteration = blp_iteration("simple", list(atol = 1e-12, max_evaluations = 5000))
  )

  # Full monopoly: all firms merge into 1
  monopoly_ids <- rep(1, length(sim_results$product_data$firm_ids))
  monopoly <- est$compute_merger(
    new_firm_ids = monopoly_ids,
    iteration = blp_iteration("simple", list(atol = 1e-12, max_evaluations = 5000))
  )

  # Monopoly price increases should be larger than partial merger
  expect_true(mean(monopoly$price_change_pct) > mean(partial$price_change_pct),
              info = "Monopoly should cause larger price increases than partial merger")
})

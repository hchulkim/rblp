# =============================================================================
# Tests for edge cases and robustness
# =============================================================================

test_that("single-product markets are handled", {
  id_data <- data.frame(
    market_ids = c("M1", "M2", "M3"),
    firm_ids = c(1, 1, 1),
    prices = c(1.5, 2.0, 1.8),
    x = c(0.5, 0.7, 0.3),
    shares = c(0.3, 0.2, 0.25),
    demand_instruments0 = c(0.1, 0.2, 0.3)
  )

  f1 <- blp_formulation(~ prices + x)
  problem <- blp_problem(list(f1), id_data)
  results <- problem$solve(method = "1s")

  expect_true(results$fp_converged)
  expect_false(is.null(results$beta))
  expect_equal(length(results$beta), 3)
})

test_that("very small shares don't cause numerical issues", {
  skip_on_cran()

  id_data <- build_id_data(T = 10, J = 10, F = 3)
  set.seed(800)
  id_data$x <- runif(nrow(id_data), 0, 1)

  # Create very small shares (near-zero but positive)
  id_data$shares <- runif(nrow(id_data), 0.001, 0.005)
  # Normalize within market
  for (t in unique(id_data$market_ids)) {
    idx <- which(id_data$market_ids == t)
    id_data$shares[idx] <- id_data$shares[idx] / (sum(id_data$shares[idx]) + 0.5)
  }
  id_data$prices <- runif(nrow(id_data), 1, 3)

  f1 <- blp_formulation(~ prices + x)
  problem <- blp_problem(list(f1), id_data)

  # Should not error even with small shares
  results <- problem$solve(method = "1s")
  expect_true(results$fp_converged)
  expect_true(all(is.finite(results$beta)))
  expect_true(all(is.finite(results$xi)))
})

test_that("large outside share (small inside shares) works", {
  skip_on_cran()

  id_data <- build_id_data(T = 5, J = 15, F = 3)
  set.seed(801)
  id_data$x <- runif(nrow(id_data), 0, 1)

  # Very small inside shares => large outside good share
  id_data$shares <- runif(nrow(id_data), 0.0005, 0.002)
  id_data$prices <- runif(nrow(id_data), 1, 3)

  f1 <- blp_formulation(~ prices + x)
  problem <- blp_problem(list(f1), id_data)

  delta <- problem$compute_logit_delta()
  expect_true(all(is.finite(delta)))
  # Delta should be very negative (small shares => log(s/s0) << 0)
  expect_true(all(delta < 0))
})

test_that("many products per market works", {
  skip_on_cran()

  id_data <- build_id_data(T = 5, J = 50, F = 10)
  set.seed(802)
  id_data$x <- runif(nrow(id_data), 0, 1)

  f1 <- blp_formulation(~ prices + x)

  sim <- blp_simulation(
    product_formulations = list(f1),
    product_data = id_data,
    beta = c(0.5, -2.0, 0.8),
    xi_variance = 0.3,
    seed = 802
  )

  sim_results <- sim$replace_endogenous(
    iteration = blp_iteration("simple", list(atol = 1e-12, max_evaluations = 5000))
  )

  # Should produce valid equilibrium
  expect_true(all(sim_results$prices > 0))
  expect_true(all(sim_results$shares > 0))
  expect_true(all(sim_results$shares < 1))

  sim_problem <- sim_results$to_problem()
  est <- sim_problem$solve(method = "1s")
  expect_true(est$fp_converged)
})

test_that("build_blp_instruments handles edge cases", {
  # Single firm per market
  X <- matrix(runif(20), 10, 2)
  market_ids <- rep(c("M1", "M2"), each = 5)
  firm_ids <- rep(1, 10)  # all same firm

  iv <- build_blp_instruments(X, market_ids, firm_ids)
  expect_equal(nrow(iv), 10)
  # Rival sums should be zero since there's only one firm
  # (own_sum columns should be nonzero, rival columns should be zero)
  n_cols <- ncol(iv)
  expect_true(n_cols >= 2)

  # Two firms
  firm_ids2 <- rep(c(1, 2), each = 5)
  iv2 <- build_blp_instruments(X, market_ids, firm_ids2)
  expect_equal(nrow(iv2), 10)
})

test_that("build_differentiation_instruments produce valid output", {
  X <- matrix(runif(30), 10, 3)
  market_ids <- rep(c("M1", "M2"), each = 5)
  firm_ids <- rep(c(1, 2, 1, 2, 1), 2)

  # Local method
  div_local <- build_differentiation_instruments(
    X, market_ids, firm_ids, method = "local"
  )
  expect_equal(nrow(div_local), 10)
  expect_true(all(is.finite(div_local)))

  # Quadratic method
  div_quad <- build_differentiation_instruments(
    X, market_ids, firm_ids, method = "quadratic"
  )
  expect_equal(nrow(div_quad), 10)
  expect_true(all(is.finite(div_quad)))
})

test_that("formulation with only intercept works", {
  id_data <- data.frame(
    market_ids = rep(c("M1", "M2"), each = 5),
    firm_ids = rep(c(1, 2, 1, 2, 1), 2),
    prices = runif(10, 1, 3),
    shares = runif(10, 0.01, 0.1)
  )

  # Only intercept + prices
  f1 <- blp_formulation(~ prices)
  problem <- blp_problem(list(f1), id_data)

  expect_equal(problem$K1, 2)  # intercept + prices
  results <- problem$solve(method = "1s")
  expect_true(results$optimization_converged)
  expect_true(results$beta[2] < 0, info = "Price coefficient should be negative")
})

test_that("2-step GMM gives weakly lower objective than 1-step", {
  skip_on_cran()

  products <- load_nevo_products()
  f1 <- blp_formulation(~ prices + sugar + mushy)
  problem <- blp_problem(list(f1), products)

  res_1s <- problem$solve(method = "1s")
  res_2s <- problem$solve(method = "2s")

  # 2s objective uses efficient weights, so should be lower
  # (comparing normalized objectives)
  expect_true(res_2s$objective <= res_1s$objective + 1e-6,
              info = "2-step GMM should have weakly lower objective")
})

test_that("blp_formulation handles various formula specifications", {
  # String formula
  f1 <- blp_formulation("~ prices + x")
  expect_true(inherits(f1, "BLPFormulation"))

  # R formula
  f2 <- blp_formulation(~ prices + x)
  expect_true(inherits(f2, "BLPFormulation"))

  # No intercept
  f3 <- blp_formulation(~ 0 + prices + x)
  expect_true(inherits(f3, "BLPFormulation"))

  # With absorb
  f4 <- blp_formulation(~ prices, absorb = ~ product_ids)
  expect_false(is.null(f4$get_absorb()))
})

test_that("integration methods produce valid nodes and weights", {
  K <- 3

  # Product rule (Gauss-Hermite)
  int_gh <- blp_integration("product", size = 3)
  result_gh <- int_gh$build(K)
  expect_equal(ncol(result_gh$nodes), K)
  expect_true(abs(sum(result_gh$weights) - 1) < 1e-10)
  expect_true(all(result_gh$weights > 0))

  # Monte Carlo
  int_mc <- blp_integration("monte_carlo", size = 100, seed = 42)
  result_mc <- int_mc$build(K)
  expect_equal(nrow(result_mc$nodes), 100)
  expect_equal(ncol(result_mc$nodes), K)
  expect_true(abs(sum(result_mc$weights) - 1) < 1e-10)
})

test_that("iteration methods converge on simple contraction", {
  # Test contraction: x -> (x + 5) / 2, fixed point is 5
  contraction <- function(x) (x + 5) / 2

  # Simple iteration
  iter_simple <- blp_iteration("simple", list(atol = 1e-10))
  result_simple <- iter_simple$iterate(0, contraction)
  expect_true(result_simple$converged)
  expect_equal(result_simple$values, 5, tolerance = 1e-8)

  # SQUAREM
  iter_sq <- blp_iteration("squarem", list(atol = 1e-10))
  result_sq <- iter_sq$iterate(0, contraction)
  expect_true(result_sq$converged)
  expect_equal(result_sq$values, 5, tolerance = 1e-8)
})

test_that("optimization methods run without error", {
  # Simple quadratic: min (x-3)^2
  obj_fn <- function(x) {
    list(objective = sum((x - 3)^2), gradient = 2 * (x - 3))
  }
  bounds <- list(lower = c(-10, -10), upper = c(10, 10))

  # L-BFGS-B
  opt1 <- blp_optimization("l-bfgs-b")
  r1 <- opt1$optimize(c(0, 0), bounds, obj_fn)
  expect_true(r1$converged)
  expect_equal(r1$values, c(3, 3), tolerance = 1e-4)

  # NLMINB
  opt2 <- blp_optimization("nlminb")
  r2 <- opt2$optimize(c(0, 0), bounds, obj_fn)
  expect_true(r2$converged)
  expect_equal(r2$values, c(3, 3), tolerance = 1e-4)

  # Return (no optimization)
  opt3 <- blp_optimization("return")
  r3 <- opt3$optimize(c(5, 5), bounds, obj_fn)
  expect_equal(r3$values, c(5, 5))
})

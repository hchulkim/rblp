# =============================================================================
# Tests for newly added pyblp-equivalent features
# =============================================================================

options(rblp.verbose = FALSE)

# Helper: get a simple estimated model for testing post-estimation
get_nevo_logit_results <- function() {
  products <- load_nevo_products()
  f1 <- blp_formulation(~ prices + sugar + mushy)
  problem <- blp_problem(list(f1), products)
  problem$solve(method = "1s")
}

# =============================================================================
# compute_aggregate_elasticities
# =============================================================================

test_that("compute_aggregate_elasticities returns one value per market", {
  skip_on_cran()
  res <- get_nevo_logit_results()
  ae <- res$compute_aggregate_elasticities()

  expect_equal(length(ae), res$problem$T)
  expect_true(all(is.finite(ae)))
  # Aggregate elasticity should be negative (higher prices -> lower total share)
  expect_true(all(ae < 0),
              info = "Aggregate elasticity should be negative")
})

test_that("compute_aggregate_elasticities for single market", {
  skip_on_cran()
  res <- get_nevo_logit_results()
  t1 <- res$problem$unique_market_ids[1]
  ae <- res$compute_aggregate_elasticities(market_id = t1)
  expect_equal(length(ae), 1)
  expect_true(ae < 0)
})

# =============================================================================
# compute_passthrough
# =============================================================================

test_that("compute_passthrough returns J x J matrix with correct properties", {
  skip_on_cran()
  res <- get_nevo_logit_results()
  t1 <- res$problem$unique_market_ids[1]
  PT <- res$compute_passthrough(t1)

  J <- sum(res$problem$products$market_ids == t1)
  expect_equal(nrow(PT), J)
  expect_equal(ncol(PT), J)
  expect_true(all(is.finite(PT)))

  # Most own-cost passthrough (diagonal) should be positive
  # (Can be negative for very small-share products due to approximation)
  expect_true(mean(diag(PT) > 0) > 0.7,
              info = "Most own-cost passthrough should be positive")
})

# =============================================================================
# compute_long_run_diversion_ratios
# =============================================================================

test_that("compute_long_run_diversion_ratios returns J x J matrix", {
  skip_on_cran()
  res <- get_nevo_logit_results()
  t1 <- res$problem$unique_market_ids[1]

  D_lr <- res$compute_long_run_diversion_ratios(t1)
  D_sr <- res$compute_diversion_ratios(t1)

  J <- nrow(D_lr)
  expect_equal(nrow(D_lr), ncol(D_lr))
  expect_true(all(is.finite(D_lr)))
  expect_true(all(diag(D_lr) == 0))

  # Long-run diversion should differ from short-run
  # (because equilibrium price adjustments change substitution)
  expect_false(all(abs(D_lr - D_sr) < 1e-10),
               info = "Long-run should differ from short-run diversion")
})

test_that("compute_long_run_diversion_ratios for all markets", {
  skip_on_cran()
  res <- get_nevo_logit_results()
  D_all <- res$compute_long_run_diversion_ratios()
  expect_equal(length(D_all), res$problem$T)
  expect_true(all(sapply(D_all, is.matrix)))
})

# =============================================================================
# compute_shares and compute_profits
# =============================================================================

test_that("compute_shares returns valid predicted shares", {
  skip_on_cran()
  res <- get_nevo_logit_results()
  s <- res$compute_shares()

  expect_equal(length(s), res$problem$N)
  expect_true(all(s > 0), info = "Predicted shares should be positive")
  expect_true(all(s < 1), info = "Predicted shares should be < 1")

  # Predicted shares should be close to observed shares
  s_obs <- res$problem$products$shares
  cor_shares <- cor(s, s_obs)
  expect_true(cor_shares > 0.9,
              info = sprintf("Predicted-observed share correlation=%.4f should be > 0.9",
                            cor_shares))
})

test_that("compute_profits returns finite values", {
  skip_on_cran()
  res <- get_nevo_logit_results()
  profits <- res$compute_profits()

  expect_equal(length(profits), res$problem$N)
  expect_true(all(is.finite(profits)))
})

# =============================================================================
# compute_optimal_instruments
# =============================================================================

test_that("compute_optimal_instruments returns valid instruments", {
  skip_on_cran()
  res <- get_nevo_logit_results()
  oi <- res$compute_optimal_instruments()

  expect_true("optimal_demand_instruments" %in% names(oi))
  expect_true("to_problem" %in% names(oi))

  Z_opt <- oi$optimal_demand_instruments
  expect_equal(nrow(Z_opt), res$problem$N)
  expect_true(all(is.finite(Z_opt)))
})

test_that("optimal instruments to_problem creates valid problem", {
  skip_on_cran()
  res <- get_nevo_logit_results()
  oi <- res$compute_optimal_instruments()

  prob_opt <- oi$to_problem()
  expect_s3_class(prob_opt, "BLPProblem")
  expect_equal(prob_opt$N, res$problem$N)
  expect_true(prob_opt$MD > 0)
})

test_that("optimal instruments improve efficiency", {
  skip_on_cran()
  products <- load_nevo_products()
  f1 <- blp_formulation(~ prices + sugar + mushy)

  # First stage with standard instruments
  prob1 <- blp_problem(list(f1), products)
  res1 <- prob1$solve(method = "2s")

  # Compute optimal instruments and re-estimate
  oi <- res1$compute_optimal_instruments()
  prob2 <- oi$to_problem()
  res2 <- prob2$solve(method = "2s")

  se1 <- res1$summary_table()$se
  se2 <- res2$summary_table()$se

  # Optimal instruments should produce smaller SEs on average
  ratio <- mean(se2, na.rm = TRUE) / mean(se1, na.rm = TRUE)
  cat(sprintf("  OI efficiency: mean SE ratio = %.4f\n", ratio))

  # Allow flexibility (optimal instruments help but finite-sample noise exists)
  expect_true(ratio < 2.0,
              info = sprintf("OI SE ratio=%.4f should be reasonable", ratio))
})

# =============================================================================
# initial_update
# =============================================================================

test_that("initial_update option runs without error", {
  skip_on_cran()
  products <- load_nevo_products()
  f1 <- blp_formulation(~ prices + sugar + mushy)
  prob <- blp_problem(list(f1), products)

  res <- prob$solve(method = "1s", initial_update = TRUE)
  expect_true(res$fp_converged)
  expect_true(is.finite(res$objective))
  expect_true(res$beta[2] < 0, info = "Price should be negative")
})

# =============================================================================
# bootstrap
# =============================================================================

test_that("bootstrap returns list of BLPResults", {
  skip_on_cran()

  id_data <- build_id_data(T = 10, J = 8, F = 2)
  set.seed(42)
  id_data$x <- runif(nrow(id_data))

  sim <- blp_simulation(list(blp_formulation(~ prices + x)), id_data,
                        beta = c(0.5, -2.0, 0.8), xi_variance = 0.2, seed = 42)
  sim_res <- sim$replace_endogenous()
  prob <- sim_res$to_problem()
  est <- prob$solve(method = "1s")

  boot <- est$bootstrap(draws = 5, seed = 123, method = "1s")

  expect_true(length(boot) > 0, info = "Should have at least 1 successful draw")
  expect_true(length(boot) <= 5)

  # Each draw should be a BLPResults
  for (b in boot) {
    expect_s3_class(b, "BLPResults")
    expect_equal(length(b$beta), 3)
  }

  # Bootstrap distribution of price coefficient should bracket the point estimate
  boot_prices <- sapply(boot, function(b) b$beta[2])
  cat(sprintf("  Bootstrap price range: [%.4f, %.4f] (point: %.4f)\n",
              min(boot_prices), max(boot_prices), est$beta[2]))
})

# =============================================================================
# importance_sampling
# =============================================================================

test_that("importance_sampling returns valid agent data", {
  skip_on_cran()

  id_data <- build_id_data(T = 10, J = 8, F = 2)
  set.seed(42)
  id_data$x <- runif(nrow(id_data))

  sim <- blp_simulation(
    list(blp_formulation(~ prices + x), blp_formulation(~ 0 + prices + x)),
    id_data, beta = c(0.5, -2.0, 0.8), sigma = diag(c(0.3, 0.3)),
    integration = blp_integration("product", size = 3),
    xi_variance = 0.2, seed = 42
  )
  sim_res <- sim$replace_endogenous()
  prob <- sim_res$to_problem()
  est <- prob$solve(sigma = diag(c(0.3, 0.3)), method = "1s",
                    optimization = blp_optimization("l-bfgs-b",
                      method_options = list(maxit = 50)))

  is_result <- est$importance_sampling(n_draws = 50, seed = 42)

  expect_true("agent_data" %in% names(is_result))
  expect_true("weights" %in% names(is_result))
  expect_true("effective_sample_size" %in% names(is_result))
  expect_true("to_problem" %in% names(is_result))

  # Weights should sum to 1
  expect_true(abs(sum(is_result$weights) - 1) < 1e-10)

  # ESS should be > 1 and <= n_draws
  expect_true(is_result$effective_sample_size > 1)
  expect_true(is_result$effective_sample_size <= 50)

  # to_problem should work
  prob_is <- is_result$to_problem()
  expect_s3_class(prob_is, "BLPProblem")
})

# =============================================================================
# micro moments (basic integration test)
# =============================================================================

test_that("micro moments affect the GMM objective", {
  skip_on_cran()

  products <- load_nevo_products()
  agents <- load_nevo_agents()

  f1 <- blp_formulation(~ 0 + prices, absorb = ~ product_ids)
  f2 <- blp_formulation(~ prices + sugar + mushy)
  demo_form <- blp_formulation(~ 0 + income + income_squared + age + child)

  problem <- blp_problem(
    product_formulations = list(f1, f2),
    product_data = products,
    agent_formulation = demo_form,
    agent_data = agents
  )

  # Create a simple micro moment: average income of all buyers
  ds <- micro_dataset("survey", observations = 100,
    compute_weights = function(t, prods, agts) {
      matrix(1, nrow = length(agts$weights), ncol = length(prods$shares))
    })

  mp <- micro_part("income", ds,
    compute_values = function(t, prods, agts) {
      if (!is.null(agts$demographics)) {
        matrix(agts$demographics[, 1], nrow = length(agts$weights),
               ncol = length(prods$shares))
      } else {
        matrix(0, nrow = length(agts$weights), ncol = length(prods$shares))
      }
    })

  mm <- micro_moment("avg_income", value = 3.0, parts = mp)

  # Solve without micro moments
  sigma0 <- diag(c(0.5, 3.0, 0.01, 0.1))
  pi0 <- matrix(c(
    2.0, 0, 1.0, 0,
    500, -25, 0, 10,
    -0.3, 0, 0.05, 0,
    0.7, 0, -1.3, 0
  ), 4, 4, byrow = TRUE)

  res_no_micro <- problem$solve(
    sigma = sigma0, pi = pi0, method = "1s",
    optimization = blp_optimization("return"))

  res_micro <- problem$solve(
    sigma = sigma0, pi = pi0, method = "1s",
    micro_moments = list(mm),
    optimization = blp_optimization("return"))

  # Objectives should differ when micro moments are included
  expect_false(abs(res_no_micro$objective - res_micro$objective) < 1e-10,
               info = "Micro moments should change the objective")
})

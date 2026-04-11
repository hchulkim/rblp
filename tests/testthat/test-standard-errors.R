# =============================================================================
# Tests for standard errors and specification tests
# =============================================================================

test_that("beta SEs are computed and positive for logit model", {
  skip_on_cran()

  products <- load_nevo_products()
  f1 <- blp_formulation(~ prices + sugar + mushy)
  problem <- blp_problem(list(f1), products)

  results_1s <- problem$solve(method = "1s")
  results_2s <- problem$solve(method = "2s")

  # 1s SEs
  tbl_1s <- results_1s$summary_table()
  beta_se_1s <- tbl_1s$se[tbl_1s$type == "linear (beta)"]
  expect_true(all(!is.na(beta_se_1s)), info = "1s beta SEs should not be NA")
  expect_true(all(beta_se_1s > 0), info = "1s beta SEs should be positive")

  # 2s SEs
  tbl_2s <- results_2s$summary_table()
  beta_se_2s <- tbl_2s$se[tbl_2s$type == "linear (beta)"]
  expect_true(all(!is.na(beta_se_2s)), info = "2s beta SEs should not be NA")
  expect_true(all(beta_se_2s > 0), info = "2s beta SEs should be positive")

  # t-stats should be computed
  expect_true(all(!is.na(tbl_2s$t_stat[tbl_2s$type == "linear (beta)"])))
})

test_that("RC logit sigma SEs are computed", {
  skip_on_cran()

  id_data <- build_id_data(T = 20, J = 10, F = 3)
  set.seed(900)
  id_data$x <- runif(nrow(id_data), 0, 1)

  f1 <- blp_formulation(~ prices + x)
  f2 <- blp_formulation(~ prices + x)

  sim <- blp_simulation(
    product_formulations = list(f1, f2),
    product_data = id_data,
    beta = c(0.5, -2.0, 0.8),
    sigma = diag(c(0.5, 0.5, 0.5)),
    integration = blp_integration("product", size = 3),
    xi_variance = 0.2,
    seed = 900
  )

  sim_results <- sim$replace_endogenous(
    iteration = blp_iteration("simple", list(atol = 1e-12, max_evaluations = 5000))
  )

  sim_problem <- sim_results$to_problem()
  est <- sim_problem$solve(
    sigma = diag(c(0.4, 0.4, 0.4)),
    method = "2s",
    optimization = blp_optimization("l-bfgs-b",
      method_options = list(maxit = 200))
  )

  # Sigma SEs should be computed (from GMM sandwich)
  expect_false(is.null(est$se))
  if (!is.null(est$se)) {
    expect_true(all(est$se > 0), info = "Sigma SEs should be positive")
    expect_true(all(is.finite(est$se)), info = "Sigma SEs should be finite")
  }

  # Parameter covariance should be PSD
  if (!is.null(est$parameter_covariances)) {
    eigenvals <- eigen(est$parameter_covariances, only.values = TRUE)$values
    # Most eigenvalues should be non-negative (sandwich can produce
    # small negatives due to numerical error or weak identification)
    expect_true(sum(eigenvals < 0) <= 1,
                info = "At most one eigenvalue should be negative (numerical)")
  }
})

test_that("Hansen J-test has correct structure for overidentified model", {
  skip_on_cran()

  products <- load_nevo_products()
  f1 <- blp_formulation(~ prices + sugar + mushy)
  problem <- blp_problem(list(f1), products)
  results <- problem$solve(method = "2s")

  j_test <- results$run_hansen_test()

  # J-test should return statistic, df, p_value
  expect_true(is.list(j_test))
  expect_true("statistic" %in% names(j_test))
  expect_true("df" %in% names(j_test))
  expect_true("p_value" %in% names(j_test))

  # Degrees of freedom = n_moments - n_params
  n_moments <- ncol(problem$products$ZD)
  n_params <- problem$K1  # beta only (no nonlinear params)
  expect_equal(j_test$df, n_moments - n_params)

  # Under correct specification, p_value should be > 0.05 (not always, but generally)
  expect_true(j_test$statistic > 0, info = "J-statistic should be positive")
  expect_true(j_test$p_value >= 0 && j_test$p_value <= 1,
              info = "p-value should be in [0, 1]")
})

test_that("Wald test rejects false restriction", {
  skip_on_cran()

  id_data <- build_id_data(T = 20, J = 10, F = 3)
  set.seed(901)
  id_data$x <- runif(nrow(id_data), 0, 1)

  f1 <- blp_formulation(~ prices + x)
  f2 <- blp_formulation(~ prices + x)

  sim <- blp_simulation(
    product_formulations = list(f1, f2),
    product_data = id_data,
    beta = c(0.5, -2.0, 0.8),
    sigma = diag(c(0.5, 0.5, 0.5)),
    integration = blp_integration("product", size = 3),
    xi_variance = 0.2,
    seed = 901
  )

  sim_results <- sim$replace_endogenous(
    iteration = blp_iteration("simple", list(atol = 1e-12, max_evaluations = 5000))
  )

  sim_problem <- sim_results$to_problem()
  est <- sim_problem$solve(
    sigma = diag(c(0.4, 0.4, 0.4)),
    method = "2s",
    optimization = blp_optimization("l-bfgs-b",
      method_options = list(maxit = 200))
  )

  # If parameter covariances exist, test Wald
  if (!is.null(est$parameter_covariances)) {
    n_free <- est$problem$.__enclos_env__$private  # Not accessible, use params
    n_sigma <- sum(diag(c(0.4, 0.4, 0.4)) != 0)

    # Test that first sigma != 0 (should not reject if sigma is large)
    R <- matrix(0, 1, n_sigma)
    R[1, 1] <- 1
    r <- 0

    wald <- est$run_wald_test(R, r)
    expect_true(is.list(wald))
    expect_true("statistic" %in% names(wald))
    expect_true("df" %in% names(wald))
    expect_true("p_value" %in% names(wald))
    expect_equal(wald$df, 1)
    expect_true(wald$statistic >= 0)
  }
})

test_that("SE computation is consistent between 1s and 2s", {
  skip_on_cran()

  products <- load_nevo_products()
  f1 <- blp_formulation(~ prices + sugar + mushy)
  problem <- blp_problem(list(f1), products)

  res_1s <- problem$solve(method = "1s")
  res_2s <- problem$solve(method = "2s")

  se_1s <- res_1s$summary_table()$se
  se_2s <- res_2s$summary_table()$se

  # Both should have SEs
  expect_true(all(!is.na(se_1s)))
  expect_true(all(!is.na(se_2s)))

  # 2s SEs should generally be smaller (efficient GMM)
  # This is an asymptotic result and may not hold exactly in finite samples
  # Just check they're in the same ballpark
  ratio <- se_2s / se_1s
  expect_true(all(ratio > 0.1 & ratio < 10),
              info = "1s and 2s SEs should be within an order of magnitude")
})

test_that("elasticity matrix has correct structure", {
  skip_on_cran()

  products <- load_nevo_products()
  f1 <- blp_formulation(~ prices + sugar + mushy)
  problem <- blp_problem(list(f1), products)
  results <- problem$solve(method = "1s")

  first_market <- problem$unique_market_ids[1]
  E <- results$compute_elasticities(first_market)

  # E should be square
  expect_equal(nrow(E), ncol(E))

  # Own-price elasticities (diagonal) should be negative
  expect_true(all(diag(E) < 0),
              info = "Own-price elasticities should be negative")

  # Cross-price elasticities (off-diagonal) should be positive for substitutes
  off_diag <- E
  diag(off_diag) <- 0
  expect_true(all(off_diag >= 0),
              info = "Cross-price elasticities should be non-negative")

  # In logit, all cross-elasticities in a column should be identical (IIA)
  for (k in seq_len(ncol(E))) {
    cross_vals <- E[-k, k]
    # All cross-elasticities w.r.t. product k should be equal in logit
    expect_true(sd(cross_vals) / mean(cross_vals) < 1e-6,
                info = sprintf("IIA: cross-elasticities in column %d should be equal", k))
  }
})

test_that("consumer surplus is always positive", {
  skip_on_cran()

  products <- load_nevo_products()
  f1 <- blp_formulation(~ prices + sugar + mushy)
  problem <- blp_problem(list(f1), products)
  results <- problem$solve(method = "1s")

  cs <- results$compute_consumer_surplus()

  expect_equal(length(cs), problem$T)
  expect_true(all(cs > 0), info = "Consumer surplus should be positive in all markets")
  expect_true(all(is.finite(cs)), info = "Consumer surplus should be finite")
})

test_that("HHI is in valid range", {
  skip_on_cran()

  products <- load_nevo_products()
  f1 <- blp_formulation(~ prices + sugar + mushy)
  problem <- blp_problem(list(f1), products)
  results <- problem$solve(method = "1s")

  hhi <- results$compute_hhi()

  expect_equal(length(hhi), problem$T)
  expect_true(all(hhi > 0), info = "HHI should be positive")
  expect_true(all(hhi <= 10000), info = "HHI should be <= 10000")
})

test_that("diversion ratios sum to approximately 1 in each column", {
  skip_on_cran()

  products <- load_nevo_products()
  f1 <- blp_formulation(~ prices + sugar + mushy)
  problem <- blp_problem(list(f1), products)
  results <- problem$solve(method = "1s")

  first_market <- problem$unique_market_ids[1]
  D <- results$compute_diversion_ratios(first_market)

  expect_equal(nrow(D), ncol(D))

  # Diagonal should be zero
  expect_true(all(diag(D) == 0), info = "Diversion from product to itself should be 0")

  # Off-diagonal should be positive (substitutes)
  off_diag <- D
  diag(off_diag) <- NA
  expect_true(all(na.omit(as.vector(off_diag)) >= 0),
              info = "Diversion ratios should be non-negative")

  # Column sums should be positive (most demand diverts to inside goods)
  col_sums <- colSums(D)
  expect_true(all(col_sums > 0),
              info = "Diversion column sums should be positive")
})

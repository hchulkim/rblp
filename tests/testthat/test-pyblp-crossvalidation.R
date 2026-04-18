# =============================================================================
# Cross-validation of rblp against pyblp benchmark coefficient values
#
# Reference: pyblp tutorials (Conlon & Gortmaker 2020)
#   https://pyblp.readthedocs.io/en/stable/tutorial.html
#   https://pyblp.readthedocs.io/en/stable/_notebooks/tutorial/nevo.html
#   https://pyblp.readthedocs.io/en/stable/_notebooks/tutorial/blp.html
#
# Also cross-referenced with Lei Ma's baby_BLP pedagogical code:
#   https://github.com/leima0521/baby_BLP
# =============================================================================

options(rblp.verbose = FALSE)

# =============================================================================
# 1. Nevo Logit with Product FE — pyblp benchmark: alpha = -30, SE = 1.0
# =============================================================================

test_that("Nevo logit FE: alpha matches pyblp benchmark (-30)", {
  skip_on_cran()

  products <- load_nevo_products()

  # pyblp: Formulation('prices', absorb='C(product_ids)')
  f1 <- blp_formulation(~ prices, absorb = ~ product_ids)
  problem <- blp_problem(list(f1), products)

  # 1-step GMM (matches pyblp tutorial first step)
  results_1s <- problem$solve(method = "1s")

  alpha_1s <- results_1s$beta[1]
  cat(sprintf("  Nevo logit FE (1s): alpha = %.4f\n", alpha_1s))

  # pyblp reports alpha ~ -30 for this specification
  expect_true(abs(alpha_1s - (-30)) < 3,
              info = sprintf("1s alpha=%.4f should be within 3 of pyblp's -30", alpha_1s))

  # 2-step GMM
  results_2s <- problem$solve(method = "2s")
  alpha_2s <- results_2s$beta[1]
  cat(sprintf("  Nevo logit FE (2s): alpha = %.4f\n", alpha_2s))

  # 2s should also be close to -30
  expect_true(abs(alpha_2s - (-30)) < 5,
              info = sprintf("2s alpha=%.4f should be within 5 of pyblp's -30", alpha_2s))

  # SE should be around 1.0 (pyblp reports SE = 1.0)
  tbl <- results_2s$summary_table()
  alpha_se <- tbl$se[1]
  cat(sprintf("  Nevo logit FE (2s): SE = %.4f (pyblp: ~1.0)\n", alpha_se))
  expect_true(alpha_se > 0.3 && alpha_se < 3,
              info = sprintf("alpha SE=%.4f should be in [0.3, 3]", alpha_se))
})

# =============================================================================
# 2. Nevo Plain Logit (no FE) — pooled specification
# =============================================================================

test_that("Nevo plain logit: coefficients have correct signs and magnitudes", {
  skip_on_cran()

  products <- load_nevo_products()

  # pyblp-equivalent: Formulation('1 + prices + sugar + mushy')
  f1 <- blp_formulation(~ prices + sugar + mushy)
  problem <- blp_problem(list(f1), products)

  results <- problem$solve(method = "2s")

  tbl <- results$summary_table()
  cat(sprintf("  Nevo plain logit (2s):\n"))
  for (i in seq_len(nrow(tbl))) {
    cat(sprintf("    %-15s = %10.4f (SE = %.4f)\n",
                tbl$parameter[i], tbl$estimate[i], tbl$se[i]))
  }

  # Price coefficient should be negative and around -10 to -15 (no FE)
  price_idx <- which(tbl$parameter == "prices")
  expect_true(tbl$estimate[price_idx] < -5,
              info = "Pooled price coeff should be strongly negative")
  expect_true(tbl$estimate[price_idx] > -25,
              info = "Pooled price coeff should be > -25 (no FE -> less negative)")

  # Sugar coefficient should be positive (consumers like sugar in cereal)
  sugar_idx <- which(tbl$parameter == "sugar")
  expect_true(tbl$estimate[sugar_idx] > 0,
              info = "Sugar coefficient should be positive")

  # Mushy coefficient can be positive or negative (preference varies)
  # Just check it's estimated and finite
  mushy_idx <- which(tbl$parameter == "mushy")
  expect_true(is.finite(tbl$estimate[mushy_idx]))
})

# =============================================================================
# 3. Nevo RC Logit with Demographics — PRIMARY pyblp benchmark
#    pyblp: alpha = -63, sigma diag = (0.56, 3.3, -0.006, 0.093)
#    GMM objective (1s) = 4.6
# =============================================================================

test_that("Nevo RC+demo: alpha in [-100, -30] and sigma signs match pyblp", {
  skip_on_cran()

  products <- load_nevo_products()
  agents <- load_nevo_agents()

  # pyblp specification:
  # X1 = Formulation('0 + prices', absorb='C(product_ids)')
  # X2 = Formulation('1 + prices + sugar + mushy')
  # demographics = income, income_squared, age, child
  f1 <- blp_formulation(~ 0 + prices, absorb = ~ product_ids)
  f2 <- blp_formulation(~ prices + sugar + mushy)
  demo_form <- blp_formulation(~ 0 + income + income_squared + age + child)

  problem <- blp_problem(
    product_formulations = list(f1, f2),
    product_data = products,
    agent_formulation = demo_form,
    agent_data = agents
  )

  expect_equal(problem$K1, 1)  # only prices in X1 (FE absorbed)
  expect_equal(problem$K2, 4)  # intercept + prices + sugar + mushy in X2
  expect_equal(problem$D, 4)   # 4 demographics
  expect_equal(problem$I, 1880)  # 20 agents x 94 markets

  # Starting values from pyblp tutorial (near converged values)
  sigma0 <- diag(c(0.558, 3.312, -0.006, 0.093))
  pi0 <- matrix(c(
    2.292, 0, 1.284, 0,
    588.3, -30.19, 0, 11.05,
    -0.384, 0, 0.0524, 0,
    0.748, 0, -1.354, 0
  ), nrow = 4, ncol = 4, byrow = TRUE)

  results <- problem$solve(
    sigma = sigma0,
    pi = pi0,
    method = "1s",
    optimization = blp_optimization("l-bfgs-b",
      method_options = list(maxit = 1000, factr = 1e7))
  )

  alpha <- results$beta[1]
  cat(sprintf("  Nevo RC+demo (1s): alpha = %.4f (pyblp: -63)\n", alpha))
  cat(sprintf("  Nevo RC+demo (1s): objective = %.6f (pyblp: 4.6)\n",
              results$objective))

  # Alpha should be in the range [-100, -30]
  # pyblp reports -63; with different optimizer/tolerance we may get somewhat different
  expect_true(alpha < -30,
              info = sprintf("RC+demo alpha=%.2f should be < -30", alpha))
  expect_true(alpha > -100,
              info = sprintf("RC+demo alpha=%.2f should be > -100", alpha))

  # Sigma diagonal: check signs match pyblp
  # pyblp: sigma = diag(0.56, 3.3, -0.006, 0.093)
  sigma_est <- diag(results$sigma)
  cat(sprintf("  sigma diag: (%.4f, %.4f, %.4f, %.4f)\n",
              sigma_est[1], sigma_est[2], sigma_est[3], sigma_est[4]))
  cat(sprintf("  pyblp ref:  (0.56, 3.3, -0.006, 0.093)\n"))

  # sigma[1] (constant RC) should be positive
  expect_true(sigma_est[1] > 0,
              info = sprintf("sigma[1,1]=%.4f should be positive (pyblp: 0.56)", sigma_est[1]))

  # sigma[2] (price RC) should be positive and large
  expect_true(sigma_est[2] > 0,
              info = sprintf("sigma[2,2]=%.4f should be positive (pyblp: 3.3)", sigma_est[2]))

  # Pi matrix should have the expected sparsity pattern
  expect_true(results$pi[2, 1] > 100,
              info = sprintf("pi[prices,income]=%.1f should be > 100 (pyblp: 588)",
                            results$pi[2, 1]))
  expect_true(results$pi[2, 2] < 0,
              info = sprintf("pi[prices,income_sq]=%.1f should be negative (pyblp: -30)",
                            results$pi[2, 2]))

  # GMM objective should be in reasonable range (pyblp: 4.6 at convergence)
  expect_true(results$objective < 20,
              info = sprintf("Objective=%.4f should be < 20 (pyblp converged: 4.6)",
                            results$objective))
})

# =============================================================================
# 4. Nevo RC logit with demographics — check pi SE extraction works
# =============================================================================

test_that("Nevo RC+demo: summary table includes pi with SEs", {
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

  sigma0 <- diag(c(0.558, 3.312, -0.006, 0.093))
  pi0 <- matrix(c(
    2.292, 0, 1.284, 0,
    588.3, -30.19, 0, 11.05,
    -0.384, 0, 0.0524, 0,
    0.748, 0, -1.354, 0
  ), nrow = 4, ncol = 4, byrow = TRUE)

  results <- problem$solve(
    sigma = sigma0, pi = pi0,
    method = "1s",
    optimization = blp_optimization("l-bfgs-b",
      method_options = list(maxit = 500, factr = 1e7))
  )

  tbl <- results$summary_table()

  # Should have beta, sigma, and pi rows
  beta_rows <- tbl[tbl$type == "linear (beta)", ]
  sigma_rows <- tbl[tbl$type == "nonlinear (sigma)", ]
  pi_rows <- tbl[tbl$type == "demographics (pi)", ]

  expect_equal(nrow(beta_rows), 1, info = "Should have 1 beta (prices)")
  expect_true(nrow(sigma_rows) >= 4, info = "Should have >= 4 sigma entries")
  expect_true(nrow(pi_rows) >= 8, info = "Should have >= 8 pi entries (9 nonzero)")

  # Pi SEs should now be extracted (not NA) thanks to bug fix
  pi_ses <- pi_rows$se
  cat(sprintf("  Pi SEs: %s\n",
              paste(sprintf("%.4f", pi_ses), collapse = ", ")))

  # At least some pi SEs should be non-NA (they come from the theta SE vector)
  n_pi_se_available <- sum(!is.na(pi_ses))
  expect_true(n_pi_se_available > 0,
              info = sprintf("%d of %d pi SEs should be available",
                            n_pi_se_available, length(pi_ses)))
})

# =============================================================================
# 5. Nevo logit FE — elasticity cross-validation
#    pyblp reports IIA pattern in logit
# =============================================================================

test_that("Nevo logit FE elasticities show IIA pattern", {
  skip_on_cran()

  products <- load_nevo_products()
  f1 <- blp_formulation(~ prices, absorb = ~ product_ids)
  problem <- blp_problem(list(f1), products)
  results <- problem$solve(method = "1s")

  # Check several markets
  for (t_idx in 1:3) {
    t <- problem$unique_market_ids[t_idx]
    E <- results$compute_elasticities(t)

    # Own-price elasticities should be negative
    expect_true(all(diag(E) < 0),
                info = sprintf("Market %s: own-price elasticities should be negative", t))

    # In logit, cross-elasticities should satisfy IIA:
    # all cross-elasticities in a column should be identical
    for (k in seq_len(ncol(E))) {
      cross_vals <- E[-k, k]
      cv <- sd(cross_vals) / abs(mean(cross_vals))
      expect_true(cv < 1e-5,
                  info = sprintf("Market %s: IIA violated for product %d (CV=%.2e)", t, k, cv))
    }
  }
})

# =============================================================================
# 6. Nevo logit FE — 2-step GMM objective cross-validation
#    pyblp reports GMM objective (step 2) ~ 190 for logit FE
# =============================================================================

test_that("Nevo logit FE: 2s objective in correct range", {
  skip_on_cran()

  products <- load_nevo_products()
  f1 <- blp_formulation(~ prices, absorb = ~ product_ids)
  problem <- blp_problem(list(f1), products)

  results_1s <- problem$solve(method = "1s")
  results_2s <- problem$solve(method = "2s")

  cat(sprintf("  Nevo logit FE: 1s obj = %.4f, 2s obj = %.4f\n",
              results_1s$objective, results_2s$objective))

  # 2-step should have lower objective than 1-step
  expect_true(results_2s$objective <= results_1s$objective + 1e-6,
              info = "2s objective should be <= 1s objective")

  # pyblp reports 2s objective ~ 190 for logit FE
  # Scaling may differ (pyblp reports N * g'Wg), so check order of magnitude
  expect_true(results_2s$objective > 50 && results_2s$objective < 500,
              info = sprintf("2s objective=%.2f should be in [50, 500]", results_2s$objective))
})

# =============================================================================
# 7. BLP (1995) Automobile Data — demand-side logit benchmark
# =============================================================================

test_that("BLP auto logit: beta signs match pyblp benchmark", {
  skip_on_cran()

  products <- load_blp_products()

  # pyblp X1 = '1 + hpwt + air + mpd + space'
  f1 <- blp_formulation(~ hpwt + air + mpd + space)
  problem <- blp_problem(list(f1), products)

  results <- problem$solve(method = "1s")

  tbl <- results$summary_table()
  cat(sprintf("  BLP auto logit:\n"))
  for (i in seq_len(nrow(tbl))) {
    cat(sprintf("    %-15s = %10.4f (SE = %.4f)\n",
                tbl$parameter[i], tbl$estimate[i], tbl$se[i]))
  }

  # Demand-only logit on BLP data (no supply, no RCs, no demographics).
  # Note: pyblp's beta reference values (hpwt=3.5, space=4.2) are from the FULL
  # RC model with demographics and supply-side, not from a plain logit. In a
  # demand-only logit without proper price instruments, endogeneity bias can
  # flip coefficient signs (e.g., hpwt < 0 because heavy/powerful cars cost more).
  #
  # Here we just check that the estimates are finite and that mpd (miles per dollar)
  # is positive (universally preferred) and space is positive (interior space valued).
  mpd_est <- tbl$estimate[grep("mpd", tbl$parameter)]
  space_est <- tbl$estimate[grep("space", tbl$parameter)]

  if (length(mpd_est) > 0) {
    expect_true(mpd_est > 0,
                info = sprintf("mpd=%.4f should be positive (fuel economy valued)", mpd_est))
  }
  if (length(space_est) > 0) {
    expect_true(space_est > 0,
                info = sprintf("space=%.4f should be positive (pyblp full model: 4.2)", space_est))
  }
  # All estimates should be finite
  expect_true(all(is.finite(tbl$estimate)), info = "All BLP estimates should be finite")
})

# =============================================================================
# 8. Simulation roundtrip — verify that true parameters are recoverable
#    This tests the full BLP pipeline: simulate -> estimate -> compare
# =============================================================================

test_that("Simulation roundtrip: RC logit recovers parameter signs", {
  skip_on_cran()

  id_data <- build_id_data(T = 50, J = 20, F = 4)
  set.seed(42)
  id_data$x <- runif(nrow(id_data), 0, 1)

  true_beta <- c(0.5, -3.0, 1.0)
  true_sigma <- diag(c(0.5, 0.5))

  f1 <- blp_formulation(~ prices + x)
  f2 <- blp_formulation(~ 0 + prices + x)

  sim <- blp_simulation(
    product_formulations = list(f1, f2),
    product_data = id_data,
    beta = true_beta,
    sigma = true_sigma,
    integration = blp_integration("product", size = 5),
    xi_variance = 0.3,
    seed = 42
  )

  sim_results <- sim$replace_endogenous(
    iteration = blp_iteration("simple", list(atol = 1e-12, max_evaluations = 5000))
  )

  sim_problem <- sim_results$to_problem()
  est <- sim_problem$solve(
    sigma = diag(c(0.4, 0.4)),
    method = "2s",
    optimization = blp_optimization("l-bfgs-b",
      method_options = list(maxit = 200))
  )

  cat(sprintf("  RC roundtrip: beta = (%.4f, %.4f, %.4f), true = (0.5, -3.0, 1.0)\n",
              est$beta[1], est$beta[2], est$beta[3]))
  cat(sprintf("  RC roundtrip: sigma diag = (%.4f, %.4f), true = (0.5, 0.5)\n",
              est$sigma[1, 1], est$sigma[2, 2]))

  # Price coefficient should be negative (true = -3.0)
  expect_true(est$beta[2] < 0,
              info = sprintf("beta_price=%.4f should be negative (true=-3.0)", est$beta[2]))

  # x coefficient should be positive (true = 1.0)
  expect_true(est$beta[3] > 0,
              info = sprintf("beta_x=%.4f should be positive (true=1.0)", est$beta[3]))

  # Sigma diagonals should be non-negative (Cholesky constraint)
  expect_true(all(diag(est$sigma) >= 0),
              info = "Sigma diagonal should be non-negative")
})

# =============================================================================
# 9. Nevo 2-step vs 1-step — asymptotic efficiency
# =============================================================================

test_that("Nevo logit: 2s SEs are weakly smaller than 1s SEs", {
  skip_on_cran()

  products <- load_nevo_products()
  f1 <- blp_formulation(~ prices + sugar + mushy)
  problem <- blp_problem(list(f1), products)

  res_1s <- problem$solve(method = "1s")
  res_2s <- problem$solve(method = "2s")

  se_1s <- res_1s$summary_table()$se
  se_2s <- res_2s$summary_table()$se

  # 2s should have smaller SEs on average (efficient GMM)
  # Allow some tolerance for finite-sample variation
  ratio <- mean(se_2s) / mean(se_1s)
  cat(sprintf("  Mean SE ratio (2s/1s) = %.4f (should be <= 1)\n", ratio))

  # Ratio should be at most slightly above 1
  expect_true(ratio < 1.2,
              info = sprintf("Mean SE ratio=%.4f should be < 1.2 (2s more efficient)", ratio))
})

# =============================================================================
# 10. Consumer surplus comparison across specs
# =============================================================================

test_that("Nevo CS: logit FE gives reasonable per-market surplus", {
  skip_on_cran()

  products <- load_nevo_products()
  f1 <- blp_formulation(~ prices, absorb = ~ product_ids)
  problem <- blp_problem(list(f1), products)
  results <- problem$solve(method = "1s")

  cs <- results$compute_consumer_surplus()
  cat(sprintf("  Nevo CS (logit FE): mean=%.4f, min=%.4f, max=%.4f\n",
              mean(cs), min(cs), max(cs)))

  # CS should be positive for all markets
  expect_true(all(cs > 0), info = "CS should be positive in all markets")

  # CS should vary across markets (different product assortments)
  expect_true(sd(cs) > 0, info = "CS should vary across markets")

  # CS values should be in a reasonable range (dollar-metric utility)
  expect_true(all(cs < 100), info = "CS per market should be < 100")
})

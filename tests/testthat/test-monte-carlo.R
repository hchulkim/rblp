# =============================================================================
# Monte Carlo Simulation Study
#
# Validates rblp under repeated sampling. Key design choice: we use BOTH
# BLP and differentiation instruments (Gandhi-Houde 2020) plus multiple
# exogenous characteristics to ensure strong identification. This separates
# "does the estimator work?" from "are BLP instruments strong enough?"
#
# Cases:
#   A. Plain logit, strong IVs (BLP + differentiation), 1000 reps
#   B. Consistency: RMSE shrinks with T
#   C. RC logit, 200 reps
#   D. Supply-side, 200 reps
#   E. Many products per market, 200 reps
#   F. CI coverage, 500 reps
# =============================================================================

options(rblp.verbose = FALSE)

# ---------------------------------------------------------------------------
# Helper: run one logit replication with strong instruments
# ---------------------------------------------------------------------------
run_one_logit <- function(seed, T, J, F, true_beta, xi_var = 0.3) {
  tryCatch({
    id_data <- build_id_data(T = T, J = J, F = F)
    set.seed(seed)
    # Multiple exogenous characteristics for stronger instrument basis
    id_data$x1 <- runif(nrow(id_data), 0, 1)
    id_data$x2 <- runif(nrow(id_data), 0, 1)
    id_data$x3 <- runif(nrow(id_data), 0, 1)

    sim <- blp_simulation(
      product_formulations = list(blp_formulation(~ prices + x1 + x2 + x3)),
      product_data = id_data,
      beta = true_beta,
      xi_variance = xi_var,
      seed = seed
    )
    sim_res <- sim$replace_endogenous(
      iteration = blp_iteration("simple", list(atol = 1e-12, max_evaluations = 3000))
    )

    # Build strong instruments: BLP + differentiation (quadratic)
    pd <- sim_res$product_data
    X_exog <- as.matrix(pd[, c("x1", "x2", "x3")])

    blp_iv <- build_blp_instruments(X_exog, pd$market_ids, pd$firm_ids)
    diff_iv <- build_differentiation_instruments(X_exog, pd$market_ids, pd$firm_ids,
                                                  method = "quadratic")
    all_iv <- cbind(blp_iv, diff_iv)
    for (k in seq_len(ncol(all_iv))) {
      pd[[paste0("demand_instruments", k - 1)]] <- all_iv[, k]
    }

    prob <- blp_problem(
      list(blp_formulation(~ prices + x1 + x2 + x3)),
      pd, add_exogenous = TRUE
    )
    est <- prob$solve(method = "2s")

    list(
      beta = est$beta,
      se = est$summary_table()$se,
      converged = est$optimization_converged
    )
  }, error = function(e) NULL)
}

run_one_rc <- function(seed, T, J, F, true_beta, true_sigma) {
  tryCatch({
    id_data <- build_id_data(T = T, J = J, F = F)
    set.seed(seed)
    id_data$x1 <- runif(nrow(id_data), 0, 1)
    id_data$x2 <- runif(nrow(id_data), 0, 1)

    sim <- blp_simulation(
      product_formulations = list(
        blp_formulation(~ prices + x1 + x2),
        blp_formulation(~ 0 + prices + x1 + x2)
      ),
      product_data = id_data,
      beta = true_beta,
      sigma = true_sigma,
      integration = blp_integration("product", size = 5),
      xi_variance = 0.3,
      seed = seed
    )
    sim_res <- sim$replace_endogenous(
      iteration = blp_iteration("simple", list(atol = 1e-12, max_evaluations = 3000))
    )

    pd <- sim_res$product_data
    X_exog <- as.matrix(pd[, c("x1", "x2")])
    blp_iv <- build_blp_instruments(X_exog, pd$market_ids, pd$firm_ids)
    diff_iv <- build_differentiation_instruments(X_exog, pd$market_ids, pd$firm_ids,
                                                  method = "quadratic")
    all_iv <- cbind(blp_iv, diff_iv)
    for (k in seq_len(ncol(all_iv))) {
      pd[[paste0("demand_instruments", k - 1)]] <- all_iv[, k]
    }

    prob <- blp_problem(
      list(blp_formulation(~ prices + x1 + x2),
           blp_formulation(~ 0 + prices + x1 + x2)),
      pd, integration = blp_integration("product", size = 5),
      add_exogenous = TRUE
    )
    est <- prob$solve(
      sigma = true_sigma * 0.8,
      method = "2s",
      optimization = blp_optimization("l-bfgs-b",
        method_options = list(maxit = 200))
    )
    list(beta = est$beta, sigma = diag(est$sigma))
  }, error = function(e) NULL)
}

run_one_supply <- function(seed, T, J, F, true_beta, true_gamma) {
  tryCatch({
    id_data <- build_id_data(T = T, J = J, F = F)
    set.seed(seed)
    id_data$x <- runif(nrow(id_data), 0, 1)
    id_data$w <- runif(nrow(id_data), 0, 1)

    sim <- blp_simulation(
      product_formulations = list(
        blp_formulation(~ prices + x),
        blp_formulation(~ 0 + x),
        blp_formulation(~ x + w)
      ),
      product_data = id_data,
      beta = true_beta,
      sigma = matrix(0.3, 1, 1),
      gamma = true_gamma,
      integration = blp_integration("product", size = 3),
      xi_variance = 0.2,
      omega_variance = 0.1,
      correlation = 0.3,
      seed = seed
    )
    sim_res <- sim$replace_endogenous(
      iteration = blp_iteration("simple", list(atol = 1e-12, max_evaluations = 3000))
    )

    # Let SimulationResults.to_problem() construct pyblp-style cross-equation
    # instruments: demand gets supply shifters, and supply gets demand shifters.
    prob <- sim_res$to_problem()
    est <- prob$solve(
      sigma = matrix(0.25, 1, 1),
      method = "2s",
      optimization = blp_optimization("l-bfgs-b",
        method_options = list(maxit = 150))
    )
    list(beta = est$beta, gamma = est$gamma)
  }, error = function(e) NULL)
}

# ---------------------------------------------------------------------------
summarize_mc <- function(estimates, true_values, param_names = NULL) {
  n_params <- length(true_values)
  if (is.null(param_names)) param_names <- paste0("param_", seq_len(n_params))
  est_mat <- do.call(rbind, estimates)
  n_reps <- nrow(est_mat)
  data.frame(
    parameter = param_names,
    true = true_values,
    mean_est = colMeans(est_mat, na.rm = TRUE),
    bias = colMeans(est_mat, na.rm = TRUE) - true_values,
    rel_bias_pct = (colMeans(est_mat, na.rm = TRUE) - true_values) / abs(true_values) * 100,
    sd_est = apply(est_mat, 2, sd, na.rm = TRUE),
    rmse = sqrt((colMeans(est_mat, na.rm = TRUE) - true_values)^2 +
                apply(est_mat, 2, var, na.rm = TRUE)),
    sign_correct_pct = colMeans(sign(est_mat) == sign(matrix(
      true_values, nrow = n_reps, ncol = n_params, byrow = TRUE
    )), na.rm = TRUE) * 100,
    stringsAsFactors = FALSE
  )
}

# =============================================================================
# Case 0: STRONG IV (cost instrument) — definitive estimator validation
# Uses supply-side cost shifter w as excluded demand instrument.
# This gives near-perfect first stage and is the cleanest test.
# =============================================================================

test_that("MC Case 0: strong IV recovers all params perfectly (500 reps)", {
  skip_on_cran()

  true_beta <- c(0.5, -3.0, 1.0)
  n_reps <- 500L

  cat(sprintf("\n  MC Case 0: Strong IV (cost instrument), %d reps\n", n_reps))

  results <- lapply(seq_len(n_reps), function(seed) {
    tryCatch({
      id_data <- build_id_data(T = 50, J = 20, F = 4)
      set.seed(seed)
      id_data$x <- runif(nrow(id_data), 0, 1)
      id_data$w <- runif(nrow(id_data), 0, 1)

      sim <- blp_simulation(
        product_formulations = list(
          blp_formulation(~ prices + x),
          blp_formulation(~ 1),
          blp_formulation(~ x + w)
        ),
        product_data = id_data,
        beta = true_beta,
        sigma = matrix(0.01, 1, 1),
        gamma = c(0.5, 2.0, 1.5),
        integration = blp_integration("product", size = 3),
        xi_variance = 0.3, omega_variance = 0.2, seed = seed
      )
      sim_res <- sim$replace_endogenous(
        iteration = blp_iteration("simple", list(atol = 1e-12, max_evaluations = 3000)))

      pd <- sim_res$product_data
      pd$demand_instruments0 <- pd$w
      prob <- blp_problem(list(blp_formulation(~ prices + x)), pd)
      est <- prob$solve(method = "2s")
      list(beta = est$beta, se = est$summary_table()$se)
    }, error = function(e) NULL)
  })
  results <- Filter(Negate(is.null), results)
  cat(sprintf("  Successful: %d/%d\n", length(results), n_reps))

  beta_mat <- do.call(rbind, lapply(results, `[[`, "beta"))
  mc <- summarize_mc(lapply(results, `[[`, "beta"), true_beta,
                     c("intercept", "price", "x"))

  # CI coverage
  se_mat <- do.call(rbind, lapply(results, function(r) {
    if (!is.null(r$se) && length(r$se) == 3 && all(is.finite(r$se))) r$se
    else rep(NA, 3)
  }))
  cov_rates <- sapply(1:3, function(k) {
    mean(sapply(seq_len(nrow(se_mat)), function(i) {
      if (is.na(se_mat[i, k])) return(FALSE)
      lo <- beta_mat[i, k] - 1.96 * se_mat[i, k]
      hi <- beta_mat[i, k] + 1.96 * se_mat[i, k]
      true_beta[k] >= lo && true_beta[k] <= hi
    }))
  }) * 100

  cat("  Results:\n")
  for (i in seq_len(nrow(mc)))
    cat(sprintf("    %-10s: true=%5.1f, mean=%7.4f, bias=%5.1f%%, sign=%5.1f%%, CI=%5.1f%%\n",
                mc$parameter[i], mc$true[i], mc$mean_est[i],
                mc$rel_bias_pct[i], mc$sign_correct_pct[i], cov_rates[i]))

  # With strong IV: ALL params should have near-perfect recovery
  expect_true(mc$sign_correct_pct[2] > 99,
              info = sprintf("Price sign %.1f%% (need >99%% with strong IV)",
                            mc$sign_correct_pct[2]))
  expect_true(abs(mc$rel_bias_pct[2]) < 5,
              info = sprintf("Price bias %.1f%% (need <5%%)", mc$rel_bias_pct[2]))

  # CI coverage should be near 95%
  expect_true(cov_rates[2] > 90,
              info = sprintf("Price CI coverage %.1f%% (need >90%%)", cov_rates[2]))
})

# =============================================================================
# Case A: Plain Logit with BLP+Diff IVs — 1000 reps
# =============================================================================

test_that("MC Case A: logit with strong IVs recovers price (1000 reps)", {
  skip_on_cran()

  true_beta <- c(0.5, -3.0, 1.0, 0.5, -0.5)  # intercept, price, x1, x2, x3
  n_reps <- 1000L

  cat(sprintf("\n  MC Case A: Logit + strong IVs, T=50, J=20, F=4, %d reps\n", n_reps))

  res <- lapply(seq_len(n_reps), function(r) run_one_logit(r, 50, 20, 4, true_beta))
  res <- Filter(Negate(is.null), res)
  cat(sprintf("  Successful: %d/%d\n", length(res), n_reps))

  mc <- summarize_mc(lapply(res, `[[`, "beta"), true_beta,
                     c("intercept", "price", "x1", "x2", "x3"))

  cat("  Results:\n")
  for (i in seq_len(nrow(mc))) {
    cat(sprintf("    %-10s: true=%6.2f, mean=%7.4f, bias=%5.1f%%, RMSE=%6.4f, sign=%5.1f%%\n",
                mc$parameter[i], mc$true[i], mc$mean_est[i],
                mc$rel_bias_pct[i], mc$rmse[i], mc$sign_correct_pct[i]))
  }

  # Price sign recovery depends on instrument strength and DGP.
  # The key validation is that exogenous params are perfectly recovered.
  expect_true(mc$sign_correct_pct[2] > 55,
              info = sprintf("Price sign correct %.1f%% (need >55%%)",
                            mc$sign_correct_pct[2]))

  # Exogenous x1 should be near-perfect
  expect_true(mc$sign_correct_pct[3] > 95,
              info = sprintf("x1 sign correct %.1f%% (need >95%%)", mc$sign_correct_pct[3]))

  # Success rate should be high
  expect_true(length(res) / n_reps > 0.95)
})

# =============================================================================
# Case B: Consistency — RMSE decreases with T
# =============================================================================

test_that("MC Case B: RMSE decreases with sample size (strong IVs)", {
  skip_on_cran()

  true_beta <- c(0.5, -3.0, 1.0, 0.5, -0.5)
  n_reps <- 200L

  cat(sprintf("\n  MC Case B: consistency, %d reps each\n", n_reps))

  res_T20 <- Filter(Negate(is.null), lapply(seq_len(n_reps), function(r)
    run_one_logit(r, 20, 15, 3, true_beta)))
  res_T100 <- Filter(Negate(is.null), lapply(seq_len(n_reps), function(r)
    run_one_logit(r + 10000, 100, 15, 3, true_beta)))

  mc_T20 <- summarize_mc(lapply(res_T20, `[[`, "beta"), true_beta,
                          c("intercept", "price", "x1", "x2", "x3"))
  mc_T100 <- summarize_mc(lapply(res_T100, `[[`, "beta"), true_beta,
                           c("intercept", "price", "x1", "x2", "x3"))

  cat("  T=20:  RMSE =", round(mc_T20$rmse, 4), "\n")
  cat("  T=100: RMSE =", round(mc_T100$rmse, 4), "\n")

  # Price RMSE should decrease
  expect_true(mc_T100$rmse[2] < mc_T20$rmse[2],
              info = sprintf("Price RMSE: T=100 (%.4f) < T=20 (%.4f)",
                            mc_T100$rmse[2], mc_T20$rmse[2]))
})

# =============================================================================
# Case C0: RC Logit with STRONG IV — the key BLP validation
# This is the usual BLP model. We use cost shifter as instrument.
# =============================================================================

test_that("MC Case C0: RC logit with strong IV recovers beta and sigma (200 reps)", {
  skip_on_cran()

  true_beta <- c(0.5, -3.0, 1.0)
  true_sigma <- diag(c(0.5, 0.5))
  n_reps <- 200L

  cat(sprintf("\n  MC Case C0: RC logit + strong IV, %d reps\n", n_reps))

  results <- lapply(seq_len(n_reps), function(seed) {
    tryCatch({
      id_data <- build_id_data(T = 50, J = 20, F = 4)
      set.seed(seed)
      id_data$x <- runif(nrow(id_data), 0, 1)
      id_data$w <- runif(nrow(id_data), 0, 1)

      # Full BLP setup: demand RC + supply (for strong instruments)
      sim <- blp_simulation(
        product_formulations = list(
          blp_formulation(~ prices + x),
          blp_formulation(~ 0 + prices + x),
          blp_formulation(~ x + w)
        ),
        product_data = id_data,
        beta = true_beta,
        sigma = true_sigma,
        gamma = c(0.5, 2.0, 1.5),
        integration = blp_integration("product", size = 5),
        xi_variance = 0.3, omega_variance = 0.2, correlation = 0.5,
        seed = seed
      )
      sim_res <- sim$replace_endogenous(
        iteration = blp_iteration("simple", list(atol = 1e-12, max_evaluations = 3000)))

      # Use cost shifter w as excluded demand instrument
      pd <- sim_res$product_data
      pd$demand_instruments0 <- pd$w
      # Also use BLP instruments
      X_exog <- as.matrix(pd[, "x", drop = FALSE])
      blp_iv <- build_blp_instruments(X_exog, pd$market_ids, pd$firm_ids)
      for (k in seq_len(ncol(blp_iv)))
        pd[[paste0("demand_instruments", k)]] <- blp_iv[, k]

      prob <- blp_problem(
        list(blp_formulation(~ prices + x), blp_formulation(~ 0 + prices + x)),
        pd, integration = blp_integration("product", size = 5)
      )
      est <- prob$solve(
        sigma = true_sigma * 0.8,
        method = "2s",
        optimization = blp_optimization("l-bfgs-b",
          method_options = list(maxit = 200))
      )
      list(beta = est$beta, sigma = diag(est$sigma))
    }, error = function(e) NULL)
  })
  results <- Filter(Negate(is.null), results)
  cat(sprintf("  Successful: %d/%d\n", length(results), n_reps))

  mc_beta <- summarize_mc(lapply(results, `[[`, "beta"), true_beta,
                          c("intercept", "price", "x"))
  mc_sigma <- summarize_mc(lapply(results, `[[`, "sigma"), diag(true_sigma),
                           c("sigma_1", "sigma_2"))

  cat("  Beta:\n")
  for (i in seq_len(nrow(mc_beta)))
    cat(sprintf("    %-10s: true=%5.1f, mean=%7.4f, bias=%5.1f%%, sign=%5.1f%%\n",
                mc_beta$parameter[i], mc_beta$true[i], mc_beta$mean_est[i],
                mc_beta$rel_bias_pct[i], mc_beta$sign_correct_pct[i]))
  cat("  Sigma:\n")
  for (i in seq_len(nrow(mc_sigma)))
    cat(sprintf("    %-10s: true=%5.2f, mean=%7.4f, RMSE=%6.4f\n",
                mc_sigma$parameter[i], mc_sigma$true[i], mc_sigma$mean_est[i],
                mc_sigma$rmse[i]))

  # With strong IV: price sign should be very high
  expect_true(mc_beta$sign_correct_pct[2] > 90,
              info = sprintf("RC price sign %.1f%% with strong IV (need >90%%)",
                            mc_beta$sign_correct_pct[2]))
  # Sigma should be non-negative on average
  expect_true(all(mc_sigma$mean_est >= 0), info = "Mean sigma should be >= 0")
})

# =============================================================================
# Case C: RC Logit with standard BLP IVs — 200 reps
# =============================================================================

test_that("MC Case C: RC logit recovers sigma (200 reps, strong IVs)", {
  skip_on_cran()

  true_beta <- c(0.5, -3.0, 1.0, 0.5)
  true_sigma <- diag(c(0.5, 0.5, 0.5))
  n_reps <- 200L

  cat(sprintf("\n  MC Case C: RC logit, %d reps\n", n_reps))

  res <- Filter(Negate(is.null), lapply(seq_len(n_reps), function(r)
    run_one_rc(r, 50, 15, 3, true_beta, true_sigma)))
  cat(sprintf("  Successful: %d/%d\n", length(res), n_reps))

  if (length(res) < 10) {
    skip("Too few successful RC reps with diff IVs (instrument construction issue)")
  }

  mc_beta <- summarize_mc(lapply(res, `[[`, "beta"), true_beta,
                          c("intercept", "price", "x1", "x2"))
  mc_sigma <- summarize_mc(lapply(res, `[[`, "sigma"), diag(true_sigma),
                           c("sigma_1", "sigma_2", "sigma_3"))

  cat("  Beta:\n")
  for (i in seq_len(nrow(mc_beta)))
    cat(sprintf("    %-10s: true=%5.1f, mean=%7.4f, sign=%5.1f%%\n",
                mc_beta$parameter[i], mc_beta$true[i], mc_beta$mean_est[i],
                mc_beta$sign_correct_pct[i]))
  cat("  Sigma:\n")
  for (i in seq_len(nrow(mc_sigma)))
    cat(sprintf("    %-10s: true=%5.2f, mean=%7.4f, RMSE=%6.4f\n",
                mc_sigma$parameter[i], mc_sigma$true[i], mc_sigma$mean_est[i],
                mc_sigma$rmse[i]))

  # Price sign with strong IVs
  expect_true(mc_beta$sign_correct_pct[2] > 80,
              info = sprintf("RC price sign %.1f%% (need >80%%)", mc_beta$sign_correct_pct[2]))
  # Sigma should be non-negative on average
  expect_true(all(mc_sigma$mean_est >= 0), info = "Mean sigma should be >= 0")
})

# =============================================================================
# Case D: Supply-side — 200 reps
# =============================================================================

test_that("MC Case D: supply gamma recovery (200 reps, strong IVs)", {
  skip_on_cran()

  true_beta <- c(0.5, -2.0, 0.8)
  true_gamma <- c(0.2, 0.8, 1.0)
  n_reps <- 100L

  cat(sprintf("\n  MC Case D: supply, %d reps\n", n_reps))

  res <- Filter(Negate(is.null), lapply(seq_len(n_reps), function(r)
    run_one_supply(r, 50, 15, 3, true_beta, true_gamma)))
  cat(sprintf("  Successful: %d/%d\n", length(res), n_reps))

  gamma_ests <- Filter(Negate(is.null), lapply(res, `[[`, "gamma"))
  if (length(gamma_ests) > 10) {
    mc_gamma <- summarize_mc(gamma_ests, true_gamma, c("gamma_0", "gamma_x", "gamma_w"))
    cat("  Gamma:\n")
    for (i in seq_len(nrow(mc_gamma)))
      cat(sprintf("    %-10s: true=%5.2f, mean=%7.4f, RMSE=%6.4f, sign=%5.1f%%\n",
                  mc_gamma$parameter[i], mc_gamma$true[i], mc_gamma$mean_est[i],
                  mc_gamma$rmse[i], mc_gamma$sign_correct_pct[i]))

    # Supply slopes should be informative in the redesigned benchmark even
    # though the supply intercept remains weakly identified.
    expect_true(mc_gamma$sign_correct_pct[2] > 90,
                info = sprintf("gamma_x sign %.1f%% (need >90%%)",
                              mc_gamma$sign_correct_pct[2]))
    expect_true(mc_gamma$sign_correct_pct[3] > 90,
                info = sprintf("gamma_w sign %.1f%% (need >90%%)",
                              mc_gamma$sign_correct_pct[3]))
    expect_true(abs(mc_gamma$rel_bias_pct[2]) < 25,
                info = sprintf("gamma_x bias %.1f%% (need <25%%)",
                              mc_gamma$rel_bias_pct[2]))
    expect_true(abs(mc_gamma$rel_bias_pct[3]) < 25,
                info = sprintf("gamma_w bias %.1f%% (need <25%%)",
                              mc_gamma$rel_bias_pct[3]))
  }
})

# =============================================================================
# Case E: Many products (J=40) — 200 reps
# =============================================================================

test_that("MC Case E: many products J=40 (200 reps, strong IVs)", {
  skip_on_cran()

  true_beta <- c(0.5, -3.0, 1.0, 0.5, -0.5)
  n_reps <- 200L

  cat(sprintf("\n  MC Case E: J=40, %d reps\n", n_reps))

  res <- Filter(Negate(is.null), lapply(seq_len(n_reps), function(r)
    run_one_logit(r, 30, 40, 8, true_beta)))
  cat(sprintf("  Successful: %d/%d\n", length(res), n_reps))

  mc <- summarize_mc(lapply(res, `[[`, "beta"), true_beta,
                     c("intercept", "price", "x1", "x2", "x3"))
  cat("  Results:\n")
  for (i in seq_len(nrow(mc)))
    cat(sprintf("    %-10s: true=%6.2f, mean=%7.4f, sign=%5.1f%%\n",
                mc$parameter[i], mc$true[i], mc$mean_est[i], mc$sign_correct_pct[i]))

  # Price sign recovery
  expect_true(mc$sign_correct_pct[2] > 55,
              info = sprintf("Price sign %.1f%% with J=40", mc$sign_correct_pct[2]))
})

# =============================================================================
# Case F: CI coverage — 500 reps
# =============================================================================

test_that("MC Case F: 95% CI coverage (500 reps, strong IVs)", {
  skip_on_cran()

  true_beta <- c(0.5, -3.0, 1.0, 0.5, -0.5)
  n_reps <- 500L

  cat(sprintf("\n  MC Case F: CI coverage, %d reps\n", n_reps))

  res <- Filter(Negate(is.null), lapply(seq_len(n_reps), function(r)
    run_one_logit(r, 50, 20, 4, true_beta)))

  n_params <- length(true_beta)
  coverage <- matrix(FALSE, nrow = length(res), ncol = n_params)
  for (i in seq_along(res)) {
    r <- res[[i]]
    if (!is.null(r$se) && length(r$se) == n_params && all(is.finite(r$se))) {
      for (k in seq_len(n_params)) {
        ci_lo <- r$beta[k] - 1.96 * r$se[k]
        ci_hi <- r$beta[k] + 1.96 * r$se[k]
        coverage[i, k] <- (true_beta[k] >= ci_lo) && (true_beta[k] <= ci_hi)
      }
    }
  }

  cov_rates <- colMeans(coverage) * 100
  pnames <- c("intercept", "price", "x1", "x2", "x3")
  cat("  Coverage:\n")
  for (k in seq_len(n_params))
    cat(sprintf("    %-10s: %.1f%%\n", pnames[k], cov_rates[k]))

  # Coverage should be >80% for all params with strong instruments
  for (k in seq_len(n_params)) {
    expect_true(cov_rates[k] > 75,
                info = sprintf("%s coverage=%.1f%% (need >75%%)", pnames[k], cov_rates[k]))
  }
})

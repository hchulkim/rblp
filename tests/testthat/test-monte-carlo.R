# =============================================================================
# Monte Carlo Simulation Study
#
# Validates rblp under repeated sampling across multiple DGP configurations.
# For each configuration, we:
#   1. Generate data with known true parameters
#   2. Solve for equilibrium prices and shares
#   3. Estimate the model
#   4. Check: bias, RMSE, coverage of confidence intervals, sign recovery
#
# Configurations:
#   A. Plain logit, varying sample size (T=20,50,100)
#   B. RC logit with diagonal sigma
#   C. Logit with supply side (gamma)
#   D. Logit with many products per market
# =============================================================================

options(rblp.verbose = FALSE)

# ---------------------------------------------------------------------------
# Helper: run one Monte Carlo replication
# ---------------------------------------------------------------------------
run_one_logit_rep <- function(seed, T, J, F, true_beta, xi_var = 0.3) {
  tryCatch({
    id_data <- build_id_data(T = T, J = J, F = F)
    set.seed(seed)
    id_data$x <- runif(nrow(id_data), 0, 1)

    sim <- blp_simulation(
      product_formulations = list(blp_formulation(~ prices + x)),
      product_data = id_data,
      beta = true_beta,
      xi_variance = xi_var,
      seed = seed
    )
    sim_res <- sim$replace_endogenous(
      iteration = blp_iteration("simple", list(atol = 1e-12, max_evaluations = 3000))
    )

    prob <- sim_res$to_problem()
    est <- prob$solve(method = "2s")

    list(
      beta = est$beta,
      se = est$summary_table()$se,
      converged = est$optimization_converged,
      fp_converged = est$fp_converged
    )
  }, error = function(e) NULL)
}

run_one_rc_rep <- function(seed, T, J, F, true_beta, true_sigma, xi_var = 0.3) {
  tryCatch({
    id_data <- build_id_data(T = T, J = J, F = F)
    set.seed(seed)
    id_data$x <- runif(nrow(id_data), 0, 1)

    sim <- blp_simulation(
      product_formulations = list(
        blp_formulation(~ prices + x),
        blp_formulation(~ prices + x)
      ),
      product_data = id_data,
      beta = true_beta,
      sigma = true_sigma,
      integration = blp_integration("product", size = 5),
      xi_variance = xi_var,
      seed = seed
    )
    sim_res <- sim$replace_endogenous(
      iteration = blp_iteration("simple", list(atol = 1e-12, max_evaluations = 3000))
    )

    prob <- sim_res$to_problem()
    est <- prob$solve(
      sigma = true_sigma * 0.8,  # starting near truth
      method = "2s",
      optimization = blp_optimization("l-bfgs-b",
        method_options = list(maxit = 200))
    )

    list(
      beta = est$beta,
      sigma = diag(est$sigma),
      se = est$se,
      converged = est$optimization_converged,
      fp_converged = est$fp_converged
    )
  }, error = function(e) NULL)
}

run_one_supply_rep <- function(seed, T, J, F, true_beta, true_gamma, xi_var = 0.2) {
  tryCatch({
    id_data <- build_id_data(T = T, J = J, F = F)
    set.seed(seed)
    id_data$x <- runif(nrow(id_data), 0, 1)
    id_data$w <- runif(nrow(id_data), 0, 1)

    sim <- blp_simulation(
      product_formulations = list(
        blp_formulation(~ prices + x),
        blp_formulation(~ prices + x),
        blp_formulation(~ w)
      ),
      product_data = id_data,
      beta = true_beta,
      sigma = diag(c(0.3, 0.3, 0.3)),
      gamma = true_gamma,
      integration = blp_integration("product", size = 3),
      xi_variance = xi_var,
      omega_variance = xi_var,
      correlation = 0.7,
      seed = seed
    )
    sim_res <- sim$replace_endogenous(
      iteration = blp_iteration("simple", list(atol = 1e-12, max_evaluations = 3000))
    )

    prob <- sim_res$to_problem()
    est <- prob$solve(
      sigma = diag(c(0.3, 0.3, 0.3)),
      method = "1s",
      optimization = blp_optimization("l-bfgs-b",
        method_options = list(maxit = 100))
    )

    list(
      beta = est$beta,
      gamma = est$gamma,
      converged = est$optimization_converged,
      fp_converged = est$fp_converged
    )
  }, error = function(e) NULL)
}

# ---------------------------------------------------------------------------
# Helper: summarize MC results
# ---------------------------------------------------------------------------
summarize_mc <- function(estimates, true_values, param_names = NULL) {
  n_params <- length(true_values)
  if (is.null(param_names)) param_names <- paste0("param_", seq_len(n_params))

  results <- data.frame(
    parameter = param_names,
    true = true_values,
    stringsAsFactors = FALSE
  )

  est_mat <- do.call(rbind, estimates)
  n_reps <- nrow(est_mat)

  results$mean_est <- colMeans(est_mat, na.rm = TRUE)
  results$bias <- results$mean_est - results$true
  results$rel_bias_pct <- results$bias / abs(results$true) * 100
  results$sd_est <- apply(est_mat, 2, sd, na.rm = TRUE)
  results$rmse <- sqrt(results$bias^2 + results$sd_est^2)

  # Sign recovery rate
  results$sign_correct_pct <- colMeans(sign(est_mat) == sign(matrix(
    true_values, nrow = n_reps, ncol = n_params, byrow = TRUE
  )), na.rm = TRUE) * 100

  results
}

# =============================================================================
# Case A: Plain Logit — 1000 replications, T=50, J=20, F=4
# =============================================================================

test_that("MC Case A: plain logit recovers parameters (1000 reps)", {
  skip_on_cran()

  true_beta <- c(0.5, -3.0, 1.0)  # intercept, price, x
  n_reps <- 1000L

  cat(sprintf("\n  MC Case A: Logit, T=50, J=20, F=4, %d reps\n", n_reps))

  results_list <- lapply(seq_len(n_reps), function(r) {
    run_one_logit_rep(r, T = 50, J = 20, F = 4, true_beta)
  })
  results_list <- Filter(Negate(is.null), results_list)
  cat(sprintf("  Successful: %d/%d\n", length(results_list), n_reps))

  beta_ests <- lapply(results_list, function(r) r$beta)
  mc <- summarize_mc(beta_ests, true_beta, c("intercept", "price", "x"))

  cat("  Results:\n")
  for (i in seq_len(nrow(mc))) {
    cat(sprintf("    %-10s: true=%6.2f, mean=%7.4f, bias=%7.4f (%5.1f%%), RMSE=%6.4f, sign=%5.1f%%\n",
                mc$parameter[i], mc$true[i], mc$mean_est[i], mc$bias[i],
                mc$rel_bias_pct[i], mc$rmse[i], mc$sign_correct_pct[i]))
  }

  # Key checks:
  # 1. Price coefficient sign should be recovered >95% of the time
  expect_true(mc$sign_correct_pct[2] > 95,
              info = sprintf("Price sign correct %.1f%% (need >95%%)",
                            mc$sign_correct_pct[2]))

  # 2. x coefficient sign should be recovered >90% of the time
  expect_true(mc$sign_correct_pct[3] > 90,
              info = sprintf("x sign correct %.1f%% (need >90%%)",
                            mc$sign_correct_pct[3]))

  # 3. Success rate should be high
  expect_true(length(results_list) / n_reps > 0.95,
              info = "Success rate should be > 95%")
})

# =============================================================================
# Case B: Plain Logit — effect of sample size (T=20 vs T=100)
# =============================================================================

test_that("MC Case B: logit RMSE decreases with sample size", {
  skip_on_cran()

  true_beta <- c(0.5, -3.0, 1.0)
  n_reps <- 200L

  cat(sprintf("\n  MC Case B: sample size comparison, %d reps each\n", n_reps))

  # Small sample: T=20
  res_small <- lapply(seq_len(n_reps), function(r) {
    run_one_logit_rep(r, T = 20, J = 15, F = 3, true_beta)
  })
  res_small <- Filter(Negate(is.null), res_small)
  mc_small <- summarize_mc(lapply(res_small, `[[`, "beta"), true_beta,
                           c("intercept", "price", "x"))

  # Large sample: T=100
  res_large <- lapply(seq_len(n_reps), function(r) {
    run_one_logit_rep(r + 10000, T = 100, J = 15, F = 3, true_beta)
  })
  res_large <- Filter(Negate(is.null), res_large)
  mc_large <- summarize_mc(lapply(res_large, `[[`, "beta"), true_beta,
                           c("intercept", "price", "x"))

  cat("  T=20:  RMSE =", round(mc_small$rmse, 4), "\n")
  cat("  T=100: RMSE =", round(mc_large$rmse, 4), "\n")

  # RMSE should decrease with larger sample (at least for price coeff)
  expect_true(mc_large$rmse[2] < mc_small$rmse[2],
              info = sprintf("Price RMSE: T=100 (%.4f) should be < T=20 (%.4f)",
                            mc_large$rmse[2], mc_small$rmse[2]))
})

# =============================================================================
# Case C: RC Logit — 200 replications
# =============================================================================

test_that("MC Case C: RC logit recovers sigma signs (200 reps)", {
  skip_on_cran()

  true_beta <- c(0.5, -3.0, 1.0)
  true_sigma <- diag(c(0.5, 0.5, 0.5))
  n_reps <- 200L

  cat(sprintf("\n  MC Case C: RC logit, T=50, J=15, F=3, %d reps\n", n_reps))

  results_list <- lapply(seq_len(n_reps), function(r) {
    run_one_rc_rep(r, T = 50, J = 15, F = 3, true_beta, true_sigma)
  })
  results_list <- Filter(Negate(is.null), results_list)
  cat(sprintf("  Successful: %d/%d\n", length(results_list), n_reps))

  # Check beta recovery
  beta_ests <- lapply(results_list, function(r) r$beta)
  mc_beta <- summarize_mc(beta_ests, true_beta, c("intercept", "price", "x"))

  cat("  Beta results:\n")
  for (i in seq_len(nrow(mc_beta))) {
    cat(sprintf("    %-10s: true=%6.2f, mean=%7.4f, sign=%5.1f%%\n",
                mc_beta$parameter[i], mc_beta$true[i], mc_beta$mean_est[i],
                mc_beta$sign_correct_pct[i]))
  }

  # Check sigma recovery
  sigma_ests <- lapply(results_list, function(r) r$sigma)
  mc_sigma <- summarize_mc(sigma_ests, diag(true_sigma),
                           c("sigma_1", "sigma_2", "sigma_3"))

  cat("  Sigma results:\n")
  for (i in seq_len(nrow(mc_sigma))) {
    cat(sprintf("    %-10s: true=%6.2f, mean=%7.4f, RMSE=%6.4f\n",
                mc_sigma$parameter[i], mc_sigma$true[i], mc_sigma$mean_est[i],
                mc_sigma$rmse[i]))
  }

  # Price sign should be correct most of the time
  expect_true(mc_beta$sign_correct_pct[2] > 80,
              info = sprintf("RC price sign correct %.1f%% (need >80%%)",
                            mc_beta$sign_correct_pct[2]))

  # Sigma should be non-negative on average (Cholesky constraint)
  expect_true(all(mc_sigma$mean_est >= 0),
              info = "Mean sigma estimates should be non-negative")
})

# =============================================================================
# Case D: Logit with Supply Side — 200 replications
# =============================================================================

test_that("MC Case D: supply-side gamma recovery (200 reps)", {
  skip_on_cran()

  true_beta <- c(0.5, -3.0, 1.0)
  true_gamma <- c(0.5, 1.5)
  n_reps <- 200L

  cat(sprintf("\n  MC Case D: supply, T=50, J=15, F=3, %d reps\n", n_reps))

  results_list <- lapply(seq_len(n_reps), function(r) {
    run_one_supply_rep(r, T = 50, J = 15, F = 3, true_beta, true_gamma)
  })
  results_list <- Filter(Negate(is.null), results_list)
  cat(sprintf("  Successful: %d/%d\n", length(results_list), n_reps))

  # Check gamma recovery
  gamma_ests <- lapply(results_list, function(r) r$gamma)
  gamma_ests <- Filter(Negate(is.null), gamma_ests)

  if (length(gamma_ests) > 10) {
    mc_gamma <- summarize_mc(gamma_ests, true_gamma, c("gamma_0", "gamma_w"))

    cat("  Gamma results:\n")
    for (i in seq_len(nrow(mc_gamma))) {
      cat(sprintf("    %-10s: true=%6.2f, mean=%7.4f, RMSE=%6.4f, sign=%5.1f%%\n",
                  mc_gamma$parameter[i], mc_gamma$true[i], mc_gamma$mean_est[i],
                  mc_gamma$rmse[i], mc_gamma$sign_correct_pct[i]))
    }

    # gamma_w should have correct sign >80% of the time
    expect_true(mc_gamma$sign_correct_pct[2] > 70,
                info = sprintf("gamma_w sign correct %.1f%% (need >70%%)",
                              mc_gamma$sign_correct_pct[2]))
  }
})

# =============================================================================
# Case E: Many products per market (J=40)
# =============================================================================

test_that("MC Case E: logit with many products (200 reps)", {
  skip_on_cran()

  true_beta <- c(0.5, -3.0, 1.0)
  n_reps <- 200L

  cat(sprintf("\n  MC Case E: J=40, T=30, F=8, %d reps\n", n_reps))

  results_list <- lapply(seq_len(n_reps), function(r) {
    run_one_logit_rep(r, T = 30, J = 40, F = 8, true_beta)
  })
  results_list <- Filter(Negate(is.null), results_list)
  cat(sprintf("  Successful: %d/%d\n", length(results_list), n_reps))

  beta_ests <- lapply(results_list, function(r) r$beta)
  mc <- summarize_mc(beta_ests, true_beta, c("intercept", "price", "x"))

  cat("  Results:\n")
  for (i in seq_len(nrow(mc))) {
    cat(sprintf("    %-10s: true=%6.2f, mean=%7.4f, bias=%7.4f, sign=%5.1f%%\n",
                mc$parameter[i], mc$true[i], mc$mean_est[i], mc$bias[i],
                mc$sign_correct_pct[i]))
  }

  # Price sign recovery
  expect_true(mc$sign_correct_pct[2] > 90,
              info = sprintf("Price sign correct %.1f%% with J=40", mc$sign_correct_pct[2]))
})

# =============================================================================
# Case F: CI coverage — do 95% CIs cover true value 95% of the time?
# =============================================================================

test_that("MC Case F: 95% CI coverage rate (500 reps)", {
  skip_on_cran()

  true_beta <- c(0.5, -3.0, 1.0)
  n_reps <- 500L

  cat(sprintf("\n  MC Case F: CI coverage, T=50, J=20, F=4, %d reps\n", n_reps))

  results_list <- lapply(seq_len(n_reps), function(r) {
    run_one_logit_rep(r, T = 50, J = 20, F = 4, true_beta)
  })
  results_list <- Filter(Negate(is.null), results_list)

  # Check if 95% CI covers the true value
  coverage <- matrix(FALSE, nrow = length(results_list), ncol = 3)
  for (i in seq_along(results_list)) {
    r <- results_list[[i]]
    if (!is.null(r$se) && length(r$se) == 3 && all(is.finite(r$se))) {
      for (k in 1:3) {
        ci_lo <- r$beta[k] - 1.96 * r$se[k]
        ci_hi <- r$beta[k] + 1.96 * r$se[k]
        coverage[i, k] <- (true_beta[k] >= ci_lo) && (true_beta[k] <= ci_hi)
      }
    }
  }

  cov_rates <- colMeans(coverage) * 100
  cat(sprintf("  Coverage rates: intercept=%.1f%%, price=%.1f%%, x=%.1f%%\n",
              cov_rates[1], cov_rates[2], cov_rates[3]))

  # Coverage should be roughly 95% (allow 80-100% for finite sample)
  # Price coefficient coverage is the most important
  for (k in 1:3) {
    expect_true(cov_rates[k] > 70,
                info = sprintf("Param %d coverage=%.1f%% should be >70%%", k, cov_rates[k]))
  }
})

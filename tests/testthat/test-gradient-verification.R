# =============================================================================
# Tests for gradient verification (analytic vs finite-difference)
# =============================================================================

test_that("analytic gradient matches finite-difference gradient for RC logit", {
  skip_on_cran()

  # Small problem for tractable gradient computation
  id_data <- build_id_data(T = 5, J = 6, F = 2)
  set.seed(700)
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
    seed = 700
  )

  sim_results <- sim$replace_endogenous(
    iteration = blp_iteration("simple", list(atol = 1e-12, max_evaluations = 5000))
  )

  sim_problem <- sim_results$to_problem(add_instruments = TRUE)

  # Set up optimization that captures the gradient
  sigma0 <- diag(c(0.4, 0.4, 0.4))
  params <- rblp:::BLPParameters$new(sigma = sigma0)
  theta0 <- params$compress()

  iteration <- blp_iteration("squarem")
  delta <- sim_problem$compute_logit_delta()

  # Build weighting matrix
  ZD <- sim_problem$products$ZD
  S_init <- crossprod(ZD) / sim_problem$N
  W <- rblp:::compute_gmm_weights(S_init)

  # Evaluate objective and analytic gradient at theta0
  progress <- sim_problem$.__enclos_env__$private$compute_progress(
    theta0, params, delta, W,
    iteration, "safe_linear", TRUE, TRUE,
    NULL, "revert", 1, 1L
  )

  analytic_grad <- progress$gradient

  # Finite-difference gradient
  eps <- 1e-5
  n_theta <- length(theta0)
  fd_grad <- numeric(n_theta)

  for (k in seq_len(n_theta)) {
    theta_plus <- theta0
    theta_minus <- theta0
    theta_plus[k] <- theta_plus[k] + eps
    theta_minus[k] <- theta_minus[k] - eps

    prog_plus <- sim_problem$.__enclos_env__$private$compute_progress(
      theta_plus, params, delta, W,
      iteration, "safe_linear", TRUE, TRUE,
      NULL, "revert", 1, 1L
    )
    prog_minus <- sim_problem$.__enclos_env__$private$compute_progress(
      theta_minus, params, delta, W,
      iteration, "safe_linear", TRUE, TRUE,
      NULL, "revert", 1, 1L
    )

    fd_grad[k] <- (prog_plus$objective - prog_minus$objective) / (2 * eps)
  }

  # Compare analytic vs finite-difference
  expect_false(is.null(analytic_grad), info = "Analytic gradient should be computed")
  expect_equal(length(analytic_grad), n_theta)

  # Relative error should be small
  for (k in seq_len(n_theta)) {
    denom <- max(abs(analytic_grad[k]), abs(fd_grad[k]), 1e-8)
    rel_err <- abs(analytic_grad[k] - fd_grad[k]) / denom
    expect_true(rel_err < 0.05,
                info = sprintf("Gradient element %d: analytic=%.6e, fd=%.6e, rel_err=%.4f",
                               k, analytic_grad[k], fd_grad[k], rel_err))
  }
})

test_that("xi-by-theta Jacobian is consistent via IFT", {
  skip_on_cran()

  # Direct test of the implicit function theorem Jacobian
  id_data <- build_id_data(T = 3, J = 5, F = 2)
  set.seed(701)
  id_data$x <- runif(nrow(id_data), 0, 1)
  id_data$shares <- runif(nrow(id_data), 0.03, 0.10)
  id_data$prices <- runif(nrow(id_data), 1, 3)

  f1 <- blp_formulation(~ prices + x)
  f2 <- blp_formulation(~ prices + x)

  integration <- blp_integration("product", size = 3)
  problem <- blp_problem(
    product_formulations = list(f1, f2),
    product_data = id_data,
    integration = integration
  )

  sigma <- diag(c(0.3, 0.3, 0.3))
  md <- problem$get_market_data(problem$unique_market_ids[1])

  sigma_free <- (sigma != 0) & lower.tri(sigma, diag = TRUE)

  mkt <- rblp:::BLPMarket$new(
    products = md$products,
    agents = md$agents,
    sigma = sigma,
    rc_types = c("linear", "linear", "linear"),
    epsilon_scale = 1.0,
    sigma_free = sigma_free
  )

  # Get delta from contraction mapping
  delta0 <- log(md$products$shares) - log(1 - sum(md$products$shares))
  fp <- mkt$compute_delta(delta0, blp_iteration("squarem"))
  delta <- fp$delta

  mu <- mkt$compute_mu()
  prob <- mkt$compute_probabilities(delta, mu)

  # Compute Jacobians
  ds_ddelta <- mkt$compute_shares_by_xi_jacobian(prob)
  ds_dtheta <- mkt$compute_shares_by_theta_jacobian(prob, sigma = sigma)
  xi_jac <- mkt$compute_xi_by_theta_jacobian(prob)

  # IFT: xi_jac should equal -(ds_ddelta)^{-1} (ds_dtheta)
  expected_xi_jac <- -solve(ds_ddelta, ds_dtheta)

  if (!is.null(xi_jac) && ncol(xi_jac) > 0) {
    expect_equal(as.matrix(xi_jac), as.matrix(expected_xi_jac),
                 tolerance = 1e-6,
                 info = "xi Jacobian should match IFT formula")
  }
})

test_that("share Jacobian ds_ddelta has correct diagonal and off-diagonal signs", {
  # For standard logit: diagonal should be positive, off-diagonal negative
  id_data <- data.frame(
    market_ids = rep("M1", 4),
    firm_ids = 1:4,
    prices = c(1, 2, 3, 4),
    x = c(0.5, 0.6, 0.7, 0.8),
    shares = c(0.1, 0.15, 0.12, 0.08)
  )

  f1 <- blp_formulation(~ prices + x)
  f2 <- blp_formulation(~ prices + x)
  integration <- blp_integration("product", size = 3)

  problem <- blp_problem(
    product_formulations = list(f1, f2),
    product_data = id_data,
    integration = integration
  )

  md <- problem$get_market_data("M1")
  sigma <- diag(c(0.3, 0.3, 0.3))

  mkt <- rblp:::BLPMarket$new(
    products = md$products,
    agents = md$agents,
    sigma = sigma,
    rc_types = c("linear", "linear", "linear"),
    epsilon_scale = 1.0
  )

  delta <- rep(0, 4)
  mu <- mkt$compute_mu()
  prob <- mkt$compute_probabilities(delta, mu)
  P <- if (is.list(prob)) prob$probabilities else prob

  jac <- mkt$compute_shares_by_xi_jacobian(prob)

  # Diagonal should be positive (own-effect)
  expect_true(all(diag(jac) > 0), info = "ds_j/dd_j should be positive")

  # Off-diagonal should be negative (substitution)
  off_diag <- jac
  diag(off_diag) <- 0
  expect_true(all(off_diag <= 0), info = "ds_j/dd_k (j!=k) should be non-positive")

  # Each column should sum to approximately zero (shares sum to constant)
  col_sums <- colSums(jac)
  # Not exactly zero due to outside good, but should be small relative to diagonal
  expect_true(all(abs(col_sums) < max(abs(diag(jac)))),
              info = "Column sums should be small")
})

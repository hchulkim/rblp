# =============================================================================
# Tests for data.table compatibility
# =============================================================================

test_that("blp_problem accepts data.table product data", {
  skip_on_cran()
  skip_if_not_installed("data.table")

  products <- data.table::as.data.table(load_nevo_products())
  expect_true(data.table::is.data.table(products))

  f1 <- blp_formulation(~ prices + sugar + mushy)
  problem <- blp_problem(list(f1), products)
  expect_equal(problem$N, 2256)

  results <- problem$solve(method = "1s")
  expect_true(results$beta[2] < 0, info = "Price should be negative")
})

test_that("blp_problem accepts data.table agent data", {
  skip_on_cran()
  skip_if_not_installed("data.table")

  products <- data.table::as.data.table(load_nevo_products())
  agents <- data.table::as.data.table(load_nevo_agents())
  expect_true(data.table::is.data.table(agents))

  f1 <- blp_formulation(~ 0 + prices, absorb = ~ product_ids)
  f2 <- blp_formulation(~ prices + sugar + mushy)
  demo_form <- blp_formulation(~ 0 + income + income_squared + age + child)

  problem <- blp_problem(list(f1, f2), products,
                         agent_formulation = demo_form,
                         agent_data = agents)
  expect_equal(problem$I, 1880)
  expect_equal(problem$D, 4)
})

test_that("blp_simulation accepts data.table product data", {
  skip_on_cran()
  skip_if_not_installed("data.table")

  id_dt <- data.table::as.data.table(build_id_data(T = 5, J = 5, F = 2))
  data.table::set(id_dt, j = "x", value = runif(nrow(id_dt)))

  sim <- blp_simulation(
    list(blp_formulation(~ prices + x)),
    id_dt,
    beta = c(0.5, -2, 0.8),
    seed = 42
  )

  sim_res <- sim$replace_endogenous()
  expect_true(all(sim_res$shares > 0))
  expect_true(all(sim_res$prices > 0))

  # to_problem should work from data.table-originated simulation
  prob <- sim_res$to_problem()
  est <- prob$solve(method = "1s")
  expect_equal(length(est$beta), 3)
})

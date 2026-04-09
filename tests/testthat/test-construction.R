test_that("build_blp_instruments creates correct size", {
  X <- matrix(c(1, 2, 3, 4, 5, 6), ncol = 2)
  market_ids <- c(1, 1, 1)
  firm_ids <- c(1, 1, 2)

  iv <- build_blp_instruments(X, market_ids, firm_ids)
  expect_equal(nrow(iv), 3)
  expect_equal(ncol(iv), 4)  # 2 own + 2 rival
})

test_that("build_blp_instruments sums correctly", {
  X <- matrix(c(10, 20, 30, 1, 2, 3), ncol = 2)
  market_ids <- c(1, 1, 1)
  firm_ids <- c(1, 1, 2)

  iv <- build_blp_instruments(X, market_ids, firm_ids)

  # Product 1 (firm 1): own_other = X[2,] = (20, 2), rival = X[3,] = (30, 3)
  expect_equal(unname(iv[1, 1]), 20)
  expect_equal(unname(iv[1, 2]), 2)
  expect_equal(unname(iv[1, 3]), 30)
  expect_equal(unname(iv[1, 4]), 3)
})

test_that("build_id_data creates balanced panel", {
  id_data <- build_id_data(T = 3, J = 6, F = 2)
  expect_equal(nrow(id_data), 18)  # 3 * 6
  expect_equal(length(unique(id_data$market_ids)), 3)

  # Each market should have equal firm representation
  for (t in unique(id_data$market_ids)) {
    firms <- id_data$firm_ids[id_data$market_ids == t]
    expect_equal(length(firms), 6)
    expect_equal(length(unique(firms)), 2)
  }
})

test_that("build_differentiation_instruments works", {
  X <- matrix(c(1, 2, 5, 0.1, 0.2, 0.5), ncol = 2)
  market_ids <- c(1, 1, 1)
  firm_ids <- c(1, 1, 2)

  iv_local <- build_differentiation_instruments(X, market_ids, firm_ids, "local")
  iv_quad <- build_differentiation_instruments(X, market_ids, firm_ids, "quadratic")

  expect_equal(nrow(iv_local), 3)
  expect_equal(ncol(iv_local), 4)
  expect_equal(nrow(iv_quad), 3)
  expect_equal(ncol(iv_quad), 4)
})

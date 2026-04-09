test_that("BLPFormulation creates design matrix from formula", {
  f <- blp_formulation(~ x1 + x2)
  df <- data.frame(x1 = c(1, 2, 3), x2 = c(4, 5, 6))
  X <- f$build_matrix(df)
  expect_equal(ncol(X), 3)  # intercept + x1 + x2
  expect_equal(nrow(X), 3)
  expect_equal(colnames(X)[1], "(Intercept)")
})

test_that("BLPFormulation detects prices after build", {
  f <- blp_formulation(~ prices + x1)
  df <- data.frame(prices = 1:3, x1 = 4:6)
  f$build_matrix(df)
  expect_true(f$has_prices())
  expect_false(f$has_shares())
})

test_that("BLPFormulation works without intercept", {
  f <- blp_formulation(~ 0 + x1 + x2)
  df <- data.frame(x1 = 1:3, x2 = 4:6)
  X <- f$build_matrix(df)
  expect_equal(ncol(X), 2)
})

test_that("BLPFormulation from string works", {
  f <- blp_formulation("~ prices + sugar + mushy")
  df <- data.frame(prices = 1:3, sugar = 4:6, mushy = c(0, 1, 0))
  f$build_matrix(df)
  expect_true(f$has_prices())
})

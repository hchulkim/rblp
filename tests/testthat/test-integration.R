test_that("BLPIntegration product rule produces correct size", {
  int <- blp_integration("product", size = 5)
  result <- int$build(2)  # 2 dimensions
  expect_equal(ncol(result$nodes), 2)
  expect_equal(nrow(result$nodes), 25)  # 5^2
  expect_equal(length(result$weights), 25)
  expect_true(abs(sum(result$weights) - 1) < 1e-10)
})

test_that("BLPIntegration monte carlo works", {
  int <- blp_integration("monte_carlo", size = 100, seed = 42)
  result <- int$build(3)
  expect_equal(ncol(result$nodes), 3)
  expect_equal(nrow(result$nodes), 100)
  expect_equal(length(result$weights), 100)
  expect_true(abs(sum(result$weights) - 1) < 1e-10)
})

test_that("BLPIntegration halton works", {
  int <- blp_integration("halton", size = 50, seed = 0)
  result <- int$build(2)
  expect_equal(ncol(result$nodes), 2)
  expect_equal(nrow(result$nodes), 50)
})

test_that("Gauss-Hermite integration is accurate", {
  int <- blp_integration("product", size = 7)
  result <- int$build(1)
  # Integration of f(x) = 1 over N(0,1) should be 1
  expect_true(abs(sum(result$weights) - 1) < 1e-10)
  # E[X] should be 0
  expect_true(abs(sum(result$weights * result$nodes[, 1])) < 1e-10)
  # E[X^2] should be 1
  expect_true(abs(sum(result$weights * result$nodes[, 1]^2) - 1) < 1e-8)
})

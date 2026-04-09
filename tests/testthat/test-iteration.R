test_that("Simple iteration converges for contraction", {
  iter <- blp_iteration("simple", list(atol = 1e-12, max_evaluations = 1000))
  contraction <- function(x) x / 2 + 1  # fixed point at x = 2
  result <- iter$iterate(0, contraction)
  expect_true(result$converged)
  expect_true(abs(result$values - 2) < 1e-10)
})

test_that("SQUAREM iteration converges faster", {
  iter_simple <- blp_iteration("simple", list(atol = 1e-12))
  iter_squarem <- blp_iteration("squarem", list(atol = 1e-12))

  contraction <- function(x) x / 2 + 1
  r1 <- iter_simple$iterate(0, contraction)
  r2 <- iter_squarem$iterate(0, contraction)

  expect_true(r1$converged)
  expect_true(r2$converged)
  expect_true(abs(r1$values - 2) < 1e-10)
  expect_true(abs(r2$values - 2) < 1e-10)
})

test_that("Return iteration returns initial values", {
  iter <- blp_iteration("return")
  result <- iter$iterate(42, function(x) x)
  expect_equal(result$values, 42)
  expect_true(result$converged)
})

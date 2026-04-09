test_that("Logit model produces reasonable estimates", {

  # Load Nevo data
  products <- load_nevo_products()
  expect_true(nrow(products) > 0)
  expect_true("market_ids" %in% names(products))
  expect_true("shares" %in% names(products))
  expect_true("prices" %in% names(products))

  # Formulation
  f1 <- blp_formulation(~ prices)

  # Create problem
  problem <- blp_problem(list(f1), products)
  expect_s3_class(problem, "BLPProblem")
  expect_equal(problem$K1, 2)  # intercept + prices
  expect_equal(problem$K2, 0)

  # Solve logit (no nonlinear params)
  results <- problem$solve()
  expect_s3_class(results, "BLPResults")

  # Beta should be estimated
  expect_equal(length(results$beta), 2)
  # Price coefficient should be negative
  expect_true(results$beta[2] < 0)

  # Objective should be finite and non-negative
  expect_true(is.finite(results$objective))
  expect_true(results$objective >= 0)

  # Xi should be residuals
  expect_equal(length(results$xi), nrow(products))
})

test_that("Logit model with multiple characteristics works", {

  products <- load_nevo_products()

  f1 <- blp_formulation(~ prices + sugar + mushy)
  problem <- blp_problem(list(f1), products)

  expect_equal(problem$K1, 4)  # intercept + 3 vars

  results <- problem$solve()
  expect_equal(length(results$beta), 4)
  expect_true(results$beta[which(colnames(problem$products$X1) == "prices")] < 0)
})

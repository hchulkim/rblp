test_that("BLPParameters compress and expand sigma", {
  sigma <- matrix(c(1, 0.5, 0, 2), 2, 2)
  params <- rblp:::BLPParameters$new(sigma = sigma)

  theta <- params$compress()
  expect_equal(length(theta), 3)  # 2 diagonal + 1 lower tri

  expanded <- params$expand(theta)
  expect_equal(expanded$sigma[1, 1], 1)
  expect_equal(expanded$sigma[2, 1], 0.5)
  expect_equal(expanded$sigma[2, 2], 2)
  expect_equal(expanded$sigma[1, 2], 0)  # upper tri stays 0
})

test_that("BLPParameters handles fixed elements", {
  sigma <- matrix(c(1, 0, 0, 2), 2, 2)
  params <- rblp:::BLPParameters$new(sigma = sigma)

  theta <- params$compress()
  expect_equal(length(theta), 2)  # only diagonal

  expanded <- params$expand(c(3, 4))
  expect_equal(expanded$sigma[1, 1], 3)
  expect_equal(expanded$sigma[2, 2], 4)
  expect_equal(expanded$sigma[2, 1], 0)
})

test_that("BLPParameters bounds are correct", {
  sigma <- matrix(c(1, 0.5, 0, 2), 2, 2)
  params <- rblp:::BLPParameters$new(sigma = sigma)
  bounds <- params$get_bounds()

  expect_equal(length(bounds$lower), 3)
  # Diagonal should have lower bound 0
  expect_equal(bounds$lower[1], 0)
  # Off-diagonal should have lower bound -Inf
  expect_equal(bounds$lower[2], -Inf)
})

test_that("BLPParameters labels are generated", {
  sigma <- matrix(c(1, 0.5, 0, 2), 2, 2)
  params <- rblp:::BLPParameters$new(sigma = sigma)
  labels <- params$get_labels()
  expect_equal(length(labels), 3)
  expect_equal(labels[1], "sigma[1,1]")
})

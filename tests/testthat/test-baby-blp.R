# =============================================================================
# Replication of Lei Ma's baby_BLP (https://github.com/leima0521/baby_BLP)
# OTC drug data: 11 products, 2 stores, 48 weeks = 96 markets
#
# Reference specifications:
#   make_plots.R: 5 logit models (OLS + IV with cost as instrument)
#   baby_blp.jl: BLP with random coefficient on constant, 500 simulation draws
#
# Key features of this data:
#   - Price is essentially exogenous (per README) — "retail scanner data"
#   - Very small inside shares (~0.1-0.6% per product), large outside option
#   - 11 OTC drug products across 4 brands (Tylenol, Advil, Bayer, Store brand)
#   - Without product FE, price coefficient is POSITIVE (omitted brand quality)
#   - With product FE, price coefficient is negative
# =============================================================================

options(rblp.verbose = FALSE)

load_baby_blp_data <- function() {
  csv_path <- file.path(testthat::test_path(), "otc_baby_blp.csv")
  if (!file.exists(csv_path)) {
    csv_path <- "tests/testthat/otc_baby_blp.csv"
  }
  if (!file.exists(csv_path)) skip("baby_BLP OTC data not found")

  otc <- read.csv(csv_path, stringsAsFactors = FALSE)

  data.frame(
    market_ids = as.character(otc$mkt),
    firm_ids = as.character(otc$product),
    shares = otc$mkt_share,
    prices = otc$price,
    promotion = otc$promotion,
    cost = otc$cost,
    product = factor(otc$product),
    store = factor(otc$store),
    ln_mkt_share_diff = otc$ln_mkt_share_diff,
    stringsAsFactors = FALSE
  )
}

# =============================================================================
# 1. OLS logit baselines — reproduce make_plots.R models m1, m2
# =============================================================================

test_that("baby_BLP: OLS logit baselines match lm()", {
  skip_on_cran()
  products <- load_baby_blp_data()

  # m1: no FE — price is POSITIVE due to omitted brand quality
  m1 <- lm(ln_mkt_share_diff ~ prices + promotion, data = products)
  cat(sprintf("\n  m1 (no FE):  price=%.4f, promo=%.4f\n",
              coef(m1)["prices"], coef(m1)["promotion"]))
  expect_true(coef(m1)["prices"] > 0,
              info = "Without FE, price is positive (omitted brand quality bias)")
  expect_true(coef(m1)["promotion"] > 0,
              info = "Promotion should increase demand")

  # m2: with product FE — price turns NEGATIVE (correct sign)
  m2 <- lm(ln_mkt_share_diff ~ prices + promotion + product, data = products)
  cat(sprintf("  m2 (prod FE): price=%.4f, promo=%.4f\n",
              coef(m2)["prices"], coef(m2)["promotion"]))
  expect_true(coef(m2)["prices"] < 0,
              info = "With product FE, price should be negative")
  expect_true(coef(m2)["promotion"] > 0,
              info = "Promotion should be positive")
})

# =============================================================================
# 2. rblp logit with product FE — cross-validate against OLS+FE
# =============================================================================

test_that("baby_BLP: rblp logit FE matches OLS FE coefficients", {
  skip_on_cran()
  products <- load_baby_blp_data()

  # OLS+FE reference
  m2 <- lm(ln_mkt_share_diff ~ prices + promotion + product, data = products)
  ols_price <- coef(m2)["prices"]
  ols_promo <- coef(m2)["promotion"]

  # rblp logit with absorbed product FE
  f1 <- blp_formulation(~ prices + promotion, absorb = ~ product)
  problem <- blp_problem(list(f1), products)
  results <- problem$solve(method = "1s")

  cat(sprintf("  rblp FE: price=%.4f, promo=%.4f\n",
              results$beta[1], results$beta[2]))
  cat(sprintf("  OLS FE:  price=%.4f, promo=%.4f\n", ols_price, ols_promo))

  # rblp should produce negative price coefficient with FE
  expect_true(results$beta[1] < 0,
              info = sprintf("rblp FE price=%.4f should be negative", results$beta[1]))
  expect_true(results$beta[2] > 0,
              info = sprintf("rblp FE promo=%.4f should be positive", results$beta[2]))

  # Should be close to OLS FE (same estimator when no excluded instruments)
  expect_true(abs(results$beta[1] - ols_price) < 0.5,
              info = sprintf("rblp price=%.4f vs OLS=%.4f, gap > 0.5",
                            results$beta[1], ols_price))
})

# =============================================================================
# 3. rblp IV logit with cost instrument + product FE
#    Matches make_plots.R model m5
# =============================================================================

test_that("baby_BLP: rblp IV logit with cost instrument", {
  skip_on_cran()
  products <- load_baby_blp_data()

  # Add cost as excluded instrument
  products$demand_instruments0 <- products$cost

  f1 <- blp_formulation(~ prices + promotion, absorb = ~ product)
  problem <- blp_problem(list(f1), products)

  cat(sprintf("  IV problem: K1=%d, MD=%d, N=%d, T=%d\n",
              problem$K1, problem$MD, problem$N, problem$T))

  results <- problem$solve(method = "2s")
  tbl <- results$summary_table()

  cat(sprintf("  IV+FE (2s): price=%.4f (SE=%.4f), promo=%.4f (SE=%.4f)\n",
              results$beta[1], tbl$se[1], results$beta[2], tbl$se[2]))

  # SEs should be finite and positive
  expect_true(all(is.finite(tbl$se) & tbl$se > 0),
              info = "IV SEs should be finite and positive")

  # Objective should be finite
  expect_true(is.finite(results$objective),
              info = "IV objective should be finite")

  # Promotion should remain positive
  expect_true(results$beta[2] > 0,
              info = "IV promo coefficient should be positive")
})

# =============================================================================
# 4. BLP RC logit — matching baby_blp.jl specification
#    Random coefficient on constant only (sigma is scalar)
# =============================================================================

test_that("baby_BLP: RC logit with random coefficient on constant", {
  skip_on_cran()
  products <- load_baby_blp_data()

  # Include cost as instrument (matching baby_blp.jl Z = [X, cost])
  products$demand_instruments0 <- products$cost

  # X1 = prices + promotion + product FE (absorbed)
  # X2 = intercept only (random coefficient on constant)
  f1 <- blp_formulation(~ prices + promotion, absorb = ~ product)
  f2 <- blp_formulation(~ 1)

  problem <- blp_problem(
    product_formulations = list(f1, f2),
    product_data = products,
    integration = blp_integration("product", size = 7)
  )

  expect_equal(problem$K2, 1, info = "K2 should be 1 (constant only)")
  cat(sprintf("  RC: K1=%d, K2=%d, MD=%d, I=%d\n",
              problem$K1, problem$K2, problem$MD, problem$I))

  # Solve with sigma starting value = 1 (matching baby_blp.jl)
  results <- problem$solve(
    sigma = matrix(1, 1, 1),
    method = "1s",
    optimization = blp_optimization("l-bfgs-b",
      method_options = list(maxit = 500))
  )

  tbl <- results$summary_table()
  cat(sprintf("  RC: price=%.4f, promo=%.4f, sigma=%.4f\n",
              results$beta[1], results$beta[2], results$sigma[1, 1]))
  cat(sprintf("  RC: objective=%.6f, converged=%s\n",
              results$objective, results$optimization_converged))

  # Sigma should be non-negative (Cholesky lower bound)
  expect_true(results$sigma[1, 1] >= 0,
              info = "sigma should be non-negative")

  # Fixed point should converge
  expect_true(results$fp_converged,
              info = "Fixed point should converge")

  # Objective should be finite
  expect_true(is.finite(results$objective),
              info = "Objective should be finite")

  # Summary table should have sigma row
  sigma_rows <- tbl[tbl$type == "nonlinear (sigma)", ]
  expect_true(nrow(sigma_rows) >= 1,
              info = "Should have at least 1 sigma in summary")
})

# =============================================================================
# 5. Elasticities on OTC data
# =============================================================================

test_that("baby_BLP: elasticities have correct structure", {
  skip_on_cran()
  products <- load_baby_blp_data()

  f1 <- blp_formulation(~ prices + promotion, absorb = ~ product)
  problem <- blp_problem(list(f1), products)
  results <- problem$solve(method = "1s")

  E <- results$compute_elasticities("1")

  expect_equal(nrow(E), 11, info = "11 products per market")
  expect_equal(ncol(E), 11)

  # Own-price elasticities should be negative
  expect_true(all(diag(E) < 0),
              info = "Own-price elasticities should be negative")

  # Cross-price should be non-negative
  off_diag <- E; diag(off_diag) <- 0
  expect_true(all(off_diag >= -1e-10),
              info = "Cross-price elasticities should be non-negative")

  cat(sprintf("  Own-price elast range: [%.3f, %.3f]\n", min(diag(E)), max(diag(E))))

  # IIA check: in logit, cross-elasticities in a column are identical
  for (k in 1:3) {
    cross <- E[-k, k]
    cv <- sd(cross) / abs(mean(cross))
    expect_true(cv < 1e-5,
                info = sprintf("IIA violated for product %d (CV=%.2e)", k, cv))
  }
})

# =============================================================================
# 6. Consumer surplus on OTC data
# =============================================================================

test_that("baby_BLP: consumer surplus is positive across all OTC markets", {
  skip_on_cran()
  products <- load_baby_blp_data()

  f1 <- blp_formulation(~ prices + promotion, absorb = ~ product)
  problem <- blp_problem(list(f1), products)
  results <- problem$solve(method = "1s")

  cs <- results$compute_consumer_surplus()

  expect_equal(length(cs), 96, info = "96 markets")
  expect_true(all(cs > 0), info = "CS should be positive in all markets")
  expect_true(all(is.finite(cs)), info = "CS should be finite")

  cat(sprintf("  CS: mean=%.6f, min=%.6f, max=%.6f\n", mean(cs), min(cs), max(cs)))
})

# =============================================================================
# 7. HHI on OTC data — 11 "firms" (products treated as firms)
# =============================================================================

test_that("baby_BLP: HHI on OTC data", {
  skip_on_cran()
  products <- load_baby_blp_data()

  f1 <- blp_formulation(~ prices + promotion, absorb = ~ product)
  problem <- blp_problem(list(f1), products)
  results <- problem$solve(method = "1s")

  hhi <- results$compute_hhi()

  expect_equal(length(hhi), 96)
  expect_true(all(hhi > 0 & hhi <= 10000))
  expect_true(all(is.finite(hhi)))

  cat(sprintf("  HHI: mean=%.0f, min=%.0f, max=%.0f\n", mean(hhi), min(hhi), max(hhi)))
})

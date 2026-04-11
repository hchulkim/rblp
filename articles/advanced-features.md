# Advanced Features: Optimal Instruments, Bootstrap, and More

## Overview

This vignette demonstrates advanced `rblp` features that go beyond basic
estimation: aggregate elasticities, cost passthrough, long-run
diversion, optimal instruments, parametric bootstrap, and importance
sampling. These tools are essential for policy analysis and robustness
checks.

## Setup: Estimate a Baseline Model

``` r
library(rblp)

products <- load_nevo_products()
f1 <- blp_formulation(~ prices + sugar + mushy)
problem <- blp_problem(list(f1), products)
results <- problem$solve(method = "2s")
results$summary_table()[, c("parameter", "estimate", "se")]
#>     parameter     estimate          se
#> 1 (Intercept)  -2.78233619 0.107316633
#> 2      prices -11.75633541 0.845317475
#> 3       sugar   0.04752765 0.004146691
#> 4       mushy   0.05934698 0.051360224
```

## Aggregate Elasticities

The aggregate elasticity measures the percentage change in **total
inside market share** from a uniform 1% increase in all prices. Unlike
the $J \times J$ product-level elasticity matrix, this is a scalar per
market — useful for quick summaries.

``` r
ae <- results$compute_aggregate_elasticities()

cat(sprintf("Aggregate elasticities across %d markets:\n", length(ae)))
#> Aggregate elasticities across 94 markets:
cat(sprintf("  Mean: %.4f\n", mean(ae)))
#>   Mean: -0.7501
cat(sprintf("  Range: [%.4f, %.4f]\n", min(ae), max(ae)))
#>   Range: [-1.2829, -0.4076]
```

All aggregate elasticities are negative: when all inside goods become
more expensive, consumers substitute toward the outside option.

For a specific market:

``` r
ae_m1 <- results$compute_aggregate_elasticities(
  market_id = problem$unique_market_ids[1]
)
cat(sprintf("Market '%s' aggregate elasticity: %.4f\n",
            problem$unique_market_ids[1], ae_m1))
#> Market 'C01Q1' aggregate elasticity: -0.8302
```

## Cost-to-Price Passthrough

The passthrough matrix $dp/dc$ measures how a \$1 increase in the
marginal cost of product $k$ changes the equilibrium price of product
$j$. Under Bertrand-Nash pricing, this is derived from the implicit
function theorem applied to the system of first-order conditions.

``` r
t1 <- problem$unique_market_ids[1]
PT <- results$compute_passthrough(t1)

cat(sprintf("Passthrough matrix: %d x %d\n", nrow(PT), ncol(PT)))
#> Passthrough matrix: 24 x 24
cat("\nOwn-cost passthrough (diagonal, first 5 products):\n")
#> 
#> Own-cost passthrough (diagonal, first 5 products):
cat(sprintf("  %s\n", paste(round(diag(PT)[1:5], 4), collapse = ", ")))
#>   1.1685, 1.1002, 1.1776, 1.0723, 1.2613
```

Under perfect competition, own-cost passthrough equals 1 (costs are
fully passed through to prices). Under oligopoly, it is typically less
than 1 — firms absorb some of the cost increase.

## Long-Run Diversion Ratios

Standard (short-run) diversion ratios hold rivals’ prices fixed.
Long-run diversion accounts for equilibrium price adjustments: when
product $j$’s cost rises and its price increases, rivals also adjust
their prices. This is more relevant for merger analysis.

``` r
D_sr <- results$compute_diversion_ratios(t1)
D_lr <- results$compute_long_run_diversion_ratios(t1)

cat("Short-run vs Long-run diversion (product 1 -> products 2-5):\n")
#> Short-run vs Long-run diversion (product 1 -> products 2-5):
for (k in 2:5) {
  cat(sprintf("  1 -> %d: SR = %.4f, LR = %.4f\n",
              k, D_sr[k, 1], D_lr[k, 1]))
}
#>   1 -> 2: SR = 0.0125, LR = 0.0142
#>   1 -> 3: SR = 0.0126, LR = 0.0143
#>   1 -> 4: SR = 0.0125, LR = 0.0142
#>   1 -> 5: SR = 0.0126, LR = 0.0144
```

Long-run diversion is typically larger because rivals raise prices when
the focal product becomes more expensive, redirecting even more demand.

## Predicted Shares and Profits

``` r
s_pred <- results$compute_shares()
profits <- results$compute_profits()

cat(sprintf("Predicted shares: range [%.6f, %.6f]\n", min(s_pred), max(s_pred)))
#> Predicted shares: range [0.000182, 0.446883]
cat(sprintf("Correlation with observed: %.6f\n",
            cor(s_pred, problem$products$shares)))
#> Correlation with observed: 1.000000
cat(sprintf("Profits: range [%.6f, %.6f]\n", min(profits), max(profits)))
#> Profits: range [0.000019, 0.084624]
```

## Optimal Instruments

Standard BLP instruments (sums of rival characteristics) may be weak.
Optimal instruments (Chamberlain 1987, BLP 1999) project the Jacobian of
the structural error onto the instrument space, yielding efficient GMM
estimates.

``` r
oi <- results$compute_optimal_instruments()
cat(sprintf("Optimal demand instruments: %d x %d\n",
            nrow(oi$optimal_demand_instruments),
            ncol(oi$optimal_demand_instruments)))
#> Optimal demand instruments: 2256 x 4

# Re-estimate with optimal instruments
prob_opt <- oi$to_problem()
res_opt <- prob_opt$solve(method = "2s")

# Compare SEs
se_standard <- results$summary_table()$se
se_optimal <- res_opt$summary_table()$se
comp <- data.frame(
  parameter = results$summary_table()$parameter,
  SE_standard = round(se_standard, 4),
  SE_optimal = round(se_optimal, 4),
  ratio = round(se_optimal / se_standard, 4)
)
print(comp, row.names = FALSE)
#>    parameter SE_standard SE_optimal  ratio
#>  (Intercept)      0.1073     0.1094 1.0197
#>       prices      0.8453     0.8584 1.0155
#>        sugar      0.0041     0.0042 1.0147
#>        mushy      0.0514     0.0527 1.0264
```

## Parametric Bootstrap

The parametric bootstrap provides standard errors for post-estimation
quantities (elasticities, CS, merger effects) by resampling the
structural errors and re-estimating:

``` r
# Use a simple simulation for speed
id_data <- build_id_data(T = 10, J = 8, F = 2)
set.seed(42)
id_data$x <- runif(nrow(id_data))

sim <- blp_simulation(list(blp_formulation(~ prices + x)), id_data,
                      beta = c(0.5, -2.0, 0.8), xi_variance = 0.2, seed = 42)
sim_res <- sim$replace_endogenous()
sim_prob <- sim_res$to_problem()
est <- sim_prob$solve(method = "1s")

# 10 bootstrap draws (use more in practice)
boot <- est$bootstrap(draws = 10, seed = 123, method = "1s")
cat(sprintf("Successful bootstrap draws: %d\n", length(boot)))
#> Successful bootstrap draws: 10

# Bootstrap distribution of price coefficient
boot_prices <- sapply(boot, function(b) b$beta[2])
cat(sprintf("Price coefficient:\n"))
#> Price coefficient:
cat(sprintf("  Point estimate: %.4f\n", est$beta[2]))
#>   Point estimate: 1.5240
cat(sprintf("  Bootstrap mean: %.4f\n", mean(boot_prices)))
#>   Bootstrap mean: 2.4238
cat(sprintf("  Bootstrap SD:   %.4f\n", sd(boot_prices)))
#>   Bootstrap SD:   2.8593
cat(sprintf("  Bootstrap 95%% CI: [%.4f, %.4f]\n",
            quantile(boot_prices, 0.025), quantile(boot_prices, 0.975)))
#>   Bootstrap 95% CI: [-1.1730, 7.5441]
```

## Importance Sampling

Importance sampling concentrates integration nodes in regions of high
probability mass, improving accuracy with fewer nodes. This is useful
for RC models where the mixing distribution has moved far from the
standard normal.

``` r
# RC model for importance sampling
sim_rc <- blp_simulation(
  list(blp_formulation(~ prices + x), blp_formulation(~ prices + x)),
  id_data, beta = c(0.5, -2.0, 0.8), sigma = diag(c(0.3, 0.3, 0.3)),
  integration = blp_integration("product", size = 3),
  xi_variance = 0.2, seed = 42
)
sim_rc_res <- sim_rc$replace_endogenous()
prob_rc <- sim_rc_res$to_problem()
est_rc <- prob_rc$solve(sigma = diag(c(0.3, 0.3, 0.3)), method = "1s",
  optimization = blp_optimization("l-bfgs-b", method_options = list(maxit = 50)))

is_result <- est_rc$importance_sampling(n_draws = 100, seed = 42)
cat(sprintf("Importance sampling:\n"))
#> Importance sampling:
cat(sprintf("  Draws: %d\n", nrow(is_result$nodes)))
#>   Draws: 100
cat(sprintf("  Effective sample size: %.1f (of %d)\n",
            is_result$effective_sample_size, nrow(is_result$nodes)))
#>   Effective sample size: 1.7 (of 100)

# Re-estimate with importance-sampled nodes
prob_is <- is_result$to_problem()
cat(sprintf("  Problem with IS: I = %d agents\n", prob_is$I))
#>   Problem with IS: I = 1000 agents
```

## Initial Update

The `initial_update` option evaluates the model at starting parameter
values and updates the weighting matrix before the first GMM step. This
is especially useful when the initial $W = (Z\prime Z/N)^{- 1}$ is a
poor approximation:

``` r
f1 <- blp_formulation(~ prices + sugar + mushy)
prob <- blp_problem(list(f1), load_nevo_products())

res_no_update <- prob$solve(method = "2s", initial_update = FALSE)
res_with_update <- prob$solve(method = "2s", initial_update = TRUE)

cat(sprintf("Without initial_update: obj = %.4f\n", res_no_update$objective))
#> Without initial_update: obj = 139.8994
cat(sprintf("With initial_update:    obj = %.4f\n", res_with_update$objective))
#> With initial_update:    obj = 139.9303
```

## Summary of Available Methods

| Method                                | Description                                 |
|---------------------------------------|---------------------------------------------|
| `compute_elasticities()`              | $J \times J$ price elasticity matrix        |
| `compute_aggregate_elasticities()`    | Scalar per-market aggregate elasticity      |
| `compute_diversion_ratios()`          | Short-run (fixed prices) diversion          |
| `compute_long_run_diversion_ratios()` | Long-run (equilibrium) diversion            |
| `compute_passthrough()`               | $J \times J$ cost-to-price passthrough      |
| `compute_costs()`                     | Recover marginal costs from Bertrand FOCs   |
| `compute_markups()`                   | Lerner index $(p - c)/p$                    |
| `compute_shares()`                    | Predicted shares at estimated parameters    |
| `compute_profits()`                   | Per-product profits $(p - c) \cdot s$       |
| `compute_consumer_surplus()`          | Small-Rosen CS per market                   |
| `compute_hhi()`                       | Herfindahl-Hirschman Index per market       |
| `compute_merger()`                    | Full merger simulation with new equilibrium |
| `compute_optimal_instruments()`       | Feasible efficient instruments              |
| `bootstrap()`                         | Parametric bootstrap for inference          |
| `importance_sampling()`               | Improved integration nodes                  |
| `run_hansen_test()`                   | Hansen J-test for overidentification        |
| `run_wald_test()`                     | Wald test for parameter restrictions        |
| `summary_table()`                     | Publication-ready parameter table           |

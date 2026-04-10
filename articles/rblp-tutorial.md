# Getting Started with rblp

## Overview

**rblp** is an R implementation of the BLP (Berry, Levinsohn, Pakes
1995) framework for estimating demand models for differentiated
products. It is translated from the
[pyblp](https://pyblp.readthedocs.io/) Python package by Conlon and
Gortmaker (2020).

The package supports:

- **Plain logit** and **nested logit** demand estimation
- **Random coefficients logit** (mixed logit) with agent-level
  heterogeneity
- **Demographics** via Pi matrix interactions
- **Product fixed effects** via the `absorb` parameter
- **GMM estimation** (1-step and 2-step)
- **Post-estimation**: elasticities, consumer surplus, HHI, diversion
  ratios
- **Simulation**: generate equilibrium data and test estimation
- **Merger simulation**: compute post-merger prices and welfare effects

## Quick Start

``` r
library(rblp)

# Load the Nevo (2000) cereal data
products <- load_nevo_products()
head(products[, c("market_ids", "product_ids", "shares", "prices", "sugar", "mushy")])
#>   market_ids product_ids      shares     prices sugar mushy
#> 1      C01Q1       F1B04 0.012417212 0.07208794     2     1
#> 2      C01Q1       F1B06 0.007809387 0.11417849    18     1
#> 3      C01Q1       F1B07 0.012994511 0.13239066     4     1
#> 4      C01Q1       F1B09 0.005769961 0.13034408     3     0
#> 5      C01Q1       F1B11 0.017934141 0.15482331    12     0
#> 6      C01Q1       F1B13 0.026601892 0.13704921    14     0
```

The data has 2256 product-market observations across 94 markets. The
`demand_instruments0` through `demand_instruments19` columns contain 20
pre-computed BLP-style excluded instruments.

## Example 1: Plain Logit

The simplest demand model. We regress the log odds ratio
$\delta_{jt} = \log\left( s_{jt} \right) - \log\left( s_{0t} \right)$ on
product characteristics using IV-GMM.

``` r
# Define the linear model: intercept + prices + sugar + mushy
f1 <- blp_formulation(~ prices + sugar + mushy)

# Create and solve the problem
logit_problem <- blp_problem(list(f1), products)
logit_results <- logit_problem$solve(method = "1s")

print(logit_results)
#> BLP Estimation Results
#>   Method: 1S GMM
#>   Objective: 2.047932e+02
#>   Optimization converged: TRUE
#>   FP converged: TRUE (94 total iterations)
#> 
#> Parameter Estimates:
#>  parameter   estimate   se       t_stat 
#>  (Intercept) -2.810674  0.109432 -25.684
#>  prices      -11.699731 0.858444 -13.629
#>  sugar       0.048381   0.004208 11.498 
#>  mushy       0.043127   0.052717 0.818
```

The price coefficient is negative as expected – higher prices reduce
utility. The exogenous characteristics (sugar, mushy) serve as
instruments for prices.

## Example 2: Logit with Product Fixed Effects

This matches the pyblp Nevo tutorial specification. Product fixed
effects absorb all time-invariant product characteristics (brand, sugar,
mushy, etc.), so we only estimate the price coefficient.

``` r
# Absorb product-level fixed effects (Frisch-Waugh-Lovell demeaning)
f1_fe <- blp_formulation(~ prices, absorb = ~ product_ids)

fe_problem <- blp_problem(list(f1_fe), products)
fe_results <- fe_problem$solve(method = "1s")

print(fe_results)
#> BLP Estimation Results
#>   Method: 1S GMM
#>   Objective: 1.797148e+02
#>   Optimization converged: TRUE
#>   FP converged: TRUE (664 total iterations)
#> 
#> Parameter Estimates:
#>  parameter estimate   se       t_stat 
#>  prices    -30.420493 1.030311 -29.526
```

The price coefficient of -30.42 matches pyblp’s reported value of
approximately -30.

## Example 3: Random Coefficients Logit

The random coefficients (RC) logit allows consumer heterogeneity in
preferences. Each consumer $i$ has individual-specific taste
coefficients:

$$u_{ijt} = \delta_{jt} + \sum\limits_{k}x_{jk}\sigma_{k}\nu_{ik} + \varepsilon_{ijt}$$

where $\nu_{ik}$ are standard normal draws and $\sigma_{k}$ captures the
dispersion of tastes.

``` r
# X1: linear demand characteristics
f1_rc <- blp_formulation(~ prices + sugar + mushy)
# X2: random coefficient characteristics (same as X1 here)
f2_rc <- blp_formulation(~ prices + sugar + mushy)

rc_problem <- blp_problem(
  product_formulations = list(f1_rc, f2_rc),
  product_data = products,
  integration = blp_integration("product", size = 3)  # Gauss-Hermite quadrature
)

# Sigma: diagonal matrix with starting values (0 = fixed, non-zero = free)
initial_sigma <- diag(c(0.5, 0.5, 0.5, 0.5))

rc_results <- rc_problem$solve(
  sigma = initial_sigma,
  method = "1s",
  optimization = blp_optimization("l-bfgs-b",
    method_options = list(maxit = 200, factr = 1e7))
)

print(rc_results)
```

**Key parameters:**

- `integration = blp_integration("product", size = 3)` uses
  Gauss-Hermite product rule with 3 nodes per dimension ($3^{4} = 81$
  integration points per market)
- `sigma` is the K2 x K2 Cholesky root. Zeros on the diagonal fix that
  parameter; non-zero values are starting values for optimization
- `method = "2s"` runs two-step efficient GMM (default); `"1s"` runs
  one-step

## Example 4: RC Logit with Demographics

Demographics allow observed consumer characteristics to interact with
product characteristics via the Pi matrix:

$$u_{ijt} = \delta_{jt} + \sum\limits_{k}x_{jk}\left( \sigma_{k}\nu_{ik} + \sum\limits_{d}\pi_{kd}D_{id} \right) + \varepsilon_{ijt}$$

``` r
# Load agent (consumer) data with demographics
agents <- load_nevo_agents()
head(agents[, c("market_ids", "weights", "income", "age", "child")])

# Demographics formulation (no intercept)
demo_form <- blp_formulation(~ 0 + income + income_squared + age + child)

# pyblp specification: X1 = '0 + prices' with product FE absorbed
# X2 = '1 + prices + sugar + mushy' (with intercept)
demo_problem <- blp_problem(
  product_formulations = list(
    blp_formulation(~ 0 + prices, absorb = ~ product_ids),  # X1
    blp_formulation(~ prices + sugar + mushy)                # X2
  ),
  product_data = products,
  agent_formulation = demo_form,
  agent_data = agents
)

# Sigma (4x4 diagonal) and Pi (4x4, K2 x D)
# Starting values near pyblp's converged estimates
sigma0 <- diag(c(0.558, 3.312, -0.006, 0.093))
pi0 <- matrix(c(
  2.292, 0, 1.284, 0,
  588.3, -30.19, 0, 11.05,
  -0.384, 0, 0.0524, 0,
  0.748, 0, -1.354, 0
), nrow = 4, ncol = 4, byrow = TRUE)

demo_results <- demo_problem$solve(
  sigma = sigma0, pi = pi0,
  method = "1s",
  optimization = blp_optimization("l-bfgs-b",
    method_options = list(maxit = 1000, factr = 1e7))
)

print(demo_results)
# Price coefficient ~ -63, matching pyblp's Nevo tutorial
```

## Post-Estimation

After estimation, compute policy-relevant quantities:

### Elasticities

``` r
# Compute elasticities for the first market
first_market <- logit_results$problem$unique_market_ids[1]
E <- logit_results$compute_elasticities(first_market)

# Own-price elasticities (diagonal)
cat("Own-price elasticities (first 5 products):\n")
#> Own-price elasticities (first 5 products):
print(round(diag(E)[1:5], 3))
#> [1] -0.833 -1.325 -1.529 -1.516 -1.779
```

### Consumer Surplus

``` r
cs <- logit_results$compute_consumer_surplus()
cat("Consumer surplus (first 5 markets):\n")
#> Consumer surplus (first 5 markets):
print(round(cs[1:5], 4))
#>  C01Q1  C03Q1  C04Q1  C05Q1  C07Q1 
#> 0.0503 0.0458 0.0944 0.0472 0.0772
```

### HHI

``` r
hhi <- logit_results$compute_hhi()
cat("HHI (first 5 markets):\n")
#> HHI (first 5 markets):
print(round(hhi[1:5], 0))
#> C01Q1 C03Q1 C04Q1 C05Q1 C07Q1 
#>   711   643  1435   566  1048
```

## Simulation

Generate synthetic equilibrium data to test the estimator:

``` r
# Create a balanced panel: 50 markets, 20 products/market, 4 firms
id_data <- build_id_data(T = 50, J = 20, F = 4)
set.seed(42)
id_data$x <- runif(nrow(id_data), 0, 1)

# True parameters: beta = (intercept=0.5, price=-2, x=0.8)
sim <- blp_simulation(
  product_formulations = list(blp_formulation(~ prices + x)),
  product_data = id_data,
  beta = c(0.5, -2, 0.8),
  xi_variance = 0.3,
  seed = 42
)

# Solve for equilibrium prices and shares
sim_results <- sim$replace_endogenous(
  iteration = blp_iteration("simple", list(atol = 1e-12, max_evaluations = 5000))
)

cat(sprintf("Price range: [%.3f, %.3f]\n", min(sim_results$prices), max(sim_results$prices)))
#> Price range: [1.550, 1.715]
cat(sprintf("Share range: [%.5f, %.5f]\n", min(sim_results$shares), max(sim_results$shares)))
#> Share range: [0.00351, 0.15537]
```

Now estimate back the true parameters:

``` r
sim_problem <- sim_results$to_problem()
sim_est <- sim_problem$solve(method = "1s")

cat("True vs Estimated:\n")
#> True vs Estimated:
cat(sprintf("  intercept: true=0.50, est=%.4f\n", sim_est$beta[1]))
#>   intercept: true=0.50, est=-1.1803
cat(sprintf("  price:     true=-2.00, est=%.4f\n", sim_est$beta[2]))
#>   price:     true=-2.00, est=-0.9726
cat(sprintf("  x:         true=0.80, est=%.4f\n", sim_est$beta[3]))
#>   x:         true=0.80, est=0.8334
```

## Merger Simulation

Simulate the effect of a merger between firms 2 and 3:

``` r
new_firm_ids <- sim_results$product_data$firm_ids
new_firm_ids[new_firm_ids == 3] <- 2  # merge firm 3 into firm 2

merger <- sim_est$compute_merger(
  new_firm_ids = new_firm_ids,
  iteration = blp_iteration("simple", list(atol = 1e-12, max_evaluations = 5000))
)

cat(sprintf("Mean price change: %.2f%%\n", mean(merger$price_change_pct)))
#> Mean price change: 7.61%
cat(sprintf("Max price change:  %.2f%%\n", max(merger$price_change_pct)))
#> Max price change:  25.48%
cat(sprintf("Mean CS change:    %.4f\n", mean(merger$delta_cs)))
#> Mean CS change:    -0.0786
```

## Key Classes and Functions

| Function                                                                                                                | Purpose                                                 |
|-------------------------------------------------------------------------------------------------------------------------|---------------------------------------------------------|
| [`blp_formulation()`](https://hchulkim.github.io/rblp/reference/blp_formulation.md)                                     | Define model formulas (with optional `absorb` for FE)   |
| [`blp_problem()`](https://hchulkim.github.io/rblp/reference/blp_problem.md)                                             | Create an estimation problem from data and formulations |
| [`blp_simulation()`](https://hchulkim.github.io/rblp/reference/blp_simulation.md)                                       | Create a simulation with known parameters               |
| [`blp_integration()`](https://hchulkim.github.io/rblp/reference/blp_integration.md)                                     | Configure integration (product rule, Monte Carlo, etc.) |
| [`blp_iteration()`](https://hchulkim.github.io/rblp/reference/blp_iteration.md)                                         | Configure fixed-point iteration (SQUAREM, simple, etc.) |
| [`blp_optimization()`](https://hchulkim.github.io/rblp/reference/blp_optimization.md)                                   | Configure optimization (L-BFGS-B, Nelder-Mead, etc.)    |
| [`build_blp_instruments()`](https://hchulkim.github.io/rblp/reference/build_blp_instruments.md)                         | Construct BLP-style instruments                         |
| [`build_differentiation_instruments()`](https://hchulkim.github.io/rblp/reference/build_differentiation_instruments.md) | Gandhi-Houde differentiation instruments                |
| [`build_id_data()`](https://hchulkim.github.io/rblp/reference/build_id_data.md)                                         | Create balanced panel identifiers                       |

## Comparison with pyblp

rblp aims to replicate pyblp’s functionality in R. Key correspondence:

| pyblp                                            | rblp                                                |
|--------------------------------------------------|-----------------------------------------------------|
| `Formulation('prices', absorb='C(product_ids)')` | `blp_formulation(~ prices, absorb = ~ product_ids)` |
| `Problem(formulation, product_data)`             | `blp_problem(list(f1), product_data)`               |
| `problem.solve()`                                | `problem$solve()`                                   |
| `results.compute_elasticities()`                 | `results$compute_elasticities()`                    |
| `Integration('product', size=5)`                 | `blp_integration("product", size = 5)`              |
| `Simulation(...)`                                | `blp_simulation(...)`                               |

## Further Reading

The following vignettes cover specific topics in greater depth:

- **[`vignette("logit-nested-logit")`](https://hchulkim.github.io/rblp/articles/logit-nested-logit.md)**:
  Derivation and estimation of the logit and nested logit demand models,
  including the logit inversion, price endogeneity, and the role of
  instruments.
- **[`vignette("random-coefficients")`](https://hchulkim.github.io/rblp/articles/random-coefficients.md)**:
  The random coefficients (mixed) logit model, integration methods, the
  sigma/pi parameterization, and convergence diagnostics.
- **[`vignette("post-estimation-mergers")`](https://hchulkim.github.io/rblp/articles/post-estimation-mergers.md)**:
  Post-estimation analysis including elasticities, diversion ratios,
  consumer surplus, markups, merger simulation, and specification tests.
- **[`vignette("mixtape-exercises")`](https://hchulkim.github.io/rblp/articles/mixtape-exercises.md)**:
  Exercises inspired by the Mixtape Sessions demand estimation course,
  progressing from OLS logit through IV logit to random coefficients
  logit.
- **[`vignette("simulation-monte-carlo")`](https://hchulkim.github.io/rblp/articles/simulation-monte-carlo.md)**:
  Using the BLPSimulation class for roundtrip validation, instrument
  construction, and Monte Carlo study design.

## References

- Berry, S., Levinsohn, J., & Pakes, A. (1995). Automobile Prices in
  Market Equilibrium. *Econometrica*, 63(4), 841-890.
- Nevo, A. (2000). A Practitioner’s Guide to Estimation of
  Random-Coefficients Logit Models of Demand. *Journal of Economics &
  Management Strategy*, 9(4), 513-548.
- Conlon, C. & Gortmaker, J. (2020). Best Practices for Differentiated
  Products Demand Estimation with pyblp. *RAND Journal of Economics*,
  51(4), 1108-1161.

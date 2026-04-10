# rblp

**BLP Demand Estimation for Differentiated Products in R**

An R implementation of the BLP (Berry, Levinsohn, Pakes 1995) framework
for estimating demand models for differentiated products, translated
from the [pyblp](https://pyblp.readthedocs.io/) Python package by Conlon
and Gortmaker (2020).

------------------------------------------------------------------------

## Features

- **Plain logit** and **nested logit** demand estimation
- **Random coefficients logit** (mixed logit) with consumer
  heterogeneity
- **Demographics** via Pi matrix interactions
- **Product fixed effects** via the `absorb` parameter
  (Frisch-Waugh-Lovell)
- **GMM estimation** (1-step and 2-step)
- **Integration**: Gauss-Hermite product rule, Monte Carlo, Halton
  sequences
- **Post-estimation**: elasticities, consumer surplus, HHI, diversion
  ratios
- **Simulation**: generate equilibrium data with known parameters
- **Merger simulation**: compute post-merger prices and welfare effects
- **BLP and differentiation instruments**:
  [`build_blp_instruments()`](https://hchulkim.github.io/rblp/reference/build_blp_instruments.md),
  [`build_differentiation_instruments()`](https://hchulkim.github.io/rblp/reference/build_differentiation_instruments.md)

------------------------------------------------------------------------

## Installation

``` r
# Install from GitHub
# install.packages("devtools")
devtools::install_github("hchulkim/rblp")
```

------------------------------------------------------------------------

## Quick Start

``` r
library(rblp)

# Load Nevo (2000) cereal data
products <- load_nevo_products()

# Plain logit: delta_jt = beta0 + alpha*p_jt + beta*X_j + xi_jt
f1 <- blp_formulation(~ prices + sugar + mushy)
problem <- blp_problem(list(f1), products)
results <- problem$solve(method = "1s")
print(results)
```

### Logit with Product Fixed Effects

``` r
# Absorb product-level FE (matches pyblp Nevo tutorial, alpha ~ -30)
f1_fe <- blp_formulation(~ prices, absorb = ~ product_ids)
fe_problem <- blp_problem(list(f1_fe), products)
fe_results <- fe_problem$solve(method = "1s")
print(fe_results)
```

### Random Coefficients Logit

``` r
f1_rc <- blp_formulation(~ prices + sugar + mushy)
f2_rc <- blp_formulation(~ prices + sugar + mushy)

rc_problem <- blp_problem(
  product_formulations = list(f1_rc, f2_rc),
  product_data = products,
  integration = blp_integration("product", size = 3)
)

rc_results <- rc_problem$solve(
  sigma = diag(c(0.5, 0.5, 0.5, 0.5)),
  method = "1s",
  optimization = blp_optimization("l-bfgs-b")
)
print(rc_results)
```

### RC Logit with Demographics

``` r
agents <- load_nevo_agents()
demo_form <- blp_formulation(~ 0 + income + income_squared + age + child)

# pyblp specification: prices with product FE in X1, full X2
demo_problem <- blp_problem(
  product_formulations = list(
    blp_formulation(~ 0 + prices, absorb = ~ product_ids),
    blp_formulation(~ prices + sugar + mushy)
  ),
  product_data = products,
  agent_formulation = demo_form,
  agent_data = agents
)

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
# Price coefficient ~ -63, matching pyblp
```

------------------------------------------------------------------------

## Post-Estimation

``` r
# Own-price elasticities
E <- results$compute_elasticities(market_id = "C01Q1")
diag(E)

# Consumer surplus
cs <- results$compute_consumer_surplus()

# HHI
hhi <- results$compute_hhi()

# Merger simulation
new_firm_ids <- products$firm_ids
new_firm_ids[new_firm_ids == 3] <- 2
merger <- results$compute_merger(new_firm_ids = new_firm_ids)
```

------------------------------------------------------------------------

## Simulation

``` r
# Create synthetic data with known parameters
id_data <- build_id_data(T = 50, J = 20, F = 4)
id_data$x <- runif(nrow(id_data))

sim <- blp_simulation(
  product_formulations = list(blp_formulation(~ prices + x)),
  product_data = id_data,
  beta = c(0.5, -2, 0.8),
  xi_variance = 0.3,
  seed = 42
)

# Solve for equilibrium
sim_data <- sim$replace_endogenous()

# Estimate back
sim_problem <- sim_data$to_problem()
est <- sim_problem$solve(method = "1s")
```

------------------------------------------------------------------------

## Comparison with pyblp

| pyblp                                            | rblp                                                |
|--------------------------------------------------|-----------------------------------------------------|
| `Formulation('prices', absorb='C(product_ids)')` | `blp_formulation(~ prices, absorb = ~ product_ids)` |
| `Problem(formulation, product_data)`             | `blp_problem(list(f1), product_data)`               |
| `problem.solve()`                                | `problem$solve()`                                   |
| `results.compute_elasticities()`                 | `results$compute_elasticities()`                    |
| `Integration('product', size=5)`                 | `blp_integration("product", size = 5)`              |
| `Simulation(...)`                                | `blp_simulation(...)`                               |

------------------------------------------------------------------------

## References

- Berry, S., Levinsohn, J., & Pakes, A. (1995). Automobile Prices in
  Market Equilibrium. *Econometrica*, 63(4), 841-890.
- Nevo, A. (2000). A Practitioner’s Guide to Estimation of
  Random-Coefficients Logit Models of Demand. *Journal of Economics &
  Management Strategy*, 9(4), 513-548.
- Conlon, C. & Gortmaker, J. (2020). Best Practices for Differentiated
  Products Demand Estimation with pyblp. *RAND Journal of Economics*,
  51(4), 1108-1161.
  [doi:10.1287/mksc.2023.1440](https://doi.org/10.1287/mksc.2023.1440)

------------------------------------------------------------------------

## License

GPL (\>= 2)

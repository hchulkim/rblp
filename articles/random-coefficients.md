# Random Coefficients Logit

## Introduction

The random coefficients (RC) logit – also known as the mixed logit or
BLP model – is the workhorse of modern demand estimation for
differentiated products. Unlike the plain logit, it generates flexible
substitution patterns that are not subject to the IIA property.

This vignette covers:

1.  The model and its parameters ($\sigma$, $\pi$)
2.  Integration methods for simulating the mixing distribution
3.  Estimation with and without demographics
4.  GMM methods and convergence diagnostics

For the plain logit and nested logit, see
[`vignette("logit-nested-logit")`](https://hchulkim.github.io/rblp/articles/logit-nested-logit.md).
For post-estimation analysis, see
[`vignette("post-estimation-mergers")`](https://hchulkim.github.io/rblp/articles/post-estimation-mergers.md).

## The Random Coefficients Model

### Utility Specification

Consumer $i$ in market $t$ receives utility from product $j$:

$$u_{ijt} = \delta_{jt} + \mu_{ijt} + \varepsilon_{ijt}$$

where:

- $\delta_{jt} = x_{jt}\prime\beta + \alpha p_{jt} + \xi_{jt}$ is the
  **mean utility** (common to all consumers)
- $\mu_{ijt} = \sum_{k}x_{jkt}\left( \sigma_{k}\nu_{ik} + \pi_{kd}D_{id} \right)$
  is the **individual deviation** from the mean
- $\varepsilon_{ijt}$ is a Type I EV idiosyncratic shock

The individual-specific component $\mu_{ijt}$ has two sources of
heterogeneity:

- **Unobserved heterogeneity** ($\sigma$): $\nu_{ik} \sim N(0,1)$ are
  standard normal draws, and $\sigma_{k}$ captures the dispersion of
  tastes for characteristic $k$
- **Observed heterogeneity** ($\pi$): $D_{id}$ are observed demographics
  (income, age, etc.), and $\pi_{kd}$ captures how demographic $d$
  shifts preferences for characteristic $k$

### The Sigma Matrix

The $\sigma$ parameter is a $K_{2} \times K_{2}$ lower-triangular matrix
(Cholesky root). In many applications, $\sigma$ is restricted to be
diagonal, meaning that unobserved taste shocks for different
characteristics are independent.

**Convention**: when specifying starting values for sigma:

- **Zero** on the diagonal means that parameter is **fixed at zero**
  (not estimated)
- **Non-zero** values are **starting values** for optimization

For example, with $K_{2} = 4$ characteristics (intercept, prices, sugar,
mushy):

``` r
# Diagonal sigma: all four RC variances are free
sigma <- diag(c(0.5, 3.0, 0.01, 0.1))

# Fix sugar and mushy RC at zero, estimate only intercept and price RC
sigma <- diag(c(0.5, 3.0, 0, 0))
```

### The Pi Matrix

When demographics are available, the $\pi$ matrix is $K_{2} \times D$,
where $D$ is the number of demographic variables. Element $\pi_{kd}$
captures how demographic $d$ shifts the taste coefficient for
characteristic $k$.

## Integration Methods

The RC logit requires numerical integration over the distribution of
unobserved heterogeneity $\nu$. The choice of integration method affects
both accuracy and computation time.

### Product Rule (Gauss-Hermite)

The product rule uses Gauss-Hermite quadrature nodes, taking the tensor
product across dimensions. With `size` nodes per dimension and $K_{2}$
random coefficients, this yields $\text{size}^{K_{2}}$ integration
points per market.

``` r
# 5 nodes per dimension, K2=4 -> 5^4 = 625 points per market
int_product <- blp_integration("product", size = 5)

# 3 nodes per dimension (faster, less accurate)
int_small <- blp_integration("product", size = 3)  # 3^4 = 81 points
```

The product rule is **deterministic** (no simulation noise) and provides
exact integration for polynomials up to degree
$2 \times \text{size} - 1$. However, it suffers from the curse of
dimensionality: the number of points grows exponentially with $K_{2}$.

### Monte Carlo

Standard Monte Carlo draws random points from the mixing distribution.
With `size` draws:

``` r
# 1000 Monte Carlo draws per market
int_mc <- blp_integration("monte_carlo", size = 1000, seed = 42)
```

Monte Carlo integration is flexible but introduces simulation noise that
can cause roughness in the GMM objective function, making optimization
difficult.

### Halton Sequences

Halton sequences are quasi-random low-discrepancy sequences that provide
more uniform coverage than pseudo-random draws:

``` r
# 100 Halton draws per market (requires randtoolbox package)
int_halton <- blp_integration("halton", size = 100, seed = 42)
```

### Recommendation

For most applications:

- **Small $K_{2}$ (2–4)**: product rule with `size = 5` or `size = 7`
- **Moderate $K_{2}$ (5–8)**: Halton sequences with 100–500 points
- **Large $K_{2}$ (9+)**: Monte Carlo with 1000+ draws

## RC Logit without Demographics

The simplest RC logit has only unobserved heterogeneity (no
demographics). This requires two formulations: one for the linear
component ($X_{1}$) and one for the nonlinear component ($X_{2}$).

### Problem Setup

``` r
library(rblp)
products <- load_nevo_products()

# X1: linear demand (intercept + prices + sugar + mushy)
f1 <- blp_formulation(~ prices + sugar + mushy)

# X2: characteristics with random coefficients (same as X1 here)
f2 <- blp_formulation(~ prices + sugar + mushy)

# Create the problem with Gauss-Hermite integration
rc_problem <- blp_problem(
  product_formulations = list(f1, f2),
  product_data = products,
  integration = blp_integration("product", size = 3)  # 3^4 = 81 points
)

cat(sprintf("K1 (linear params): %d\n", rc_problem$K1))
#> K1 (linear params): 4
cat(sprintf("K2 (random coefficients): %d\n", rc_problem$K2))
#> K2 (random coefficients): 4
cat(sprintf("Total integration nodes: %d (%d per market)\n",
            rc_problem$I, rc_problem$I / rc_problem$T))
#> Total integration nodes: 7614 (81 per market)
cat(sprintf("Total observations: %d\n", rc_problem$N))
#> Total observations: 2256
cat(sprintf("Markets: %d\n", rc_problem$T))
#> Markets: 94
```

The problem has 4 random coefficients (intercept, prices, sugar, mushy),
so with `size = 3` Gauss-Hermite nodes per dimension, we get
$3^{4} = 81$ integration points per market.

### Estimation

``` r
# Sigma starting values (4x4 diagonal)
initial_sigma <- diag(c(0.5, 0.5, 0.5, 0.5))

# Solve with one-step GMM
rc_results <- rc_problem$solve(
  sigma = initial_sigma,
  method = "1s",
  optimization = blp_optimization("l-bfgs-b",
    method_options = list(maxit = 200, factr = 1e7))
)

print(rc_results)
```

**Interpretation**:

- The estimated ${\widehat{\sigma}}_{k}$ values capture the standard
  deviation of taste heterogeneity for each characteristic
- Large ${\widehat{\sigma}}_{\text{price}}$ means consumers differ
  substantially in price sensitivity
- These taste differences generate realistic substitution: consumers who
  are less price-sensitive substitute toward premium products

## RC Logit with Demographics (Full Nevo Specification)

The full Nevo (2000) specification includes observed demographics
(income, age, child presence) interacting with product characteristics.
This is the canonical BLP specification that matches the pyblp Nevo
tutorial.

### Problem Setup

``` r
agents <- load_nevo_agents()

cat("Agent data columns:\n")
#> Agent data columns:
cat(paste(names(agents), collapse = ", "), "\n\n")
#> market_ids, city_ids, quarter, weights, nodes0, nodes1, nodes2, nodes3, income, income_squared, age, child
cat(sprintf("Agents: %d (%d per market)\n", nrow(agents),
            nrow(agents) / length(unique(agents$market_ids))))
#> Agents: 1880 (20 per market)

# Inspect the agent data
head(agents[, c("market_ids", "weights", "income", "income_squared", "age", "child")])
#>   market_ids weights     income income_squared          age      child
#> 1      C01Q1    0.05  0.4951235       8.331304 -0.230109009 -0.2308511
#> 2      C01Q1    0.05  0.3787622       6.121865 -2.532694102  0.7691489
#> 3      C01Q1    0.05  0.1050146       1.030803 -0.006965458 -0.2308511
#> 4      C01Q1    0.05 -1.4854809     -25.583605 -0.827946010  0.7691489
#> 5      C01Q1    0.05 -0.3165969      -6.517009 -0.230109009 -0.2308511
#> 6      C01Q1    0.05 -0.3372762      -6.878070  0.851696161 -0.2308511
```

``` r
# X1: linear demand with product fixed effects (absorb product dummies)
# Only price is estimated; intercept, sugar, mushy are absorbed
f1 <- blp_formulation(~ 0 + prices, absorb = ~ product_ids)

# X2: nonlinear demand (intercept + prices + sugar + mushy get RC)
f2 <- blp_formulation(~ prices + sugar + mushy)

# Agent formulation: demographics that interact with X2
f_demo <- blp_formulation(~ 0 + income + income_squared + age + child)

# Create the problem (agent data provides nodes, weights, and demographics)
demo_problem <- blp_problem(
  product_formulations = list(f1, f2),
  product_data = products,
  agent_formulation = f_demo,
  agent_data = agents
)

cat(sprintf("K1 (linear params): %d (only prices, FE absorbed)\n", demo_problem$K1))
#> K1 (linear params): 1 (only prices, FE absorbed)
cat(sprintf("K2 (random coefficients): %d\n", demo_problem$K2))
#> K2 (random coefficients): 4
cat(sprintf("D (demographics): %d\n", demo_problem$D))
#> D (demographics): 4
cat(sprintf("Agents per market (from data): %d\n", demo_problem$I))
#> Agents per market (from data): 1880
cat(sprintf("Sigma matrix: %d x %d\n", demo_problem$K2, demo_problem$K2))
#> Sigma matrix: 4 x 4
cat(sprintf("Pi matrix: %d x %d\n", demo_problem$K2, demo_problem$D))
#> Pi matrix: 4 x 4
```

The sigma and pi matrices encode the parameter structure:

``` r
# Sigma (4x4 diagonal): RC standard deviations
# Starting values near pyblp's converged estimates
sigma0 <- diag(c(0.558, 3.312, -0.006, 0.093))

# Pi (4x4): interactions between X2 characteristics and demographics
# Rows = X2 vars (intercept, prices, sugar, mushy)
# Cols = demographics (income, income_squared, age, child)
pi0 <- matrix(c(
  2.292, 0, 1.284, 0,        # intercept interactions
  588.3, -30.19, 0, 11.05,   # price interactions
  -0.384, 0, 0.0524, 0,      # sugar interactions
  0.748, 0, -1.354, 0        # mushy interactions
), nrow = 4, ncol = 4, byrow = TRUE)

# Solve
demo_results <- demo_problem$solve(
  sigma = sigma0,
  pi = pi0,
  method = "1s",
  optimization = blp_optimization("l-bfgs-b",
    method_options = list(maxit = 1000, factr = 1e7))
)

print(demo_results)
# Price coefficient ~ -63, matching pyblp's Nevo tutorial
```

**Key points about the Nevo specification**:

- The agent data includes 20 agents per market with pre-computed nodes
  (columns `nodes0`–`nodes3`), integration weights, and demographic
  variables
- Product fixed effects in $X_{1}$ absorb brand-level unobserved
  quality, so only the price coefficient is estimated in the linear part
- The large number of moment conditions (20 excluded instruments +
  exogenous characteristics) helps identify the many nonlinear
  parameters
- The Pi matrix zeros encode economic restrictions (e.g.,
  `income_squared` only interacts with the price coefficient)

## Two-Step vs. One-Step GMM

### One-Step GMM (`method = "1s"`)

Uses the initial weighting matrix $W = (Z\prime Z/N)^{- 1}$ (2SLS
weight) throughout. This produces **consistent** estimates but is not
asymptotically efficient.

### Two-Step GMM (`method = "2s"`)

1.  **Step 1**: Estimate with the 2SLS weight
    $W_{1} = (Z\prime Z/N)^{- 1}$
2.  **Update**: Compute the optimal weight $W^{*} = {\widehat{S}}^{- 1}$
    from Step 1 residuals
3.  **Step 2**: Re-estimate with $W^{*}$

Two-step GMM is **efficient** (achieves the smallest asymptotic variance
among GMM estimators) but can exhibit finite-sample bias, especially
with many moment conditions.

``` r
# One-step (faster, consistent)
results_1s <- problem$solve(sigma = sigma0, method = "1s")

# Two-step (efficient, slower)
results_2s <- problem$solve(sigma = sigma0, method = "2s")
```

## Convergence Diagnostics

The RC logit objective function can be non-convex, so checking
convergence is important.

### Optimizer Convergence

``` r
# Check if the optimizer converged
results$optimization_converged

# Number of iterations and function evaluations
results$optimization_iterations
results$optimization_evaluations

# Final objective value
results$objective

# Gradient at solution (should be near zero)
max(abs(results$gradient))
```

### Fixed-Point Convergence

The BLP contraction mapping must converge in every market at every
evaluation:

``` r
# Did all fixed points converge?
results$fp_converged

# Total contraction iterations (across all markets and evaluations)
results$fp_iterations
```

### Best Practices

1.  **Try multiple starting values**: the GMM objective may have local
    minima
2.  **Start from the logit solution**: use logit estimates as a sanity
    check
3.  **Use tight contraction tolerances**: `blp_iteration("squarem")` is
    the default and is faster than the simple contraction
4.  **Increase integration accuracy**: more nodes reduce simulation
    error
5.  **Check the gradient norm**: it should be small at convergence

## Optimization Options

rblp supports several optimizers:

``` r
# L-BFGS-B (default): quasi-Newton with bounds, uses gradient
opt_lbfgsb <- blp_optimization("l-bfgs-b",
  method_options = list(maxit = 1000, factr = 1e7))

# Nelder-Mead: derivative-free simplex method
opt_nm <- blp_optimization("nelder-mead",
  method_options = list(maxit = 5000))

# nlminb: PORT optimization, good with bounds
opt_nlminb <- blp_optimization("nlminb",
  method_options = list(iter.max = 1000, eval.max = 2000))
```

For most applications, **L-BFGS-B** is recommended because it uses
analytical gradients (computed via the implicit function theorem) and
respects parameter bounds.

## Iteration Options

The BLP contraction mapping can be accelerated with SQUAREM:

``` r
# SQUAREM (default): accelerated fixed-point iteration
iter_squarem <- blp_iteration("squarem")

# Simple contraction: safer but slower
iter_simple <- blp_iteration("simple",
  list(atol = 1e-14, max_evaluations = 5000))
```

SQUAREM typically converges in 10–50% fewer iterations than the simple
contraction, with no loss of accuracy.

## Summary of Key Formulas

| Component            | Formula                                                                                                               |
|----------------------|-----------------------------------------------------------------------------------------------------------------------|
| Mean utility         | $\delta_{jt} = x_{jt}\prime\beta + \alpha p_{jt} + \xi_{jt}$                                                          |
| Individual deviation | $\mu_{ijt} = \sum_{k}x_{jkt}\left( \sigma_{k}\nu_{ik} + \pi_{kd}D_{id} \right)$                                       |
| Choice probability   | $s_{ijt} = \frac{\exp\left( \delta_{jt} + \mu_{ijt} \right)}{1 + \sum_{k}\exp\left( \delta_{kt} + \mu_{ikt} \right)}$ |
| Predicted share      | $s_{jt}(\delta,\theta) = \int s_{ijt}\, dF(\nu,D)$                                                                    |
| Contraction          | $\delta^{r + 1} = \delta^{r} + \log s^{\text{obs}} - \log s\left( \delta^{r},\theta \right)$                          |
| Moment condition     | $E\left\lbrack Z_{jt}\prime\xi_{jt}(\theta) \right\rbrack = 0$                                                        |
| GMM objective        | $Q(\theta) = g(\theta)\prime Wg(\theta)$, where $g = Z\prime\xi/N$                                                    |

## References

- Berry, S., Levinsohn, J., & Pakes, A. (1995). Automobile Prices in
  Market Equilibrium. *Econometrica*, 63(4), 841-890.
- Nevo, A. (2000). A Practitioner’s Guide to Estimation of
  Random-Coefficients Logit Models of Demand. *Journal of Economics &
  Management Strategy*, 9(4), 513-548.
- Conlon, C. & Gortmaker, J. (2020). Best Practices for Differentiated
  Products Demand Estimation with pyblp. *RAND Journal of Economics*,
  51(4), 1108-1161.
- Varadhan, R. & Roland, C. (2008). Simple and Globally Convergent
  Methods for Accelerating the Convergence of Any EM Algorithm.
  *Scandinavian Journal of Statistics*, 35(2), 335-353.

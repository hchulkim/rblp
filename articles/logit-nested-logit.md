# Logit and Nested Logit Models

## Introduction

This vignette covers the simplest demand models in the BLP framework:
the plain logit and the logit with product fixed effects. These models
serve as useful starting points before moving to random coefficients
specifications (see
[`vignette("random-coefficients")`](https://hchulkim.github.io/rblp/articles/random-coefficients.md)).

## The Logit Demand Model

### Utility Specification

In the logit demand model, consumer $i$ in market $t$ receives utility
from product $j$:

$$u_{ijt} = \delta_{jt} + \varepsilon_{ijt}$$

where $\delta_{jt} = x_{jt}\prime\beta + \alpha p_{jt} + \xi_{jt}$ is
the mean utility that is common across all consumers, and
$\varepsilon_{ijt}$ is an i.i.d. Type I Extreme Value (Gumbel) taste
shock. The term $\xi_{jt}$ is an unobserved (to the econometrician)
product-market quality.

The consumer chooses the outside option (good 0) if no inside good
provides higher utility than $u_{i0t} = \varepsilon_{i0t}$.

### The Logit Inversion

Under the Type I EV assumption, the market share of product $j$ takes
the well-known logit form:

$$s_{jt} = \frac{\exp\left( \delta_{jt} \right)}{1 + \sum\limits_{k = 1}^{J_{t}}\exp\left( \delta_{kt} \right)}$$

Berry (1994) showed that this can be analytically inverted to recover
mean utilities from observed shares:

$$\delta_{jt} = \log\left( s_{jt} \right) - \log\left( s_{0t} \right)$$

where $s_{0t} = 1 - \sum_{j}s_{jt}$ is the outside good share. This
inversion is exact and requires no iteration – a major computational
advantage of the plain logit.

### Why Instruments Are Needed

Substituting the linear specification for $\delta_{jt}$, the estimating
equation is:

$$\log\left( s_{jt} \right) - \log\left( s_{0t} \right) = x_{jt}\prime\beta + \alpha p_{jt} + \xi_{jt}$$

Prices $p_{jt}$ are set by firms that observe $\xi_{jt}$ (unobserved
quality), creating endogeneity:
$\text{Cov}\left( p_{jt},\xi_{jt} \right) \neq 0$. OLS estimation of
$\alpha$ is therefore biased (typically toward zero, understating price
sensitivity).

Instrumental variables (IV) are needed to consistently estimate the
price coefficient. The moment condition is:

$$E\left\lbrack Z_{jt}\prime\xi_{jt} \right\rbrack = 0$$

where $Z_{jt}$ are instruments that are correlated with prices but
uncorrelated with unobserved quality. Common choices include:

- **BLP instruments**: sums of characteristics of other own-firm and
  rival products
- **Differentiation instruments** (Gandhi & Houde 2020): measures of
  product isolation in characteristic space
- **Cost shifters**: input prices, exchange rates, etc.

The Nevo cereal dataset comes with 20 pre-computed excluded instruments
(`demand_instruments0` through `demand_instruments19`).

## Example: Plain Logit with Nevo Data

``` r
library(rblp)

# Load the Nevo (2000) cereal data
products <- load_nevo_products()

cat(sprintf("Observations: %d\n", nrow(products)))
#> Observations: 2256
cat(sprintf("Markets: %d\n", length(unique(products$market_ids))))
#> Markets: 94
cat(sprintf("Products per market: %d\n", length(unique(products$product_ids))))
#> Products per market: 24
```

The plain logit includes an intercept and observable product
characteristics in the linear index:

``` r
# Define the linear demand formulation: intercept + prices + sugar + mushy
f1 <- blp_formulation(~ prices + sugar + mushy)

# Create the problem (single formulation = logit)
logit_problem <- blp_problem(list(f1), products)

# Solve with one-step GMM
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

Key observations:

- The **price coefficient** ($\alpha$) is negative, as expected: higher
  prices reduce utility.
- The exogenous product characteristics (intercept, sugar, mushy) and
  the 20 pre-computed instruments serve as excluded instruments for the
  price variable.
- This is a IV-GMM estimator, not OLS, so we are addressing price
  endogeneity.

### Post-Estimation Checks

``` r
# Elasticities for the first market
first_market <- logit_results$problem$unique_market_ids[1]
E <- logit_results$compute_elasticities(first_market)

cat("Own-price elasticities (first 5 products):\n")
#> Own-price elasticities (first 5 products):
print(round(diag(E)[1:5], 3))
#> [1] -0.833 -1.325 -1.529 -1.516 -1.779

# Consumer surplus
cs <- logit_results$compute_consumer_surplus()
cat("\nConsumer surplus (first 5 markets):\n")
#> 
#> Consumer surplus (first 5 markets):
print(round(cs[1:5], 4))
#>  C01Q1  C03Q1  C04Q1  C05Q1  C07Q1 
#> 0.0503 0.0458 0.0944 0.0472 0.0772
```

A well-known limitation of the plain logit is the **Independence of
Irrelevant Alternatives (IIA)** property: the ratio of any two products’
market shares depends only on their own mean utilities, not on what
other products are available. This generates unrealistic substitution
patterns. Product fixed effects (below) and random coefficients
([`vignette("random-coefficients")`](https://hchulkim.github.io/rblp/articles/random-coefficients.md))
address this.

## Example: Logit with Product Fixed Effects

When panel data is available (the same products observed across multiple
markets/time periods), product fixed effects can absorb all
time-invariant product characteristics. This is the **within estimator**
applied to the logit share inversion.

The specification absorbs $\alpha_{j}$ (product dummies), so that the
only remaining coefficient is the price coefficient $\alpha$:

$$\log\left( s_{jt} \right) - \log\left( s_{0t} \right) = \alpha_{j} + \alpha p_{jt} + \xi_{jt}$$

In rblp, the `absorb` argument implements the Frisch-Waugh-Lovell (FWL)
theorem, demeaning all variables within the groups defined by the absorb
formula.

``` r
# Absorb product-level fixed effects
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

The price coefficient of -30.42 is substantially more negative than the
plain logit estimate. This is because the product fixed effects control
for time-invariant unobserved quality, reducing the upward bias from the
positive correlation between prices and quality.

This result matches the pyblp Nevo tutorial specification, where the
analogous code is:

``` python
# pyblp equivalent
formulation = pyblp.Formulation('0 + prices', absorb='C(product_ids)')
```

### Comparing Elasticities

``` r
first_market <- fe_results$problem$unique_market_ids[1]
E_fe <- fe_results$compute_elasticities(first_market)

cat("Own-price elasticities with FE (first 5 products):\n")
#> Own-price elasticities with FE (first 5 products):
print(round(diag(E_fe)[1:5], 3))
#> [1] -2.166 -3.446 -3.975 -3.942 -4.625
```

The fixed-effects estimates yield much larger own-price elasticities (in
absolute value), reflecting the more negative price coefficient. This is
a common finding: controlling for product quality reveals greater
consumer price sensitivity.

## Nested Logit

The nested logit relaxes IIA by grouping products into **nests**
(categories). Substitution within a nest is stronger than substitution
across nests. The utility specification is:

$$u_{ijt} = \delta_{jt} + \zeta_{igt} + (1 - \rho)\varepsilon_{ijt}$$

where $g$ indexes the nest that product $j$ belongs to, $\zeta_{igt}$ is
a nest-specific shock, and $\rho \in \lbrack 0,1)$ is the **nesting
parameter**. When $\rho = 0$, the nested logit reduces to the plain
logit; as $\left. \rho\rightarrow 1 \right.$, products within a nest
become perfect substitutes.

The share equation becomes:

$$\log\left( s_{jt} \right) - \log\left( s_{0t} \right) = x_{jt}\prime\beta + \alpha p_{jt} + \rho\log\left( s_{j|g,t} \right) + \xi_{jt}$$

where $s_{j|g,t}$ is the within-nest share of product $j$.

### Nested Logit with rblp

To estimate a nested logit in rblp, include a `nesting_ids` column in
the product data and pass the `rho` parameter to
[`solve()`](https://rdrr.io/r/base/solve.html):

``` r
# Suppose products have a 'nesting_ids' column (e.g., product category)
# products$nesting_ids <- products$brand_category  # hypothetical

# Same formulation as logit
f1 <- blp_formulation(~ prices + sugar + mushy)
nested_problem <- blp_problem(list(f1), products)

# Solve with nesting parameter (starting value)
nested_results <- nested_problem$solve(
  rho = 0.5,  # starting value for the nesting parameter
  method = "1s"
)

# The estimated rho captures within-nest correlation
print(nested_results)
```

The nested logit offers a middle ground between the plain logit (too
restrictive substitution) and the full random coefficients model
(computationally demanding). It is particularly useful when products
have a natural categorical structure (e.g., cereal brands by
manufacturer, cars by segment).

## Summary

| Model               | IIA                  | Computation                        | Parameters                         |
|---------------------|----------------------|------------------------------------|------------------------------------|
| Plain logit         | Yes                  | Analytical inversion, very fast    | $\beta$, $\alpha$                  |
| Logit + FE          | Yes                  | Analytical + demeaning, fast       | $\alpha$ (FE absorbed)             |
| Nested logit        | Relaxed within nests | Analytical + optimize $\rho$       | $\beta$, $\alpha$, $\rho$          |
| Random coefficients | No                   | Contraction mapping + optimization | $\beta$, $\alpha$, $\Sigma$, $\Pi$ |

For richer substitution patterns that do not rely on pre-specified
nesting structures, see
[`vignette("random-coefficients")`](https://hchulkim.github.io/rblp/articles/random-coefficients.md).

## References

- Berry, S. (1994). Estimating Discrete-Choice Models of Product
  Differentiation. *RAND Journal of Economics*, 25(2), 242-262.
- Berry, S., Levinsohn, J., & Pakes, A. (1995). Automobile Prices in
  Market Equilibrium. *Econometrica*, 63(4), 841-890.
- Nevo, A. (2000). A Practitioner’s Guide to Estimation of
  Random-Coefficients Logit Models of Demand. *Journal of Economics &
  Management Strategy*, 9(4), 513-548.
- Gandhi, A. & Houde, J.-F. (2020). Measuring Substitution Patterns in
  Differentiated-Products Industries. NBER Working Paper.

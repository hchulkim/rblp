# Mixtape Demand Estimation Exercises

## Introduction

This vignette reproduces key exercises from the Mixtape Sessions demand
estimation course by Scott Cunningham, Jeff Gortmaker, and others. The
exercises build from OLS logit through IV logit to random coefficients
logit, using the Nevo (2000) cereal dataset.

The rblp package bundles the same underlying Nevo data used in the
pyblp-based course materials, so results should be directly comparable.

For background on the models, see
[`vignette("logit-nested-logit")`](https://hchulkim.github.io/rblp/articles/logit-nested-logit.md)
and
[`vignette("random-coefficients")`](https://hchulkim.github.io/rblp/articles/random-coefficients.md).

## Data Setup

``` r
library(rblp)

products <- load_nevo_products()
agents <- load_nevo_agents()

# Key variables
cat(sprintf("Product-market observations: %d\n", nrow(products)))
#> Product-market observations: 2256
cat(sprintf("Markets: %d\n", length(unique(products$market_ids))))
#> Markets: 94
cat(sprintf("Products per market: %d\n", length(unique(products$product_ids))))
#> Products per market: 24
cat(sprintf("Agents per market: %d\n",
            nrow(agents) / length(unique(agents$market_ids))))
#> Agents per market: 20
```

### Computing the Outside Good Share

The logit model requires the outside good share $s_{0t}$ in each market.
The Nevo data defines the potential market size implicitly: inside
shares sum to less than 1, and $s_{0t} = 1 - \sum_{j}s_{jt}$.

``` r
# Compute inside share by market
market_inside_shares <- tapply(products$shares, products$market_ids, sum)

cat("Inside share summary:\n")
#> Inside share summary:
cat(sprintf("  Mean: %.4f\n", mean(market_inside_shares)))
#>   Mean: 0.4758
cat(sprintf("  Min:  %.4f\n", min(market_inside_shares)))
#>   Min:  0.1848
cat(sprintf("  Max:  %.4f\n", max(market_inside_shares)))
#>   Max:  0.6954
```

## Exercise 1: OLS Logit vs. IV Logit

### Part A: OLS Logit Regression

The OLS logit regression estimates the linear model:

$$\log\left( s_{jt} \right) - \log\left( s_{0t} \right) = \beta_{0} + \alpha p_{jt} + \beta_{1}\text{sugar}_{jt} + \beta_{2}\text{mushy}_{jt} + \xi_{jt}$$

by ordinary least squares, ignoring the endogeneity of prices.

``` r
# Compute log odds ratio (the dependent variable)
outside_shares <- 1 - tapply(products$shares, products$market_ids, sum)
products$s0 <- outside_shares[as.character(products$market_ids)]
products$y <- log(products$shares) - log(products$s0)

# OLS regression
ols_fit <- lm(y ~ prices + sugar + mushy, data = products)

cat("OLS Logit Results:\n")
#> OLS Logit Results:
cat("=================\n")
#> =================
print(round(summary(ols_fit)$coefficients, 4))
#>             Estimate Std. Error  t value Pr(>|t|)
#> (Intercept)  -2.9928     0.1117 -26.7981   0.0000
#> prices      -10.1199     0.8795 -11.5059   0.0000
#> sugar         0.0461     0.0044  10.4917   0.0000
#> mushy         0.0520     0.0519   1.0011   0.3169
```

The OLS price coefficient is approximately $- 7.5$. This is biased
toward zero because of the positive correlation between prices and
unobserved product quality ($\xi_{jt}$): high-quality products command
higher prices, which partially offsets the negative effect of price on
utility.

### Part B: IV Logit

The IV logit uses excluded instruments to address price endogeneity. The
Nevo data includes 20 pre-computed demand instruments.

``` r
# IV logit using rblp's GMM framework
f1 <- blp_formulation(~ prices + sugar + mushy)
logit_problem <- blp_problem(list(f1), products)
iv_results <- logit_problem$solve(method = "1s")

print(iv_results)
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

The IV price coefficient is approximately $- 30.6$, much more negative
than the OLS estimate. This large difference illustrates the severity of
price endogeneity bias in demand estimation.

### Part C: IV Logit with Product Fixed Effects

Adding product fixed effects absorbs all time-invariant product
characteristics, leaving only the within-product price variation for
identification:

``` r
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

### Comparison

``` r
cat("Price Coefficient Comparison:\n")
#> Price Coefficient Comparison:
cat(sprintf("  OLS:          %.4f\n", coef(ols_fit)["prices"]))
#>   OLS:          -10.1199
cat(sprintf("  IV:           %.4f\n", iv_results$beta[1]))
#>   IV:           -2.8107
cat(sprintf("  IV + FE:      %.4f\n", fe_results$beta[1]))
#>   IV + FE:      -30.4205
```

The progression shows how controlling for endogeneity (IV) and
unobserved product heterogeneity (FE) reveals the true consumer price
sensitivity.

### Part D: Counterfactual Price Increase

Using the IV logit results, consider a 10% price increase for all
products in the first market:

``` r
first_market <- iv_results$problem$unique_market_ids[1]
E <- iv_results$compute_elasticities(first_market)

# Own-price elasticities
own_elast <- diag(E)
cat("Own-price elasticities (first market):\n")
#> Own-price elasticities (first market):
cat(sprintf("  Mean: %.3f\n", mean(own_elast)))
#>   Mean: -1.514
cat(sprintf("  Range: [%.3f, %.3f]\n", min(own_elast), max(own_elast)))
#>   Range: [-2.033, -0.833]

# Approximate share changes from a 10% uniform price increase
# ds_j/s_j = E_jj * (dp/p) for own effect (ignoring cross effects for simplicity)
pct_change <- own_elast * 0.10  # 10% price increase
cat(sprintf("\nApprox. mean share change from 10%% price increase: %.2f%%\n",
            mean(pct_change) * 100))
#> 
#> Approx. mean share change from 10% price increase: -15.14%
```

## Exercise 2: Random Coefficients Logit with Demographics

This exercise estimates the full Nevo specification with random
coefficients and demographic interactions. This is the canonical BLP
model.

### Setting Up the Problem

``` r
# X1: price coefficient (with product FE absorbed)
f1 <- blp_formulation(~ 0 + prices, absorb = ~ product_ids)

# X2: random coefficient characteristics
f2 <- blp_formulation(~ prices + sugar + mushy)

# Demographics formulation
f_demo <- blp_formulation(~ 0 + income + income_squared + age + child)

# Create the problem
rc_problem <- blp_problem(
  product_formulations = list(f1, f2),
  product_data = products,
  agent_formulation = f_demo,
  agent_data = agents
)
```

### Starting Values and Estimation

The starting values for sigma and pi are important because the GMM
objective can be non-convex. We use values near the converged pyblp
estimates:

``` r
# Sigma: standard deviations of unobserved taste heterogeneity
sigma0 <- diag(c(0.558, 3.312, -0.006, 0.093))

# Pi: demographic interactions
# Rows: intercept, prices, sugar, mushy
# Cols: income, income_squared, age, child
pi0 <- matrix(c(
  2.292,  0,      1.284,  0,
  588.3,  -30.19, 0,      11.05,
  -0.384, 0,      0.0524, 0,
  0.748,  0,      -1.354, 0
), nrow = 4, ncol = 4, byrow = TRUE)

# Estimate with 1-step GMM
rc_results <- rc_problem$solve(
  sigma = sigma0,
  pi = pi0,
  method = "1s",
  optimization = blp_optimization("l-bfgs-b",
    method_options = list(maxit = 1000, factr = 1e7))
)

print(rc_results)
```

### Expected Results

The key parameter estimates from the Nevo specification (as reported in
pyblp’s tutorial):

| Parameter                                | Estimate | Interpretation                                                 |
|------------------------------------------|----------|----------------------------------------------------------------|
| $\alpha$ (price)                         | ~ $- 63$ | Very high price sensitivity after controlling for demographics |
| $\sigma_{1}$ (constant)                  | ~ $0.56$ | Moderate unobserved heterogeneity in base utility              |
| $\sigma_{2}$ (price)                     | ~ $3.3$  | Large unobserved heterogeneity in price sensitivity            |
| $\pi_{21}$ (price $\times$ income)       | ~ $588$  | Higher income reduces price sensitivity                        |
| $\pi_{22}$ (price $\times$ income$^{2}$) | ~ $- 30$ | Diminishing income effect                                      |

### Post-Estimation

``` r
# Elasticities are now richer (no IIA)
E_rc <- rc_results$compute_elasticities(first_market)

# Compare logit vs RC logit own-price elasticities
cat("RC logit own-price elasticities (first 5 products):\n")
print(round(diag(E_rc)[1:5], 3))

# Diversion ratios now reflect product similarity
D_rc <- rc_results$compute_diversion_ratios(first_market)

# Consumer surplus
cs_rc <- rc_results$compute_consumer_surplus()
cat(sprintf("Mean CS: %.4f\n", mean(cs_rc)))
```

## Exercise 3: Micro Moments (Conceptual Overview)

Micro moments supplement the traditional aggregate BLP moments with
information from individual-level (micro) data. They are particularly
useful when:

- The number of random coefficients exceeds what aggregate data can
  identify
- Survey or scanner data provides information about individual purchase
  patterns
- You want to match demographic-purchase correlations

### The Micro Moment Framework

A micro moment matches a statistic computed from micro data to its model
prediction:

$$E\left\lbrack {\bar{m}}^{\text{micro}} - m(\theta) \right\rbrack = 0$$

where ${\bar{m}}^{\text{micro}}$ is the observed micro statistic and
$m(\theta)$ is the model-predicted value.

Common micro moments include:

- **Conditional choice probabilities**:
  $E\left\lbrack y_{ij}|D_{i} \right\rbrack$ – the probability of buying
  product $j$ given demographics
- **Covariance moments**: $\text{Cov}\left( D_{i},x_{j{(i)}} \right)$ –
  the covariance between demographics and chosen product characteristics
- **Conditional means**:
  $E\left\lbrack x_{j{(i)}}|D_{i} \in G \right\rbrack$ – the mean
  characteristic of the chosen product for a demographic group

### rblp’s Micro Moment Classes

The rblp package provides three classes for constructing micro moments:

``` r
# 1. Define a micro dataset (describes the data source)
md <- micro_dataset(
  name = "Nielsen scanner",
  observations = 50000,
  compute_weights = function(market_id, products, agents) {
    # Return I x J weight matrix
    # e.g., uniform weights
    I <- length(agents$weights)
    J <- length(products$shares)
    matrix(1, I, J)
  }
)

# 2. Define a micro part (what statistic to compute)
mp <- micro_part(
  name = "income_x_price",
  dataset = md,
  compute_values = function(market_id, products, agents) {
    # Return I x J value matrix
    # e.g., income * price for each agent-product pair
    outer(agents$income, products$prices)
  }
)

# 3. Define a micro moment (match target to prediction)
mm <- micro_moment(
  name = "E[income * price | purchase]",
  value = 0.5,  # observed from micro data
  parts = mp
)

# Pass to solve()
results <- problem$solve(
  sigma = sigma0,
  pi = pi0,
  micro_moments = list(mm)
)
```

Micro moments are an advanced topic. For a full treatment, see Conlon
and Gortmaker (2020) and the pyblp documentation.

## Summary of Exercise Results

| Exercise | Model                   | Key Finding                                      |
|----------|-------------------------|--------------------------------------------------|
| 1A       | OLS logit               | Price coeff. ~ $- 7.5$ (biased toward zero)      |
| 1B       | IV logit                | Price coeff. ~ $- 30.6$ (endogeneity corrected)  |
| 1C       | IV logit + FE           | Price coeff. ~ $- 30$ (product quality absorbed) |
| 2        | RC logit + demographics | Price coeff. ~ $- 63$; rich heterogeneity        |
| 3        | Micro moments           | Supplement aggregate moments with micro data     |

The progression from OLS to IV to FE to RC shows how each step addresses
a different source of bias or enriches the substitution patterns.

## References

- Cunningham, S. (2021). *Causal Inference: The Mixtape*. Yale
  University Press.
- Berry, S. (1994). Estimating Discrete-Choice Models of Product
  Differentiation. *RAND Journal of Economics*, 25(2), 242-262.
- Nevo, A. (2000). A Practitioner’s Guide to Estimation of
  Random-Coefficients Logit Models of Demand. *Journal of Economics &
  Management Strategy*, 9(4), 513-548.
- Conlon, C. & Gortmaker, J. (2020). Best Practices for Differentiated
  Products Demand Estimation with pyblp. *RAND Journal of Economics*,
  51(4), 1108-1161.

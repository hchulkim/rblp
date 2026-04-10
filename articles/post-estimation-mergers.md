# Post-Estimation and Merger Simulation

## Introduction

After estimating a demand model, the structural parameters can be used
to compute policy-relevant quantities: elasticities, diversion ratios,
consumer surplus, markups, and merger simulations. This vignette
demonstrates these post-estimation tools.

We start with a logit model for fast, evaluable examples, then discuss
how results extend to the random coefficients case. For estimation
details, see
[`vignette("logit-nested-logit")`](https://hchulkim.github.io/rblp/articles/logit-nested-logit.md)
and
[`vignette("random-coefficients")`](https://hchulkim.github.io/rblp/articles/random-coefficients.md).

## Setup: Estimate a Logit Model

``` r
library(rblp)

products <- load_nevo_products()

# Logit with product fixed effects
f1 <- blp_formulation(~ prices, absorb = ~ product_ids)
problem <- blp_problem(list(f1), products)
results <- problem$solve(method = "1s")

cat(sprintf("Price coefficient: %.4f\n", results$beta[1]))
#> Price coefficient: -30.4205
```

## Elasticities

### Own-Price and Cross-Price Elasticities

The price elasticity matrix $E$ for market $t$ has entries:

$$E_{jk} = \frac{\partial s_{j}}{\partial p_{k}} \cdot \frac{p_{k}}{s_{j}}$$

- **Diagonal** entries ($E_{jj}$) are own-price elasticities: the
  percent change in $j$’s share for a 1% increase in $j$’s price. These
  should be negative.
- **Off-diagonal** entries ($E_{jk}$, $j \neq k$) are cross-price
  elasticities: the percent change in $j$’s share for a 1% increase in
  $k$’s price. In a substitutes market, these are positive.

``` r
# Compute elasticities for the first market
first_market <- results$problem$unique_market_ids[1]
E <- results$compute_elasticities(first_market)

cat("Elasticity matrix dimensions:", dim(E), "\n")
#> Elasticity matrix dimensions: 24 24

# Own-price elasticities (diagonal)
own_elast <- diag(E)
cat("\nOwn-price elasticities summary:\n")
#> 
#> Own-price elasticities summary:
cat(sprintf("  Mean: %.3f\n", mean(own_elast)))
#>   Mean: -3.935
cat(sprintf("  Min:  %.3f\n", min(own_elast)))
#>   Min:  -5.286
cat(sprintf("  Max:  %.3f\n", max(own_elast)))
#>   Max:  -2.166

# Cross-price elasticities (first row, first 5 columns)
cat("\nCross-price elasticities (product 1 w.r.t. products 1-5):\n")
#> 
#> Cross-price elasticities (product 1 w.r.t. products 1-5):
print(round(E[1, 1:min(5, ncol(E))], 4))
#> [1] -2.1657  0.0271  0.0523  0.0229  0.0845
```

### Logit vs. RC Logit Elasticities

In the plain logit, cross-price elasticities follow the IIA pattern: if
product $k$’s price increases, all other products gain share
proportionally to their current share. This means a luxury cereal and a
budget cereal respond identically to a third product’s price change,
which is unrealistic.

In the RC logit, cross-price elasticities depend on how similar products
are in the random coefficient characteristic space. Products that appeal
to similar consumer segments have higher cross-price elasticities,
generating realistic substitution patterns.

## Diversion Ratios

The diversion ratio $D_{jk}$ measures: if product $j$ is removed (or its
price increases), what fraction of $j$’s lost demand is captured by
product $k$?

$$D_{jk} = - \frac{\partial s_{k}/\partial p_{j}}{\partial s_{j}/\partial p_{j}}$$

Diversion ratios are central to merger analysis. The **upward pricing
pressure** (UPP) test for a merger between firms owning products $j$ and
$k$ is proportional to $D_{jk} \times m_{k}$, where $m_{k}$ is product
$k$’s margin.

``` r
# Compute diversion ratios for the first market
D <- results$compute_diversion_ratios(first_market)

cat("Diversion ratio matrix dimensions:", dim(D), "\n")
#> Diversion ratio matrix dimensions: 24 24

# Where does demand from product 1 go?
cat("\nDiversion from product 1 (first 5 products):\n")
#> 
#> Diversion from product 1 (first 5 products):
print(round(D[1, 1:min(5, ncol(D))], 4))
#> [1] 0.0000 0.0079 0.0132 0.0058 0.0182

# What fraction goes to the outside good?
cat(sprintf("\nDiversion to outside good from product 1: %.4f\n",
            1 - sum(D[1, ])))
#> 
#> Diversion to outside good from product 1: 0.5622
```

In the logit model, diversion is proportional to market shares (another
manifestation of IIA). In the RC logit, diversion reflects the actual
similarity of products to heterogeneous consumers.

## Consumer Surplus

Consumer surplus in the logit family follows the Small-Rosen (1981)
log-sum formula:

$$CS_{t} = \frac{1}{\alpha} \cdot E_{i}\left\lbrack \log\left( 1 + \sum\limits_{j = 1}^{J_{t}}\exp\left( \delta_{jt} + \mu_{ijt} \right) \right) \right\rbrack$$

where $\alpha$ is the (negative) price coefficient (marginal utility of
income). In the plain logit, $\mu_{ijt} = 0$ and the expectation is
trivial. In the RC logit, the expectation is taken over the simulated
consumer population.

``` r
cs <- results$compute_consumer_surplus()

cat("Consumer surplus per market (first 10 markets):\n")
#> Consumer surplus per market (first 10 markets):
print(round(cs[1:10], 4))
#>  C01Q1  C03Q1  C04Q1  C05Q1  C07Q1  C08Q1  C11Q1  C12Q1  C13Q1  C14Q1 
#> 0.0193 0.0176 0.0363 0.0182 0.0297 0.0180 0.0213 0.0204 0.0197 0.0249

cat(sprintf("\nMean CS across markets: %.4f\n", mean(cs)))
#> 
#> Mean CS across markets: 0.0220
cat(sprintf("Std. dev. of CS: %.4f\n", sd(cs)))
#> Std. dev. of CS: 0.0071
```

Changes in consumer surplus across policy scenarios (e.g., before and
after a merger) provide a monetary measure of consumer welfare impact.

## HHI (Herfindahl-Hirschman Index)

The HHI measures market concentration:

$$HHI_{t} = \sum\limits_{f}\left( \sum\limits_{j \in f}s_{jt} \right)^{2} \times 10000$$

where the inner sum aggregates shares within each firm $f$. The HHI
ranges from near zero (perfect competition) to 10,000 (monopoly). The
U.S. DOJ considers markets with HHI above 2,500 to be highly
concentrated.

``` r
hhi <- results$compute_hhi()

cat("HHI per market (first 10 markets):\n")
#> HHI per market (first 10 markets):
print(round(hhi[1:10], 0))
#> C01Q1 C03Q1 C04Q1 C05Q1 C07Q1 C08Q1 C11Q1 C12Q1 C13Q1 C14Q1 
#>   711   643  1435   566  1048   540   679   712   793   817

cat(sprintf("\nMean HHI: %.0f\n", mean(hhi)))
#> 
#> Mean HHI: 834
```

## Costs and Markups

### Recovering Marginal Costs

Given estimated demand, marginal costs can be backed out from the
Nash-Bertrand first-order conditions. In equilibrium, each firm sets
prices to satisfy:

$$p_{j} - mc_{j} = - \left( \Omega^{- 1}s \right)_{j}$$

where $\Omega$ is the ownership-weighted matrix of share derivatives.
Rearranging:

$$mc_{j} = p_{j} + \left( \Omega^{- 1}s \right)_{j}$$

``` r
costs <- results$compute_costs()

cat("Marginal costs (first 10 products):\n")
#> Marginal costs (first 10 products):
print(round(costs[1:10], 4))
#>  [1] 0.0348 0.0769 0.0951 0.0930 0.1175 0.0997 0.1069 0.0909 0.1123 0.0657

cat(sprintf("\nPrice range:  [%.4f, %.4f]\n",
            min(products$prices), max(products$prices)))
#> 
#> Price range:  [0.0455, 0.2257]
cat(sprintf("Cost range:   [%.4f, %.4f]\n", min(costs), max(costs)))
#> Cost range:   [-0.0001, 0.1878]
```

### Markups (Lerner Index)

The Lerner index measures the fraction of price that is markup:

$$L_{j} = \frac{p_{j} - mc_{j}}{p_{j}}$$

``` r
markups <- results$compute_markups()

cat("Markups (first 10 products):\n")
#> Markups (first 10 products):
print(round(markups[1:10], 4))
#>  [1] 0.5176 0.3268 0.2818 0.2862 0.2410 0.2722 0.2587 0.2910 0.2494 0.3941

cat(sprintf("\nMean markup: %.4f\n", mean(markups)))
#> 
#> Mean markup: 0.3292
cat(sprintf("Median markup: %.4f\n", median(markups)))
#> Median markup: 0.3116
```

## Merger Simulation

Merger simulation is the central counterfactual exercise in structural
IO. The idea: when two firms merge, the merged entity internalizes the
cross-price effects between the previously competing products, leading
to higher prices.

### Steps

1.  **Estimate demand** and recover marginal costs $mc$ under the
    pre-merger ownership structure
2.  **Update the ownership matrix** to reflect the merger (products of
    the acquired firm now belong to the acquirer)
3.  **Solve for new equilibrium prices** under the post-merger
    ownership, holding costs fixed
4.  **Compare** pre- and post-merger prices, shares, and consumer
    surplus

### Simulation Example

We use the simulation framework for a clean merger example:

``` r
# Create a balanced panel: 50 markets, 20 products, 4 firms
id_data <- build_id_data(T = 50, J = 20, F = 4)
set.seed(42)
id_data$x <- runif(nrow(id_data), 0, 1)

# Simulate with known parameters
sim <- blp_simulation(
  product_formulations = list(blp_formulation(~ prices + x)),
  product_data = id_data,
  beta = c(0.5, -2, 0.8),
  xi_variance = 0.3,
  seed = 42
)

sim_data <- sim$replace_endogenous(
  iteration = blp_iteration("simple", list(atol = 1e-12, max_evaluations = 5000))
)

# Estimate
sim_problem <- sim_data$to_problem()
sim_results <- sim_problem$solve(method = "1s")
```

Now simulate a merger between firms 2 and 3:

``` r
# Create post-merger firm IDs: merge firm 3 into firm 2
new_firm_ids <- sim_data$product_data$firm_ids
new_firm_ids[new_firm_ids == 3] <- 2

cat(sprintf("Pre-merger firms:  %s\n", paste(sort(unique(sim_data$product_data$firm_ids)), collapse = ", ")))
#> Pre-merger firms:  1, 2, 3, 4
cat(sprintf("Post-merger firms: %s\n", paste(sort(unique(new_firm_ids)), collapse = ", ")))
#> Post-merger firms: 1, 2, 4

# Run the merger simulation
merger <- sim_results$compute_merger(
  new_firm_ids = new_firm_ids,
  iteration = blp_iteration("simple", list(atol = 1e-12, max_evaluations = 5000))
)

cat(sprintf("\nMean price change:   %.2f%%\n", mean(merger$price_change_pct)))
#> 
#> Mean price change:   7.61%
cat(sprintf("Max price change:    %.2f%%\n", max(merger$price_change_pct)))
#> Max price change:    25.48%
cat(sprintf("Mean CS change:      %.4f\n", mean(merger$delta_cs)))
#> Mean CS change:      -0.0786
```

### Interpreting Merger Results

``` r
# Which products see the biggest price increases?
# Products of the merging firms should see larger increases
merging_products <- sim_data$product_data$firm_ids %in% c(2, 3)
non_merging <- !merging_products

cat("Price changes by group:\n")
#> Price changes by group:
cat(sprintf("  Merging firms (2 & 3):     mean = %.2f%%\n",
            mean(merger$price_change_pct[merging_products])))
#>   Merging firms (2 & 3):     mean = 14.22%
cat(sprintf("  Non-merging firms (1 & 4): mean = %.2f%%\n",
            mean(merger$price_change_pct[non_merging])))
#>   Non-merging firms (1 & 4): mean = 1.00%

# Consumer surplus change
cat(sprintf("\nTotal CS change: %.4f\n", sum(merger$delta_cs)))
#> 
#> Total CS change: -3.9287
cat(sprintf("CS change is %s\n",
            ifelse(sum(merger$delta_cs) < 0, "negative (consumers harmed)",
                   "positive (consumers benefit)")))
#> CS change is negative (consumers harmed)
```

### Pre-Merger vs. Post-Merger HHI

``` r
# Compute pre-merger and post-merger HHI
hhi_pre <- sim_results$compute_hhi()

cat(sprintf("Pre-merger mean HHI:  %.0f\n", mean(hhi_pre)))
#> Pre-merger mean HHI:  1250
cat(sprintf("Delta HHI estimate: merging two of four equal firms\n"))
#> Delta HHI estimate: merging two of four equal firms
```

## Specification Tests

### Hansen J-Test

The Hansen J-test checks the validity of the overidentifying
restrictions. Under the null that all instruments are valid, the
J-statistic is asymptotically chi-squared with degrees of freedom equal
to the number of moment conditions minus the number of parameters.

``` r
# Only meaningful for two-step GMM
results_2s <- problem$solve(method = "2s")
j_test <- results_2s$run_hansen_test()

cat(sprintf("J-statistic: %.4f\n", j_test$statistic))
cat(sprintf("Degrees of freedom: %d\n", j_test$df))
cat(sprintf("p-value: %.4f\n", j_test$p_value))
```

A small p-value suggests that at least some instruments may be invalid
(correlated with the structural error $\xi$).

### Wald Test

The Wald test can be used to test restrictions on the nonlinear
parameters (sigma, pi). For example, to test whether all random
coefficient variances are jointly zero (i.e., the logit is adequate):

``` r
# Test H0: all sigma parameters = 0
n_theta <- length(results$se)
R <- diag(n_theta)
r <- rep(0, n_theta)

wald <- results$run_wald_test(R, r)
cat(sprintf("Wald statistic: %.4f, p-value: %.4f\n",
            wald$statistic, wald$p_value))
```

## Summary of Post-Estimation Methods

| Method                                | What It Computes                     | Output              |
|---------------------------------------|--------------------------------------|---------------------|
| `compute_elasticities(market_id)`     | $J \times J$ price elasticity matrix | Matrix              |
| `compute_diversion_ratios(market_id)` | $J \times J$ diversion ratio matrix  | Matrix              |
| `compute_consumer_surplus()`          | Small-Rosen CS per market            | Named vector        |
| `compute_hhi()`                       | HHI per market                       | Named vector        |
| `compute_costs()`                     | Marginal costs from Bertrand FOC     | Vector (length $N$) |
| `compute_markups()`                   | Lerner index $(p - mc)/p$            | Vector (length $N$) |
| `compute_merger(new_firm_ids)`        | Post-merger equilibrium              | List                |
| `run_hansen_test()`                   | Overidentification J-test            | List                |
| `run_wald_test(R, r)`                 | Wald test on parameters              | List                |
| `summary_table()`                     | All parameter estimates with SEs     | Data frame          |

## References

- Small, K. & Rosen, H. (1981). Applied Welfare Economics with Discrete
  Choice Models. *Econometrica*, 49(1), 105-130.
- Werden, G. (1996). A Robust Test for Consumer Welfare Enhancing
  Mergers Among Sellers of Differentiated Products. *Journal of
  Industrial Economics*, 44(4), 409-413.
- Nevo, A. (2000). Mergers with Differentiated Products: The Case of the
  Ready-to-Eat Cereal Industry. *RAND Journal of Economics*, 31(3),
  395-421.
- Conlon, C. & Gortmaker, J. (2020). Best Practices for Differentiated
  Products Demand Estimation with pyblp. *RAND Journal of Economics*,
  51(4), 1108-1161.

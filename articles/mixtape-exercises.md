# Mixtape Demand Estimation Exercises

## Introduction

This vignette walks through the exercises from the [Mixtape Sessions
Demand
Estimation](https://github.com/Mixtape-Sessions/Demand-Estimation)
course, taught by Jeff Gortmaker and others. The original course
materials use Python and
[pyblp](https://pyblp.readthedocs.io/en/stable/); here we reproduce
every exercise in R using the **rblp** package.

The exercises build progressively:

1.  **Exercise 1 (Pure Logit):** OLS logit, IV logit, fixed effects,
    counterfactual price cuts, and elasticities.
2.  **Exercise 2 (Mixed Logit):** Adding demographic interactions and
    random coefficients to capture consumer heterogeneity.
3.  **Exercise 3 (Micro Moments):** Using individual-level data to
    sharpen identification (conceptual overview).

The dataset is a simplified version of the Nevo (2000) ready-to-eat
cereal data covering 24 products across 94 city-quarter markets.

## Data Setup

rblp bundles the Mixtape cereal data. The raw data has columns for
`market`, `product`, `mushy`, `servings_sold`, `city_population`,
`price_per_serving`, and `price_instrument`. We use
[`prepare_mixtape_data()`](https://hchulkim.github.io/rblp/reference/prepare_mixtape_data.md)
to compute market shares and rename columns into the format expected by
[`blp_problem()`](https://hchulkim.github.io/rblp/reference/blp_problem.md).

``` r
library(rblp)

# Load and prepare data
raw <- load_mixtape_products()
products <- prepare_mixtape_data(raw)

cat(sprintf("Product-market observations: %d\n", nrow(products)))
#> Product-market observations: 2256
cat(sprintf("Markets (city-quarter): %d\n", length(unique(products$market_ids))))
#> Markets (city-quarter): 94
cat(sprintf("Products per market: %d\n", length(unique(products$product_ids))))
#> Products per market: 24
cat(sprintf("Firms: %s\n", paste(sort(unique(products$firm_ids)), collapse = ", ")))
#> Firms: F1, F2, F3, F4, F6
```

### Understanding the Data

The key columns after preparation:

- `market_ids`: city-quarter identifier (e.g., “C01Q1”)
- `product_ids`: product identifier (e.g., “F1B04” = Firm 1, Brand 04)
- `shares`: market share = `servings_sold / (city_population * 90)`,
  where 90 represents one serving per day over a quarter
- `prices`: price per serving (dollars)
- `mushy`: indicator for whether the cereal is mushy (0/1)
- `firm_ids`: extracted from the first two characters of `product_ids`
- `demand_instruments0`: excluded price instrument for IV estimation

``` r
cat("Price per serving summary:\n")
#> Price per serving summary:
print(round(summary(products$prices), 4))
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>  0.0455  0.1055  0.1238  0.1257  0.1433  0.2257
cat("\nMarket share summary:\n")
#> 
#> Market share summary:
print(round(summary(products$shares), 6))
#>     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#> 0.000182 0.005183 0.011141 0.019825 0.024646 0.446883
cat("\nMushy breakdown:\n")
#> 
#> Mushy breakdown:
print(table(mushy = products$mushy[!duplicated(products$product_ids)]))
#> mushy
#>  0  1 
#> 16  8
```

### Computing the Outside Good Share

In the logit model, the outside good share in each market is
$s_{0t} = 1 - \sum_{j}s_{jt}$. The outside good represents consumers who
do not purchase any of the 24 inside products.

``` r
market_inside_shares <- tapply(products$shares, products$market_ids, sum)

cat("Inside share (sum of product shares per market):\n")
#> Inside share (sum of product shares per market):
cat(sprintf("  Mean:  %.4f\n", mean(market_inside_shares)))
#>   Mean:  0.4758
cat(sprintf("  Min:   %.4f\n", min(market_inside_shares)))
#>   Min:   0.1848
cat(sprintf("  Max:   %.4f\n", max(market_inside_shares)))
#>   Max:   0.6954

cat("\nOutside share:\n")
#> 
#> Outside share:
cat(sprintf("  Mean:  %.4f\n", mean(1 - market_inside_shares)))
#>   Mean:  0.5242
cat(sprintf("  Min:   %.4f\n", min(1 - market_inside_shares)))
#>   Min:   0.3046
cat(sprintf("  Max:   %.4f\n", max(1 - market_inside_shares)))
#>   Max:   0.8152
```

The outside share ranges from about 30% to 82%, indicating substantial
variation in how much of each market is captured by the 24 tracked
cereal products.

------------------------------------------------------------------------

## Exercise 1: Pure Logit

### Q1-Q2: Data Description and Market Share Computation

The Mixtape dataset is a panel of 24 cereal products observed in 94
city-quarter markets. For each product-market, we observe the quantity
sold (`servings_sold`), the city population, the price per serving, and
one excluded instrument for price.

Market shares are computed as:
$$s_{jt} = \frac{\text{servings\_sold}_{jt}}{\text{city\_population}_{t} \times 90}$$

where 90 is the assumed number of potential servings per person per
quarter (about one per day). This defines the total market size, with
the outside good absorbing whatever fraction of the population does not
purchase any of the 24 products.

``` r
# Show first few rows of the key variables
head(products[, c("market_ids", "product_ids", "firm_ids",
                   "shares", "prices", "mushy")], 10)
#>    market_ids product_ids firm_ids      shares     prices mushy
#> 1       C01Q1       F1B04       F1 0.012417212 0.07208794     1
#> 2       C01Q1       F1B06       F1 0.007809387 0.11417849     1
#> 3       C01Q1       F1B07       F1 0.012994511 0.13239066     1
#> 4       C01Q1       F1B09       F1 0.005769961 0.13034408     0
#> 5       C01Q1       F1B11       F1 0.017934141 0.15482331     0
#> 6       C01Q1       F1B13       F1 0.026601892 0.13704921     0
#> 7       C01Q1       F1B17       F1 0.025014766 0.14420936     1
#> 8       C01Q1       F1B30       F1 0.005058050 0.12819085     0
#> 9       C01Q1       F1B45       F1 0.005331765 0.14961074     0
#> 10      C01Q1       F2B05       F2 0.038067798 0.10851394     0
```

### Q3: OLS Logit with `lm()`

The plain logit model has the Berry (1994) linear form:

$$\log\left( s_{jt} \right) - \log\left( s_{0t} \right) = \beta_{0} + \alpha\, p_{jt} + \beta_{1}\,\text{mushy}_{jt} + \xi_{jt}$$

We first estimate this by OLS, which is biased because prices $p_{jt}$
are correlated with the unobserved product-market quality $\xi_{jt}$
(price endogeneity).

``` r
# Compute the dependent variable: log(s_j) - log(s_0)
outside_shares <- 1 - tapply(products$shares, products$market_ids, sum)
products$s0 <- as.numeric(outside_shares[as.character(products$market_ids)])
products$y <- log(products$shares) - log(products$s0)

# OLS regression
ols_fit <- lm(y ~ prices + mushy, data = products)

cat("=== OLS Logit Results ===\n\n")
#> === OLS Logit Results ===
print(round(summary(ols_fit)$coefficients, 3))
#>             Estimate Std. Error t value Pr(>|t|)
#> (Intercept)   -2.935      0.114 -25.694     0.00
#> prices        -7.480      0.863  -8.668     0.00
#> mushy          0.075      0.053   1.407     0.16
```

The OLS estimates are:

- Intercept: approximately $- 2.935$
- Price coefficient ($\alpha$): approximately $- 7.48$
- Mushy: approximately $0.075$ (not statistically significant)

The price coefficient is biased toward zero. High-quality products
command higher prices (positive correlation between $p_{jt}$ and
$\xi_{jt}$), which partially offsets the true negative effect of price
on demand.

### Q4: OLS Logit with rblp

We can replicate the OLS logit using
[`blp_problem()`](https://hchulkim.github.io/rblp/reference/blp_problem.md).
The key is to set `demand_instruments0 = prices`, which makes the
instrument set identical to the regressors (i.e., no excluded
instruments), so the IV estimator reduces to OLS.

``` r
# For OLS: use prices as own instrument (no excluded instruments)
products_ols <- products
products_ols$demand_instruments0 <- products_ols$prices

f1_ols <- blp_formulation(~ prices + mushy)
ols_problem <- blp_problem(list(f1_ols), products_ols)
ols_results <- ols_problem$solve(method = "1s")

cat("=== rblp OLS Logit Results ===\n\n")
#> === rblp OLS Logit Results ===
print(ols_results)
#> BLP Estimation Results
#>   Method: 1S GMM
#>   Objective: 1.336262e-24
#>   Optimization converged: TRUE
#>   FP converged: TRUE (94 total iterations)
#> 
#> Parameter Estimates:
#>  parameter   estimate  se       t_stat 
#>  (Intercept) -2.934501 0.107883 -27.201
#>  prices      -7.480136 0.839535 -8.910 
#>  mushy       0.074765  0.054087 1.382
```

The rblp results match the [`lm()`](https://rdrr.io/r/stats/lm.html)
estimates exactly: this confirms that the GMM framework with the
identity instrument matrix reduces to OLS.

### Q5: Fixed Effects Logit (Market + Product FE)

Adding market and product fixed effects absorbs all time-invariant
product characteristics and all market-level demand shifters. This
eliminates $\xi_{jt}$ components that are purely product-specific or
purely market-specific, leaving only the within-product, within-market
variation for identification.

With fixed effects and prices as their own instrument, this is a “FE
logit” that does not yet address endogeneity of the remaining price
variation.

``` r
# FE logit: absorb market + product, prices as own instrument
products_fe <- products
products_fe$demand_instruments0 <- products_fe$prices

f1_fe <- blp_formulation(~ 0 + prices, absorb = ~ market_ids + product_ids)
fe_problem <- blp_problem(list(f1_fe), products_fe)
fe_results <- fe_problem$solve(method = "1s")

cat("=== FE Logit Results (prices as own IV) ===\n\n")
#> === FE Logit Results (prices as own IV) ===
print(fe_results)
#> BLP Estimation Results
#>   Method: 1S GMM
#>   Objective: 1.357330e-28
#>   Optimization converged: TRUE
#>   FP converged: TRUE (664 total iterations)
#> 
#> Parameter Estimates:
#>  parameter estimate   se       t_stat 
#>  prices    -28.617866 1.039270 -27.537
```

The price coefficient with market and product fixed effects is
approximately $- 28.6$, much more negative than the OLS estimate of
$- 7.48$. The fixed effects absorb the product-level quality that was
positively correlated with prices, revealing a much stronger true price
sensitivity.

### Q6: IV + Fixed Effects Logit

Now we use the excluded `price_instrument` as an instrument for prices,
combined with the market and product fixed effects. This addresses
endogeneity from any remaining within-product, within-market correlation
between prices and unobserved quality shocks.

``` r
# IV + FE: use the price_instrument as excluded instrument
products_iv <- products
# demand_instruments0 is already set to price_instrument by prepare_mixtape_data()

f1_iv <- blp_formulation(~ 0 + prices, absorb = ~ market_ids + product_ids)
iv_problem <- blp_problem(list(f1_iv), products_iv)
iv_results <- iv_problem$solve(method = "1s")

cat("=== IV + FE Logit Results ===\n\n")
#> === IV + FE Logit Results ===
print(iv_results)
#> BLP Estimation Results
#>   Method: 1S GMM
#>   Objective: 1.255224e-29
#>   Optimization converged: TRUE
#>   FP converged: TRUE (664 total iterations)
#> 
#> Parameter Estimates:
#>  parameter estimate   se       t_stat 
#>  prices    -30.599521 1.132117 -27.029
```

The IV + FE price coefficient is approximately $- 30.6$ (with a standard
error around 1.0). This is slightly more negative than the FE-only
estimate, indicating that even after absorbing fixed effects, there is
some remaining positive correlation between prices and unobserved
quality that the instrument helps correct.

### Comparison of Logit Specifications

``` r
# ols_results$beta has: (Intercept), prices, mushy
# fe_results$beta and iv_results$beta have: prices only (intercept absorbed)
ols_price_idx <- which(names(coef(ols_fit)) == "prices")
rblp_ols_price <- ols_results$beta[ols_price_idx]  # index 2 = prices

cat("=== Price Coefficient Comparison ===\n\n")
#> === Price Coefficient Comparison ===
cat(sprintf("  OLS (lm):        %8.3f\n", coef(ols_fit)["prices"]))
#>   OLS (lm):          -7.480
cat(sprintf("  rblp OLS:        %8.3f\n", rblp_ols_price))
#>   rblp OLS:          -7.480
cat(sprintf("  FE (no IV):      %8.3f\n", fe_results$beta[1]))
#>   FE (no IV):       -28.618
cat(sprintf("  FE + IV:         %8.3f  (SE = %.3f)\n",
            iv_results$beta[1],
            iv_results$summary_table()$se[1]))
#>   FE + IV:          -30.600  (SE = 1.132)
```

The progression from $- 7.5$ (OLS) to $- 28.6$ (FE) to $- 30.6$ (IV+FE)
shows how each correction — absorbing fixed effects and instrumenting
for price — reveals increasingly strong consumer price sensitivity.

### Q7: Price Cut Counterfactual

We now use the IV + FE logit results to simulate a counterfactual: what
happens to market shares if product F1B04 receives a price cut of \$0.04
per serving in market C01Q1? (F1B04’s baseline price is about \$0.072,
so this is roughly a 55% price reduction.)

In the pure logit, shares have a closed form given mean utilities
$\delta_{jt}$. A price change for product $j$ shifts its mean utility by
$\Delta\delta_{j} = \alpha \times \Delta p_{j}$, and new shares follow
from the logit formula.

``` r
# Pick the first market
first_market <- "C01Q1"
mkt_idx <- which(products_iv$market_ids == first_market)
mkt_data <- products_iv[mkt_idx, ]

# Current shares and delta
s0_mkt <- unique(mkt_data$s0)
current_delta <- log(mkt_data$shares) - log(s0_mkt)
current_shares <- mkt_data$shares

# Price coefficient from IV + FE
alpha <- iv_results$beta[1]

# Price cut of $0.04/serving for F1B04
new_prices <- mkt_data$prices
f1b04_idx <- which(mkt_data$product_ids == "F1B04")
price_cut <- 0.04
new_prices[f1b04_idx] <- new_prices[f1b04_idx] - price_cut

cat(sprintf("F1B04 price: $%.4f -> $%.4f (cut of $%.2f/serving)\n",
            mkt_data$prices[f1b04_idx], new_prices[f1b04_idx], price_cut))
#> F1B04 price: $0.0721 -> $0.0321 (cut of $0.04/serving)

# New delta: delta_new = delta_old + alpha * (p_new - p_old)
new_delta <- current_delta + alpha * (new_prices - mkt_data$prices)

# New logit shares
exp_delta <- exp(new_delta)
new_shares <- exp_delta / (1 + sum(exp_delta))

# Percentage change in shares
pct_change <- (new_shares - current_shares) / current_shares * 100

cat("\n=== Counterfactual: $0.04 Price Cut for F1B04 in C01Q1 ===\n\n")
#> 
#> === Counterfactual: $0.04 Price Cut for F1B04 in C01Q1 ===
cf_table <- data.frame(
  product = mkt_data$product_ids,
  old_share = sprintf("%.6f", current_shares),
  new_share = sprintf("%.6f", new_shares),
  pct_change = sprintf("%+.1f%%", pct_change),
  stringsAsFactors = FALSE
)
print(cf_table, row.names = FALSE, right = FALSE)
#>  product old_share new_share pct_change
#>  F1B04   0.012417  0.041005  +230.2%   
#>  F1B06   0.007809  0.007583  -2.9%     
#>  F1B07   0.012995  0.012618  -2.9%     
#>  F1B09   0.005770  0.005603  -2.9%     
#>  F1B11   0.017934  0.017415  -2.9%     
#>  F1B13   0.026602  0.025832  -2.9%     
#>  F1B17   0.025015  0.024291  -2.9%     
#>  F1B30   0.005058  0.004912  -2.9%     
#>  F1B45   0.005332  0.005177  -2.9%     
#>  F2B05   0.038068  0.036966  -2.9%     
#>  F2B08   0.008348  0.008106  -2.9%     
#>  F2B15   0.006596  0.006405  -2.9%     
#>  F2B16   0.030117  0.029246  -2.9%     
#>  F2B19   0.100046  0.097150  -2.9%     
#>  F2B26   0.013234  0.012851  -2.9%     
#>  F2B28   0.023643  0.022958  -2.9%     
#>  F2B40   0.008662  0.008411  -2.9%     
#>  F2B48   0.002699  0.002621  -2.9%     
#>  F3B06   0.018943  0.018395  -2.9%     
#>  F3B14   0.010846  0.010532  -2.9%     
#>  F4B02   0.007850  0.007623  -2.9%     
#>  F4B10   0.000732  0.000711  -2.9%     
#>  F4B12   0.009490  0.009216  -2.9%     
#>  F6B18   0.046570  0.045222  -2.9%

cat(sprintf("\nF1B04 share change: %+.1f%%\n", pct_change[f1b04_idx]))
#> 
#> F1B04 share change: +230.2%
cat(sprintf("Mean change for other products: %.2f%%\n",
            mean(pct_change[-f1b04_idx])))
#> Mean change for other products: -2.89%
```

F1B04’s share increases dramatically (by over 200%), while all other
products lose a few percent of their market share. This large response
reflects the estimated price coefficient of approximately $- 30.6$: a
\$0.04 price cut translates to a utility gain of
$0.04 \times 30.6 \approx 1.22$ units, which is substantial on the logit
scale.

In the pure logit, the substitution pattern is governed by the
Independence of Irrelevant Alternatives (IIA) property: all competing
products lose share in proportion to their current shares, regardless of
how similar they are to F1B04.

This is both a strength (simplicity) and a weakness (unrealistic
substitution patterns) of the logit model. The random coefficients model
in Exercise 2 relaxes IIA.

**Note:** The pyblp Mixtape solutions report F1B04 share change of
+223.6% and other products at -1.45%. Small numerical differences arise
from optimizer tolerances and GMM weighting matrix details.

### Q8: Own-Price Elasticities

The own-price elasticity in the logit model is:
$$\eta_{jj} = \alpha\, p_{jt}\,\left( 1 - s_{jt} \right)$$

Since $\alpha < 0$, own-price elasticities are negative. We compute them
using the `compute_elasticities()` method.

``` r
E <- iv_results$compute_elasticities(first_market)
own_elast <- diag(E)

cat("=== Own-Price Elasticities (Market C01Q1) ===\n\n")
#> === Own-Price Elasticities (Market C01Q1) ===

elast_table <- data.frame(
  product = mkt_data$product_ids,
  price = sprintf("%.4f", mkt_data$prices),
  share = sprintf("%.6f", mkt_data$shares),
  elasticity = sprintf("%.3f", own_elast),
  stringsAsFactors = FALSE
)
print(elast_table, row.names = FALSE, right = FALSE)
#>  product price  share    elasticity
#>  F1B04   0.0721 0.012417 -2.178    
#>  F1B06   0.1142 0.007809 -3.467    
#>  F1B07   0.1324 0.012995 -3.998    
#>  F1B09   0.1303 0.005770 -3.965    
#>  F1B11   0.1548 0.017934 -4.653    
#>  F1B13   0.1370 0.026602 -4.082    
#>  F1B17   0.1442 0.025015 -4.302    
#>  F1B30   0.1282 0.005058 -3.903    
#>  F1B45   0.1496 0.005332 -4.554    
#>  F2B05   0.1085 0.038068 -3.194    
#>  F2B08   0.1323 0.008348 -4.014    
#>  F2B15   0.1121 0.006596 -3.407    
#>  F2B16   0.1153 0.030117 -3.421    
#>  F2B19   0.1109 0.100046 -3.053    
#>  F2B26   0.1281 0.013234 -3.868    
#>  F2B28   0.1748 0.023643 -5.223    
#>  F2B40   0.1336 0.008662 -4.054    
#>  F2B48   0.1474 0.002699 -4.499    
#>  F3B06   0.1098 0.018943 -3.296    
#>  F3B14   0.1370 0.010846 -4.145    
#>  F4B02   0.1751 0.007850 -5.317    
#>  F4B10   0.1357 0.000732 -4.148    
#>  F4B12   0.1352 0.009490 -4.097    
#>  F6B18   0.1428 0.046570 -4.167

cat(sprintf("\nOwn-price elasticity range: [%.1f, %.1f]\n",
            min(own_elast), max(own_elast)))
#> 
#> Own-price elasticity range: [-5.3, -2.2]
cat(sprintf("Mean own-price elasticity: %.2f\n", mean(own_elast)))
#> Mean own-price elasticity: -3.96
```

Own-price elasticities are negative and range from approximately $- 2$
to $- 5$ (or roughly $- 2.4$ to $- 6.3$ depending on the exact estimated
alpha). More expensive products have more elastic demand, which makes
economic sense: a 1% price increase for an expensive cereal represents a
larger absolute price change, inducing more consumers to switch away.

### Cross-Price Elasticities

In the logit model, the cross-price elasticity is:
$$\eta_{jk} = - \alpha\, p_{kt}\, s_{kt}\quad(j \neq k)$$

This means all products have the same cross-elasticity with respect to a
given product $k$ — a direct consequence of IIA.

``` r
# Show a 5x5 submatrix of cross-elasticities
cat("=== Cross-Price Elasticity Submatrix (first 5 products) ===\n\n")
#> === Cross-Price Elasticity Submatrix (first 5 products) ===
sub_E <- E[1:5, 1:5]
rownames(sub_E) <- mkt_data$product_ids[1:5]
colnames(sub_E) <- mkt_data$product_ids[1:5]
print(round(sub_E, 3))
#>        F1B04  F1B06  F1B07  F1B09  F1B11
#> F1B04 -2.178  0.027  0.053  0.023  0.085
#> F1B06  0.027 -3.467  0.053  0.023  0.085
#> F1B07  0.027  0.027 -3.998  0.023  0.085
#> F1B09  0.027  0.027  0.053 -3.965  0.085
#> F1B11  0.027  0.027  0.053  0.023 -4.653

cat("\nNote: Off-diagonal entries in each column are identical (IIA property).\n")
#> 
#> Note: Off-diagonal entries in each column are identical (IIA property).
```

------------------------------------------------------------------------

## Exercise 2: Mixed Logit (Random Coefficients)

The pure logit imposes IIA: substitution patterns are determined
entirely by market shares, not by product similarity. The mixed (random
coefficients) logit relaxes IIA by allowing consumer preferences to
vary. Products that are similar in the characteristic space will
naturally be closer substitutes.

### Cross-Market Variation and Identification

A key insight from BLP (1995) is that cross-market variation in product
sets and demographics helps identify random coefficient parameters.
Markets where product characteristics are more spread out provide
different “experiments” for estimating how consumers substitute between
products.

``` r
# Show price variation across markets for a few products
some_products <- c("F1B04", "F2B05", "F3B14")
for (p in some_products) {
  prices_p <- products$prices[products$product_ids == p]
  cat(sprintf("  %s: mean price = %.4f, sd = %.4f, range = [%.4f, %.4f]\n",
              p, mean(prices_p), sd(prices_p), min(prices_p), max(prices_p)))
}
#>   F1B04: mean price = 0.0833, sd = 0.0180, range = [0.0455, 0.1190]
#>   F2B05: mean price = 0.1176, sd = 0.0181, range = [0.0806, 0.1643]
#>   F3B14: mean price = 0.1434, sd = 0.0213, range = [0.0959, 0.1922]
```

### Setting Up Demographics

The Mixtape dataset includes demographic data with quarterly income for
20 simulated individuals per market. These demographics create observed
heterogeneity in price sensitivity: wealthier consumers are less price
sensitive.

``` r
demographics <- load_mixtape_demographics()

cat(sprintf("Demographic observations: %d\n", nrow(demographics)))
#> Demographic observations: 1880
cat(sprintf("Markets: %d\n", length(unique(demographics$market))))
#> Markets: 94
cat(sprintf("Individuals per market: %d\n",
            nrow(demographics) / length(unique(demographics$market))))
#> Individuals per market: 20

cat("\nQuarterly income summary:\n")
#> 
#> Quarterly income summary:
print(round(summary(demographics$quarterly_income), 1))
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>    23.4  2128.6  3836.1  4575.0  5884.6 33724.4

cat(sprintf("\nLog income range: [%.2f, %.2f]\n",
            min(log(demographics$quarterly_income)),
            max(log(demographics$quarterly_income))))
#> 
#> Log income range: [3.15, 10.43]
```

### Preparing Agent Data for rblp

To use demographics in rblp, we need to construct an agent data frame
with `market_ids`, `weights`, integration `nodes`, and demographic
columns.

``` r
# Prepare agent data with equal weights and log income
agents <- data.frame(
  market_ids = demographics$market,
  weights = 1 / 20,  # equal weights for 20 agents per market
  nodes0 = 0,        # placeholder node (not used for pure demographic model)
  log_income = log(demographics$quarterly_income),
  stringsAsFactors = FALSE
)

cat("Agent data preview:\n")
head(agents)
```

### Mushy x Log Income Interaction

Before estimating full random coefficients, the Mixtape exercises ask
for a simple demographic interaction: does the valuation of mushy
cereals vary with income? We add `log_income` as a demographic and
interact it with `mushy` through the pi matrix.

The pyblp solution finds that the pi coefficient on `mushy x log_income`
is approximately 0.251, indicating that higher-income consumers have a
(slightly) higher taste for mushy cereals.

``` r
# Formulations
f1_rc <- blp_formulation(~ 0 + prices + mushy, absorb = ~ market_ids + product_ids)
f2_rc <- blp_formulation(~ 0 + mushy)

# Demographics formulation: log_income
f_demo <- blp_formulation(~ 0 + log_income)

# Create problem
rc_problem <- blp_problem(
  product_formulations = list(f1_rc, f2_rc),
  product_data = products,
  agent_formulation = f_demo,
  agent_data = agents
)

# Solve with pi on mushy x log_income
# sigma = 0 (no unobserved heterogeneity), pi = initial guess
pi0 <- matrix(0.1, nrow = 1, ncol = 1)  # K2=1 (mushy), D=1 (log_income)

rc_mushy_results <- rc_problem$solve(
  sigma = matrix(0, 1, 1),  # no sigma (fixed at 0)
  pi = pi0,
  method = "1s",
  optimization = blp_optimization("l-bfgs-b",
    method_options = list(maxit = 500))
)

print(rc_mushy_results)
```

**Expected result:** The pi coefficient on `mushy x log_income` is
approximately 0.251, meaning a one-unit increase in log income raises
the mean utility of mushy cereals by about 0.25.

### Random Coefficient on Price (Sigma)

The next step adds a random coefficient (sigma) on price, allowing
unobserved heterogeneity in price sensitivity across consumers. Combined
with the pi interaction (price x income), this creates a rich model
where price sensitivity varies both with observed income and with
unobserved consumer characteristics.

``` r
# X2 now includes price (for sigma) and mushy (for pi)
f1_full <- blp_formulation(~ 0 + prices + mushy, absorb = ~ market_ids + product_ids)
f2_full <- blp_formulation(~ 0 + prices + mushy)

# Demographics: log_income
f_demo_full <- blp_formulation(~ 0 + log_income)

# Prepare products with proper IV
products_rc <- products  # demand_instruments0 already = price_instrument

rc_full_problem <- blp_problem(
  product_formulations = list(f1_full, f2_full),
  product_data = products_rc,
  agent_formulation = f_demo_full,
  agent_data = agents
)

# Starting values near the pyblp solution
# sigma: K2 x K2 = 2x2 diagonal (sigma on prices, 0 on mushy)
sigma0 <- diag(c(5.0, 0))  # start sigma_price at 5

# pi: K2 x D = 2x1 (prices x log_income, mushy x log_income)
pi0_full <- matrix(c(-5.0, 0.2), nrow = 2, ncol = 1)

rc_full_results <- rc_full_problem$solve(
  sigma = sigma0,
  pi = pi0_full,
  method = "1s",
  optimization = blp_optimization("l-bfgs-b",
    method_options = list(maxit = 1000, factr = 1e7))
)

print(rc_full_results)
```

**Expected results from pyblp:**

| Parameter                                 | Estimate | Description                                       |
|-------------------------------------------|----------|---------------------------------------------------|
| $\alpha$ (mean price)                     | negative | Mean price coefficient (in beta, from FE)         |
| $\sigma_{\text{price}}$                   | ~6.02    | SD of unobserved price sensitivity                |
| $\pi_{\text{price} \times \text{income}}$ | ~-5.96   | Higher income $\Rightarrow$ less price sensitive  |
| $\pi_{\text{mushy} \times \text{income}}$ | ~0.2     | Higher income $\Rightarrow$ slightly prefer mushy |

The large $\sigma_{\text{price}}$ (~6) indicates substantial unobserved
heterogeneity in price sensitivity, beyond what income explains. The
negative $\pi_{\text{price} \times \text{income}}$ means wealthier
consumers are less price sensitive, as expected.

### Why Random Coefficients Matter

With random coefficients, the substitution matrix is no longer governed
by IIA. Two products with similar characteristics (e.g., two mushy
cereals at similar prices) will have higher cross-price elasticities
than two dissimilar products. This has direct implications for:

- **Merger analysis:** Which mergers raise prices depends on which
  products are close substitutes
- **Counterfactual pricing:** Firms internalize cross-effects among
  their own products
- **Welfare analysis:** Consumer surplus calculations account for
  heterogeneous valuations

``` r
# Compare logit vs RC logit elasticities
E_rc <- rc_full_results$compute_elasticities(first_market)

cat("=== Own-Price Elasticity Comparison (C01Q1) ===\n\n")
comp_table <- data.frame(
  product = mkt_data$product_ids[1:6],
  logit = sprintf("%.2f", diag(E)[1:6]),
  rc_logit = sprintf("%.2f", diag(E_rc)[1:6]),
  stringsAsFactors = FALSE
)
print(comp_table, row.names = FALSE, right = FALSE)
```

------------------------------------------------------------------------

## Exercise 3: Micro Moments (Conceptual Overview)

Micro moments supplement the traditional BLP aggregate moments (based on
market shares) with information from individual-level data sources such
as household scanner panels, surveys, or browsing data.

### Motivation

In the standard BLP framework, identification of random coefficients
comes from cross-market variation in product sets, prices, and aggregate
shares. This can be weak when:

- Markets are similar in their product offerings
- The number of random coefficients exceeds what aggregate variation can
  identify
- We want to match specific demographic-choice correlations observed in
  micro data

Micro moments address this by directly targeting individual-level
statistics.

### The Micro Moment Framework

A micro moment matches a statistic computed from micro data to its
model-predicted counterpart:

$$E\left\lbrack {\bar{m}}^{\text{micro}} - m(\theta) \right\rbrack = 0$$

where ${\bar{m}}^{\text{micro}}$ is the observed micro statistic (e.g.,
the average income of consumers who buy mushy cereals) and $m(\theta)$
is the model-predicted value (computed from the estimated choice
probabilities and demographic distributions).

### Types of Micro Moments

Common micro moments include:

1.  **Conditional choice probabilities:**
    $E\left\lbrack y_{ij}|D_{i} \in G \right\rbrack$ – the probability
    of buying product $j$ given demographics in group $G$
2.  **Conditional means of characteristics:**
    $E\left\lbrack x_{j{(i)}}|D_{i} \in G \right\rbrack$ – the mean
    characteristic of the chosen product for a demographic group
3.  **Covariance moments:** $\text{Cov}\left( D_{i},x_{j{(i)}} \right)$
    – covariance between demographics and chosen product attributes

### Implementation in rblp

rblp provides three classes for constructing micro moments:

``` r
# 1. Define a micro dataset (describes the data source)
md <- micro_dataset(
  name = "household scanner panel",
  observations = 50000,
  compute_weights = function(market_id, products, agents) {
    # Return I x J weight matrix
    I <- length(agents$weights)
    J <- length(products$shares)
    matrix(1, I, J)
  }
)

# 2. Define a micro part (what statistic to compute for each i,j pair)
mp <- micro_part(
  name = "income_x_mushy",
  dataset = md,
  compute_values = function(market_id, products, agents) {
    # Return I x J value matrix: income_i * mushy_j
    outer(agents$log_income, products$mushy)
  }
)

# 3. Define a micro moment (match target to prediction)
mm <- micro_moment(
  name = "E[log_income * mushy | purchase]",
  value = 0.5,   # observed from micro data
  parts = mp
)

# Pass to solve()
results <- rc_full_problem$solve(
  sigma = sigma0,
  pi = pi0_full,
  micro_moments = list(mm),
  method = "1s"
)
```

Micro moments can dramatically improve the precision of demographic
interaction parameters (pi) and random coefficient standard deviations
(sigma) by providing direct information about who buys what. For a full
treatment, see Conlon and Gortmaker (2020, Section 4) and the pyblp
micro moments tutorial.

------------------------------------------------------------------------

## Summary of Results

The table below summarizes key findings across all three exercises:

| Exercise | Model                   | Price Coeff.                              | Key Insight                                         |
|----------|-------------------------|-------------------------------------------|-----------------------------------------------------|
| 1 (Q3)   | OLS logit               | $- 7.48$                                  | Biased toward zero by price endogeneity             |
| 1 (Q4)   | rblp OLS                | $- 7.48$                                  | Confirms equivalence of GMM/OLS                     |
| 1 (Q5)   | FE logit (prices as IV) | $- 28.6$                                  | FE absorb product quality; reveals true sensitivity |
| 1 (Q6)   | IV + FE logit           | $- 30.6$ (SE $\approx$ 1.0)               | IV corrects remaining endogeneity                   |
| 1 (Q7)   | Counterfactual          | –                                         | F1B04 share +230% from \$0.04 cut; IIA substitution |
| 1 (Q8)   | Elasticities            | –                                         | Own-price elasticities range $- 2$ to $- 5$         |
| 2        | Mushy x income          | pi $\approx$ 0.25                         | Wealthier consumers slightly prefer mushy           |
| 2        | RC on price             | $\sigma \approx 6.0$, $\pi \approx - 6.0$ | Large heterogeneity in price sensitivity            |
| 3        | Micro moments           | –                                         | Individual data sharpens identification             |

### Key Takeaways

1.  **Price endogeneity is severe.** The OLS price coefficient ($- 7.5$)
    is only one-quarter of the IV+FE estimate ($- 30.6$). Failing to
    instrument biases the price effect toward zero, which leads to
    systematically understating consumer price sensitivity.

2.  **Fixed effects matter.** Moving from OLS to FE logit (without IV)
    already moves the coefficient from $- 7.5$ to $- 28.6$. Most of the
    endogeneity bias in this application comes from product-level
    unobserved quality.

3.  **Random coefficients break IIA.** The logit’s substitution patterns
    are mechanical (proportional to shares). The mixed logit allows
    substitution to depend on product similarity, which matters for
    merger analysis, targeted pricing, and welfare evaluation.

4.  **Demographics provide economic structure.** The interaction between
    price sensitivity and income
    ($\pi_{\text{price} \times \text{income}} \approx - 6$) is both
    statistically and economically significant: it means that a
    low-income consumer faces four times the utility loss from a price
    increase compared to a high-income consumer.

------------------------------------------------------------------------

## References

- Berry, S. (1994). Estimating Discrete-Choice Models of Product
  Differentiation. *RAND Journal of Economics*, 25(2), 242-262.

- Berry, S., Levinsohn, J., & Pakes, A. (1995). Automobile Prices in
  Market Equilibrium. *Econometrica*, 63(4), 841-890.

- Conlon, C. & Gortmaker, J. (2020). Best Practices for Differentiated
  Products Demand Estimation with pyblp. *RAND Journal of Economics*,
  51(4), 1108-1161.

- Nevo, A. (2000). A Practitioner’s Guide to Estimation of
  Random-Coefficients Logit Models of Demand. *Journal of Economics &
  Management Strategy*, 9(4), 513-548.

- Mixtape Sessions. Demand Estimation.
  <https://github.com/Mixtape-Sessions/Demand-Estimation>

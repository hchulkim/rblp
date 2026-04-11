# Preparing Data for rblp

## Overview

`rblp` expects product and agent data in a specific format. This guide
explains how to prepare your data for BLP estimation. Both `data.frame`
and `data.table` inputs are supported — `rblp` converts internally.

## Product Data

The product data must be a data frame (or data.table) with one row per
product-market observation. The following columns are recognized:

### Required Columns

| Column       | Type      | Description                                                                               |
|--------------|-----------|-------------------------------------------------------------------------------------------|
| `market_ids` | character | Market identifier (e.g., city-quarter, store-week)                                        |
| `shares`     | numeric   | Market share $s_{jt}$ (must be in $(0,1)$; sum across products in a market must be $< 1$) |
| `prices`     | numeric   | Product price $p_{jt}$                                                                    |

### Commonly Used Columns

| Column           | Type              | Description                                                                          |
|------------------|-------------------|--------------------------------------------------------------------------------------|
| `firm_ids`       | character/numeric | Firm identifier (required for supply-side, merger simulation, BLP instruments)       |
| `product_ids`    | character/factor  | Product identifier (for fixed effects via `absorb = ~ product_ids`)                  |
| `nesting_ids`    | character         | Nest identifier (for nested logit; products in the same nest are closer substitutes) |
| `clustering_ids` | character         | Cluster identifier (for clustered standard errors)                                   |

### Instrument Columns

Excluded demand and supply instruments are detected by name pattern:

| Pattern                                         | Description                                                                        |
|-------------------------------------------------|------------------------------------------------------------------------------------|
| `demand_instruments0`, `demand_instruments1`, … | Excluded demand-side instruments                                                   |
| `supply_instruments0`, `supply_instruments1`, … | Excluded supply-side instruments (only used when a supply formulation is provided) |

Exogenous product characteristics from $X_{1}$ (everything except prices
and shares) are automatically added to the instrument set. You only need
to supply **excluded** instruments — those that shift prices but do not
enter utility directly.

### Product Characteristics

Any additional numeric columns can be used as product characteristics in
formulations (e.g., `sugar`, `mushy`, `hpwt`, `mpd`, `space`).

## Example: Building Product Data from Scratch

``` r
library(rblp)

# Suppose you have raw sales data
set.seed(42)
raw <- data.frame(
  city     = rep(c("NYC", "LA", "CHI"), each = 30),
  quarter  = rep(rep(1:2, each = 15), 3),
  brand    = rep(paste0("brand_", 1:5), 18),
  firm     = rep(c("A", "A", "B", "B", "C"), 18),
  units    = rpois(90, lambda = 500),
  pop      = rep(c(1e6, 8e5, 5e5), each = 30),
  price    = runif(90, 2, 8),
  sugar    = rep(runif(5, 1, 15), 18),
  cost     = runif(90, 1, 5)
)

# Step 1: Create market identifier
raw$market_ids <- paste(raw$city, raw$quarter, sep = "_")

# Step 2: Compute market shares
# Market share = units / market size (population)
# The outside share s_0 = 1 - sum(s_j) must be positive!
raw$shares <- raw$units / raw$pop

# Verify: inside shares sum to < 1 in each market
inside_totals <- tapply(raw$shares, raw$market_ids, sum)
stopifnot(all(inside_totals < 1))
cat("Inside share range:", range(inside_totals), "\n")
#> Inside share range: 0.007537 0.015024

# Step 3: Rename columns to match rblp conventions
products <- data.frame(
  market_ids = raw$market_ids,
  firm_ids   = raw$firm,
  product_ids = raw$brand,
  shares     = raw$shares,
  prices     = raw$price,
  sugar      = raw$sugar,
  stringsAsFactors = FALSE
)

# Step 4: Add instruments
# BLP instruments: sums of rival/own-firm characteristics
X_exog <- as.matrix(products[, "sugar", drop = FALSE])
blp_iv <- build_blp_instruments(X_exog, products$market_ids, products$firm_ids)
for (k in seq_len(ncol(blp_iv))) {
  products[[paste0("demand_instruments", k - 1)]] <- blp_iv[, k]
}

cat("Product data ready:", nrow(products), "observations\n")
#> Product data ready: 90 observations
cat("Columns:", paste(names(products), collapse = ", "), "\n")
#> Columns: market_ids, firm_ids, product_ids, shares, prices, sugar, demand_instruments0, demand_instruments1
```

## Using data.table

`rblp` accepts `data.table` objects directly. No conversion needed:

``` r
library(data.table)

# Read data as data.table
products <- fread("my_data.csv")

# Add columns using data.table syntax
products[, market_ids := paste(city, quarter, sep = "_")]
products[, shares := units / pop]

# Build instruments
X_exog <- as.matrix(products[, .(sugar)])
iv <- build_blp_instruments(X_exog, products$market_ids, products$firm_ids)
products[, paste0("demand_instruments", 0:(ncol(iv)-1)) := as.data.frame(iv)]

# Pass directly to blp_problem (no as.data.frame() needed)
problem <- blp_problem(list(blp_formulation(~ prices + sugar)), products)
```

## Agent Data (for Random Coefficients with Demographics)

If you have consumer demographics (income, age, etc.) and want to
estimate how preferences vary with observables, you need agent data:

### Required Columns

| Column       | Type      | Description                                      |
|--------------|-----------|--------------------------------------------------|
| `market_ids` | character | Must match product data market identifiers       |
| `weights`    | numeric   | Agent weights (must sum to 1 within each market) |

### Optional Columns

| Pattern               | Description                                                                                                                                                                                                                 |
|-----------------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `nodes0`, `nodes1`, … | Integration nodes (draws from the mixing distribution). If using [`blp_integration()`](https://hchulkim.github.io/rblp/reference/blp_integration.md), these are generated automatically and you do not need to supply them. |
| Demographic columns   | Any numeric columns referenced in the agent formulation (e.g., `income`, `age`, `child`)                                                                                                                                    |

### Example: Agent Data with Demographics

``` r
agents <- load_nevo_agents()
cat("Agent data columns:", paste(names(agents), collapse = ", "), "\n")
#> Agent data columns: market_ids, city_ids, quarter, weights, nodes0, nodes1, nodes2, nodes3, income, income_squared, age, child
cat("Agents per market:", nrow(agents) / length(unique(agents$market_ids)), "\n")
#> Agents per market: 20
cat("Weights sum per market:", tapply(agents$weights, agents$market_ids, sum)[1], "\n")
#> Weights sum per market: 1
```

## Integration Nodes (No Agent Data Required)

For random coefficients without demographics, use
[`blp_integration()`](https://hchulkim.github.io/rblp/reference/blp_integration.md)
to generate integration nodes automatically:

``` r
# Gauss-Hermite product rule: 5^K2 nodes per market
int_gh <- blp_integration("product", size = 5)

# Monte Carlo: 200 draws per market
int_mc <- blp_integration("monte_carlo", size = 200, seed = 42)

# Halton quasi-random: 100 draws per market
int_halton <- blp_integration("halton", size = 100)
```

## Formulation Reference

### Linear Demand (X1)

``` r
# With intercept (default)
blp_formulation(~ prices + sugar + mushy)

# Without intercept
blp_formulation(~ 0 + prices + sugar + mushy)

# With absorbed fixed effects
blp_formulation(~ prices, absorb = ~ product_ids)

# String formula also works
blp_formulation("~ prices + sugar + mushy")
```

### Nonlinear Demand (X2) — Random Coefficients

``` r
# Variables whose coefficients vary across consumers
blp_formulation(~ prices + sugar + mushy)

# Intercept only (random coefficient on constant)
blp_formulation(~ 1)
```

### Supply Side (X3) — Cost Equation

``` r
# Cost shifters: mc_j = X3_j' gamma + omega_j
blp_formulation(~ log_hpwt + air + log_mpg + log_space + trend)
```

### Demographics (Agent Formulation)

``` r
# No intercept (demographics interact with X2 characteristics)
blp_formulation(~ 0 + income + income_squared + age + child)
```

## Common Pitfalls

1.  **Shares must be strictly between 0 and 1.** Zero or negative shares
    cause `log(s)` to fail. Filter out zero-sales observations or add a
    small positive constant.

2.  **Inside shares must sum to less than 1.** The outside good share
    $s_{0} = 1 - \sum_{j}s_{j}$ must be positive. If your shares sum to
    1, you need to redefine the market size to include non-purchasers.

3.  **Price column must be named `prices`** (with an “s”). This is how
    `rblp` identifies the endogenous regressor for instrument
    construction and elasticity computation.

4.  **Instrument columns must follow the naming convention:**
    `demand_instruments0`, `demand_instruments1`, etc. Arbitrary column
    names will not be detected as instruments.

5.  **Agent weights must sum to 1** within each market. If using
    [`blp_integration()`](https://hchulkim.github.io/rblp/reference/blp_integration.md),
    this is handled automatically.

6.  **With `absorb = ~ product_ids`**, the `product_ids` column must
    exist in the data. The formula variable name must match the column
    name.

## pyblp Correspondence

| pyblp                                                                | rblp                                                                                       |
|----------------------------------------------------------------------|--------------------------------------------------------------------------------------------|
| `Formulation('prices + sugar + mushy')`                              | `blp_formulation(~ prices + sugar + mushy)`                                                |
| `Formulation('0 + prices', absorb='C(product_ids)')`                 | `blp_formulation(~ 0 + prices, absorb = ~ product_ids)`                                    |
| `Problem([f1, f2, f3], product_data, agent_formulation, agent_data)` | `blp_problem(list(f1, f2, f3), products, agent_formulation = f_demo, agent_data = agents)` |
| `Integration('product', size=5)`                                     | `blp_integration("product", size = 5)`                                                     |
| `Iteration('squarem')`                                               | `blp_iteration("squarem")`                                                                 |
| `Optimization('l-bfgs-b')`                                           | `blp_optimization("l-bfgs-b")`                                                             |
| `results.compute_elasticities()`                                     | `results$compute_elasticities()`                                                           |
| `results.compute_consumer_surplus()`                                 | `results$compute_consumer_surplus()`                                                       |

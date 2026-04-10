# Create a BLP Problem

Main entry point for BLP demand estimation.

## Usage

``` r
blp_problem(
  product_formulations,
  product_data,
  agent_formulation = NULL,
  agent_data = NULL,
  integration = NULL,
  rc_types = NULL,
  epsilon_scale = 1,
  costs_type = "linear",
  add_exogenous = TRUE
)
```

## Arguments

- product_formulations:

  List of BLPFormulation objects (1-3)

- product_data:

  Data frame with market_ids, shares, prices, etc.

- agent_formulation:

  Optional demographics formulation

- agent_data:

  Optional agent data frame

- integration:

  Optional BLPIntegration object

- rc_types:

  Character vector of random coefficient types

- epsilon_scale:

  Epsilon scaling (default 1)

- costs_type:

  "linear" or "log"

- add_exogenous:

  Whether to add exogenous regressors to instruments

## Value

A BLPProblem object

## Examples

``` r
if (FALSE) { # \dontrun{
# Logit model
f1 <- blp_formulation(~ prices + sugar + mushy)
problem <- blp_problem(list(f1), product_data)
results <- problem$solve()

# Random coefficients
f2 <- blp_formulation(~ prices + sugar + mushy)
problem <- blp_problem(list(f1, f2), product_data,
                        integration = blp_integration("product", 5))
results <- problem$solve(sigma = diag(3))
} # }
```

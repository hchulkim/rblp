# Create a BLP Simulation

Create a BLP Simulation

## Usage

``` r
blp_simulation(
  product_formulations,
  product_data,
  beta,
  sigma = NULL,
  pi = NULL,
  gamma = NULL,
  rho = NULL,
  agent_formulation = NULL,
  agent_data = NULL,
  integration = NULL,
  xi = NULL,
  omega = NULL,
  xi_variance = 1,
  omega_variance = 1,
  correlation = 0.9,
  rc_types = NULL,
  costs_type = "linear",
  seed = NULL
)
```

## Arguments

- product_formulations:

  List of BLPFormulation objects

- product_data:

  Data frame with market_ids, firm_ids, and characteristics

- beta:

  Demand linear coefficients

- sigma:

  Optional Cholesky root of RC covariance

- pi:

  Optional demographics interaction

- gamma:

  Optional supply coefficients

- rho:

  Optional nesting parameters

- agent_formulation:

  Optional demographics formulation

- agent_data:

  Optional agent data

- integration:

  Optional BLPIntegration

- xi:

  Optional demand errors

- omega:

  Optional supply errors

- xi_variance:

  Variance of xi

- omega_variance:

  Variance of omega

- correlation:

  Correlation between xi and omega

- rc_types:

  Character vector of RC types

- costs_type:

  "linear" or "log"

- seed:

  Random seed

## Value

A BLPSimulation object

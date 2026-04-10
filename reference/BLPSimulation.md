# BLP Simulation

Simulate BLP-type equilibrium data for testing and Monte Carlo studies.

## Super class

[`rblp::BLPEconomy`](https://hchulkim.github.io/rblp/reference/BLPEconomy.md)
-\> `BLPSimulation`

## Public fields

- `beta`:

  True demand linear parameters

- `sigma`:

  True Cholesky root of RC covariance

- `pi`:

  True demographics interaction

- `gamma`:

  True supply parameters

- `rho`:

  True nesting parameters

- `xi`:

  Demand structural errors

- `omega`:

  Supply structural errors

## Methods

### Public methods

- [`BLPSimulation$new()`](#method-BLPSimulation-new)

- [`BLPSimulation$replace_endogenous()`](#method-BLPSimulation-replace_endogenous)

- [`BLPSimulation$print()`](#method-BLPSimulation-print)

- [`BLPSimulation$clone()`](#method-BLPSimulation-clone)

Inherited methods

- [`rblp::BLPEconomy$compute_logit_delta()`](https://hchulkim.github.io/rblp/reference/BLPEconomy.html#method-compute_logit_delta)
- [`rblp::BLPEconomy$get_market_data()`](https://hchulkim.github.io/rblp/reference/BLPEconomy.html#method-get_market_data)

------------------------------------------------------------------------

### Method `new()`

Create a BLP simulation

#### Usage

    BLPSimulation$new(
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

#### Arguments

- `product_formulations`:

  List of BLPFormulation objects (1-3)

- `product_data`:

  Data frame with market_ids, firm_ids, and characteristics

- `beta`:

  Demand linear coefficients (required)

- `sigma`:

  Optional K2 x K2 Cholesky root

- `pi`:

  Optional K2 x D demographics interaction

- `gamma`:

  Optional supply coefficients

- `rho`:

  Optional nesting parameters

- `agent_formulation`:

  Optional demographics formulation

- `agent_data`:

  Optional agent data

- `integration`:

  Optional BLPIntegration

- `xi`:

  Optional demand errors (drawn if NULL)

- `omega`:

  Optional supply errors (drawn if NULL)

- `xi_variance`:

  Variance of xi (default 1)

- `omega_variance`:

  Variance of omega (default 1)

- `correlation`:

  Correlation between xi and omega (default 0.9)

- `rc_types`:

  Character vector of RC types

- `costs_type`:

  "linear" or "log"

- `seed`:

  Random seed

------------------------------------------------------------------------

### Method `replace_endogenous()`

Solve for equilibrium prices and shares

#### Usage

    BLPSimulation$replace_endogenous(iteration = NULL, constant_costs = TRUE)

#### Arguments

- `iteration`:

  BLPIteration for fixed-point iteration

- `constant_costs`:

  Whether costs are independent of shares

#### Returns

A BLPSimulationResults object

------------------------------------------------------------------------

### Method [`print()`](https://rdrr.io/r/base/print.html)

Print the simulation

#### Usage

    BLPSimulation$print(...)

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    BLPSimulation$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

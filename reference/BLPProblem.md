# BLP Problem

The main user-facing class for BLP demand estimation. Inherits from
BLPEconomy and adds the solve() interface.

## Super class

[`rblp::BLPEconomy`](https://hchulkim.github.io/rblp/reference/BLPEconomy.md)
-\> `BLPProblem`

## Methods

### Public methods

- [`BLPProblem$new()`](#method-BLPProblem-new)

- [`BLPProblem$solve()`](#method-BLPProblem-solve)

- [`BLPProblem$print()`](#method-BLPProblem-print)

- [`BLPProblem$clone()`](#method-BLPProblem-clone)

Inherited methods

- [`rblp::BLPEconomy$compute_logit_delta()`](https://hchulkim.github.io/rblp/reference/BLPEconomy.html#method-compute_logit_delta)
- [`rblp::BLPEconomy$get_market_data()`](https://hchulkim.github.io/rblp/reference/BLPEconomy.html#method-get_market_data)

------------------------------------------------------------------------

### Method `new()`

Create a BLP Problem

#### Usage

    BLPProblem$new(
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

#### Arguments

- `product_formulations`:

  List of BLPFormulation objects (1-3). First: linear demand (X1).
  Second: nonlinear demand (X2). Third: supply (X3).

- `product_data`:

  Data frame with columns: market_ids, shares, prices, and optionally
  firm_ids, nesting_ids, demand_instruments\*, supply_instruments\*

- `agent_formulation`:

  Optional BLPFormulation for demographics

- `agent_data`:

  Optional data frame with columns: market_ids, weights, nodes\*

- `integration`:

  Optional BLPIntegration object

- `rc_types`:

  Character vector of random coefficient types ("linear", "log",
  "logit")

- `epsilon_scale`:

  Scaling for epsilon (default 1)

- `costs_type`:

  "linear" or "log"

- `add_exogenous`:

  Whether to add exogenous X1/X3 columns to instruments

------------------------------------------------------------------------

### Method [`solve()`](https://rdrr.io/r/base/solve.html)

Solve the BLP estimation problem

#### Usage

    BLPProblem$solve(
      sigma = NULL,
      pi = NULL,
      rho = NULL,
      beta = NULL,
      gamma = NULL,
      sigma_bounds = NULL,
      pi_bounds = NULL,
      rho_bounds = NULL,
      method = "2s",
      optimization = NULL,
      iteration = NULL,
      fp_type = "safe_linear",
      W_type = "robust",
      se_type = "robust",
      initial_W = NULL,
      initial_update = FALSE,
      scale_objective = TRUE,
      center_moments = TRUE,
      delta_behavior = "first",
      micro_moments = NULL,
      error_behavior = "revert",
      error_punishment = 1,
      processes = 1L
    )

#### Arguments

- `sigma`:

  Initial K2 x K2 Cholesky root (0 = fixed, non-zero = free)

- `pi`:

  Initial K2 x D demographics interaction (0 = fixed)

- `rho`:

  Initial nesting parameters (0 = fixed)

- `beta`:

  Demand linear coefficients (NA = concentrated out)

- `gamma`:

  Supply linear coefficients (NA = concentrated out)

- `sigma_bounds`:

  List of lower, upper bound matrices for sigma

- `pi_bounds`:

  List of lower, upper bound matrices for pi

- `rho_bounds`:

  List of lower, upper bound vectors for rho

- `method`:

  '1s' or '2s' for one-step or two-step GMM

- `optimization`:

  BLPOptimization object

- `iteration`:

  BLPIteration object (for contraction mapping)

- `fp_type`:

  Fixed point type: "safe_linear", "linear", "nonlinear"

- `W_type`:

  Weighting matrix type: "robust", "clustered", "unadjusted"

- `se_type`:

  Standard error type: "robust", "clustered", "unadjusted"

- `initial_W`:

  Optional initial weighting matrix

- `scale_objective`:

  Whether to scale objective by N

- `center_moments`:

  Whether to center moment conditions

- `delta_behavior`:

  How to initialize delta: "first", "logit", "last"

- `micro_moments`:

  Optional list of MicroMoment objects

- `error_behavior`:

  "revert", "punish", or "raise"

- `error_punishment`:

  Punishment scale for objective on error

- `processes`:

  Number of parallel processes (1 = no parallel)

#### Returns

A BLPResults object

------------------------------------------------------------------------

### Method [`print()`](https://rdrr.io/r/base/print.html)

Print the problem

#### Usage

    BLPProblem$print(...)

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    BLPProblem$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

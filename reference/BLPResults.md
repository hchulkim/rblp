# BLP Estimation Results

Stores estimation results and provides post-estimation methods.

## Public fields

- `problem`:

  The originating BLPProblem

- `sigma`:

  Estimated Cholesky root of random coefficient covariance

- `pi`:

  Estimated demographics interaction matrix

- `rho`:

  Estimated nesting parameters

- `beta`:

  Estimated demand linear parameters

- `gamma`:

  Estimated supply linear parameters

- `delta`:

  Estimated mean utilities

- `xi`:

  Demand-side structural error

- `omega`:

  Supply-side structural error

- `objective`:

  GMM objective value

- `gradient`:

  GMM gradient at solution

- `hessian`:

  GMM Hessian at solution

- `se`:

  Standard errors for nonlinear parameters

- `parameter_covariances`:

  Covariance matrix for nonlinear parameters

- `W`:

  Final weighting matrix

- `optimization_converged`:

  Whether optimization converged

- `optimization_iterations`:

  Number of optimization iterations

- `optimization_evaluations`:

  Number of function evaluations

- `fp_converged`:

  Whether all fixed points converged

- `fp_iterations`:

  Total fixed-point iterations

- `method`:

  GMM method ("1s" or "2s")

- `se_type`:

  Standard error type

## Methods

### Public methods

- [`BLPResults$new()`](#method-BLPResults-new)

- [`BLPResults$compute_elasticities()`](#method-BLPResults-compute_elasticities)

- [`BLPResults$compute_diversion_ratios()`](#method-BLPResults-compute_diversion_ratios)

- [`BLPResults$compute_costs()`](#method-BLPResults-compute_costs)

- [`BLPResults$compute_markups()`](#method-BLPResults-compute_markups)

- [`BLPResults$compute_consumer_surplus()`](#method-BLPResults-compute_consumer_surplus)

- [`BLPResults$compute_hhi()`](#method-BLPResults-compute_hhi)

- [`BLPResults$compute_merger()`](#method-BLPResults-compute_merger)

- [`BLPResults$run_hansen_test()`](#method-BLPResults-run_hansen_test)

- [`BLPResults$run_wald_test()`](#method-BLPResults-run_wald_test)

- [`BLPResults$sigma_squared()`](#method-BLPResults-sigma_squared)

- [`BLPResults$summary_table()`](#method-BLPResults-summary_table)

- [`BLPResults$print()`](#method-BLPResults-print)

- [`BLPResults$clone()`](#method-BLPResults-clone)

------------------------------------------------------------------------

### Method `new()`

Create results object

#### Usage

    BLPResults$new(
      problem,
      params,
      sigma,
      pi,
      rho,
      beta,
      gamma,
      delta,
      xi,
      omega,
      objective,
      gradient,
      hessian,
      se,
      parameter_covariances,
      W,
      step_results,
      optimization_converged,
      optimization_iterations,
      optimization_evaluations,
      fp_converged,
      fp_iterations,
      method,
      se_type
    )

------------------------------------------------------------------------

### Method `compute_elasticities()`

Compute own-price elasticities for a specific market

#### Usage

    BLPResults$compute_elasticities(market_id = NULL)

#### Arguments

- `market_id`:

  Market identifier

#### Returns

J x J elasticity matrix

------------------------------------------------------------------------

### Method `compute_diversion_ratios()`

Compute diversion ratios for a specific market

#### Usage

    BLPResults$compute_diversion_ratios(market_id = NULL)

#### Arguments

- `market_id`:

  Market identifier

#### Returns

J x J diversion ratio matrix

------------------------------------------------------------------------

### Method `compute_costs()`

Extract marginal costs

#### Usage

    BLPResults$compute_costs()

#### Returns

Named list with costs per market, or full vector

------------------------------------------------------------------------

### Method `compute_markups()`

Compute markups (p - c) / p

#### Usage

    BLPResults$compute_markups()

#### Returns

Markup vector

------------------------------------------------------------------------

### Method `compute_consumer_surplus()`

Compute consumer surplus per market

#### Usage

    BLPResults$compute_consumer_surplus()

#### Returns

Named numeric vector of CS per market

------------------------------------------------------------------------

### Method `compute_hhi()`

Compute HHI per market

#### Usage

    BLPResults$compute_hhi()

#### Returns

Named numeric vector of HHI per market

------------------------------------------------------------------------

### Method `compute_merger()`

Simulate merger and compute new equilibrium prices

#### Usage

    BLPResults$compute_merger(new_firm_ids, iteration = NULL, costs = NULL)

#### Arguments

- `new_firm_ids`:

  Updated firm ownership (length N)

- `iteration`:

  BLPIteration for price equilibrium

- `costs`:

  Optional pre-computed costs; computed if NULL

#### Returns

List with new_prices, new_shares, new_costs, delta_cs

------------------------------------------------------------------------

### Method `run_hansen_test()`

Run Hansen's J-test for overidentifying restrictions

#### Usage

    BLPResults$run_hansen_test()

#### Returns

List with statistic, df, p_value

------------------------------------------------------------------------

### Method `run_wald_test()`

Run Wald test for parameter restrictions

#### Usage

    BLPResults$run_wald_test(R, r = NULL)

#### Arguments

- `R`:

  Restriction matrix (q x p)

- `r`:

  Restriction values (q x 1), default 0

#### Returns

List with statistic, df, p_value

------------------------------------------------------------------------

### Method `sigma_squared()`

Extract sigma squared (Sigma %\*% Sigma')

#### Usage

    BLPResults$sigma_squared()

#### Returns

K2 x K2 covariance matrix

------------------------------------------------------------------------

### Method `summary_table()`

Get a summary table of all estimated parameters

#### Usage

    BLPResults$summary_table()

#### Returns

Data frame with estimates and standard errors

------------------------------------------------------------------------

### Method [`print()`](https://rdrr.io/r/base/print.html)

Print estimation results

#### Usage

    BLPResults$print(...)

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    BLPResults$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

# BLP Simulation Results

Results from a BLP simulation equilibrium computation.

## Public fields

- `simulation`:

  The originating BLPSimulation

- `product_data`:

  Updated product data with equilibrium prices/shares

- `delta`:

  Mean utilities at equilibrium

- `costs`:

  Marginal costs

- `prices`:

  Equilibrium prices

- `shares`:

  Equilibrium shares

## Methods

### Public methods

- [`BLPSimulationResults$new()`](#method-BLPSimulationResults-new)

- [`BLPSimulationResults$to_problem()`](#method-BLPSimulationResults-to_problem)

- [`BLPSimulationResults$print()`](#method-BLPSimulationResults-print)

- [`BLPSimulationResults$clone()`](#method-BLPSimulationResults-clone)

------------------------------------------------------------------------

### Method `new()`

Create simulation results

#### Usage

    BLPSimulationResults$new(
      simulation,
      product_data,
      delta,
      costs,
      prices,
      shares
    )

------------------------------------------------------------------------

### Method `to_problem()`

Convert to a BLPProblem for estimation

#### Usage

    BLPSimulationResults$to_problem(
      product_formulations = NULL,
      add_instruments = TRUE
    )

#### Arguments

- `product_formulations`:

  Optional override formulations

- `add_instruments`:

  Whether to auto-add BLP instruments

#### Returns

A BLPProblem object

------------------------------------------------------------------------

### Method [`print()`](https://rdrr.io/r/base/print.html)

Print simulation results

#### Usage

    BLPSimulationResults$print(...)

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    BLPSimulationResults$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

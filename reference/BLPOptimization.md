# Optimization Configuration

Configuration for the GMM outer-loop optimization.

## Methods

### Public methods

- [`BLPOptimization$new()`](#method-BLPOptimization-new)

- [`BLPOptimization$optimize()`](#method-BLPOptimization-optimize)

- [`BLPOptimization$print()`](#method-BLPOptimization-print)

- [`BLPOptimization$clone()`](#method-BLPOptimization-clone)

------------------------------------------------------------------------

### Method `new()`

Create optimization configuration

#### Usage

    BLPOptimization$new(
      method = "l-bfgs-b",
      method_options = list(),
      compute_gradient = TRUE
    )

#### Arguments

- `method`:

  Character: optimizer name

- `method_options`:

  Named list passed to the optimizer

- `compute_gradient`:

  Logical: use analytic gradients?

------------------------------------------------------------------------

### Method [`optimize()`](https://rdrr.io/r/stats/optimize.html)

Run optimization

#### Usage

    BLPOptimization$optimize(initial, bounds, objective_function)

#### Arguments

- `initial`:

  Starting parameter values

- `bounds`:

  List with lower and upper bound vectors

- `objective_function`:

  Function: theta -\> list(objective, gradient)

#### Returns

List with values, converged, iterations, evaluations

------------------------------------------------------------------------

### Method [`print()`](https://rdrr.io/r/base/print.html)

Print optimization configuration

#### Usage

    BLPOptimization$print(...)

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    BLPOptimization$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

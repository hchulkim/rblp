# Micro Moment

Matches a micro data target to a model prediction.

## Public fields

- `name`:

  Moment name

- `value`:

  Observed (target) value from micro data

- `parts`:

  List of MicroPart objects

- `compute_value_fn`:

  Function of part values -\> scalar

- `compute_gradient_fn`:

  Function of part values -\> gradient vector

## Methods

### Public methods

- [`MicroMoment$new()`](#method-MicroMoment-new)

- [`MicroMoment$compute_simulated_value()`](#method-MicroMoment-compute_simulated_value)

- [`MicroMoment$print()`](#method-MicroMoment-print)

- [`MicroMoment$clone()`](#method-MicroMoment-clone)

------------------------------------------------------------------------

### Method `new()`

Create micro moment

#### Usage

    MicroMoment$new(
      name,
      value,
      parts,
      compute_value = NULL,
      compute_gradient = NULL
    )

#### Arguments

- `name`:

  Moment name

- `value`:

  Observed scalar target

- `parts`:

  Single MicroPart or list of MicroPart objects

- `compute_value`:

  Function of part values -\> scalar (default: identity)

- `compute_gradient`:

  Function of part values -\> gradient (default: 1)

------------------------------------------------------------------------

### Method `compute_simulated_value()`

Compute simulated moment value

#### Usage

    MicroMoment$compute_simulated_value(
      economy,
      delta,
      sigma = NULL,
      pi = NULL,
      rho = NULL
    )

#### Arguments

- `economy`:

  BLPEconomy object

- `delta`:

  Mean utilities

- `sigma`:

  Sigma matrix

- `pi`:

  Pi matrix

- `rho`:

  Rho vector

#### Returns

Scalar simulated value

------------------------------------------------------------------------

### Method [`print()`](https://rdrr.io/r/base/print.html)

Print the moment

#### Usage

    MicroMoment$print(...)

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    MicroMoment$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

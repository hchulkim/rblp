# BLP Parameter Manager

Manages BLP model parameter compression, expansion, bounds, and labels.

## Methods

### Public methods

- [`BLPParameters$new()`](#method-BLPParameters-new)

- [`BLPParameters$compress()`](#method-BLPParameters-compress)

- [`BLPParameters$expand()`](#method-BLPParameters-expand)

- [`BLPParameters$get_bounds()`](#method-BLPParameters-get_bounds)

- [`BLPParameters$get_labels()`](#method-BLPParameters-get_labels)

- [`BLPParameters$n_free()`](#method-BLPParameters-n_free)

- [`BLPParameters$print()`](#method-BLPParameters-print)

- [`BLPParameters$clone()`](#method-BLPParameters-clone)

------------------------------------------------------------------------

### Method `new()`

#### Usage

    BLPParameters$new(
      sigma = NULL,
      pi = NULL,
      rho = NULL,
      beta = NULL,
      gamma = NULL,
      sigma_bounds = NULL,
      pi_bounds = NULL,
      rho_bounds = NULL,
      beta_bounds = NULL,
      gamma_bounds = NULL,
      rc_types = NULL
    )

------------------------------------------------------------------------

### Method `compress()`

#### Usage

    BLPParameters$compress(
      sigma = NULL,
      pi = NULL,
      rho = NULL,
      beta = NULL,
      gamma = NULL
    )

------------------------------------------------------------------------

### Method `expand()`

#### Usage

    BLPParameters$expand(theta)

------------------------------------------------------------------------

### Method `get_bounds()`

#### Usage

    BLPParameters$get_bounds()

------------------------------------------------------------------------

### Method `get_labels()`

#### Usage

    BLPParameters$get_labels()

------------------------------------------------------------------------

### Method `n_free()`

#### Usage

    BLPParameters$n_free()

------------------------------------------------------------------------

### Method [`print()`](https://rdrr.io/r/base/print.html)

#### Usage

    BLPParameters$print(...)

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    BLPParameters$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

# Integration Configuration

Configuration for building integration nodes and weights for
approximating integrals over the distribution of random coefficients.

## Methods

### Public methods

- [`BLPIntegration$new()`](#method-BLPIntegration-new)

- [`BLPIntegration$build()`](#method-BLPIntegration-build)

- [`BLPIntegration$print()`](#method-BLPIntegration-print)

- [`BLPIntegration$clone()`](#method-BLPIntegration-clone)

------------------------------------------------------------------------

### Method `new()`

Create integration configuration

#### Usage

    BLPIntegration$new(specification, size, seed = NULL)

#### Arguments

- `specification`:

  Character: method name

- `size`:

  Integer: number of draws or quadrature level

- `seed`:

  Optional RNG seed

------------------------------------------------------------------------

### Method `build()`

Build nodes and weights

#### Usage

    BLPIntegration$build(dimensions)

#### Arguments

- `dimensions`:

  Number of integration dimensions (K2)

#### Returns

List with `nodes` (N x K2) and `weights` (N x 1)

------------------------------------------------------------------------

### Method [`print()`](https://rdrr.io/r/base/print.html)

Print integration configuration

#### Usage

    BLPIntegration$print(...)

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    BLPIntegration$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

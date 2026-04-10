# Micro Part

A single component of a micro moment calculation.

## Public fields

- `name`:

  Part name

- `dataset`:

  MicroDataset object

- `compute_values`:

  Function(market_id, products, agents) -\> I x J values

## Methods

### Public methods

- [`MicroPart$new()`](#method-MicroPart-new)

- [`MicroPart$compute()`](#method-MicroPart-compute)

- [`MicroPart$print()`](#method-MicroPart-print)

- [`MicroPart$clone()`](#method-MicroPart-clone)

------------------------------------------------------------------------

### Method `new()`

Create a micro part

#### Usage

    MicroPart$new(name, dataset, compute_values)

#### Arguments

- `name`:

  Part name

- `dataset`:

  MicroDataset object

- `compute_values`:

  Function returning agent x product value matrix

------------------------------------------------------------------------

### Method `compute()`

Compute aggregated part value across markets

#### Usage

    MicroPart$compute(economy, delta, sigma = NULL, pi = NULL, rho = NULL)

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

Scalar aggregated value

------------------------------------------------------------------------

### Method [`print()`](https://rdrr.io/r/base/print.html)

Print the part

#### Usage

    MicroPart$print(...)

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    MicroPart$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

# BLP Model Formulation

Configuration for designing model matrices from R formulas.

## Methods

### Public methods

- [`BLPFormulation$new()`](#method-BLPFormulation-new)

- [`BLPFormulation$build_matrix()`](#method-BLPFormulation-build_matrix)

- [`BLPFormulation$get_names()`](#method-BLPFormulation-get_names)

- [`BLPFormulation$has_prices()`](#method-BLPFormulation-has_prices)

- [`BLPFormulation$has_shares()`](#method-BLPFormulation-has_shares)

- [`BLPFormulation$get_absorb()`](#method-BLPFormulation-get_absorb)

- [`BLPFormulation$print()`](#method-BLPFormulation-print)

- [`BLPFormulation$clone()`](#method-BLPFormulation-clone)

------------------------------------------------------------------------

### Method `new()`

Create a new formulation

#### Usage

    BLPFormulation$new(formula, absorb = NULL)

#### Arguments

- `formula`:

  Formula string or R formula

- `absorb`:

  Optional formula for absorbed fixed effects

------------------------------------------------------------------------

### Method `build_matrix()`

Build design matrix from data

#### Usage

    BLPFormulation$build_matrix(data)

#### Arguments

- `data`:

  Data frame

#### Returns

Matrix with named columns

------------------------------------------------------------------------

### Method `get_names()`

Get column names

#### Usage

    BLPFormulation$get_names()

#### Returns

Character vector

------------------------------------------------------------------------

### Method `has_prices()`

Check if prices are in the formulation

#### Usage

    BLPFormulation$has_prices()

------------------------------------------------------------------------

### Method `has_shares()`

Check if shares are in the formulation

#### Usage

    BLPFormulation$has_shares()

------------------------------------------------------------------------

### Method `get_absorb()`

Get absorb formula (if any)

#### Usage

    BLPFormulation$get_absorb()

------------------------------------------------------------------------

### Method [`print()`](https://rdrr.io/r/base/print.html)

Print formulation

#### Usage

    BLPFormulation$print(...)

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    BLPFormulation$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

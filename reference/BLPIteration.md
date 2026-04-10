# Iteration Configuration

Configuration for the BLP contraction mapping fixed-point iteration.

## Methods

### Public methods

- [`BLPIteration$new()`](#method-BLPIteration-new)

- [`BLPIteration$iterate()`](#method-BLPIteration-iterate)

- [`BLPIteration$print()`](#method-BLPIteration-print)

- [`BLPIteration$clone()`](#method-BLPIteration-clone)

------------------------------------------------------------------------

### Method `new()`

Create iteration configuration

#### Usage

    BLPIteration$new(method = "squarem", method_options = list())

#### Arguments

- `method`:

  Character: "simple", "squarem", or "return"

- `method_options`:

  Named list of options

------------------------------------------------------------------------

### Method `iterate()`

Run fixed-point iteration

#### Usage

    BLPIteration$iterate(initial, contraction)

#### Arguments

- `initial`:

  Starting values

- `contraction`:

  Function: x -\> list(x_new, weights) or just x_new

#### Returns

List with values, converged, iterations, evaluations

------------------------------------------------------------------------

### Method [`print()`](https://rdrr.io/r/base/print.html)

Print iteration configuration

#### Usage

    BLPIteration$print(...)

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    BLPIteration$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

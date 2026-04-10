# Micro Dataset

Defines a micro-level dataset for micro moments.

## Public fields

- `name`:

  Dataset name

- `observations`:

  Number of observations in micro data

- `compute_weights`:

  Function(market_id, products, agents) -\> I x J weights

- `market_ids`:

  Optional subset of markets

## Methods

### Public methods

- [`MicroDataset$new()`](#method-MicroDataset-new)

- [`MicroDataset$print()`](#method-MicroDataset-print)

- [`MicroDataset$clone()`](#method-MicroDataset-clone)

------------------------------------------------------------------------

### Method `new()`

Create micro dataset

#### Usage

    MicroDataset$new(name, observations, compute_weights, market_ids = NULL)

#### Arguments

- `name`:

  Dataset name

- `observations`:

  Number of micro observations

- `compute_weights`:

  Function returning agent x product weight matrix

- `market_ids`:

  Optional market subset

------------------------------------------------------------------------

### Method [`print()`](https://rdrr.io/r/base/print.html)

Print the dataset

#### Usage

    MicroDataset$print(...)

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    MicroDataset$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

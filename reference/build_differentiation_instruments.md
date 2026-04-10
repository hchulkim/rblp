# Build Differentiation Instruments (Gandhi & Houde 2020)

Constructs instruments based on product proximity in characteristic
space.

## Usage

``` r
build_differentiation_instruments(
  X,
  market_ids,
  firm_ids,
  method = "local",
  interact = NULL
)
```

## Arguments

- X:

  Matrix of product characteristics (N x K)

- market_ids:

  Market identifiers (length N)

- firm_ids:

  Firm identifiers (length N)

- method:

  "local" (count nearby products) or "quadratic" (sum squared distances)

- interact:

  Optional matrix of interaction characteristics (N x L)

## Value

Matrix of instruments

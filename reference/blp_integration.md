# Create integration configuration

Create integration configuration

## Usage

``` r
blp_integration(specification, size, seed = NULL)
```

## Arguments

- specification:

  Character: "monte_carlo", "halton", "lhs", "mlhs", "product", or
  "grid"

- size:

  Integer: number of draws or quadrature level

- seed:

  Optional RNG seed for reproducibility

## Value

A BLPIntegration object

## Examples

``` r
int_mc <- blp_integration("monte_carlo", 50, seed = 42)
int_prod <- blp_integration("product", 7)
```

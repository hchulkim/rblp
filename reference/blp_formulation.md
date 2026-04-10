# Create a BLP formulation

Create a BLP formulation

## Usage

``` r
blp_formulation(formula, absorb = NULL)
```

## Arguments

- formula:

  Formula string (e.g., "1 + prices + sugar") or R formula

- absorb:

  Optional formula for absorbed fixed effects (e.g., "~ product_id")

## Value

A BLPFormulation object

## Examples

``` r
f1 <- blp_formulation("1 + prices + sugar + mushy")
f2 <- blp_formulation("0 + prices", absorb = "~ product_ids")
```

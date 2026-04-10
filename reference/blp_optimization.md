# Create optimization configuration

Create optimization configuration

## Usage

``` r
blp_optimization(
  method = "l-bfgs-b",
  method_options = list(),
  compute_gradient = TRUE
)
```

## Arguments

- method:

  Character: "l-bfgs-b" (default), "bfgs", "nelder-mead", "nlminb", or
  "return"

- method_options:

  Named list of options passed to the optimizer

- compute_gradient:

  Logical: use analytic gradients? (default TRUE)

## Value

A BLPOptimization object

## Examples

``` r
opt <- blp_optimization("l-bfgs-b")
```

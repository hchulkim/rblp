# Create iteration configuration

Create iteration configuration

## Usage

``` r
blp_iteration(method = "squarem", method_options = list())
```

## Arguments

- method:

  Character: "simple", "squarem", or "return"

- method_options:

  Named list. Key options:

  max_evaluations

  :   Maximum contraction evaluations (default 5000)

  atol

  :   Absolute convergence tolerance (default 1e-14)

  rtol

  :   Relative convergence tolerance (default 0)

## Value

A BLPIteration object

## Examples

``` r
iter <- blp_iteration("squarem", list(atol = 1e-14))
```

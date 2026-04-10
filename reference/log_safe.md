# Safe logarithm with underflow protection

Safe logarithm with underflow protection

## Usage

``` r
log_safe(x, min_val = .Machine$double.xmin)
```

## Arguments

- x:

  Numeric vector

- min_val:

  Minimum value to clamp to

## Value

log(pmax(x, min_val))

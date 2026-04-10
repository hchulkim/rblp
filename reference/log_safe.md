# Safe logarithm with underflow protection

Safe logarithm with underflow protection

## Usage

``` r
log_safe(x, min_val = 9.99999999999999e-301)
```

## Arguments

- x:

  Numeric vector

- min_val:

  Minimum value to clamp to (default 1e-300)

## Value

log(pmax(x, min_val))

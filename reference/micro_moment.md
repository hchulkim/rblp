# Create a Micro Moment

Create a Micro Moment

## Usage

``` r
micro_moment(name, value, parts, compute_value = NULL, compute_gradient = NULL)
```

## Arguments

- name:

  Moment name

- value:

  Observed target value

- parts:

  Single MicroPart or list of MicroPart objects

- compute_value:

  Function of part values -\> scalar

- compute_gradient:

  Function of part values -\> gradient

## Value

A MicroMoment object

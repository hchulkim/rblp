# Get or Set rblp Package Options

Get or Set rblp Package Options

## Usage

``` r
rblp_options(...)
```

## Arguments

- ...:

  Named arguments to set options. If empty, returns all current options.

## Value

Invisibly returns the previous values of changed options, or all options
if none changed.

## Examples

``` r
rblp_options(verbose = FALSE)
rblp_options()
#> $digits
#> [1] 7
#> 
#> $verbose
#> [1] FALSE
#> 
#> $verbose_output
#> [1] ""
#> 
#> $pseudo_inverses
#> [1] TRUE
#> 
#> $collinear_atol
#> [1] 1e-10
#> 
#> $collinear_rtol
#> [1] 1e-10
#> 
#> $psd_atol
#> [1] 1e-08
#> 
#> $psd_rtol
#> [1] 1e-08
#> 
#> $finite_differences_epsilon
#> [1] 1.490116e-08
#> 
#> $weights_tol
#> [1] 1e-10
#> 
#> $micro_computation_chunks
#> [1] 1
#> 
#> $num_processes
#> [1] 1
#> 
```

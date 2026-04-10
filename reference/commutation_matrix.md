# Commutation matrix

Returns K_m,n such that K %\*% vec(A) = vec(t(A)) for m x n matrix A.

## Usage

``` r
commutation_matrix(m, n = m)
```

## Arguments

- m:

  Number of rows

- n:

  Number of columns (default = m)

## Value

Commutation matrix (m*n x m*n)

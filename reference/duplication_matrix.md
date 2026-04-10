# Duplication matrix

Returns D_n such that D_n %\*% vech(A) = vec(A) for symmetric n x n
matrix A.

## Usage

``` r
duplication_matrix(n)
```

## Arguments

- n:

  Matrix dimension

## Value

Duplication matrix (n^2 x n\*(n+1)/2)

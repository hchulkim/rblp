# Elimination matrix

Returns L_n such that L_n %\*% vec(A) = vech(A) for symmetric n x n
matrix A.

## Usage

``` r
elimination_matrix(n)
```

## Arguments

- n:

  Matrix dimension

## Value

Elimination matrix (n\*(n+1)/2 x n^2)

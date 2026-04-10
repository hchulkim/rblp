#' Matrix inversion with cascading fallbacks
#' @param x Square matrix
#' @return List with \code{inverse} and \code{replacement} (NULL if exact, string if fallback used)
#' @keywords internal
# Cascading inversion with graceful degradation. In BLP estimation, matrices
# like (X'Z W Z'X) or (G'WG) can become near-singular when instruments are
# weak, parameters are at boundary values, or the model is poorly identified.
# Rather than failing outright, this function tries progressively coarser
# inversions: (1) exact solve or Moore-Penrose pseudo-inverse (ginv), which
# handles rank-deficient matrices by zeroing out small singular values;
# (2) diagonal approximation (inverting only the diagonal entries) as a last
# resort, which preserves the scale of each parameter's variance but ignores
# all cross-parameter correlations.
approximately_invert <- function(x) {
  if (isTRUE(getOption("rblp.pseudo_inverses", TRUE))) {
    # When pseudo_inverses mode is enabled, go directly to ginv (Moore-Penrose),
    # which is robust to singularity at the cost of slightly more computation.
    inv <- tryCatch(MASS::ginv(x), error = function(e) NULL)
    if (!is.null(inv)) return(list(inverse = inv, replacement = NULL))
    return(list(inverse = diag(1 / diag(x)), replacement = "diagonal"))
  }
  # Standard mode: try exact inversion first (fastest when the matrix is
  # well-conditioned), then fall back to ginv, then to diagonal.
  inv <- tryCatch(solve(x), error = function(e) NULL)
  if (!is.null(inv)) return(list(inverse = inv, replacement = NULL))
  inv <- tryCatch(MASS::ginv(x), error = function(e) NULL)
  if (!is.null(inv)) return(list(inverse = inv, replacement = "pseudo-inverse"))
  list(inverse = diag(1 / diag(x)), replacement = "diagonal")
}

#' Solve linear system with fallback
#' @param a Square matrix
#' @param b Right-hand side vector or matrix
#' @return Solution x such that a %*% x = b
#' @keywords internal
# Solve a linear system a %*% x = b with fallback. Prefers the numerically
# stable solve(a, b) (which uses LU decomposition without explicitly forming
# the inverse), but falls back to approximately_invert if the system is
# singular or ill-conditioned.
approximately_solve <- function(a, b) {
  x <- tryCatch(solve(a, b), error = function(e) NULL)
  if (!is.null(x)) return(x)
  inv <- approximately_invert(a)
  inv$inverse %*% b
}

#' Detect near-singularity
#' @param x Square matrix
#' @return Logical: TRUE if matrix appears singular
#' @keywords internal
# Check whether a matrix is effectively singular by computing its reciprocal
# condition number (rcond). A value near zero indicates the matrix is close
# to rank-deficient, which would make IV estimates unreliable due to weak
# identification or multicollinearity among instruments.
detect_singularity <- function(x) {
  cond <- tryCatch(rcond(x), error = function(e) 0)
  cond < getOption("rblp.collinear_atol", 1e-10)
}

#' Detect collinearity via QR decomposition
#' @param x Matrix
#' @return List with \code{detected} (logical) and \code{indices} (which columns to drop)
#' @keywords internal
# Detect collinearity among columns (e.g., instruments or regressors) using
# QR decomposition with column pivoting. The pivoted QR reorders columns by
# their contribution to the column space; columns beyond the numerical rank
# are redundant and should be dropped to avoid singularity in the IV/GMM
# estimation. Returns the indices of the linearly dependent columns.
detect_collinearity <- function(x) {
  qr_obj <- qr(x)
  rank <- qr_obj$rank
  p <- ncol(x)
  if (rank == p) return(list(detected = FALSE, indices = integer(0)))
  keep <- qr_obj$pivot[seq_len(rank)]
  drop <- setdiff(seq_len(p), keep)
  list(detected = TRUE, indices = drop)
}

#' Duplication matrix
#'
#' Returns D_n such that D_n %*% vech(A) = vec(A) for symmetric n x n matrix A.
#' @param n Matrix dimension
#' @return Duplication matrix (n^2 x n*(n+1)/2)
#' @keywords internal
# Duplication matrix D_n: maps the half-vectorization vech(A) of a symmetric
# n x n matrix (which stores only the n(n+1)/2 unique lower-triangular entries)
# back to the full vectorization vec(A) of size n^2. This is used when
# parameterizing the Cholesky factor of the random coefficients covariance
# matrix Sigma: the lower-triangular elements are the free parameters, and
# D_n reconstructs the full symmetric matrix for computing derivatives.
duplication_matrix <- function(n) {
  m <- n * (n + 1) / 2
  D <- matrix(0, n^2, m)
  col <- 0
  for (j in seq_len(n)) {
    for (i in j:n) {
      col <- col + 1
      D[(j - 1) * n + i, col] <- 1
      if (i != j) {
        D[(i - 1) * n + j, col] <- 1
      }
    }
  }
  D
}

#' Elimination matrix
#'
#' Returns L_n such that L_n %*% vec(A) = vech(A) for symmetric n x n matrix A.
#' @param n Matrix dimension
#' @return Elimination matrix (n*(n+1)/2 x n^2)
#' @keywords internal
# Elimination matrix L_n: the left-inverse of the duplication matrix.
# Extracts the unique lower-triangular elements from vec(A), so that
# L_n %*% vec(A) = vech(A). Used in derivative computations to go from
# the full n^2 vectorized derivative to the n(n+1)/2 free parameters.
elimination_matrix <- function(n) {
  m <- n * (n + 1) / 2
  L <- matrix(0, m, n^2)
  row <- 0
  for (j in seq_len(n)) {
    for (i in j:n) {
      row <- row + 1
      L[row, (j - 1) * n + i] <- 1
    }
  }
  L
}

#' Commutation matrix
#'
#' Returns K such that K \%*\% vec(A) = vec(t(A)) for m x n matrix A.
#' @param m Number of rows
#' @param n Number of columns (default = m)
#' @return Commutation matrix (m*n x m*n)
#' @keywords internal
# Commutation matrix K_{m,n}: a permutation matrix such that K %*% vec(A) = vec(A')
# for an m x n matrix A. This arises in matrix calculus when transposing
# vectorized expressions, e.g., in computing the derivative of a symmetric
# matrix function with respect to its Cholesky factor.
commutation_matrix <- function(m, n = m) {
  K <- matrix(0, m * n, m * n)
  for (i in seq_len(m)) {
    for (j in seq_len(n)) {
      K[(i - 1) * n + j, (j - 1) * m + i] <- 1
    }
  }
  K
}

#' Vectorize a matrix (column-major)
#' @param x Matrix
#' @return Numeric vector
#' @keywords internal
mat_vec <- function(x) {
  as.numeric(x)
}

#' Half-vectorize a symmetric matrix (lower triangle, column-major)
#' @param x Symmetric matrix
#' @return Numeric vector
#' @keywords internal
mat_vech <- function(x) {
  x[lower.tri(x, diag = TRUE)]
}

#' Reconstruct lower triangular matrix from vech
#' @param v Numeric vector (half-vectorization)
#' @param n Matrix dimension
#' @return Lower triangular matrix
#' @keywords internal
vech_to_lower <- function(v, n) {
  L <- matrix(0, n, n)
  L[lower.tri(L, diag = TRUE)] <- v
  L
}

#' Reconstruct full symmetric matrix from vech
#' @param v Numeric vector (half-vectorization)
#' @param n Matrix dimension
#' @return Symmetric matrix
#' @keywords internal
vech_to_full <- function(v, n) {
  L <- vech_to_lower(v, n)
  L + t(L) - diag(diag(L))
}

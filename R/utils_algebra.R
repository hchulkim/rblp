#' Matrix inversion with cascading fallbacks
#' @param x Square matrix
#' @return List with \code{inverse} and \code{replacement} (NULL if exact, string if fallback used)
#' @keywords internal
approximately_invert <- function(x) {
  if (isTRUE(getOption("rblp.pseudo_inverses", TRUE))) {
    inv <- tryCatch(MASS::ginv(x), error = function(e) NULL)
    if (!is.null(inv)) return(list(inverse = inv, replacement = NULL))
    return(list(inverse = diag(1 / diag(x)), replacement = "diagonal"))
  }
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
detect_singularity <- function(x) {
  cond <- tryCatch(rcond(x), error = function(e) 0)
  cond < getOption("rblp.collinear_atol", 1e-10)
}

#' Detect collinearity via QR decomposition
#' @param x Matrix
#' @return List with \code{detected} (logical) and \code{indices} (which columns to drop)
#' @keywords internal
detect_collinearity <- function(x) {
  qr_obj <- qr(x)
  rank <- qr_obj$rank
  p <- ncol(x)
  if (rank == p) return(list(detected = FALSE, indices = integer(0)))
  # Columns to drop: those not in the pivot
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
duplication_matrix <- function(n) {
  m <- n * (n + 1) / 2
  D <- matrix(0, n^2, m)
  col <- 0
  for (j in seq_len(n)) {
    for (i in j:n) {
      col <- col + 1
      # (i,j) element -> row (j-1)*n + i
      D[(j - 1) * n + i, col] <- 1
      if (i != j) {
        # (j,i) element -> row (i-1)*n + j
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
#' Returns K_{m,n} such that K %*% vec(A) = vec(t(A)) for m x n matrix A.
#' @param m Number of rows
#' @param n Number of columns (default = m)
#' @return Commutation matrix (m*n x m*n)
#' @keywords internal
commutation_matrix <- function(m, n = m) {
  K <- matrix(0, m * n, m * n)
  for (i in seq_len(m)) {
    for (j in seq_len(n)) {
      # vec(A) element at position (j-1)*m + i maps to vec(A') at (i-1)*n + j
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

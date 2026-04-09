#' @keywords internal
rblp_message <- function(...) {
  if (isTRUE(getOption("rblp.verbose", TRUE))) {
    message(...)
  }
}

#' Sum within groups
#' @param x Numeric vector or matrix
#' @param ids Factor or character vector of group identifiers
#' @return Named vector or matrix of group sums
#' @keywords internal
groups_sum <- function(x, ids) {
  if (is.matrix(x)) {
    uid <- unique(ids)
    out <- matrix(0, length(uid), ncol(x))
    for (i in seq_along(uid)) {
      idx <- which(ids == uid[i])
      out[i, ] <- colSums(x[idx, , drop = FALSE])
    }
    rownames(out) <- uid
    return(out)
  }
  tapply(x, ids, sum)
}

#' Mean within groups
#' @param x Numeric vector or matrix
#' @param ids Group identifiers
#' @return Named vector or matrix of group means
#' @keywords internal
groups_mean <- function(x, ids) {
  if (is.matrix(x)) {
    uid <- unique(ids)
    out <- matrix(0, length(uid), ncol(x))
    for (i in seq_along(uid)) {
      idx <- which(ids == uid[i])
      out[i, ] <- colMeans(x[idx, , drop = FALSE])
    }
    rownames(out) <- uid
    return(out)
  }
  tapply(x, ids, mean)
}

#' Expand group-level values to observation level
#' @param x_group Named vector of group-level values
#' @param ids Group identifiers for each observation
#' @return Vector expanded to observation level
#' @keywords internal
groups_expand <- function(x_group, ids) {
  uid <- names(x_group)
  if (is.null(uid)) uid <- unique(ids)
  out <- numeric(length(ids))
  for (i in seq_along(uid)) {
    idx <- which(ids == uid[i])
    out[idx] <- x_group[i]
  }
  out
}

#' Extract numbered columns from data frame into matrix
#' @param data Data frame
#' @param prefix Column name prefix (e.g., "demand_instruments")
#' @return Matrix of matched columns, or NULL if no match
#' @keywords internal
extract_columns <- function(data, prefix) {
  pattern <- paste0("^", prefix, "\\d+$")
  cols <- grep(pattern, names(data), value = TRUE)
  if (length(cols) == 0) return(NULL)
  # Sort by numeric suffix
  nums <- as.integer(sub(paste0("^", prefix), "", cols))
  cols <- cols[order(nums)]
  as.matrix(data[, cols, drop = FALSE])
}

#' Numerically stable exponentiation
#' @param x Numeric vector or matrix
#' @return exp(x) computed with overflow protection
#' @keywords internal
safe_exp <- function(x) {
  if (is.matrix(x)) {
    col_max <- apply(x, 2, max)
    return(exp(sweep(x, 2, col_max)) )
  }
  x_max <- max(x)
  exp(x - x_max)
}

#' Safe logarithm with underflow protection
#' @param x Numeric vector
#' @param min_val Minimum value to clamp to (default 1e-300)
#' @return log(pmax(x, min_val))
#' @keywords internal
log_safe <- function(x, min_val = 1e-300) {
  log(pmax(x, min_val))
}

#' Compute central finite differences
#' @param f Function taking numeric vector
#' @param x Point at which to evaluate
#' @param epsilon Step size
#' @return Numeric vector of partial derivatives
#' @keywords internal
compute_finite_differences <- function(f, x, epsilon = NULL) {
  if (is.null(epsilon)) epsilon <- getOption("rblp.finite_differences_epsilon",
                                              sqrt(.Machine$double.eps))
  n <- length(x)
  grad <- numeric(n)
  for (i in seq_len(n)) {
    x_plus <- x_minus <- x
    x_plus[i] <- x[i] + epsilon
    x_minus[i] <- x[i] - epsilon
    grad[i] <- (f(x_plus) - f(x_minus)) / (2 * epsilon)
  }
  grad
}

#' Format seconds into human-readable string
#' @param seconds Numeric
#' @return Character string
#' @keywords internal
format_seconds <- function(seconds) {
  if (seconds < 1) return(sprintf("%.0fms", seconds * 1000))
  if (seconds < 60) return(sprintf("%.1fs", seconds))
  if (seconds < 3600) return(sprintf("%.1fm", seconds / 60))
  sprintf("%.1fh", seconds / 3600)
}

#' Format number for display
#' @param x Numeric
#' @param digits Number of significant digits
#' @return Character string
#' @keywords internal
format_number <- function(x, digits = NULL) {
  if (is.null(digits)) digits <- getOption("rblp.digits", 7L)
  formatC(x, digits = digits, format = "g")
}

#' Build block-diagonal ownership matrix
#' @param firm_ids Firm IDs for each product
#' @param market_ids Market IDs for each product
#' @return Block-diagonal ownership matrix (N x N)
#' @keywords internal
build_ownership_matrix <- function(firm_ids, market_ids) {
  N <- length(firm_ids)
  ownership <- matrix(0, N, N)
  for (t in unique(market_ids)) {
    idx <- which(market_ids == t)
    fids <- firm_ids[idx]
    ownership[idx, idx] <- outer(fids, fids, "==") * 1.0
  }
  ownership
}

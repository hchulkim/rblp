# Conditional messaging: only emits diagnostic output when the rblp.verbose
# option is enabled. This controls progress reporting during estimation
# (e.g., contraction iterations, objective values) without cluttering output
# in production or testing contexts.
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
# Group-level summation: aggregates observations within each market (or other
# grouping variable). In BLP, many quantities are computed per-market --
# e.g., summing individual choice probabilities across consumers within a
# market to obtain predicted market shares, or summing moment contributions
# within clusters for clustered standard errors.
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
# Group-level means: computes within-group averages, used primarily in the
# FWL demeaning step (absorbing fixed effects) and for computing market-level
# averages of product characteristics.
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
# Expand group-level values back to observation level: the inverse of
# aggregation. For example, after computing a market-level statistic
# (like predicted outside-good share), this broadcasts it back to every
# product row belonging to that market.
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
# Extract numbered instrument columns from a data frame into a matrix.
# BLP/pyblp convention names instruments as "demand_instruments0",
# "demand_instruments1", ..., "demand_instrumentsN". This function matches
# all columns with the given prefix followed by a numeric suffix, sorts
# them in numeric order, and stacks them into the instrument matrix Z
# used in the IV/GMM estimation.
extract_columns <- function(data, prefix) {
  pattern <- paste0("^", prefix, "\\d+$")
  cols <- grep(pattern, names(data), value = TRUE)
  if (length(cols) == 0) return(NULL)
  nums <- as.integer(sub(paste0("^", prefix), "", cols))
  cols <- cols[order(nums)]
  as.matrix(data[, cols, drop = FALSE])
}

#' Numerically stable exponentiation
#' @param x Numeric vector or matrix
#' @return exp(x) computed with overflow protection
#' @keywords internal
# Numerically stable exponentiation using the log-sum-exp trick: subtract
# the column-wise maximum before exponentiating to prevent overflow (exp(700+)
# = Inf in double precision). The subtracted constant cancels in the logit
# choice probability formula s_ij = exp(V_ij) / sum_k exp(V_kj), so this
# transformation does not change the economic content. Without this, large
# mean utilities (common with product fixed effects) would cause overflow.
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
#' @param min_val Minimum value to clamp to
#' @return log(pmax(x, min_val))
#' @keywords internal
# Safe logarithm: clamps the argument away from zero to avoid log(0) = -Inf.
# This arises in the BLP contraction mapping where we compute log(s_observed)
# and log(s_predicted); if a predicted share is numerically zero (e.g., a
# product with very low utility), the log would produce -Inf and corrupt
# the iteration. The floor of 1e-300 is far below any meaningful share but
# keeps the computation finite.
log_safe <- function(x, min_val = .Machine$double.xmin) {
  log(pmax(x, min_val))
}

#' Compute central finite differences
#' @param f Function taking numeric vector
#' @param x Point at which to evaluate
#' @param epsilon Step size
#' @return Numeric vector of partial derivatives
#' @keywords internal
# Central finite differences for gradient verification. Approximates each
# partial derivative as [f(x + e_i*h) - f(x - e_i*h)] / (2h), which has
# O(h^2) truncation error (vs. O(h) for forward differences). The default
# step size h = sqrt(machine epsilon) ~ 1.5e-8 balances truncation error
# against floating-point cancellation error. Used to validate the analytic
# gradient of the BLP GMM objective during development and debugging.
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
# Build the block-diagonal ownership matrix O, where O[j,k] = 1 if products
# j and k are produced by the same firm in the same market, and 0 otherwise.
# This matrix enters the Bertrand-Nash first-order conditions for pricing:
#   p = mc - (O * dS/dp)^{-1} s
# where * is element-wise multiplication. Under a merger, firm_ids change
# to reflect the merged entity, which changes O and hence equilibrium prices.
# The block-diagonal structure ensures products in different markets do not
# interact.
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

#' Build BLP Instruments (Sums of Characteristics)
#'
#' Constructs excluded instruments using the traditional BLP approach:
#' sums of characteristics of other own-firm products and rival products.
#'
#' @param X Matrix of exogenous product characteristics (N x K)
#' @param market_ids Market identifiers (length N)
#' @param firm_ids Firm identifiers (length N)
#' @return Matrix of instruments (N x 2K): own-firm and rival columns
#' @export
#' @examples
#' \dontrun{
#' iv <- build_blp_instruments(X_exog, product_data$market_ids, product_data$firm_ids)
#' }
# Constructs the classical BLP (1995) excluded instruments for demand
# estimation. The identification argument: in a Nash-Bertrand equilibrium,
# product j's price depends on the characteristics of competing products
# (which determine the elasticity of demand j faces and hence the markup).
# But rival characteristics are exogenous to j's unobserved quality xi_j,
# providing valid instruments for the endogenous price.
#
# For each product j in market t, two types of instruments are built:
#   (1) Own-firm: sum of characteristics of OTHER products sold by j's firm.
#       These shift j's price because multi-product firms internalize
#       cannibalization -- more own products nearby -> lower markup.
#   (2) Rival: sum of characteristics of products sold by competing firms.
#       More/closer rivals -> more elastic demand -> lower markup -> lower price.
#
# Returns an N x 2K matrix: K own-firm columns + K rival columns.
build_blp_instruments <- function(X, market_ids, firm_ids) {
  X <- as.matrix(X)
  N <- nrow(X)
  K <- ncol(X)

  own_other <- matrix(0, N, K)
  rival <- matrix(0, N, K)

  markets <- unique(market_ids)

  for (t in markets) {
    idx <- which(market_ids == t)
    firms_t <- firm_ids[idx]
    X_t <- X[idx, , drop = FALSE]
    J_t <- length(idx)

    for (j in seq_len(J_t)) {
      f <- firms_t[j]
      # Sum characteristics of other products by the SAME firm (excluding j).
      # Relevance: more own-firm products in the market increase competitive
      # pressure within the firm's portfolio, affecting equilibrium pricing.
      own_idx <- which(firms_t == f)
      own_idx <- own_idx[own_idx != j]
      if (length(own_idx) > 0) {
        own_other[idx[j], ] <- colSums(X_t[own_idx, , drop = FALSE])
      }

      # Sum characteristics of all products by RIVAL firms.
      # Relevance: rival characteristics shift market-level competition,
      # affecting j's residual demand elasticity and markup.
      rival_idx <- which(firms_t != f)
      if (length(rival_idx) > 0) {
        rival[idx[j], ] <- colSums(X_t[rival_idx, , drop = FALSE])
      }
    }
  }

  x_names <- colnames(X) %||% paste0("x", seq_len(K))
  colnames(own_other) <- paste0("blp_own_", x_names)
  colnames(rival) <- paste0("blp_rival_", x_names)

  cbind(own_other, rival)
}

#' Build Differentiation Instruments (Gandhi & Houde 2020)
#'
#' Constructs instruments based on product proximity in characteristic space.
#'
#' @param X Matrix of product characteristics (N x K)
#' @param market_ids Market identifiers (length N)
#' @param firm_ids Firm identifiers (length N)
#' @param method "local" (count nearby products) or "quadratic" (sum squared distances)
#' @param interact Optional matrix of interaction characteristics (N x L)
#' @return Matrix of instruments
#' @export
# Constructs differentiation instruments following Gandhi & Houde (2020).
# Unlike BLP instruments (which sum characteristics), these measure how
# isolated or crowded a product is in characteristic space. The intuition:
# a product surrounded by close competitors faces more elastic demand and
# sets a lower markup. This proximity is a valid instrument because the
# location of rival products in characteristic space is exogenous to j's
# unobserved quality xi_j.
#
# Two methods are available:
#   "local": count of products within one SD of product j along each
#     dimension. Captures the NUMBER of close competitors (extensive margin).
#   "quadratic": sum of squared characteristic distances to all other
#     products. Captures both number and closeness (intensive margin).
#
# Both are computed separately for own-firm and rival products.
build_differentiation_instruments <- function(X, market_ids, firm_ids,
                                               method = "local",
                                               interact = NULL) {
  X <- as.matrix(X)
  N <- nrow(X)
  K <- ncol(X)
  markets <- unique(market_ids)

  # For the "local" method, compute the standard deviation of all pairwise
  # distances across markets for each characteristic. This SD serves as a
  # bandwidth: products within 1 SD are "nearby" in that dimension. Using
  # a global (not market-specific) SD prevents the instrument from being
  # contaminated by within-market composition, which could be endogenous.
  sds <- numeric(K)
  if (method == "local") {
    for (k in seq_len(K)) {
      all_dists <- numeric(0)
      for (t in markets) {
        idx <- which(market_ids == t)
        x_t <- X[idx, k]
        J_t <- length(x_t)
        if (J_t < 2) next
        for (j in seq_len(J_t)) {
          for (l in seq_len(J_t)) {
            if (j != l) all_dists <- c(all_dists, abs(x_t[j] - x_t[l]))
          }
        }
      }
      sds[k] <- stats::sd(all_dists)
    }
  }

  n_interact <- if (!is.null(interact)) ncol(as.matrix(interact)) else 0L
  n_cols <- if (n_interact > 0) K * (1 + n_interact) * 2 else K * 2
  own_iv <- matrix(0, N, K)
  rival_iv <- matrix(0, N, K)

  for (t in markets) {
    idx <- which(market_ids == t)
    firms_t <- firm_ids[idx]
    X_t <- X[idx, , drop = FALSE]
    J_t <- length(idx)

    for (j in seq_len(J_t)) {
      f <- firms_t[j]
      for (k in seq_len(K)) {
        # Own firm, other products
        own_idx <- which(firms_t == f)
        own_idx <- own_idx[own_idx != j]
        rival_idx <- which(firms_t != f)

        if (method == "local") {
          # Local: count how many other products are within 1 SD along
          # characteristic k. A high count means product j faces many
          # nearby substitutes along this dimension -> more competition.
          if (length(own_idx) > 0) {
            dists <- abs(X_t[own_idx, k] - X_t[j, k])
            own_iv[idx[j], k] <- sum(dists < sds[k])
          }
          if (length(rival_idx) > 0) {
            dists <- abs(X_t[rival_idx, k] - X_t[j, k])
            rival_iv[idx[j], k] <- sum(dists < sds[k])
          }
        } else {
          # Quadratic: sum of squared distances measures how differentiated
          # product j is from competitors. A LARGE value means j is far
          # from others -> more market power -> higher markup. A small
          # value means j is in a crowded part of characteristic space.
          if (length(own_idx) > 0) {
            dists <- (X_t[own_idx, k] - X_t[j, k])^2
            own_iv[idx[j], k] <- sum(dists)
          }
          if (length(rival_idx) > 0) {
            dists <- (X_t[rival_idx, k] - X_t[j, k])^2
            rival_iv[idx[j], k] <- sum(dists)
          }
        }
      }
    }
  }

  x_names <- colnames(X) %||% paste0("x", seq_len(K))
  prefix <- if (method == "local") "diff_local" else "diff_quad"
  colnames(own_iv) <- paste0(prefix, "_own_", x_names)
  colnames(rival_iv) <- paste0(prefix, "_rival_", x_names)

  cbind(own_iv, rival_iv)
}

#' Build Balanced ID Data
#'
#' Creates a balanced panel of market and firm IDs for simulation.
#'
#' @param T Number of markets
#' @param J Number of products per market
#' @param F Number of firms
#' @return Data frame with market_ids and firm_ids
#' @export
#
# Helper for Monte Carlo simulations: generates a balanced panel where
# every market has exactly J products distributed as evenly as possible
# across F firms. Products are assigned round-robin so that firm 1 gets
# ceil(J/F) products, firm 2 gets ceil(J/F), etc., until the remainder
# is exhausted. The resulting data frame is the skeleton onto which
# simulated characteristics, prices, and shares are attached.
build_id_data <- function(T, J, F) {
  products_per_firm <- rep(J %/% F, F)
  remainder <- J %% F
  if (remainder > 0) products_per_firm[seq_len(remainder)] <-
    products_per_firm[seq_len(remainder)] + 1L

  firm_ids_per_market <- rep(seq_len(F), times = products_per_firm)
  N <- T * J

  data.frame(
    market_ids = rep(seq_len(T), each = J),
    firm_ids = rep(firm_ids_per_market, T),
    stringsAsFactors = FALSE
  )
}

#' Build Custom Ownership Matrix
#'
#' Creates ownership matrices with custom specifications.
#'
#' @param product_data Data frame with market_ids and firm_ids
#' @param kappa "standard" (same firm = 1), "monopoly" (all 1), or "single" (identity)
#' @return Block-diagonal ownership matrix (N x N)
#' @export
#
# The ownership matrix O enters the Bertrand pricing FOC:
#   p - mc = -(O * dS/dp)^{-1} s
# where * is element-wise multiplication. O(j,k) = 1 means firms j and k
# internalize each other's profits when setting prices. Three scenarios:
#   "standard": O(j,k) = 1 iff j and k belong to same firm (multi-product).
#   "monopoly": O(j,k) = 1 for all j,k in the market (full collusion).
#   "single": O = I, each product priced independently (single-product firms).
# Comparing markups across these ownership structures quantifies the
# competitive effects of market structure changes.
build_custom_ownership <- function(product_data, kappa = "standard") {
  market_ids <- product_data$market_ids
  firm_ids <- product_data$firm_ids
  N <- length(market_ids)

  if (kappa == "standard") {
    return(build_ownership_matrix(firm_ids, market_ids))
  }

  O <- matrix(0, N, N)
  markets <- unique(market_ids)
  for (t in markets) {
    idx <- which(market_ids == t)
    J_t <- length(idx)
    if (kappa == "monopoly") {
      O[idx, idx] <- 1
    } else if (kappa == "single") {
      O[idx, idx] <- diag(J_t)
    }
  }
  O
}

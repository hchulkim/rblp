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
      # Own firm, other products
      own_idx <- which(firms_t == f)
      own_idx <- own_idx[own_idx != j]
      if (length(own_idx) > 0) {
        own_other[idx[j], ] <- colSums(X_t[own_idx, , drop = FALSE])
      }

      # Rival products
      rival_idx <- which(firms_t != f)
      if (length(rival_idx) > 0) {
        rival[idx[j], ] <- colSums(X_t[rival_idx, , drop = FALSE])
      }
    }
  }

  # Name columns
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
build_differentiation_instruments <- function(X, market_ids, firm_ids,
                                               method = "local",
                                               interact = NULL) {
  X <- as.matrix(X)
  N <- nrow(X)
  K <- ncol(X)
  markets <- unique(market_ids)

  # Compute global SDs for each characteristic (for local method)
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
          # Count nearby products
          if (length(own_idx) > 0) {
            dists <- abs(X_t[own_idx, k] - X_t[j, k])
            own_iv[idx[j], k] <- sum(dists < sds[k])
          }
          if (length(rival_idx) > 0) {
            dists <- abs(X_t[rival_idx, k] - X_t[j, k])
            rival_iv[idx[j], k] <- sum(dists < sds[k])
          }
        } else {
          # Quadratic: sum of squared distances
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

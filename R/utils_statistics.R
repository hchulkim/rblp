#' IV/GMM Estimation of Linear Parameters
#'
#' Concentrates out linear parameters via instrumental variables.
#'
#' @param X Regressor matrix (N x K)
#' @param Z Instrument matrix (N x M)
#' @param W Weighting matrix (M x M)
#' @param y Dependent variable (N x 1)
#' @param jacobian Optional Jacobian of y w.r.t. nonlinear parameters (N x P)
#' @return List with parameters, residuals, covariances, and residual_jacobian
#' @keywords internal
iv_estimate <- function(X, Z, W, y, jacobian = NULL) {
  XZ <- crossprod(X, Z)    # K x M
  ZX <- t(XZ)              # M x K
  ZY <- crossprod(Z, y)    # M x 1

  # (X'Z W Z'X)^{-1}
  bread <- XZ %*% W %*% ZX
  bread_inv <- approximately_invert(bread)$inverse

  # beta = (X'Z W Z'X)^{-1} X'Z W Z'y
  params <- as.numeric(bread_inv %*% XZ %*% W %*% ZY)
  residuals <- as.numeric(y - X %*% params)

  residual_jacobian <- NULL
  if (!is.null(jacobian)) {
    # d_residuals/d_theta = d_y/d_theta - X (X'Z W Z'X)^{-1} X'Z W Z' d_y/d_theta
    ZJ <- crossprod(Z, jacobian)  # M x P
    residual_jacobian <- jacobian - X %*% (bread_inv %*% XZ %*% W %*% ZJ)
  }

  list(
    parameters = params,
    residuals = residuals,
    covariances = bread_inv,
    residual_jacobian = residual_jacobian
  )
}

#' Compute GMM Weighting Matrix
#' @param S Moment covariance matrix
#' @return Weighting matrix W = inverse of S, symmetrized
#' @keywords internal
compute_gmm_weights <- function(S) {
  W <- approximately_invert(S)$inverse
  # Symmetrize
  (W + t(W)) / 2
}

#' Compute GMM Moment Covariances
#'
#' @param u_list List of residual vectors (one per equation)
#' @param Z_list List of instrument matrices (one per equation)
#' @param type "robust", "clustered", or "unadjusted"
#' @param clustering_ids Optional clustering identifiers
#' @return Moment covariance matrix S
#' @keywords internal
compute_gmm_moment_covariances <- function(u_list, Z_list, type = "robust",
                                            clustering_ids = NULL) {
  N <- length(u_list[[1]])
  n_eq <- length(u_list)

  if (type == "unadjusted") {
    # S = block matrix of cov(u_a, u_b) * Z_a' Z_b / N
    blocks <- list()
    offset_r <- 0
    for (a in seq_len(n_eq)) {
      offset_c <- 0
      M_a <- ncol(Z_list[[a]])
      for (b in seq_len(n_eq)) {
        M_b <- ncol(Z_list[[b]])
        cov_ab <- sum(u_list[[a]] * u_list[[b]]) / N
        block <- cov_ab * crossprod(Z_list[[a]], Z_list[[b]]) / N
        if (a == 1 && b == 1) {
          total_dim <- sum(sapply(Z_list, ncol))
          S <- matrix(0, total_dim, total_dim)
        }
        S[(offset_r + 1):(offset_r + M_a), (offset_c + 1):(offset_c + M_b)] <- block
        offset_c <- offset_c + M_b
      }
      offset_r <- offset_r + M_a
    }
    return((S + t(S)) / 2)
  }

  # Form stacked moment contributions: g_i = [u1_i * Z1_i, u2_i * Z2_i, ...]
  g_parts <- lapply(seq_len(n_eq), function(a) {
    sweep(Z_list[[a]], 1, u_list[[a]], "*")
  })
  g <- do.call(cbind, g_parts)  # N x total_moments

  if (type == "clustered" && !is.null(clustering_ids)) {
    # Sum within clusters
    uid <- unique(clustering_ids)
    g_clustered <- matrix(0, length(uid), ncol(g))
    for (i in seq_along(uid)) {
      idx <- which(clustering_ids == uid[i])
      g_clustered[i, ] <- colSums(g[idx, , drop = FALSE])
    }
    S <- crossprod(g_clustered) / N
  } else {
    S <- crossprod(g) / N
  }

  (S + t(S)) / 2
}

#' Compute GMM Parameter Covariances
#'
#' @param W Weighting matrix
#' @param S Moment covariance matrix
#' @param G Moment Jacobian (d_moments / d_theta)
#' @param se_type "robust", "clustered", or "unadjusted"
#' @return Parameter covariance matrix V
#' @keywords internal
compute_gmm_parameter_covariances <- function(W, S, G, se_type = "robust") {
  GWG <- crossprod(G, W %*% G)
  GWG_inv <- approximately_invert(GWG)$inverse

  if (se_type == "unadjusted") {
    return(GWG_inv)
  }

  # Sandwich: (G'WG)^{-1} G'W S W G (G'WG)^{-1}
  GW <- crossprod(G, W)
  V <- GWG_inv %*% GW %*% S %*% t(GW) %*% GWG_inv
  (V + t(V)) / 2
}

#' Compute Parameter Sensitivity (Andrews, Gentzkow, Shapiro 2017)
#' @param W Weighting matrix
#' @param G Moment Jacobian
#' @return Sensitivity matrix Lambda
#' @keywords internal
compute_gmm_sensitivity <- function(W, G) {
  GWG_inv <- approximately_invert(crossprod(G, W %*% G))$inverse
  -GWG_inv %*% crossprod(G, W)
}

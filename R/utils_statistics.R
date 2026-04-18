#' IV/GMM Estimation of Linear Parameters
#'
#' Concentrates out linear parameters via instrumental variables.
#'
#' @param X Regressor matrix (N x K)
#' @param Z Instrument matrix (N x M)
#' @param W Weighting matrix (M x M)
#' @param y Dependent variable (N x 1)
#' @param jacobian Optional Jacobian of y w.r.t. nonlinear parameters (N x P)
#' @param precomputed Optional list with reusable linear algebra pieces:
#'   `XZ`, `bread_inv`, and `XZW`
#' @return List with parameters, residuals, covariances, and residual_jacobian
#' @keywords internal
iv_estimate <- function(X, Z, W, y, jacobian = NULL, precomputed = NULL) {
  # Instrumental variables (2SLS/GMM) estimation of linear parameters beta,
  # concentrating them out of the GMM objective. In BLP, the "dependent
  # variable" y is the vector of mean utilities delta (which depends on the
  # nonlinear parameters theta), X contains the product characteristics,
  # and Z contains the excluded instruments plus exogenous regressors.
  # The GMM estimator for linear parameters given weighting matrix W is:
  #   beta = (X'Z W Z'X)^{-1} X'Z W Z'y
  # This is equivalent to 2SLS when W = (Z'Z)^{-1}.
  XZ <- precomputed$XZ %||% crossprod(X, Z)    # K x M
  ZX <- t(XZ)                                   # M x K
  XZW <- precomputed$XZW %||% (XZ %*% W)       # K x M
  ZY <- crossprod(Z, y)                         # M x 1

  # "Bread" matrix of the IV sandwich -- the precision of the IV estimator.
  # Inverting this gives the leading term of the parameter covariance.
  bread_inv <- precomputed$bread_inv
  if (is.null(bread_inv)) {
    bread <- XZW %*% ZX
    bread_inv <- approximately_invert(bread)$inverse
  }

  # Closed-form IV/GMM estimate of the linear parameters
  params <- as.numeric(bread_inv %*% XZW %*% ZY)
  residuals <- as.numeric(y - X %*% params)

  residual_jacobian <- NULL
  if (!is.null(jacobian)) {
    # The Jacobian of the structural residuals xi with respect to the
    # nonlinear parameters theta. By the implicit function theorem,
    # d(xi)/d(theta) = d(delta)/d(theta) - X * d(beta)/d(theta).
    # Since beta is a function of delta (which depends on theta), the chain
    # rule gives the expression below. This Jacobian is needed for the
    # analytic gradient of the GMM objective.
    ZJ <- crossprod(Z, jacobian)  # M x P
    residual_jacobian <- jacobian - X %*% (bread_inv %*% XZW %*% ZJ)
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
  # Efficient GMM uses W = S^{-1} as the weighting matrix, where S is the
  # variance of the sample moment conditions. This is the optimal weighting
  # matrix that minimizes the asymptotic variance of the GMM estimator
  # (Hansen 1982). In the first step of two-step GMM, S is typically
  # estimated under a homoskedasticity assumption (or W = (Z'Z)^{-1});
  # in the second step, S is re-estimated from first-step residuals.
  # Symmetrization corrects any numerical asymmetry from the inversion.
  W <- approximately_invert(S)$inverse
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

  # Unadjusted (homoskedastic) moment covariance: assumes E[u_i u_j | Z] is
  # constant across observations. The covariance of the moment conditions
  # simplifies to S_{ab} = sigma_{ab} * (Z_a' Z_b / N), where sigma_{ab} is
  # the cross-equation error covariance. This is the classical 2SLS assumption.
  if (type == "unadjusted") {
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

  # Robust (heteroskedasticity-consistent) moment covariance: form the
  # observation-level moment contributions g_i = [u1_i * Z1_i, u2_i * Z2_i, ...],
  # then estimate S = (1/N) sum_i g_i g_i'. This is the Eicker-Huber-White
  # "sandwich" estimator that is valid under arbitrary heteroskedasticity.
  g_parts <- lapply(seq_len(n_eq), function(a) {
    sweep(Z_list[[a]], 1, u_list[[a]], "*")
  })
  g <- do.call(cbind, g_parts)  # N x total_moments

  if (type == "clustered" && !is.null(clustering_ids)) {
    # Clustered moment covariance (Cameron, Gelbach, Miller 2011): sum the
    # moment contributions within each cluster before forming the outer product.
    # This allows for arbitrary within-cluster correlation of the errors,
    # which is important when markets or firms define natural clusters.
    uid <- unique(clustering_ids)
    g_clustered <- matrix(0, length(uid), ncol(g))
    for (i in seq_along(uid)) {
      idx <- which(clustering_ids == uid[i])
      g_clustered[i, ] <- colSums(g[idx, , drop = FALSE])
    }
    S <- crossprod(g_clustered) / N
  } else {
    # Robust (White) covariance: S = g'g / N
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
  # Parameter covariance via the GMM sandwich formula. G is the Jacobian of
  # the moment conditions with respect to all parameters (both linear beta
  # and nonlinear theta). The sandwich form:
  #   V = (G'WG)^{-1} G'W S WG (G'WG)^{-1}
  # is valid regardless of whether W equals the optimal weighting matrix S^{-1}.
  # When W = S^{-1} (efficient GMM), this simplifies to V = (G'WG)^{-1}.
  # The "bread" (G'WG)^{-1} captures the curvature of the GMM objective,
  # while the "meat" G'W S WG captures the sampling variability of the moments.
  GWG <- crossprod(G, W %*% G)
  GWG_inv <- approximately_invert(GWG)$inverse

  if (se_type == "unadjusted") {
    # Under correct specification and optimal weighting, the bread alone
    # gives the asymptotic variance (no need for the sandwich).
    return(GWG_inv)
  }

  # Full sandwich covariance: robust to misspecification of the weighting
  # matrix and to heteroskedasticity/clustering in the moment conditions.
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
  # Sensitivity of the GMM estimates to individual moment conditions
  # (Andrews, Gentzkow, Shapiro 2017). Lambda = -(G'WG)^{-1} G'W measures
  # how much each parameter would change if a particular moment condition
  # were perturbed. Large entries indicate that the estimate relies heavily
  # on specific moments, which is informative for assessing instrument quality
  # and model robustness.
  GWG_inv <- approximately_invert(crossprod(G, W %*% G))$inverse
  -GWG_inv %*% crossprod(G, W)
}

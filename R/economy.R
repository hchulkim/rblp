#' @title BLP Economy Base Class
#' @description Base class storing product/agent data and shared computations.
#' @keywords internal
BLPEconomy <- R6::R6Class("BLPEconomy",
  public = list(
    T = 0L, N = 0L, F = 0L, I = 0L,
    K1 = 0L, K2 = 0L, K3 = 0L, D = 0L,
    MD = 0L, MS = 0L, H = 0L,
    products = NULL,
    agents = NULL,
    product_formulations = NULL,
    agent_formulation = NULL,
    unique_market_ids = NULL,
    unique_firm_ids = NULL,
    unique_nesting_ids = NULL,
    rc_types = NULL,
    epsilon_scale = 1.0,
    costs_type = "linear",

    initialize = function(product_formulations, product_data,
                          agent_formulation = NULL, agent_data = NULL,
                          integration = NULL, rc_types = NULL,
                          epsilon_scale = 1.0, costs_type = "linear",
                          add_exogenous = TRUE) {
      # Normalize formulations to list
      if (inherits(product_formulations, "BLPFormulation")) {
        product_formulations <- list(product_formulations)
      }
      n_form <- length(product_formulations)
      self$product_formulations <- product_formulations
      self$agent_formulation <- agent_formulation

      # Validate product data
      stopifnot("market_ids" %in% names(product_data))
      stopifnot("shares" %in% names(product_data))
      stopifnot("prices" %in% names(product_data))

      product_data <- as.data.frame(product_data)
      self$N <- nrow(product_data)
      self$unique_market_ids <- unique(product_data$market_ids)
      self$T <- length(self$unique_market_ids)

      # Build market indices
      private$market_indices_ <- split(seq_len(self$N), product_data$market_ids)
      private$product_data_ <- product_data

      # Build the design matrices from BLPFormulation objects.
      # X1 (N x K1): linear demand characteristics that enter mean utility delta.
      #   These include the constant, product characteristics, and prices.
      #   The coefficients on X1 (beta) are concentrated out of the GMM objective.
      X1 <- product_formulations[[1]]$build_matrix(product_data)
      self$K1 <- ncol(X1)

      # X2 (N x K2): nonlinear demand characteristics whose coefficients vary
      #   across consumers via random coefficients. X2 often overlaps with X1
      #   (e.g., price appears in both). X2 enters mu_ij = X2_j' * (sigma*nu_i + pi*D_i).
      X2 <- NULL
      if (n_form >= 2) {
        X2 <- product_formulations[[2]]$build_matrix(product_data)
        self$K2 <- ncol(X2)
      }

      # X3 (N x K3): supply-side cost characteristics. If specified, the model
      #   jointly estimates demand and supply, where marginal cost is
      #   mc_j = X3_j' * gamma + omega_j. The supply-side moment condition
      #   is E[Z_s' omega] = 0.
      X3 <- NULL
      if (n_form >= 3) {
        X3 <- product_formulations[[3]]$build_matrix(product_data)
        self$K3 <- ncol(X3)
        stopifnot("firm_ids" %in% names(product_data))
      }

      # Firm IDs
      firm_ids <- product_data$firm_ids
      if (!is.null(firm_ids)) {
        self$unique_firm_ids <- unique(firm_ids)
        self$F <- length(self$unique_firm_ids)
      }

      # Nesting IDs
      nesting_ids <- product_data$nesting_ids
      if (!is.null(nesting_ids)) {
        self$unique_nesting_ids <- unique(nesting_ids)
        self$H <- length(self$unique_nesting_ids)
      }

      # Extract "excluded instruments" -- variables that appear in Z but not in X.
      # These are the instruments that provide identification for endogenous
      # regressors (e.g., BLP instruments, Hausman instruments, cost shifters).
      # Columns named demand_instruments0, demand_instruments1, ... are collected.
      ZD_excluded <- extract_columns(product_data, "demand_instruments")
      ZS_excluded <- extract_columns(product_data, "supply_instruments")

      # Assemble the full instrument matrix Z = [X_exog, Z_excluded].
      # The instrument set must include all exogenous regressors from X1
      # (the "included instruments") plus the excluded instruments. Prices and
      # shares are excluded from the exogenous columns because prices are
      # endogenous (correlated with xi) and shares are the dependent variable.
      # The intercept is included as an instrument since it is exogenous.
      # This is the standard IV/GMM convention: Z must span the exogenous
      # part of X so that the IV projection P_Z X has the right column space.
      if (add_exogenous && !is.null(X1)) {
        x1_names <- colnames(X1)
        exog_cols <- which(!grepl("prices|shares", x1_names, ignore.case = TRUE) &
                            x1_names != "(Intercept)")
        ic <- which(x1_names == "(Intercept)")
        exog_cols <- sort(c(ic, exog_cols))

        if (length(exog_cols) > 0) {
          X1_exog <- X1[, exog_cols, drop = FALSE]
          ZD <- if (!is.null(ZD_excluded)) cbind(X1_exog, ZD_excluded) else X1_exog
        } else {
          ZD <- ZD_excluded
        }
      } else {
        ZD <- ZD_excluded
      }

      # Supply-side instruments: exogenous cost shifters from X3 plus excluded
      # supply instruments. Only used when a supply formulation is present (K3 > 0).
      # Otherwise, supply_instruments* columns in the data are ignored.
      if (!is.null(X3) && self$K3 > 0) {
        if (add_exogenous) {
          x3_names <- colnames(X3)
          exog_cols3 <- which(!grepl("prices|shares", x3_names, ignore.case = TRUE))
          if (length(exog_cols3) > 0) {
            X3_exog <- X3[, exog_cols3, drop = FALSE]
            ZS <- if (!is.null(ZS_excluded)) cbind(X3_exog, ZS_excluded) else X3_exog
          } else {
            ZS <- ZS_excluded
          }
        } else {
          ZS <- ZS_excluded
        }
      } else {
        ZS <- NULL
      }

      # Frisch-Waugh-Lovell (FWL) demeaning for absorbed fixed effects.
      # When the formulation includes absorb = ~ product_ids (or market_ids),
      # we subtract group means from X1, ZD, delta, and the xi-Jacobian
      # rather than including a large dummy matrix. The FWL theorem guarantees
      # that OLS/IV on the demeaned data yields identical coefficient estimates
      # and residuals as the full dummy-variable regression, but with much
      # lower computational cost (O(N) vs O(N*G) where G is the number of groups).
      # X1 demeaning happens inside build_matrix via the formulation's absorb;
      # here we demean the instruments ZD so that the IV projection is consistent.
      absorb_form <- product_formulations[[1]]$get_absorb()
      if (!is.null(absorb_form)) {
        fe_vars <- all.vars(absorb_form)
        for (fv in fe_vars) {
          if (fv %in% names(product_data)) {
            grp <- product_data[[fv]]
            if (!is.null(ZD)) {
              for (j in seq_len(ncol(ZD))) {
                gm <- tapply(ZD[, j], grp, mean)
                ZD[, j] <- ZD[, j] - gm[match(grp, names(gm))]
              }
            }
            private$absorb_groups_ <- grp
            private$absorb_group_indices_ <- split(seq_len(self$N), grp)
          }
        }
      }

      self$MD <- if (!is.null(ZD)) ncol(ZD) else 0L
      self$MS <- if (!is.null(ZS)) ncol(ZS) else 0L

      # Build the ownership matrix O (N x N block-diagonal). Entry O_{jk} = 1
      # if products j and k are produced by the same firm in the same market,
      # 0 otherwise. This matrix encodes multi-product firm structure and is
      # used in the Bertrand-Nash pricing FOCs: the firm internalizes cross-
      # product demand effects among its own products when setting prices.
      # For merger simulation, one simply changes firm_ids to reflect the
      # post-merger ownership and re-solves for equilibrium prices.
      ownership <- NULL
      if (!is.null(firm_ids)) {
        ownership <- build_ownership_matrix(firm_ids, product_data$market_ids)
      }

      # Detect price column indices in X1 and X2
      price_col_x1 <- which(colnames(X1) == "prices")
      if (length(price_col_x1) == 0) {
        price_col_x1 <- grep("prices", colnames(X1))[1]
      }
      price_col_x2 <- NULL
      if (!is.null(X2)) {
        price_col_x2 <- which(colnames(X2) == "prices")
        if (length(price_col_x2) == 0) {
          price_col_x2 <- grep("prices", colnames(X2))[1]
        }
        if (length(price_col_x2) == 0) price_col_x2 <- NULL
      }

      # Store products
      self$products <- list(
        X1 = X1, X2 = X2, X3 = X3,
        ZD = ZD, ZS = ZS,
        shares = product_data$shares,
        prices = product_data$prices,
        firm_ids = firm_ids,
        nesting_ids = nesting_ids,
        clustering_ids = product_data$clustering_ids,
        ownership = ownership,
        price_col_x1 = if (length(price_col_x1) > 0) price_col_x1[1] else NULL,
        price_col_x2 = if (!is.null(price_col_x2)) price_col_x2[1] else NULL,
        market_ids = product_data$market_ids,
        original_data = product_data
      )

      # Build agent data
      self$rc_types <- rc_types %||% rep("linear", self$K2)
      self$epsilon_scale <- epsilon_scale
      self$costs_type <- costs_type

      if (!is.null(integration) && self$K2 > 0) {
        # Build quadrature nodes and weights to approximate the integral over
        # the mixing distribution of random coefficients. For "product" (Gauss-
        # Hermite) integration, each node represents a draw from the standard
        # normal distribution with an associated probability weight. The share
        # integral s_j = integral P_ij f(nu) dnu is approximated as
        # s_j ~ sum_i w_i * P_ij(nu_i). More nodes give higher accuracy but
        # increase computation; the number of nodes grows as size^K2 (tensor
        # product), so K2 > 3-4 typically requires Monte Carlo or sparse grids.
        int_result <- integration$build(self$K2)
        n_agents <- nrow(int_result$nodes)

        # Replicate the same integration nodes for every market. This assumes
        # the distribution of unobserved tastes is identical across markets
        # (the standard BLP assumption). Each market gets its own copy of
        # nodes/weights so that market-level operations can be parallelized.
        total_I <- n_agents * self$T
        all_nodes <- matrix(0, total_I, self$K2)
        all_weights <- numeric(total_I)
        all_market_ids <- character(total_I)

        for (i in seq_along(self$unique_market_ids)) {
          idx <- ((i - 1) * n_agents + 1):(i * n_agents)
          all_nodes[idx, ] <- int_result$nodes
          all_weights[idx] <- int_result$weights
          all_market_ids[idx] <- rep(self$unique_market_ids[i], n_agents)
        }

        demos <- NULL
        if (!is.null(agent_formulation) && !is.null(agent_data)) {
          demos <- agent_formulation$build_matrix(agent_data)
        }

        self$agents <- list(
          nodes = all_nodes,
          weights = all_weights,
          demographics = demos,
          market_ids = all_market_ids,
          original_data = agent_data
        )
        self$I <- total_I
        self$D <- if (!is.null(demos)) ncol(demos) else 0L

      } else if (!is.null(agent_data)) {
        agent_data <- as.data.frame(agent_data)
        nodes <- extract_columns(agent_data, "nodes")
        weights <- agent_data$weights
        demos <- NULL
        if (!is.null(agent_formulation)) {
          demos <- agent_formulation$build_matrix(agent_data)
        }

        self$agents <- list(
          nodes = nodes,
          weights = weights,
          demographics = demos,
          market_ids = agent_data$market_ids,
          original_data = agent_data
        )
        self$I <- nrow(agent_data)
        self$D <- if (!is.null(demos)) ncol(demos) else 0L
      } else {
        self$agents <- list(nodes = NULL, weights = NULL, demographics = NULL,
                            market_ids = NULL, original_data = NULL)
        self$I <- 0L
      }

      private$agent_data_ <- agent_data
      private$integration_ <- integration
    },

    get_market_data = function(market_id) {
      idx <- private$market_indices_[[as.character(market_id)]]
      if (is.null(idx)) stop(sprintf("Market '%s' not found", market_id))

      prods <- list(
        X1 = self$products$X1[idx, , drop = FALSE],
        X2 = if (!is.null(self$products$X2)) self$products$X2[idx, , drop = FALSE] else NULL,
        X3 = if (!is.null(self$products$X3)) self$products$X3[idx, , drop = FALSE] else NULL,
        ZD = if (!is.null(self$products$ZD)) self$products$ZD[idx, , drop = FALSE] else NULL,
        ZS = if (!is.null(self$products$ZS)) self$products$ZS[idx, , drop = FALSE] else NULL,
        shares = self$products$shares[idx],
        prices = self$products$prices[idx],
        firm_ids = if (!is.null(self$products$firm_ids)) self$products$firm_ids[idx] else NULL,
        nesting_ids = if (!is.null(self$products$nesting_ids)) self$products$nesting_ids[idx] else NULL,
        ownership = if (!is.null(self$products$ownership)) self$products$ownership[idx, idx, drop = FALSE] else NULL,
        price_col_x1 = self$products$price_col_x1,
        price_col_x2 = self$products$price_col_x2,
        beta_price = NULL  # set during estimation
      )

      # Agent data for this market
      if (!is.null(self$agents$market_ids)) {
        a_idx <- which(self$agents$market_ids == market_id)
        agts <- list(
          nodes = if (!is.null(self$agents$nodes)) self$agents$nodes[a_idx, , drop = FALSE] else NULL,
          weights = self$agents$weights[a_idx],
          demographics = if (!is.null(self$agents$demographics)) self$agents$demographics[a_idx, , drop = FALSE] else NULL
        )
      } else {
        agts <- list(nodes = NULL, weights = 1, demographics = NULL)
      }

      list(products = prods, agents = agts, indices = idx)
    },

    compute_logit_delta = function(rho = NULL) {
      # Compute the closed-form logit delta as a starting value for the BLP
      # contraction mapping. For the plain logit (no random coefficients),
      # Berry (1994) showed that the share equation s_j = exp(delta_j)/(1+sum_k exp(delta_k))
      # can be analytically inverted to:
      #   delta_j = log(s_j) - log(s_0)
      # where s_0 = 1 - sum_j s_j is the outside good share. This is exact
      # for logit and serves as a warm start for the random coefficients model.
      shares <- self$products$shares
      market_ids <- self$products$market_ids
      delta <- numeric(self$N)

      for (t in self$unique_market_ids) {
        idx <- private$market_indices_[[as.character(t)]]
        s <- shares[idx]
        s0 <- 1 - sum(s)  # outside good share

        delta[idx] <- log(s) - log(s0)

        if (!is.null(rho) && any(rho != 0)) {
          # Nested logit inversion (Berry 1994, eq. 3):
          #   delta_j = log(s_j) - log(s_0) - rho * log(s_j/s_g)
          # where s_g = sum_{k in g} s_k is the nest-level share. The extra
          # term -rho*log(s_j/s_g) accounts for the within-nest correlation:
          # higher rho means more within-nest substitution is captured by
          # the nesting structure rather than by delta.
          nids <- self$products$nesting_ids[idx]
          if (!is.null(nids)) {
            rho_val <- rho[1]
            s_g <- tapply(s, nids, sum)
            s_g_expanded <- s_g[match(nids, names(s_g))]
            delta[idx] <- delta[idx] - rho_val * (log(s) - log(s_g_expanded))
          }
        }
      }

      # If fixed effects are absorbed, demean delta for FWL consistency.
      # This ensures the initial delta lives in the same demeaned space as
      # the X1 and ZD matrices used in the IV regression.
      if (!is.null(private$absorb_groups_)) {
        grp <- private$absorb_groups_
        gm <- tapply(delta, grp, mean)
        delta <- as.numeric(delta - gm[match(grp, names(gm))])
      }

      delta
    },

    print = function(...) {
      cat(sprintf("BLPEconomy: %d markets, %d products\n", self$T, self$N))
      invisible(self)
    }
  ),
  private = list(
    market_indices_ = NULL,
    product_data_ = NULL,
    agent_data_ = NULL,
    integration_ = NULL,
    absorb_groups_ = NULL,
    absorb_group_indices_ = NULL
  )
)

#' Linked Multiple Factor Analysis (Linked MFA)
#'
#' @description
#' Linked MFA generalizes Multiple Factor Analysis (MFA) to the case where a single
#' reference block `Y` (with `N` rows) is linked to multiple blocks `X_k` that may
#' have different numbers of rows. Each row of `X_k` is mapped to a row of `Y` via
#' an index vector `row_index[[k]]`. The model estimates a shared score matrix
#' `S` for the rows of `Y` and block-specific loading matrices for `Y` and each `X_k`.
#'
#' Optionally, a feature prior can be provided to encourage corresponding or
#' similar features across different `X_k` blocks to have similar loading vectors.
#'
#' @details
#' ## Model
#' The fitted model has the form:
#' \deqn{Y \approx S B^\top}
#' \deqn{X_k \approx S[\mathrm{idx}_k,] V_k^\top}
#' where `S` is `N × ncomp`, `B` is `q × ncomp`, and each `V_k` is `p_k × ncomp`.
#'
#' ## Feature similarity prior (v1)
#' When `feature_lambda > 0` and `feature_groups` is supplied, Linked MFA applies
#' a group-shrinkage penalty that pulls the loading vectors of features in the same
#' group toward a shared group center.
#'
#' @param Y Numeric matrix/data.frame (`N × q`) serving as the reference block.
#' @param X A list of numeric matrices/data.frames. Each element `X[[k]]` is `n_k × p_k`.
#' @param row_index A list of integer vectors. `row_index[[k]]` has length `n_k`
#'   and maps rows of `X[[k]]` to rows of `Y` (values in `1..N`).
#' @param preproc A `multivarious` preprocessing pipeline (a `pre_processor`/`prepper`)
#'   or a list of them. If a list, it must have length `1 + length(X)` and will be
#'   applied to `c(list(Y), X)` in that order.
#' @param ncomp Integer number of components to extract.
#' @param normalization Block weighting scheme. `"MFA"` uses inverse squared first
#'   singular value per block; `"None"` uses uniform weights; `"custom"` uses `alpha`.
#' @param alpha Optional numeric vector of per-block weights (length `1 + length(X)`),
#'   used when `normalization = "custom"`. The first weight corresponds to `Y`.
#' @param feature_groups Feature prior specification. One of:
#'   * `NULL` (no feature prior),
#'   * `"colnames"` to group X-features with identical column names across blocks,
#'   * a `data.frame` with columns `block`, `feature`, `group` and optional `weight`.
#'     `block` refers to a name or index in `X` (not including `Y`), and `feature` is
#'     a column name or index within that block.
#' @param feature_lambda Non-negative scalar controlling strength of the feature prior.
#' @param max_iter Maximum number of alternating least-squares iterations.
#' @param tol Relative tolerance on the objective for convergence.
#' @param ridge Non-negative ridge stabilization added to normal equations.
#' @param verbose Logical; if `TRUE`, prints iteration progress.
#' @param ... Unused (reserved for future extensions).
#'
#' @return An object inheriting from `multivarious::multiblock_biprojector` with
#'   additional class `"linked_mfa"`. The object contains global scores in `s`,
#'   concatenated loadings in `v`, and block mappings in `block_indices`. Additional
#'   fields include `V_list`, `B`, `row_index`, `alpha_blocks`, and `objective_trace`.
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' N <- 30
#' Y <- matrix(rnorm(N * 5), N, 5)
#' X1 <- matrix(rnorm(20 * 10), 20, 10)
#' X2 <- matrix(rnorm(15 * 8), 15, 8)
#' idx1 <- sample.int(N, nrow(X1), replace = FALSE)
#' idx2 <- sample.int(N, nrow(X2), replace = FALSE)
#'
#' fit <- linked_mfa(Y, list(X1 = X1, X2 = X2), list(X1 = idx1, X2 = idx2), ncomp = 2)
#' stopifnot(nrow(multivarious::scores(fit)) == N)
#' }
#' @export
linked_mfa <- function(Y,
                       X,
                       row_index,
                       preproc = multivarious::center(),
                       ncomp = 2,
                       normalization = c("MFA", "None", "custom"),
                       alpha = NULL,
                       feature_groups = NULL,
                       feature_lambda = 0,
                       max_iter = 50,
                       tol = 1e-6,
                       ridge = 1e-8,
                       verbose = FALSE,
                       ...) {
  normalization <- match.arg(normalization)

  # ---------------------------------------------------------------------------
  # 0. Basic checks / coercions
  # ---------------------------------------------------------------------------
  chk::chk_true(is.matrix(Y) || is.data.frame(Y))
  Y <- as.matrix(Y)
  chk::chk_list(X)
  chk::chk_true(length(X) >= 1)
  chk::chk_list(row_index)
  chk::chk_equal(length(X), length(row_index))
  ncomp <- as.integer(ncomp)
  chk::chk_integer(ncomp)
  chk::chk_gte(ncomp, 1)
  chk::chk_numeric(feature_lambda)
  chk::chk_gte(feature_lambda, 0)
  max_iter <- as.integer(max_iter)
  chk::chk_integer(max_iter)
  chk::chk_gte(max_iter, 1)
  chk::chk_numeric(tol)
  chk::chk_gte(tol, 0)
  chk::chk_numeric(ridge)
  chk::chk_gte(ridge, 0)
  chk::chk_flag(verbose)

  X <- lapply(X, function(x) {
    chk::chk_true(is.matrix(x) || is.data.frame(x))
    as.matrix(x)
  })

  N <- nrow(Y)
  for (k in seq_along(X)) {
    idx <- row_index[[k]]
    chk::chk_true(is.numeric(idx) || is.integer(idx))
    idx <- as.integer(idx)
    chk::chk_equal(length(idx), nrow(X[[k]]))
    if (anyNA(idx)) stop("row_index cannot contain NA.", call. = FALSE)
    if (any(idx < 1L) || any(idx > N)) {
      stop("row_index values must be in 1..nrow(Y).", call. = FALSE)
    }
    row_index[[k]] <- idx
  }

  if (is.null(names(X))) {
    names(X) <- paste0("X", seq_along(X))
  }
  if (is.null(names(row_index))) {
    names(row_index) <- names(X)
  }

  # ---------------------------------------------------------------------------
  # 1. Preprocess blocks
  # ---------------------------------------------------------------------------
  blocks <- c(list(Y = Y), X)
  if (is.list(preproc) && !inherits(preproc, "pre_processor") && !inherits(preproc, "prepper")) {
    chk::chk_equal(length(preproc), length(blocks))
  }
  prep_res <- prepare_block_preprocessors(blocks, preproc, check_consistent_ncol = FALSE)
  Yp <- prep_res$Xp[[1]]
  Xp <- prep_res$Xp[-1]
  proclist <- prep_res$proclist

  # Ensure N refers to preprocessed Y rows (should match input)
  N <- nrow(Yp)

  # ---------------------------------------------------------------------------
  # 2. Block weights (alpha)
  # ---------------------------------------------------------------------------
  alpha_blocks <- if (normalization == "custom") {
    if (is.null(alpha)) {
      stop("When normalization = 'custom', 'alpha' must be provided.", call. = FALSE)
    }
    alpha <- as.numeric(alpha)
    chk::chk_equal(length(alpha), length(blocks))
    if (any(!is.finite(alpha)) || any(alpha < 0)) {
      stop("'alpha' must be finite and non-negative.", call. = FALSE)
    }
    alpha
  } else if (normalization == "None") {
    rep(1, length(blocks))
  } else { # MFA-style
    first_sv <- function(mat) {
      mdim <- min(dim(mat))
      method <- if (mdim < 3) "base" else "svds"
      tryCatch(
        multivarious::svd_wrapper(mat, ncomp = 1, method = method)$sdev[1],
        error = function(e) multivarious::svd_wrapper(mat, ncomp = 1, method = "base")$sdev[1]
      )
    }
    sapply(prep_res$Xp, function(M) 1 / (first_sv(M)^2 + ridge))
  }
  names(alpha_blocks) <- c("Y", names(X))

  # ---------------------------------------------------------------------------
  # 3. Feature groups (X-only)
  # ---------------------------------------------------------------------------
  fg <- .lmfa_parse_feature_groups(feature_groups, Xp, feature_lambda)
  feature_lambda_eff <- fg$feature_lambda

  # ---------------------------------------------------------------------------
  # 4. Initialization from Y
  # ---------------------------------------------------------------------------
  min_p <- min(c(ncol(Yp), vapply(Xp, ncol, integer(1))))
  K <- min(as.integer(ncomp), nrow(Yp), min_p)
  if (K < 1L) stop("ncomp must be at least 1 and <= min(nrow(Y), ncol(Y)).", call. = FALSE)

  init <- .lmfa_init_from_Y(Yp, K, ridge = ridge)
  S <- init$S # N x K
  B <- init$B # q x K

  # Initialize V_k by ridge regression of X_k on S[idx_k,]
  V_list <- lapply(seq_along(Xp), function(k) {
    Sk <- S[row_index[[k]], , drop = FALSE]
    .lmfa_update_loadings(Sk, Xp[[k]], ridge = ridge)$V
  })
  names(V_list) <- names(Xp)

  objective_trace <- numeric(0)
  prev_obj <- Inf

  # ---------------------------------------------------------------------------
  # 5. ALS loop
  # ---------------------------------------------------------------------------
  for (iter in seq_len(max_iter)) {
    # (a) Update B (Y loadings)
    B <- .lmfa_update_loadings(S, Yp, ridge = ridge)$V

    # (b) Update feature group centers from current V (if enabled)
    centers <- .lmfa_group_centers(V_list, fg)

    # (c) Update each V_k (X loadings), optionally with group-shrinkage prior
    for (k in seq_along(Xp)) {
      Sk <- S[row_index[[k]], , drop = FALSE]
      V_list[[k]] <- .lmfa_update_V_block(
        Sk = Sk,
        Xk = Xp[[k]],
        fg_block = fg$by_block[[k]],
        weights_block = fg$weights_by_block[[k]],
        centers = centers,
        feature_lambda = feature_lambda_eff,
        ridge = ridge
      )
    }

    # (d) Recompute centers after V update for objective + next iteration
    centers <- .lmfa_group_centers(V_list, fg)

    # (e) Update S row-by-row (shared scores for Y rows)
    S <- .lmfa_update_scores(
      Y = Yp,
      B = B,
      X_list = Xp,
      V_list = V_list,
      row_index = row_index,
      alpha_y = alpha_blocks[1],
      alpha_blocks = alpha_blocks[-1],
      ridge = ridge
    )

    # (f) Orthonormalize S and rotate all loadings to preserve fitted values
    ortho <- .lmfa_orthonormalize_scores(S)
    S <- ortho$S
    rot <- ortho$R
    B <- B %*% t(rot)
    V_list <- lapply(V_list, function(V) V %*% t(rot))

    # (g) Objective / convergence
    obj <- .lmfa_objective(
      Y = Yp,
      S = S,
      B = B,
      X_list = Xp,
      V_list = V_list,
      row_index = row_index,
      alpha_y = alpha_blocks[1],
      alpha_blocks = alpha_blocks[-1],
      fg = fg,
      centers = centers,
      feature_lambda = feature_lambda_eff
    )
    objective_trace <- c(objective_trace, obj)

    rel_change <- abs(obj - prev_obj) / (abs(prev_obj) + ridge)
    if (isTRUE(verbose)) {
      message(
        sprintf(
          "linked_mfa iter %d: obj=%.6g, rel_change=%.3g",
          iter, obj, rel_change
        )
      )
    }
    if (is.finite(prev_obj) && rel_change < tol) break
    prev_obj <- obj
  }

  # ---------------------------------------------------------------------------
  # 6. Build multivarious-compatible return object
  # ---------------------------------------------------------------------------
  block_indices <- list()
  current <- 1L
  block_indices[[1]] <- current:(current + ncol(Yp) - 1L)
  current <- current + ncol(Yp)
  for (k in seq_along(Xp)) {
    block_indices[[k + 1L]] <- current:(current + ncol(Xp[[k]]) - 1L)
    current <- current + ncol(Xp[[k]])
  }
  names(block_indices) <- c("Y", names(Xp))

  proc <- multivarious::concat_pre_processors(proclist, block_indices)

  v_concat <- do.call(rbind, c(list(B), V_list))
  sdev <- apply(S, 2, stats::sd)

  multivarious::multiblock_biprojector(
    v = v_concat,
    s = S,
    sdev = sdev,
    preproc = proc,
    block_indices = block_indices,
    # linked_mfa metadata
    B = B,
    V_list = V_list,
    row_index = row_index,
    alpha_blocks = alpha_blocks,
    normalization = normalization,
    feature_groups = feature_groups,
    feature_lambda = feature_lambda_eff,
    objective_trace = objective_trace,
    classes = "linked_mfa"
  )
}

# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

.lmfa_init_from_Y <- function(Y, ncomp, ridge = 1e-8) {
  # Use a deterministic initialization: SVD on Y
  sv <- multivarious::svd_wrapper(Y, ncomp = ncomp, method = if (min(dim(Y)) < 3) "base" else "svds")
  # sv$u is N x ncomp, sv$v is q x ncomp, sv$sdev is length ncomp
  # Use scores as u * sdev (PCA-style), loadings as v
  S <- sv$u[, seq_len(ncomp), drop = FALSE] %*% diag(sv$sdev[seq_len(ncomp)], ncomp, ncomp)
  B <- sv$v[, seq_len(ncomp), drop = FALSE]

  # Stabilize / standardize to orthonormal score columns for downstream updates
  ortho <- .lmfa_orthonormalize_scores(S)
  S <- ortho$S
  B <- B %*% t(ortho$R)

  list(S = S, B = B)
}

.lmfa_update_loadings <- function(S, X, ridge = 1e-8) {
  # Solve min ||X - S V^T||_F^2 + ridge ||V||_F^2 for V (p x k)
  k <- ncol(S)
  StS <- crossprod(S) + ridge * diag(k)
  rhs <- crossprod(S, X) # k x p
  Vt <- solve(StS, rhs)  # k x p
  list(V = t(Vt))
}

.lmfa_orthonormalize_scores <- function(S) {
  qrS <- qr(S)
  Q <- qr.Q(qrS)
  R <- qr.R(qrS)
  list(S = Q, R = R)
}

.lmfa_parse_feature_groups <- function(feature_groups, X_list, feature_lambda) {
  # X_list is preprocessed X blocks (list of matrices), not including Y
  if (is.null(feature_groups) || isTRUE(feature_lambda == 0)) {
    return(list(
      enabled = FALSE,
      feature_lambda = 0,
      by_block = lapply(X_list, function(X) rep.int(NA_character_, ncol(X))),
      weights_by_block = lapply(X_list, function(X) rep.int(1, ncol(X))),
      members = list()
    ))
  }

  if (identical(feature_groups, "colnames")) {
    # Group features by identical column names across X blocks.
    cn_all <- lapply(X_list, colnames)
    if (any(vapply(cn_all, is.null, logical(1)))) {
      stop("feature_groups='colnames' requires non-NULL colnames for every X block.", call. = FALSE)
    }
    name_to_members <- new.env(parent = emptyenv())
    for (k in seq_along(X_list)) {
      for (j in seq_len(ncol(X_list[[k]]))) {
        nm <- cn_all[[k]][j]
        cur <- if (exists(nm, envir = name_to_members, inherits = FALSE)) {
          get(nm, envir = name_to_members, inherits = FALSE)
        } else {
          list()
        }
        cur[[length(cur) + 1L]] <- list(block = k, feature = j, weight = 1)
        assign(nm, cur, envir = name_to_members)
      }
    }
    groups <- ls(envir = name_to_members)
    members <- lapply(groups, function(g) get(g, envir = name_to_members, inherits = FALSE))
    names(members) <- groups
    # Keep only groups that link at least two features.
    keep <- vapply(members, length, integer(1)) >= 2L
    members <- members[keep]
    groups <- names(members)

    by_block <- lapply(X_list, function(X) rep.int(NA_character_, ncol(X)))
    for (g in groups) {
      for (m in members[[g]]) {
        by_block[[m$block]][m$feature] <- g
      }
    }
    return(list(
      enabled = TRUE,
      feature_lambda = feature_lambda,
      by_block = by_block,
      weights_by_block = lapply(X_list, function(X) rep.int(1, ncol(X))),
      members = members
    ))
  }

  if (!is.data.frame(feature_groups)) {
    stop("feature_groups must be NULL, 'colnames', or a data.frame.", call. = FALSE)
  }

  req <- c("block", "feature", "group")
  if (!all(req %in% names(feature_groups))) {
    stop("feature_groups data.frame must contain columns: block, feature, group.", call. = FALSE)
  }
  if (!("weight" %in% names(feature_groups))) {
    feature_groups$weight <- 1
  }

  by_block <- lapply(X_list, function(X) rep.int(NA_character_, ncol(X)))
  weights_by_block <- lapply(X_list, function(X) rep.int(1, ncol(X)))

  # Helpers to resolve block + feature indices
  x_names <- names(X_list)
  resolve_block <- function(b) {
    if (is.numeric(b) || is.integer(b)) {
      bi <- as.integer(b)
      if (bi < 1L || bi > length(X_list)) stop("feature_groups$block out of range.", call. = FALSE)
      return(bi)
    }
    b <- as.character(b)
    bi <- match(b, x_names)
    if (is.na(bi)) stop("feature_groups$block not found in names(X).", call. = FALSE)
    bi
  }
  resolve_feature <- function(k, f) {
    if (is.numeric(f) || is.integer(f)) {
      fi <- as.integer(f)
      if (fi < 1L || fi > ncol(X_list[[k]])) stop("feature_groups$feature out of range.", call. = FALSE)
      return(fi)
    }
    f <- as.character(f)
    cn <- colnames(X_list[[k]])
    if (is.null(cn)) stop("feature_groups uses feature names but X block has NULL colnames.", call. = FALSE)
    fi <- match(f, cn)
    if (is.na(fi)) stop("feature_groups$feature name not found in X block colnames.", call. = FALSE)
    fi
  }

  seen <- new.env(parent = emptyenv())
  members <- list()
  for (i in seq_len(nrow(feature_groups))) {
    k <- resolve_block(feature_groups$block[[i]])
    j <- resolve_feature(k, feature_groups$feature[[i]])
    g <- as.character(feature_groups$group[[i]])
    w <- as.numeric(feature_groups$weight[[i]])
    if (!is.finite(w) || w < 0) stop("feature_groups$weight must be finite and non-negative.", call. = FALSE)

    key <- paste(k, j, sep = ":")
    if (exists(key, envir = seen, inherits = FALSE)) {
      stop("Duplicate (block, feature) entries in feature_groups.", call. = FALSE)
    }
    assign(key, TRUE, envir = seen)

    by_block[[k]][j] <- g
    weights_by_block[[k]][j] <- w
    if (is.null(members[[g]])) members[[g]] <- list()
    members[[g]][[length(members[[g]]) + 1L]] <- list(block = k, feature = j, weight = w)
  }

  # Drop singleton groups (they don't link across features/blocks).
  keep <- vapply(members, length, integer(1)) >= 2L
  members <- members[keep]
  if (length(members) == 0) {
    return(list(
      enabled = FALSE,
      feature_lambda = 0,
      by_block = lapply(X_list, function(X) rep.int(NA_character_, ncol(X))),
      weights_by_block = lapply(X_list, function(X) rep.int(1, ncol(X))),
      members = list()
    ))
  }

  # Clear group labels for any features that ended up in dropped groups.
  valid_groups <- names(members)
  for (k in seq_along(by_block)) {
    by_block[[k]][!(by_block[[k]] %in% valid_groups)] <- NA_character_
  }

  list(
    enabled = TRUE,
    feature_lambda = feature_lambda,
    by_block = by_block,
    weights_by_block = weights_by_block,
    members = members
  )
}

.lmfa_group_centers <- function(V_list, fg) {
  if (!isTRUE(fg$enabled) || length(fg$members) == 0) return(NULL)
  groups <- names(fg$members)
  if (is.null(groups)) groups <- seq_along(fg$members)
  kdim <- ncol(V_list[[1]])
  centers <- matrix(0, nrow = length(fg$members), ncol = kdim)
  denom <- numeric(length(fg$members))
  rownames(centers) <- groups
  for (gi in seq_along(fg$members)) {
    mem <- fg$members[[gi]]
    for (m in mem) {
      v <- V_list[[m$block]][m$feature, , drop = TRUE]
      w <- m$weight
      centers[gi, ] <- centers[gi, ] + w * v
      denom[gi] <- denom[gi] + w
    }
    if (denom[gi] > 0) centers[gi, ] <- centers[gi, ] / denom[gi]
  }
  list(centers = centers, denom = denom)
}

.lmfa_update_V_block <- function(Sk,
                                Xk,
                                fg_block,
                                weights_block,
                                centers,
                                feature_lambda,
                                ridge = 1e-8) {
  # Update V_k (p x K). For ungrouped features, use ridge LS.
  # For grouped features, use ridge LS toward group center.
  p <- ncol(Xk)
  K <- ncol(Sk)

  G0 <- crossprod(Sk) + ridge * diag(K)
  rhs_all <- crossprod(Sk, Xk) # K x p

  Vt0 <- solve(G0, rhs_all) # K x p, unpenalized
  V <- t(Vt0)               # p x K

  if (isTRUE(feature_lambda > 0) && !is.null(fg_block) && !all(is.na(fg_block)) && !is.null(centers)) {
    # Cache cholesky factors per lambda value to avoid repeated decompositions.
    chol_cache <- new.env(parent = emptyenv())
    get_chol <- function(lam) {
      key <- format(lam, scientific = TRUE, digits = 10)
      if (exists(key, envir = chol_cache, inherits = FALSE)) {
        return(get(key, envir = chol_cache, inherits = FALSE))
      }
      R <- chol(G0 + lam * diag(K))
      assign(key, R, envir = chol_cache)
      R
    }
    solve_spd <- function(R, b) {
      backsolve(R, forwardsolve(t(R), b))
    }

    for (j in seq_len(p)) {
      g <- fg_block[j]
      if (is.na(g)) next
      gi <- match(g, rownames(centers$centers))
      if (is.na(gi)) next
      lam <- feature_lambda * weights_block[j]
      if (!is.finite(lam) || lam <= 0) next
      center <- centers$centers[gi, , drop = TRUE]
      rhs <- rhs_all[, j, drop = FALSE] + lam * matrix(center, nrow = K)
      R <- get_chol(lam)
      vj <- solve_spd(R, rhs)
      V[j, ] <- as.numeric(vj)
    }
  }

  V
}

.lmfa_update_scores <- function(Y,
                               B,
                               X_list,
                               V_list,
                               row_index,
                               alpha_y,
                               alpha_blocks,
                               ridge = 1e-8) {
  N <- nrow(Y)
  K <- ncol(B)

  BtB <- crossprod(B)  # K x K
  YB <- Y %*% B        # N x K

  # Precompute per-block pieces for S update.
  VtV_list <- lapply(V_list, crossprod)           # each K x K
  XV_list <- lapply(seq_along(X_list), function(k) X_list[[k]] %*% V_list[[k]]) # n_k x K
  sums_by_i <- lapply(seq_along(X_list), function(k) {
    idx <- row_index[[k]]
    # rowsum returns groups present; build full N x K with zeros.
    rs <- rowsum(XV_list[[k]], idx, reorder = FALSE)
    out <- matrix(0, nrow = N, ncol = K)
    out[as.integer(rownames(rs)), ] <- rs
    out
  })
  counts_by_i <- lapply(seq_along(X_list), function(k) {
    tabulate(row_index[[k]], nbins = N)
  })

  S_new <- matrix(0, nrow = N, ncol = K)
  I_K <- diag(K)
  for (i in seq_len(N)) {
    A <- alpha_y * BtB
    b <- alpha_y * YB[i, ]
    for (k in seq_along(X_list)) {
      cnt <- counts_by_i[[k]][i]
      if (cnt == 0L) next
      a_k <- alpha_blocks[k]
      A <- A + a_k * cnt * VtV_list[[k]]
      b <- b + a_k * sums_by_i[[k]][i, ]
    }
    A <- A + ridge * I_K
    S_new[i, ] <- solve(A, b)
  }
  S_new
}

.lmfa_objective <- function(Y,
                           S,
                           B,
                           X_list,
                           V_list,
                           row_index,
                           alpha_y,
                           alpha_blocks,
                           fg,
                           centers,
                           feature_lambda) {
  # Reconstruction errors
  err_y <- alpha_y * sum((Y - S %*% t(B))^2)
  err_x <- 0
  for (k in seq_along(X_list)) {
    Sk <- S[row_index[[k]], , drop = FALSE]
    resid <- X_list[[k]] - Sk %*% t(V_list[[k]])
    err_x <- err_x + alpha_blocks[k] * sum(resid^2)
  }

  pen <- 0
  if (isTRUE(feature_lambda > 0) && isTRUE(fg$enabled) && !is.null(centers) && length(fg$members) > 0) {
    # group penalty: Σ w ||v - center||^2
    groups <- names(fg$members)
    if (is.null(groups)) groups <- seq_along(fg$members)
    for (gi in seq_along(fg$members)) {
      gname <- groups[gi]
      cvec <- centers$centers[match(gname, rownames(centers$centers)), , drop = TRUE]
      for (m in fg$members[[gi]]) {
        v <- V_list[[m$block]][m$feature, , drop = TRUE]
        w <- m$weight
        pen <- pen + w * sum((v - cvec)^2)
      }
    }
    pen <- feature_lambda * pen
  }

  err_y + err_x + pen
}

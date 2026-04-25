#' @md
#' @rdname mcca
#' @details
#' The `mcca.list` method applies MCCA to a list of data matrices or data frames.
#' This method first converts the list to a `multiblock` object and then calls
#' `mcca.multiblock()`.
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' X <- replicate(3, matrix(rnorm(50 * 20), 50, 20), simplify = FALSE)
#' fit <- mcca(X, ncomp = 3)
#' }
#' @export
mcca.list <- function(data, preproc = multivarious::center(), ncomp = 2,
                      ridge = 1e-6, block_weights = NULL, use_future = FALSE, ...) {
  data <- multidesign::multiblock(data)
  mcca.multiblock(
    data = data,
    preproc = preproc,
    ncomp = ncomp,
    ridge = ridge,
    block_weights = block_weights,
    use_future = use_future,
    ...
  )
}

# -------------------------------------------------------------------------
# Internal helpers
# -------------------------------------------------------------------------

.mcca_as_dense <- function(x) {
  if (inherits(x, "Matrix")) return(as.matrix(x))
  x
}

.mcca_chol_solve <- function(R, B) {
  # Solve (t(R) %*% R) X = B where R is upper-triangular chol factor.
  y <- forwardsolve(t(R), B, upper.tri = FALSE, transpose = FALSE)
  backsolve(R, y, upper.tri = TRUE, transpose = FALSE)
}

.mcca_fit_one_block <- function(X, ridge_mult, block_weight, block_name, n) {
  K <- .mcca_as_dense(Matrix::tcrossprod(X))
  if (!all(is.finite(K))) stop("Non-finite values encountered in block Gram matrix.")

  scale_k <- mean(diag(K))
  if (!is.finite(scale_k) || scale_k <= 0) scale_k <- 1

  kappa <- ridge_mult * scale_k
  if (!is.finite(kappa) || kappa < 0) {
    stop("ridge must be non-negative and finite.")
  }

  # If ridge is zero, try the unregularized system first (may be ill-posed after
  # centering); fall back to a tiny ridge for numerical invertibility.
  K_reg <- NULL
  R <- NULL
  if (kappa == 0) {
    R0 <- tryCatch(chol(K), error = function(e) NULL)
    if (!is.null(R0)) {
      K_reg <- K
      R <- R0
    } else {
      kappa <- 1e-8 * scale_k
      warning(
        sprintf(
          "mcca: ridge=0 yields a singular system for block '%s'; using ridge=1e-8 (scaled).",
          block_name
        ),
        call. = FALSE
      )
      K_reg <- K + kappa * diag(n)
    }
  } else {
    K_reg <- K + kappa * diag(n)
  }

  if (is.null(R)) {
    R <- tryCatch(
      chol(K_reg),
      error = function(e) {
        # Retry with increased ridge if needed
        kappa2 <- max(kappa, 1e-6 * scale_k) * 100
        K_reg2 <- K + kappa2 * diag(n)
        R2 <- chol(K_reg2)
        warning(
          sprintf(
            "mcca: increased ridge for block '%s' from %.3g to %.3g for numerical stability.",
            block_name, kappa, kappa2
          ),
          call. = FALSE
        )
        kappa <<- kappa2
        R2
      }
    )
  }

  # P = (K + kappa I)^{-1} K; K commutes with K + kappa I, so this equals K (K + kappa I)^{-1}.
  P <- .mcca_chol_solve(R, K)
  P <- (P + t(P)) / 2

  list(
    P = block_weight * P,
    P_raw = P,
    chol = R,
    kappa = kappa,
    ridge_mult = ridge_mult,
    scale_k = scale_k,
    block_weight = block_weight
  )
}

# Lift a block-local n_k x n_k matrix into N x N reference-row space via
# row-index grouping. Used by aligned_mcca / anchored_mcca.
.mcca_lift_to_N <- function(M_local, idx, N) {
  rs1 <- rowsum(M_local, idx, reorder = FALSE)
  rs2 <- rowsum(t(rs1), idx, reorder = FALSE)
  out <- matrix(0, N, N)
  g_rows <- as.integer(rownames(rs2))
  g_cols <- as.integer(rownames(rs1))
  out[g_rows, g_cols] <- rs2
  out
}

# Top-ncomp symmetric eigendecomposition with RSpectra fallback to base eigen.
.mcca_top_eigen <- function(H, ncomp, N) {
  ncomp_eff <- min(as.integer(ncomp), N)
  eig <- NULL
  if (ncomp_eff < N) {
    eig <- tryCatch(
      RSpectra::eigs_sym(H, k = ncomp_eff, which = "LA"),
      error = function(e) NULL
    )
  }
  if (is.null(eig)) {
    full <- eigen(H, symmetric = TRUE)
    eig <- list(
      values = full$values[seq_len(ncomp_eff)],
      vectors = full$vectors[, seq_len(ncomp_eff), drop = FALSE]
    )
  }
  ord <- order(eig$values, decreasing = TRUE)
  list(values = eig$values[ord], vectors = eig$vectors[, ord, drop = FALSE])
}

# Per-component balanced-MCCA solver: projected gradient ascent with Armijo
# line search on the geometric-mean objective, with orthogonal deflation
# across components.
#
# Objective (per component):
#   phi(g) = Σ_b β_b * log(g^T M_b g)        subject to ||g||^2 = 1,
#                                             g ⊥ g_1..g_{k-1}
# Gradient:
#   ∇phi(g) = Σ_b β_b · 2 M_b g / (g^T M_b g)
#           = 2 H_w g,                         w_b = β_b / (g^T M_b g)
# Stationarity says w_b · (g^T M_b g) = β_b (per-block contribution ∝ β).
# Uniform β → equal contributions across blocks (no single block dominates).
#
# Unlike a pure IRLS fixed-point iteration — which tends to limit-cycle when
# different weight regimes prefer different eigenvectors — gradient ascent
# with Armijo line search guarantees monotone improvement in phi and
# converges to a local maximum.
#
# Returns per-component eigenpairs, per-component weights, a summary weight
# vector (geometric mean across components), and iteration diagnostics.
.mcca_irls_balance <- function(M_list,
                               ncomp,
                               target = NULL,
                               w_init = NULL,
                               max_iter = 50L,
                               tol = 1e-6,
                               damping = 0.5,       # kept for signature stability; unused by GA
                               step_init = 0.5,
                               step_shrink = 0.5,
                               armijo_c = 1e-4,
                               max_line_search = 25L,
                               verbose = FALSE) {
  B <- length(M_list)
  N <- nrow(M_list[[1]])

  if (is.null(target)) target <- rep(1, B)
  target <- as.numeric(target)
  chk::chk_true(length(target) == B)
  chk::chk_true(all(is.finite(target)))
  chk::chk_true(all(target >= 0))
  if (!any(target > 0)) stop("Balanced target must have at least one positive entry.", call. = FALSE)
  target <- target / mean(target)

  if (is.null(w_init)) w_init <- rep(1, B)
  w_init <- as.numeric(w_init)
  chk::chk_true(length(w_init) == B)
  if (any(!is.finite(w_init)) || any(w_init < 0)) stop("Balanced w_init must be finite and non-negative.", call. = FALSE)
  if (!any(w_init > 0)) stop("Balanced w_init must have at least one positive entry.", call. = FALSE)
  w_init <- w_init / mean(w_init)

  ncomp_eff <- min(as.integer(ncomp), N)

  # Objective phi(g) and its Euclidean gradient.
  contribs <- function(g) {
    vapply(M_list, function(Mb) as.numeric(crossprod(g, Mb %*% g)), numeric(1))
  }
  objective <- function(g) {
    c_b <- contribs(g)
    if (any(c_b <= 0) || any(!is.finite(c_b))) return(-Inf)
    sum(target * log(c_b))
  }
  gradient <- function(g, c_b = NULL) {
    if (is.null(c_b)) c_b <- contribs(g)
    c_safe <- pmax(c_b, 1e-12)
    coeff <- target / c_safe
    grad <- rep(0, N)
    for (b in seq_len(B)) grad <- grad + coeff[b] * (2 * as.numeric(M_list[[b]] %*% g))
    grad
  }

  G <- matrix(0, N, ncomp_eff)
  lambda <- numeric(ncomp_eff)
  weights_per_comp <- matrix(0, nrow = ncomp_eff, ncol = B)
  colnames(weights_per_comp) <- names(M_list)
  converged_per_comp <- logical(ncomp_eff)
  iters_per_comp <- integer(ncomp_eff)
  trace <- vector("list", ncomp_eff)

  Q <- matrix(0, N, 0)  # orthonormal basis of previously found components

  # Project vector v onto the tangent space of the sphere at g, intersected
  # with the orthogonal complement of Q (so deflation is built in).
  project_tangent <- function(g, v) {
    if (ncol(Q) > 0) v <- v - Q %*% (crossprod(Q, v))
    v - as.numeric(crossprod(g, v)) * g
  }

  # Retract from tangent space back onto the unit sphere ∩ Q⊥.
  retract <- function(g_new) {
    if (ncol(Q) > 0) g_new <- g_new - Q %*% (crossprod(Q, g_new))
    nrm <- sqrt(sum(g_new^2))
    if (nrm < 1e-12) stop("Retraction to unit sphere failed (zero-norm iterate).", call. = FALSE)
    g_new / nrm
  }

  for (k in seq_len(ncomp_eff)) {
    # Warm start: leading eigenvector of H_{w_init} restricted to span(Q)^⊥.
    H_init <- matrix(0, N, N)
    for (b in seq_len(B)) H_init <- H_init + w_init[b] * M_list[[b]]
    H_init <- (H_init + t(H_init)) / 2
    if (ncol(Q) > 0) {
      AQ <- H_init %*% Q
      QtAQ <- crossprod(Q, AQ)
      H_init <- H_init - tcrossprod(Q, AQ) - tcrossprod(AQ, Q) + Q %*% QtAQ %*% t(Q)
      H_init <- (H_init + t(H_init)) / 2
    }
    g <- .mcca_top_eigen(H_init, 1L, N)$vectors[, 1L, drop = TRUE]
    g <- retract(g)

    obj <- objective(g)
    if (!is.finite(obj)) {
      # Degenerate init; fall back to a random orthogonal unit vector in Q⊥.
      g <- retract(rnorm(N))
      obj <- objective(g)
    }

    step <- as.numeric(step_init)
    comp_trace <- vector("list", 0)
    conv_k <- FALSE
    grad_norm_init <- NA_real_

    for (iter in seq_len(max_iter)) {
      c_b <- contribs(g)
      grad <- gradient(g, c_b = c_b)
      pgrad <- project_tangent(g, grad)
      grad_norm <- sqrt(sum(pgrad^2))
      if (iter == 1L) grad_norm_init <- grad_norm

      imbalance <- {
        c_pos <- c_b[c_b > 0]
        if (length(c_pos) > 1) max(c_pos) / min(c_pos) else 1
      }

      comp_trace[[iter]] <- list(
        comp = k, iter = iter, g = g, c_b = c_b,
        objective = obj, grad_norm = grad_norm, imbalance = imbalance, step = step
      )

      if (isTRUE(verbose)) {
        message(sprintf(
          "anchored_mcca balanced comp %d iter %d: phi = %.4g, grad_norm = %.3g, imbalance = %.3g",
          k, iter, obj, grad_norm, imbalance
        ))
      }

      # Relative gradient-norm stopping rule — more stable than absolute tol
      # for an objective whose scale depends on the block contributions.
      if (!is.finite(grad_norm)) {
        break
      }
      rel_grad <- grad_norm / max(grad_norm_init, 1e-12)
      if (grad_norm < tol || rel_grad < tol) {
        conv_k <- TRUE
        break
      }

      # Armijo backtracking line search: find step with
      #   phi(retract(g + step·pgrad)) >= phi(g) + armijo_c · step · ||pgrad||^2
      accepted <- FALSE
      step_try <- step
      obj_prev <- obj
      for (ls in seq_len(max_line_search)) {
        g_trial <- retract(g + step_try * pgrad)
        obj_trial <- objective(g_trial)
        if (is.finite(obj_trial) &&
            obj_trial >= obj + armijo_c * step_try * grad_norm^2) {
          g <- g_trial
          obj <- obj_trial
          step <- min(step_try / step_shrink, 1.0)  # tentative increase for next iter
          accepted <- TRUE
          break
        }
        step_try <- step_try * step_shrink
      }

      if (!accepted) {
        # Cannot find an improving step — treat as a local max.
        conv_k <- TRUE
        break
      }

      # Also accept convergence when relative objective change is tiny.
      if (is.finite(obj_prev) && abs(obj - obj_prev) <= tol * max(1, abs(obj_prev))) {
        conv_k <- TRUE
        break
      }
    }

    c_final <- contribs(g)
    c_final[!is.finite(c_final) | c_final < 0] <- 0
    c_safe <- pmax(c_final, 1e-12)

    w_star <- target / c_safe
    w_star <- w_star / mean(w_star)

    # λ_k = g^T H_{w*} g = Σ_b w*_b · (g^T M_b g) = Σ_b target_b = B (const at fixed point
    # when target mean is 1). Report the raw Rayleigh quotient on the anchor H
    # (uniform-weight sum of contributions) as a more interpretable eigenvalue.
    lam <- sum(c_final)

    G[, k] <- g
    lambda[k] <- lam
    weights_per_comp[k, ] <- w_star
    converged_per_comp[k] <- conv_k
    iters_per_comp[k] <- length(comp_trace)
    trace[[k]] <- comp_trace

    Q <- qr.Q(qr(cbind(Q, g)))
  }

  # Summary weights: geometric mean across components (per block).
  weights_summary <- exp(colMeans(log(pmax(weights_per_comp, 1e-12))))
  weights_summary <- weights_summary / mean(weights_summary)

  list(
    G = G,
    lambda = lambda,
    weights = weights_summary,
    weights_per_comp = weights_per_comp,
    trace = trace,
    converged = all(converged_per_comp),
    converged_per_comp = converged_per_comp,
    iters = sum(iters_per_comp),
    iters_per_comp = iters_per_comp
  )
}

# MFA-style per-block weight: 1 / (first singular value)^2 of the preprocessed block.
.mcca_mfa_weights <- function(Xp, ridge_floor = 1e-8) {
  first_sv <- function(mat) {
    method <- .muscal_svd_method(mat, ncomp = 1)
    tryCatch(
      multivarious::svd_wrapper(mat, ncomp = 1, method = method)$sdev[1],
      error = function(e) multivarious::svd_wrapper(mat, ncomp = 1, method = "base")$sdev[1]
    )
  }
  sapply(Xp, function(M) {
    sv <- first_sv(M)
    1 / (sv^2 + ridge_floor)
  })
}

.mcca_map_fun <- function(use_future) {
  if (isTRUE(use_future)) {
    if (!requireNamespace("furrr", quietly = TRUE)) {
      stop("use_future = TRUE requires the 'furrr' package.", call. = FALSE)
    }
    function(.x, .f) {
      furrr::future_map(.x, .f, .options = furrr::furrr_options(seed = TRUE))
    }
  } else {
    function(.x, .f) lapply(.x, .f)
  }
}

.mcca_make_anchored_mcca_refit_fn <- function(preproc,
                                              ncomp,
                                              normalization,
                                              alpha,
                                              ridge,
                                              max_iter,
                                              tol,
                                              use_future,
                                              fit_dots) {
  force(preproc)
  force(ncomp)
  force(normalization)
  force(alpha)
  force(ridge)
  force(max_iter)
  force(tol)
  force(use_future)
  force(fit_dots)

  function(data) {
    do.call(
      anchored_mcca,
      c(
        list(
          Y = data$Y,
          X = data$X,
          row_index = data$row_index,
          preproc = preproc,
          ncomp = ncomp,
          normalization = normalization,
          alpha = alpha,
          ridge = ridge,
          max_iter = max_iter,
          tol = tol,
          verbose = FALSE,
          use_future = use_future
        ),
        fit_dots
      )
    )
  }
}

#' @md
#' @rdname mcca
#' @details
#' The `mcca.multiblock` method implements a ridge-stabilized MAXVAR generalized
#' CCA. Internally it computes a compromise score space by eigendecomposition of
#' a weighted sum of block projection operators:
#' \deqn{H = \sum_s w_s \, X_s (X_s^T X_s + \kappa_s I)^{-1} X_s^T.}
#'
#' The ridge parameter is applied per block as \eqn{\kappa_s = \text{ridge}_s \cdot
#' \bar{d}_s}, where \eqn{\bar{d}_s} is the mean diagonal of the block Gram matrix
#' \eqn{X_s X_s^T}. This scaling makes the default ridge behave reasonably across
#' different variable scalings.
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' blocks <- replicate(3, matrix(rnorm(20 * 50), 20, 50), simplify = FALSE)
#' fit <- mcca(blocks, ncomp = 2)
#' stopifnot(ncol(multivarious::scores(fit)) == 2)
#' }
#' @export
mcca.multiblock <- function(data, preproc = multivarious::center(), ncomp = 2,
                            ridge = 1e-6, block_weights = NULL,
                            use_future = FALSE, ...) {
  fit_call <- match.call(expand.dots = FALSE)
  fit_dots <- list(...)
  chk::chk_true(length(data) > 1)
  for (i in seq_along(data)) {
    chk::chkor_vld(chk::vld_matrix(data[[i]]), chk::vld_s4_class(data[[i]], "Matrix"))
  }

  nrs <- sapply(data, nrow)
  chk::chk_true(all(nrs == nrs[1]))
  n <- nrs[1]

  S <- length(data)
  if (is.null(names(data))) names(data) <- paste0("B", seq_len(S))
  data_refit <- lapply(data, function(x) as.matrix(x))

  if (is.null(block_weights)) block_weights <- rep(1, S)
  chk::chk_numeric(block_weights)
  chk::chk_true(length(block_weights) == S)
  chk::chk_true(all(is.finite(block_weights)))
  chk::chk_true(all(block_weights >= 0))

  # Preprocess blocks independently
  prep <- prepare_block_preprocessors(data, preproc, check_consistent_ncol = FALSE)
  strata <- prep$Xp
  proclist <- prep$proclist
  fitted_proclist <- .muscal_materialize_block_preprocessors(data, proclist)

  # Block indices for concatenated loadings space
  block_indices <- list()
  ind <- 1
  for (i in seq_along(strata)) {
    block_indices[[i]] <- seq(ind, ind + ncol(strata[[i]]) - 1)
    ind <- ind + ncol(strata[[i]])
  }
  names(block_indices) <- names(data)

  proc <- multivarious::concat_pre_processors(fitted_proclist, block_indices)

  # Ridge vector (scaled per block)
  if (length(ridge) == 1) ridge <- rep(ridge, S)
  chk::chk_numeric(ridge)
  chk::chk_true(length(ridge) == S)
  chk::chk_true(all(is.finite(ridge)))
  chk::chk_true(all(ridge >= 0))

  idx <- seq_along(strata)
  map_fun <- if (isTRUE(use_future)) {
    function(.x, .f) {
      furrr::future_map(.x, .f, .options = furrr::furrr_options(seed = TRUE))
    }
  } else {
    function(.x, .f) lapply(.x, .f)
  }

  block_fits <- map_fun(idx, function(i) {
    .mcca_fit_one_block(
      X = strata[[i]],
      ridge_mult = ridge[i],
      block_weight = block_weights[i],
      block_name = names(data)[i],
      n = n
    )
  })

  H <- Reduce(`+`, lapply(block_fits, `[[`, "P"))
  H <- (H + t(H)) / 2

  # Cap ncomp
  chk::chk_numeric(ncomp)
  chk::chk_true(length(ncomp) == 1)
  chk::chk_true(is.finite(ncomp))
  chk::chk_gte(ncomp, 1)
  if (abs(ncomp - round(ncomp)) > 1e-8) {
    stop("ncomp must be an integer-like value.", call. = FALSE)
  }
  ncomp_eff <- min(as.integer(round(ncomp)), n)

  eig <- NULL
  if (ncomp_eff < n) {
    eig <- tryCatch(
      RSpectra::eigs_sym(H, k = ncomp_eff, which = "LA"),
      error = function(e) NULL
    )
  }
  if (is.null(eig)) {
    full <- eigen(H, symmetric = TRUE)
    eig <- list(
      values = full$values[seq_len(ncomp_eff)],
      vectors = full$vectors[, seq_len(ncomp_eff), drop = FALSE]
    )
  }

  ord <- order(eig$values, decreasing = TRUE)
  lambda <- eig$values[ord]
  G <- eig$vectors[, ord, drop = FALSE]

  lambda[!is.finite(lambda)] <- 0
  lambda[lambda < 0] <- 0
  sdev <- sqrt(lambda)

  # Scale scores like PCA: S = U * sdev so that crossprod(S) = diag(lambda).
  S_scores <- G %*% diag(sdev, nrow = length(sdev), ncol = length(sdev))
  colnames(S_scores) <- paste0("Comp", seq_len(ncol(S_scores)))

  # Build concatenated loading matrix v (so that X_concat %*% v == scores)
  p_tot <- sum(vapply(strata, ncol, integer(1)))
  v_concat <- matrix(0, p_tot, ncomp_eff)

  partial_scores <- vector("list", S)
  names(partial_scores) <- names(data)

  W_list <- vector("list", S)
  names(W_list) <- names(data)

  for (i in seq_along(strata)) {
    X <- strata[[i]]
    R <- block_fits[[i]]$chol
    A <- .mcca_chol_solve(R, G) # (K + kappa I)^{-1} G
    W_raw <- crossprod(X, A)    # p_i x k (canonical weights in feature space)
    W_list[[i]] <- W_raw

    # Scale so that concatenated scores match S_scores.
    inv_sdev <- ifelse(sdev > 1e-12, 1 / sdev, 0)
    W_weighted <- block_fits[[i]]$block_weight * W_raw
    V_block <- sweep(W_weighted, 2, inv_sdev, `*`)
    v_concat[block_indices[[i]], ] <- V_block
    partial_scores[[i]] <- X %*% V_block
  }

  cor_loadings <- tryCatch(
    {
      X_concat <- do.call(cbind, strata)
      C <- suppressWarnings(stats::cor(X_concat, S_scores))
      C[!is.finite(C)] <- 0
      C
    },
    error = function(e) NULL
  )

  rv <- tryCatch(
    {
      M <- compute_sim_mat(strata, function(x1, x2) MatrixCorrelation::RV(x1, x2))
      diag(M) <- 1
      rownames(M) <- colnames(M) <- names(data)
      M
    },
    error = function(e) NULL
  )

  rv2 <- tryCatch(
    {
      M <- compute_sim_mat(strata, function(x1, x2) MatrixCorrelation::RV2(x1, x2))
      diag(M) <- 1
      rownames(M) <- colnames(M) <- names(data)
      M
    },
    error = function(e) NULL
  )

  out <- multivarious::multiblock_biprojector(
    v = v_concat,
    s = S_scores,
    sdev = sdev,
    preproc = proc,
    block_indices = block_indices,
    # MCCA metadata / diagnostics
    lambda = lambda,
    ridge = ridge,
    kappa = vapply(block_fits, `[[`, numeric(1), "kappa"),
    block_weights = block_weights,
    partial_scores = partial_scores,
    canonical_weights = W_list,
    cor_loadings = cor_loadings,
    rv = rv,
    rv2 = rv2,
    names = names(data),
    classes = "mcca"
  )

  .muscal_attach_fit_contract(
    out,
    method = "mcca",
    task = "reconstruction",
    oos_types = c("scores", "reconstruction"),
    fit_call = fit_call,
    refit_supported = TRUE,
    prediction_target = "blocks",
    refit = .muscal_make_refit_spec(
      data = data_refit,
      fit_fn = function(data) {
        do.call(
          mcca,
          c(
            list(
              data = data,
              preproc = preproc,
              ncomp = ncomp,
              ridge = ridge,
              block_weights = block_weights,
              use_future = use_future
            ),
            fit_dots
          )
        )
      },
      bootstrap_fn = function(data) {
        n <- nrow(data[[1L]])
        idx <- sample.int(n, size = n, replace = TRUE)
        lapply(data, function(block) block[idx, , drop = FALSE])
      },
      permutation_fn = function(data) {
        lapply(data, function(block) block[sample.int(nrow(block)), , drop = FALSE])
      },
      resample_unit = "rows"
    )
  )
}

#' Predict from an MCCA Fit
#'
#' @param object A fitted `mcca` object.
#' @param new_data Optional matrix/data.frame or list of blocks with the same
#'   structure as the training data.
#' @param type One of `"scores"` or `"reconstruction"`.
#' @param ... Additional arguments passed to the underlying projection helpers.
#'
#' @return A numeric matrix of projected scores or reconstructed observations.
#' @export
predict.mcca <- function(object, new_data = NULL,
                         type = c("scores", "reconstruction"), ...) {
  type <- match.arg(type)

  if (is.null(new_data)) {
    if (type == "scores") {
      return(multivarious::scores(object))
    }
    stop("`new_data` must be supplied when type = 'reconstruction'.", call. = FALSE)
  }

  if (type == "scores") {
    return(multivarious::project(object, new_data, ...))
  }

  multivarious::reconstruct_new(object, new_data, ...)
}

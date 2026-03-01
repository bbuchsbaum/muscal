#' Recommend Practical Defaults for Bilinear Mixed
#'
#' @description
#' Provides data-adaptive defaults for [bilinear_mixed()] so users can
#' avoid manual specification of many tuning parameters.
#'
#' @param data A list of numeric connectivity matrices.
#' @param subject Subject identifier with length \code{length(data)}.
#' @param z Optional repeat-level design.
#' @param y Optional subject-level traits.
#' @param row_design Optional row-level covariates (same number of rows as each
#'   connectivity matrix).
#' @param mode One of \code{"auto"}, \code{"seed_axis"}, \code{"seed_repeat"},
#'   or \code{"both"}.
#' @param connectivity_type One of \code{"auto"}, \code{"cross"}, or
#'   \code{"symmetric"}.
#' @param profile One of \code{"fast"}, \code{"balanced"}, or \code{"adaptive"}.
#' @param max_r_seed Maximum suggested row rank.
#' @param max_r_roi Maximum suggested column rank.
#' @param max_k_subject Maximum suggested subject latent dimension.
#' @param sym_tol Numeric tolerance for symmetry detection when
#'   \code{connectivity_type = "auto"}.
#'
#' @return Named list of recommended arguments for [bilinear_mixed()].
#' @export
bilinear_mixed_recommend <- function(data,
                                            subject,
                                            z = NULL,
                                            y = NULL,
                                            row_design = NULL,
                                            mode = c("auto", "seed_axis", "seed_repeat", "both"),
                                            connectivity_type = c("auto", "cross", "symmetric"),
                                            profile = c("fast", "balanced", "adaptive"),
                                            max_r_seed = 6,
                                            max_r_roi = 20,
                                            max_k_subject = 5,
                                            sym_tol = 1e-8) {
  mode <- match.arg(mode)
  connectivity_type <- match.arg(connectivity_type)
  profile <- match.arg(profile)

  chk::chk_list(data)
  chk::chk_true(length(data) > 1)
  data <- lapply(data, as.matrix)
  s <- nrow(data[[1]])
  p <- ncol(data[[1]])
  same_dim <- vapply(data, function(x) nrow(x) == s && ncol(x) == p, logical(1))
  if (!all(same_dim)) stop("All matrices in `data` must have identical dimensions.", call. = FALSE)

  if (length(subject) != length(data)) {
    stop("`subject` must have length equal to length(data).", call. = FALSE)
  }
  subject <- factor(subject)
  n_subject <- nlevels(subject)

  is_square_all <- all(vapply(data, function(x) nrow(x) == ncol(x), logical(1)))
  is_symmetric_all <- FALSE
  if (is_square_all) {
    is_symmetric_all <- all(vapply(data, function(x) max(abs(x - t(x))) <= sym_tol, logical(1)))
  }
  conn_eff <- if (connectivity_type == "auto") {
    if (is_square_all && is_symmetric_all) "symmetric" else "cross"
  } else {
    connectivity_type
  }

  mode_eff <- if (mode == "auto") {
    if (!is.null(row_design)) "both" else "seed_axis"
  } else {
    mode
  }

  var_target <- switch(
    profile,
    fast = 0.75,
    balanced = 0.85,
    adaptive = 0.92
  )

  C_row <- Reduce(`+`, lapply(data, tcrossprod)) / length(data)
  C_col <- Reduce(`+`, lapply(data, crossprod)) / length(data)

  r_row <- .bc_rank_from_cov(C_row, target = var_target, cap = min(max_r_seed, s))
  r_col <- .bc_rank_from_cov(C_col, target = var_target, cap = min(max_r_roi, p))
  if (conn_eff == "symmetric") {
    rr <- min(r_row, r_col)
    r_row <- rr
    r_col <- rr
  }

  # Subject-level shrinkage proxy: within/between variance ratio.
  vec <- do.call(rbind, lapply(data, as.vector))
  subj_idx <- as.integer(subject)
  subj_mean <- rowsum(vec, subj_idx, reorder = FALSE) / as.vector(table(subject))
  between_var <- mean(apply(subj_mean, 2, stats::var), na.rm = TRUE)
  within_var <- mean((vec - subj_mean[subj_idx, , drop = FALSE])^2)
  ratio <- within_var / (between_var + 1e-8)
  ratio <- .bc_clamp(ratio, 1e-3, 1e2)

  Ysubj <- .bc_prepare_traits(y, levels(subject))
  k_default <- if (!is.null(Ysubj)) {
    min(max_k_subject, n_subject - 1L, ncol(Ysubj))
  } else {
    switch(profile, fast = 1L, balanced = 2L, adaptive = 3L)
  }
  k_default <- max(1L, as.integer(min(k_default, n_subject - 1L, r_col)))

  list(
    mode = mode_eff,
    connectivity_type = conn_eff,
    sym_tol = sym_tol,
    r_seed = r_row,
    r_roi = r_col,
    k_subject = k_default,
    include_seed_interactions = TRUE,
    lambda_w = .bc_clamp(0.10 * ratio, 1e-3, 1),
    lambda_t = .bc_clamp(ratio, 1e-3, 10),
    lambda_b = .bc_clamp(0.01 * ratio, 1e-4, 0.1),
    lambda_a = 1e-3,
    lambda_y = if (!is.null(Ysubj)) 1 else 0,
    lambda_repeat = 0,
    lambda_seed = 0,
    lambda_roi = 0,
    lambda_design_smooth = 0
  )
}

#' Easy Front-End for Bilinear Mixed
#'
#' @description
#' Lightweight wrapper around [bilinear_mixed()] that uses
#' [bilinear_mixed_recommend()] defaults. Optionally performs compact
#' subject-block CV tuning via [bilinear_mixed_tune()].
#'
#' @param data A list of numeric connectivity matrices.
#' @param subject Subject identifier.
#' @param z Optional repeat-level design.
#' @param y Optional subject-level traits.
#' @param row_design Optional row-level covariates.
#' @param mode One of \code{"auto"}, \code{"seed_axis"}, \code{"seed_repeat"},
#'   or \code{"both"}.
#' @param connectivity_type One of \code{"auto"}, \code{"cross"}, or
#'   \code{"symmetric"}.
#' @param profile One of \code{"fast"}, \code{"balanced"}, or \code{"adaptive"}.
#' @param tune Logical; if \code{TRUE}, run [bilinear_mixed_tune()].
#' @param tune_grid Optional grid for tuning; see [bilinear_mixed_tune()].
#' @param n_folds Number of subject-block CV folds for tuning.
#' @param metric One of \code{"auto"}, \code{"reconstruction"}, or
#'   \code{"trait_r2"}.
#' @param seed Integer random seed for fold assignment.
#' @param verbose Logical verbosity.
#' @param ... Optional overrides passed to fitting/tuning.
#'
#' @return Fitted \code{"bilinear_mixed"} object if \code{tune = FALSE},
#'   otherwise a \code{"bilinear_mixed_tuning"} object.
#' @export
bilinear_mixed_easy <- function(data,
                                       subject,
                                       z = NULL,
                                       y = NULL,
                                       row_design = NULL,
                                       mode = c("auto", "seed_axis", "seed_repeat", "both"),
                                       connectivity_type = c("auto", "cross", "symmetric"),
                                       profile = c("fast", "balanced", "adaptive"),
                                       tune = FALSE,
                                       tune_grid = NULL,
                                       n_folds = 3,
                                       metric = c("auto", "reconstruction", "trait_r2"),
                                       seed = 1,
                                       verbose = FALSE,
                                       ...) {
  mode <- match.arg(mode)
  connectivity_type <- match.arg(connectivity_type)
  profile <- match.arg(profile)
  metric <- match.arg(metric)
  chk::chk_flag(tune)
  chk::chk_flag(verbose)

  if (isTRUE(tune)) {
    return(
      bilinear_mixed_tune(
        data = data,
        subject = subject,
        z = z,
        y = y,
        row_design = row_design,
        mode = mode,
        connectivity_type = connectivity_type,
        profile = profile,
        grid = tune_grid,
        n_folds = n_folds,
        metric = metric,
        seed = seed,
        verbose = verbose,
        ...
      )
    )
  }

  rec <- bilinear_mixed_recommend(
    data = data,
    subject = subject,
    z = z,
    y = y,
    row_design = row_design,
    mode = mode,
    connectivity_type = connectivity_type,
    profile = profile
  )
  rec <- utils::modifyList(rec, list(...))

  do.call(
    bilinear_mixed,
    c(
      list(
        data = data,
        subject = subject,
        z = z,
        y = y,
        row_design = row_design,
        verbose = verbose
      ),
      rec
    )
  )
}

#' Tune Bilinear Mixed Hyperparameters via Subject-Block CV
#'
#' @description
#' Tunes a compact set of key parameters (`r_seed`, `r_roi`, `k_subject`,
#' optional `lambda_y`) using subject-blocked cross-validation.
#'
#' @param data A list of numeric connectivity matrices.
#' @param subject Subject identifier.
#' @param z Optional repeat-level design.
#' @param y Optional subject-level traits.
#' @param row_design Optional row-level covariates.
#' @param mode One of \code{"auto"}, \code{"seed_axis"}, \code{"seed_repeat"},
#'   or \code{"both"}.
#' @param connectivity_type One of \code{"auto"}, \code{"cross"}, or
#'   \code{"symmetric"}.
#' @param profile One of \code{"fast"}, \code{"balanced"}, or \code{"adaptive"}.
#' @param grid Optional candidate grid. Either:
#'   \itemize{
#'     \item a \code{data.frame} with columns among
#'       \code{r_seed}, \code{r_roi}, \code{k_subject}, \code{lambda_y}
#'     \item or a named list of vectors for those fields.
#'   }
#' @param n_folds Number of subject-block CV folds.
#' @param metric One of \code{"auto"}, \code{"reconstruction"}, or
#'   \code{"trait_r2"}.
#' @param seed Random seed for fold assignment.
#' @param verbose Logical verbosity.
#' @param ... Additional fixed arguments forwarded to [bilinear_mixed()].
#'
#' @return A list of class \code{"bilinear_mixed_tuning"} with
#'   \code{best_params}, \code{results}, and refit \code{fit}.
#' @export
bilinear_mixed_tune <- function(data,
                                       subject,
                                       z = NULL,
                                       y = NULL,
                                       row_design = NULL,
                                       mode = c("auto", "seed_axis", "seed_repeat", "both"),
                                       connectivity_type = c("auto", "cross", "symmetric"),
                                       profile = c("fast", "balanced", "adaptive"),
                                       grid = NULL,
                                       n_folds = 3,
                                       metric = c("auto", "reconstruction", "trait_r2"),
                                       seed = 1,
                                       verbose = FALSE,
                                       ...) {
  mode <- match.arg(mode)
  connectivity_type <- match.arg(connectivity_type)
  profile <- match.arg(profile)
  metric <- match.arg(metric)
  chk::chk_number(n_folds)
  chk::chk_gte(n_folds, 2)
  n_folds <- as.integer(round(n_folds))
  chk::chk_flag(verbose)

  rec <- bilinear_mixed_recommend(
    data = data,
    subject = subject,
    z = z,
    y = y,
    row_design = row_design,
    mode = mode,
    connectivity_type = connectivity_type,
    profile = profile
  )
  fixed <- utils::modifyList(rec, list(...))

  metric_eff <- if (metric == "auto") {
    if (is.null(y)) "reconstruction" else "trait_r2"
  } else {
    metric
  }
  if (metric_eff == "trait_r2" && is.null(y)) {
    stop("metric='trait_r2' requires subject-level traits `y`.", call. = FALSE)
  }

  grid_df <- .bc_make_tune_grid(grid, fixed, has_y = !is.null(y))
  subj <- factor(subject)
  subj_levels <- levels(subj)
  n_subject <- length(subj_levels)
  n_folds_eff <- min(n_folds, n_subject)
  folds <- .bc_subject_folds(subj_levels, n_folds = n_folds_eff, seed = seed)

  results <- data.frame(
    candidate = seq_len(nrow(grid_df)),
    r_seed = grid_df$r_seed,
    r_roi = grid_df$r_roi,
    k_subject = grid_df$k_subject,
    lambda_y = grid_df$lambda_y,
    score = if (metric_eff == "trait_r2") -Inf else Inf,
    n_success = 0L,
    stringsAsFactors = FALSE
  )

  for (gi in seq_len(nrow(grid_df))) {
    cand <- as.list(grid_df[gi, , drop = FALSE])
    params <- utils::modifyList(fixed, cand)
    fold_scores <- numeric(0)

    for (fi in seq_along(folds)) {
      test_subj <- folds[[fi]]
      idx_test <- which(subj %in% test_subj)
      idx_train <- setdiff(seq_along(data), idx_test)
      if (length(idx_train) < 2 || length(idx_test) < 1) next

      z_train <- .bc_subset_rows(z, idx_train)
      z_test <- .bc_subset_rows(z, idx_test)
      y_train <- .bc_subset_traits_input(y, subject_all = subject, keep_subjects = unique(subject[idx_train]))
      y_test <- .bc_subset_traits_input(y, subject_all = subject, keep_subjects = unique(subject[idx_test]))

      fit_i <- tryCatch(
        do.call(
          bilinear_mixed,
          c(
            list(
              data = data[idx_train],
              subject = subject[idx_train],
              z = z_train,
              y = y_train,
              row_design = row_design,
              verbose = FALSE
            ),
            params
          )
        ),
        error = function(e) e
      )
      if (inherits(fit_i, "error")) {
        if (isTRUE(verbose)) {
          message(sprintf("candidate %d fold %d failed: %s", gi, fi, conditionMessage(fit_i)))
        }
        next
      }

      score_i <- .bc_eval_fit(
        fit = fit_i,
        data = data[idx_test],
        subject = subject[idx_test],
        z = z_test,
        y = y_test,
        row_design = row_design,
        metric = metric_eff
      )
      if (is.finite(score_i)) {
        fold_scores <- c(fold_scores, score_i)
      }
    }

    if (length(fold_scores) > 0) {
      results$n_success[gi] <- length(fold_scores)
      results$score[gi] <- mean(fold_scores)
    }
  }

  valid <- results$n_success > 0 & is.finite(results$score)
  if (!any(valid)) {
    stop("All tuning candidates failed.", call. = FALSE)
  }

  best_idx <- if (metric_eff == "trait_r2") {
    which.max(results$score)
  } else {
    which.min(results$score)
  }
  best_params <- as.list(grid_df[best_idx, , drop = FALSE])
  final_params <- utils::modifyList(fixed, best_params)

  fit_full <- do.call(
    bilinear_mixed,
    c(
      list(
        data = data,
        subject = subject,
        z = z,
        y = y,
        row_design = row_design,
        verbose = verbose
      ),
      final_params
    )
  )

  out <- list(
    metric = metric_eff,
    defaults = rec,
    best_params = final_params,
    best_score = results$score[best_idx],
    results = results,
    grid = grid_df,
    fit = fit_full
  )
  class(out) <- "bilinear_mixed_tuning"
  out
}

#' @export
print.bilinear_mixed_tuning <- function(x, n = 10, ...) {
  n <- max(1L, as.integer(round(n)))
  cat("Bilinear Mixed Tuning\n")
  cat("  metric:", x$metric, "\n")
  cat("  best score:", signif(x$best_score, 4), "\n")

  if (!is.null(x$best_params) && length(x$best_params) > 0) {
    cat("  best params:\n")
    bp <- x$best_params[c("r_seed", "r_roi", "k_subject", "lambda_y")]
    bp <- bp[!vapply(bp, is.null, logical(1))]
    for (nm in names(bp)) {
      cat("   -", nm, "=", format(bp[[nm]], scientific = FALSE), "\n")
    }
  }

  if (is.data.frame(x$results) && nrow(x$results) > 0) {
    ord <- if (identical(x$metric, "trait_r2")) {
      order(-x$results$score, -x$results$n_success)
    } else {
      order(x$results$score, -x$results$n_success)
    }
    keep_cols <- intersect(
      c("candidate", "r_seed", "r_roi", "k_subject", "lambda_y", "score", "n_success"),
      colnames(x$results)
    )
    show <- utils::head(x$results[ord, keep_cols, drop = FALSE], n = n)
    cat("\nTop candidates:\n")
    print(show, row.names = FALSE)
  }

  invisible(x)
}

# -----------------------------------------------------------------------------
# Internal helpers
# -----------------------------------------------------------------------------

.bc_rank_from_cov <- function(C, target = 0.85, cap = nrow(C)) {
  cap <- max(1L, as.integer(cap))
  C <- as.matrix(C)
  C <- (C + t(C)) / 2
  C[!is.finite(C)] <- 0

  ev <- tryCatch(
    eigen(C, symmetric = TRUE, only.values = TRUE)$values,
    error = function(e) numeric(0)
  )
  if (length(ev) == 0) return(1L)
  ev <- pmax(ev, 0)
  if (sum(ev) <= 0 || !is.finite(sum(ev))) return(1L)
  cum <- cumsum(ev) / sum(ev)
  k <- which(cum >= target)[1]
  if (is.na(k)) k <- length(ev)
  as.integer(min(max(1L, k), cap))
}

.bc_clamp <- function(x, lo, hi) {
  pmin(pmax(x, lo), hi)
}

.bc_make_tune_grid <- function(grid, base, has_y = FALSE) {
  if (is.null(grid)) {
    rs <- unique(pmax(1L, c(base$r_seed - 1L, base$r_seed, base$r_seed + 1L)))
    rr <- unique(pmax(1L, c(base$r_roi - 1L, base$r_roi, base$r_roi + 1L)))
    ks <- unique(pmax(1L, c(base$k_subject - 1L, base$k_subject, base$k_subject + 1L)))
    ly <- if (has_y) unique(c(max(base$lambda_y / 3, 1e-3), base$lambda_y, base$lambda_y * 3)) else 0

    if (identical(base$connectivity_type, "symmetric")) {
      r_shared <- unique(pmax(1L, c(base$r_roi - 1L, base$r_roi, base$r_roi + 1L)))
      g <- expand.grid(r_shared = r_shared, k_subject = ks, lambda_y = ly, KEEP.OUT.ATTRS = FALSE)
      g$r_seed <- g$r_shared
      g$r_roi <- g$r_shared
      g$r_shared <- NULL
      return(g[, c("r_seed", "r_roi", "k_subject", "lambda_y"), drop = FALSE])
    }

    return(expand.grid(
      r_seed = rs,
      r_roi = rr,
      k_subject = ks,
      lambda_y = ly,
      KEEP.OUT.ATTRS = FALSE
    ))
  }

  if (is.data.frame(grid)) {
    g <- grid
  } else if (is.list(grid)) {
    g <- expand.grid(grid, KEEP.OUT.ATTRS = FALSE)
  } else {
    stop("`grid` must be NULL, a data.frame, or a named list.", call. = FALSE)
  }

  needed <- c("r_seed", "r_roi", "k_subject", "lambda_y")
  for (nm in needed) {
    if (!nm %in% colnames(g)) g[[nm]] <- base[[nm]]
  }
  g <- g[, needed, drop = FALSE]
  g$r_seed <- pmax(1L, as.integer(round(g$r_seed)))
  g$r_roi <- pmax(1L, as.integer(round(g$r_roi)))
  g$k_subject <- pmax(1L, as.integer(round(g$k_subject)))
  g$lambda_y <- as.numeric(g$lambda_y)
  g
}

.bc_subject_folds <- function(subj_levels, n_folds = 3, seed = 1) {
  set.seed(seed)
  s <- sample(subj_levels, length(subj_levels), replace = FALSE)
  split(s, rep(seq_len(n_folds), length.out = length(s)))
}

.bc_subset_rows <- function(x, idx) {
  if (is.null(x)) return(NULL)
  if (is.vector(x) && !is.list(x)) return(x[idx])
  x[idx, , drop = FALSE]
}

.bc_subset_traits_input <- function(y, subject_all, keep_subjects) {
  if (is.null(y)) return(NULL)
  keep_subjects <- unique(as.character(keep_subjects))
  subj_levels <- unique(as.character(subject_all))

  if (is.vector(y) && !is.list(y)) {
    if (!is.null(names(y))) {
      return(y[keep_subjects])
    }
    idx <- match(keep_subjects, subj_levels)
    return(as.vector(y[idx]))
  }

  Y <- as.matrix(y)
  if (!is.null(rownames(Y))) {
    return(Y[keep_subjects, , drop = FALSE])
  }
  idx <- match(keep_subjects, subj_levels)
  Y[idx, , drop = FALSE]
}

.bc_eval_fit <- function(fit, data, subject, z, y, row_design, metric = c("reconstruction", "trait_r2")) {
  metric <- match.arg(metric)
  if (metric == "reconstruction") {
    return(.bc_reconstruction_mse(fit, data, subject, z, row_design))
  }
  .bc_trait_r2(fit, data, subject, z, y, row_design)
}

.bc_reconstruction_mse <- function(fit, data, subject, z, row_design) {
  if (!is.null(fit$axis)) {
    return(.bc_reconstruction_mse_axis(fit, data, subject, z))
  }
  .bc_reconstruction_mse_repeat(fit, data, subject, z, row_design)
}

.bc_reconstruction_mse_axis <- function(fit, data, subject, z) {
  axis <- fit$axis
  Z <- .bc_as_design(z, length(data))
  t_hat <- .bc_infer_subject_scores_axis(axis, data, subject, Z, lambda_t = fit$settings$lambda_t)

  L <- axis$L
  R <- axis$R
  B <- axis$B
  beta0 <- axis$beta0
  W <- axis$W
  d <- length(beta0)
  n_r <- ncol(L)
  n_c <- ncol(R)

  subj_char <- as.character(subject)
  sse <- 0
  den <- 0
  for (i in seq_along(data)) {
    t_i <- t_hat[subj_char[i], , drop = TRUE]
    mhat <- beta0 + as.vector(W %*% t_i)
    if (ncol(Z) > 0) mhat <- mhat + as.vector(B %*% Z[i, ])
    Mhat <- matrix(mhat, nrow = n_r, ncol = n_c)
    Xhat <- L %*% Mhat %*% t(R)
    if (identical(fit$connectivity_type, "symmetric")) Xhat <- (Xhat + t(Xhat)) / 2
    diff <- data[[i]] - Xhat
    sse <- sse + sum(diff^2)
    den <- den + length(diff)
  }
  sse / max(den, 1)
}

.bc_reconstruction_mse_repeat <- function(fit, data, subject, z, row_design) {
  repf <- fit$repeat_head
  if (is.null(repf)) return(Inf)
  t_hat <- .bc_infer_subject_scores_repeat(
    repf = repf,
    data = data,
    subject = subject,
    z = z,
    row_design = row_design,
    lambda_t = fit$settings$lambda_t,
    include_interactions = fit$settings$include_seed_interactions
  )

  R <- repf$R
  B <- repf$B
  beta0 <- repf$beta0
  W <- repf$W
  n_lat <- length(beta0)

  subj_char <- as.character(subject)
  sse <- 0
  den <- 0
  for (i in seq_along(data)) {
    t_i <- t_hat[subj_char[i], , drop = TRUE]
    qmat <- .bc_row_design(row_design, nrow(data[[i]]))
    zmat <- .bc_as_design(.bc_subset_rows(z, i), 1)
    Zrep <- if (ncol(zmat) > 0) zmat[rep(1, nrow(data[[i]])), , drop = FALSE] else matrix(numeric(0), nrow(data[[i]]), 0)
    Drows <- .bc_build_repeat_design(Z = Zrep, Q = qmat, include_interactions = fit$settings$include_seed_interactions)

    base <- as.vector(beta0 + W %*% t_i)
    Uhat <- matrix(rep(base, each = nrow(data[[i]])), nrow = nrow(data[[i]]), ncol = n_lat)
    if (ncol(Drows) > 0) Uhat <- Uhat + Drows %*% t(B)
    Xhat <- Uhat %*% t(R)
    if (identical(fit$connectivity_type, "symmetric")) Xhat <- (Xhat + t(Xhat)) / 2

    diff <- data[[i]] - Xhat
    sse <- sse + sum(diff^2)
    den <- den + length(diff)
  }
  sse / max(den, 1)
}

.bc_trait_r2 <- function(fit, data, subject, z, y, row_design) {
  if (is.null(y)) return(NA_real_)
  head <- if (!is.null(fit$axis)) fit$axis else fit$repeat_head
  if (is.null(head$A)) return(NA_real_)

  t_hat <- if (!is.null(fit$axis)) {
    Z <- .bc_as_design(z, length(data))
    .bc_infer_subject_scores_axis(fit$axis, data, subject, Z, lambda_t = fit$settings$lambda_t)
  } else {
    .bc_infer_subject_scores_repeat(
      repf = fit$repeat_head,
      data = data,
      subject = subject,
      z = z,
      row_design = row_design,
      lambda_t = fit$settings$lambda_t,
      include_interactions = fit$settings$include_seed_interactions
    )
  }

  subj_levels <- rownames(t_hat)
  Ysubj <- .bc_prepare_traits(y, subj_levels)
  Yhat <- t_hat %*% head$A

  r2 <- vapply(seq_len(ncol(Ysubj)), function(j) {
    yj <- Ysubj[, j]
    yhat_j <- Yhat[, j]
    sst <- sum((yj - mean(yj))^2)
    if (sst <= 0) return(NA_real_)
    1 - sum((yj - yhat_j)^2) / sst
  }, numeric(1))
  mean(r2, na.rm = TRUE)
}

.bc_infer_subject_scores_axis <- function(axis, data, subject, Z, lambda_t = 1e-2) {
  subj <- factor(subject)
  subj_levels <- levels(subj)
  L <- axis$L
  R <- axis$R
  B <- axis$B
  beta0 <- axis$beta0
  W <- axis$W
  K <- ncol(W)

  Mobs <- t(vapply(data, function(X) as.vector(crossprod(L, X) %*% R), numeric(length(beta0))))
  t_hat <- matrix(0, nrow = length(subj_levels), ncol = K)
  rownames(t_hat) <- subj_levels
  WW <- crossprod(W)

  for (si in seq_along(subj_levels)) {
    idx <- which(subj == subj_levels[si])
    Ri <- length(idx)
    fixed <- matrix(beta0, nrow = Ri, ncol = length(beta0), byrow = TRUE)
    if (ncol(Z) > 0) fixed <- fixed + Z[idx, , drop = FALSE] %*% t(B)
    e_sum <- colSums(Mobs[idx, , drop = FALSE] - fixed)
    lhs <- Ri * WW + (lambda_t + 1e-8) * diag(K)
    rhs <- as.vector(crossprod(W, e_sum))
    t_hat[si, ] <- tryCatch(
      solve(lhs, rhs),
      error = function(e) as.vector(MASS::ginv(lhs) %*% rhs)
    )
  }
  t_hat
}

.bc_infer_subject_scores_repeat <- function(repf,
                                            data,
                                            subject,
                                            z,
                                            row_design,
                                            lambda_t = 1e-2,
                                            include_interactions = TRUE) {
  subj <- factor(subject)
  subj_levels <- levels(subj)
  R <- repf$R
  B <- repf$B
  beta0 <- repf$beta0
  W <- repf$W
  K <- ncol(W)

  row_cov <- .bc_row_design(row_design, nrow(data[[1]]))
  Z <- .bc_as_design(z, length(data))
  U_long <- matrix(0, nrow = length(data) * nrow(data[[1]]), ncol = length(beta0))
  D_long <- matrix(0, nrow = nrow(U_long), ncol = ncol(B))
  subj_long <- character(nrow(U_long))

  idx <- 1L
  for (i in seq_along(data)) {
    XiR <- data[[i]] %*% R
    nr <- nrow(XiR)
    rng <- idx:(idx + nr - 1L)
    U_long[rng, ] <- XiR
    Zrep <- if (ncol(Z) > 0) Z[rep(i, nr), , drop = FALSE] else matrix(numeric(0), nr, 0)
    D_long[rng, ] <- .bc_build_repeat_design(Z = Zrep, Q = row_cov, include_interactions = include_interactions)
    subj_long[rng] <- as.character(subj[i])
    idx <- idx + nr
  }

  subj_long <- factor(subj_long, levels = subj_levels)
  t_hat <- matrix(0, nrow = length(subj_levels), ncol = K)
  rownames(t_hat) <- subj_levels
  WW <- crossprod(W)
  for (si in seq_along(subj_levels)) {
    idx <- which(subj_long == subj_levels[si])
    Ri <- length(idx)
    fixed <- matrix(beta0, nrow = Ri, ncol = length(beta0), byrow = TRUE)
    if (ncol(D_long) > 0) fixed <- fixed + D_long[idx, , drop = FALSE] %*% t(B)
    e_sum <- colSums(U_long[idx, , drop = FALSE] - fixed)
    lhs <- Ri * WW + (lambda_t + 1e-8) * diag(K)
    rhs <- as.vector(crossprod(W, e_sum))
    t_hat[si, ] <- tryCatch(
      solve(lhs, rhs),
      error = function(e) as.vector(MASS::ginv(lhs) %*% rhs)
    )
  }
  t_hat
}

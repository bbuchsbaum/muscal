#' Bilinear Repeated-Measures Mixed Model
#'
#' @description
#' Fits a matrix-aware latent model for repeated connectivity matrices with
#' optional subject-level supervision. The method supports two complementary
#' representations:
#' \itemize{
#'   \item \code{"seed_axis"}: treat row entities as the row axis of each matrix (\eqn{s \times p}).
#'   \item \code{"seed_repeat"}: treat row profiles as repeated observations nested
#'   within subject/repeat.
#'   \item \code{"both"}: fit both heads with a shared ROI basis.
#' }
#'
#' @param data A list of numeric matrices. Each element is one observed
#'   connectivity matrix (typically row-entity \eqn{\times} ROI, or symmetric
#'   ROI \eqn{\times} ROI).
#' @param subject Subject identifier for each matrix in \code{data}.
#' @param z Optional repeat-level design matrix/data frame/vector (same number of
#'   rows as \code{length(data)}).
#' @param y Optional subject-level traits. Can be a named numeric vector, a
#'   matrix/data frame with one row per subject, or a vector ordered by subject
#'   level.
#' @param row_design Optional row-entity design (vector/matrix with one row per
#'   matrix row). Used by \code{mode = "seed_repeat"}.
#' @param seed_design Deprecated alias of \code{row_design}, retained for
#'   backward compatibility.
#' @param mode One of \code{"seed_axis"}, \code{"seed_repeat"}, or \code{"both"}.
#' @param connectivity_type One of \code{"auto"}, \code{"cross"}, or
#'   \code{"symmetric"}.
#' @param symmetric Deprecated logical alias for \code{connectivity_type =
#'   "symmetric"}.
#' @param sym_tol Numeric tolerance for automatic symmetry detection in
#'   \code{connectivity_type = "auto"}.
#' @param complexity One of \code{"manual"}, \code{"fast"},
#'   \code{"balanced"}, or \code{"adaptive"}. If not \code{"manual"},
#'   data-adaptive defaults are inferred for unspecified tuning parameters.
#' @param r_seed Rank for seed mode basis (row basis).
#' @param r_roi Rank for ROI mode basis (column basis).
#' @param k_subject Subject latent dimensionality.
#' @param include_seed_interactions Logical; if \code{TRUE}, includes
#'   \eqn{z \otimes row\_design}{z x row_design} interactions in the seed-repeat head.
#' @param repeat_id Optional repeat/condition index per matrix. Required when
#'   using \code{laplacian_repeat}.
#' @param laplacian_repeat Optional repeat Laplacian matrix.
#' @param laplacian_seed Optional seed-axis Laplacian.
#' @param laplacian_roi Optional ROI-axis Laplacian.
#' @param laplacian_design Optional design Laplacian on columns of \code{z}.
#' @param lambda_w Ridge penalty for \code{W}.
#' @param lambda_t Ridge penalty for subject scores \code{t}.
#' @param lambda_b Ridge penalty for design effects.
#' @param lambda_a Ridge penalty for supervision map \code{A}.
#' @param lambda_y Weight for supervised trait term.
#' @param lambda_repeat Smoothing strength for repeat smoothing.
#' @param lambda_seed Basis smoothing strength for seed basis.
#' @param lambda_roi Basis smoothing strength for ROI basis.
#' @param lambda_design_smooth Laplacian smoothing for design effects.
#' @param max_iter Maximum ALS iterations for the core mixed model.
#' @param tol Relative convergence tolerance for the core mixed model.
#' @param verbose Logical; if \code{TRUE}, prints iteration diagnostics.
#' @param ... Reserved for future extensions.
#'
#' @return A list of class \code{"bilinear_mixed"} containing fitted
#'   head(s), bases, core-model parameters, and backprojected effect maps.
#'
#' @export
bilinear_mixed <- function(data,
                                  subject,
                                  z = NULL,
                                  y = NULL,
                                  row_design = NULL,
                                  seed_design = NULL,
                                  mode = c("seed_axis", "seed_repeat", "both"),
                                  connectivity_type = c("auto", "cross", "symmetric"),
                                  symmetric = FALSE,
                                  sym_tol = 1e-8,
                                  complexity = c("manual", "fast", "balanced", "adaptive"),
                                  r_seed = 2,
                                  r_roi = 8,
                                  k_subject = 2,
                                  include_seed_interactions = TRUE,
                                  repeat_id = NULL,
                                  laplacian_repeat = NULL,
                                  laplacian_seed = NULL,
                                  laplacian_roi = NULL,
                                  laplacian_design = NULL,
                                  lambda_w = 1e-2,
                                  lambda_t = 1e-2,
                                  lambda_b = 1e-3,
                                  lambda_a = 1e-3,
                                  lambda_y = 0,
                                  lambda_repeat = 0,
                                  lambda_seed = 0,
                                  lambda_roi = 0,
                                  lambda_design_smooth = 0,
                                  max_iter = 50,
                                  tol = 1e-6,
                                  verbose = FALSE,
                                  ...) {
  UseMethod("bilinear_mixed")
}

#' @rdname bilinear_mixed
#' @export
bilinear_mixed.default <- function(data, ...) {
  stop(
    "Unsupported data type for bilinear_mixed(). ",
    "Use a list of numeric matrices.",
    call. = FALSE
  )
}

#' @rdname bilinear_mixed
#' @export
bilinear_mixed.list <- function(data,
                                       subject,
                                       z = NULL,
                                       y = NULL,
                                       row_design = NULL,
                                       seed_design = NULL,
                                       mode = c("seed_axis", "seed_repeat", "both"),
                                       connectivity_type = c("auto", "cross", "symmetric"),
                                       symmetric = FALSE,
                                       sym_tol = 1e-8,
                                       complexity = c("manual", "fast", "balanced", "adaptive"),
                                       r_seed = 2,
                                       r_roi = 8,
                                       k_subject = 2,
                                       include_seed_interactions = TRUE,
                                       repeat_id = NULL,
                                       laplacian_repeat = NULL,
                                       laplacian_seed = NULL,
                                       laplacian_roi = NULL,
                                       laplacian_design = NULL,
                                       lambda_w = 1e-2,
                                       lambda_t = 1e-2,
                                       lambda_b = 1e-3,
                                       lambda_a = 1e-3,
                                       lambda_y = 0,
                                       lambda_repeat = 0,
                                       lambda_seed = 0,
                                       lambda_roi = 0,
                                       lambda_design_smooth = 0,
                                       max_iter = 50,
                                       tol = 1e-6,
                                       verbose = FALSE,
                                       ...) {
  mode_missing <- missing(mode)
  connectivity_type_missing <- missing(connectivity_type)
  r_seed_missing <- missing(r_seed)
  r_roi_missing <- missing(r_roi)
  k_subject_missing <- missing(k_subject)
  include_seed_interactions_missing <- missing(include_seed_interactions)
  lambda_w_missing <- missing(lambda_w)
  lambda_t_missing <- missing(lambda_t)
  lambda_b_missing <- missing(lambda_b)
  lambda_a_missing <- missing(lambda_a)
  lambda_y_missing <- missing(lambda_y)
  lambda_repeat_missing <- missing(lambda_repeat)
  lambda_seed_missing <- missing(lambda_seed)
  lambda_roi_missing <- missing(lambda_roi)
  lambda_design_smooth_missing <- missing(lambda_design_smooth)
  complexity <- match.arg(complexity)

  chk::chk_list(data)
  chk::chk_true(length(data) > 1)

  dims <- lapply(data, dim)
  if (any(vapply(dims, is.null, logical(1)))) {
    stop("Each element of `data` must be a matrix/data.frame.", call. = FALSE)
  }
  data <- lapply(data, as.matrix)
  s <- nrow(data[[1]])
  p <- ncol(data[[1]])
  same_dim <- vapply(data, function(x) nrow(x) == s && ncol(x) == p, logical(1))
  if (!all(same_dim)) {
    stop("All matrices in `data` must have identical dimensions.", call. = FALSE)
  }

  n_obs <- length(data)
  if (length(subject) != n_obs) {
    stop("`subject` must have length equal to length(data).", call. = FALSE)
  }
  subject <- factor(subject)
  subj_levels <- levels(subject)
  n_subject <- nlevels(subject)

  if (!is.null(seed_design)) {
    if (!is.null(row_design)) {
      stop("Provide only one of `row_design` or `seed_design`.", call. = FALSE)
    }
    row_design <- seed_design
  }

  if (complexity != "manual") {
    mode_for_rec <- if (mode_missing) "auto" else match.arg(mode, c("seed_axis", "seed_repeat", "both"))
    connectivity_for_rec <- if (isTRUE(symmetric)) {
      "symmetric"
    } else if (connectivity_type_missing) {
      "auto"
    } else {
      match.arg(connectivity_type, c("auto", "cross", "symmetric"))
    }

    rec <- bilinear_mixed_recommend(
      data = data,
      subject = subject,
      z = z,
      y = y,
      row_design = row_design,
      mode = mode_for_rec,
      connectivity_type = connectivity_for_rec,
      profile = complexity,
      sym_tol = sym_tol
    )

    if (mode_missing) mode <- rec$mode
    if (connectivity_type_missing && !isTRUE(symmetric)) connectivity_type <- rec$connectivity_type
    if (r_seed_missing) r_seed <- rec$r_seed
    if (r_roi_missing) r_roi <- rec$r_roi
    if (k_subject_missing) k_subject <- rec$k_subject
    if (include_seed_interactions_missing) include_seed_interactions <- rec$include_seed_interactions
    if (lambda_w_missing) lambda_w <- rec$lambda_w
    if (lambda_t_missing) lambda_t <- rec$lambda_t
    if (lambda_b_missing) lambda_b <- rec$lambda_b
    if (lambda_a_missing) lambda_a <- rec$lambda_a
    if (lambda_y_missing) lambda_y <- rec$lambda_y
    if (lambda_repeat_missing) lambda_repeat <- rec$lambda_repeat
    if (lambda_seed_missing) lambda_seed <- rec$lambda_seed
    if (lambda_roi_missing) lambda_roi <- rec$lambda_roi
    if (lambda_design_smooth_missing) lambda_design_smooth <- rec$lambda_design_smooth
  }

  mode <- match.arg(mode, c("seed_axis", "seed_repeat", "both"))
  connectivity_type <- match.arg(connectivity_type, c("auto", "cross", "symmetric"))

  Z <- .bc_as_design(z, n_obs)
  Pz <- ncol(Z)

  if (!is.null(repeat_id)) {
    if (length(repeat_id) != n_obs) {
      stop("`repeat_id` must have length equal to length(data).", call. = FALSE)
    }
    repeat_id <- as.integer(repeat_id)
  }

  Ysubj <- .bc_prepare_traits(y, subj_levels)
  if (!is.null(Ysubj) && lambda_y > 0 && k_subject < 1) {
    stop("Supervised fitting requires k_subject >= 1.", call. = FALSE)
  }

  is_square_all <- all(vapply(data, function(x) nrow(x) == ncol(x), logical(1)))
  is_symmetric_all <- FALSE
  if (is_square_all) {
    is_symmetric_all <- all(vapply(
      data,
      function(x) max(abs(x - t(x))) <= sym_tol,
      logical(1)
    ))
  }

  if (isTRUE(symmetric)) connectivity_type <- "symmetric"
  if (connectivity_type == "auto") {
    connectivity_type <- if (is_square_all && is_symmetric_all) "symmetric" else "cross"
  }

  if (connectivity_type == "symmetric") {
    if (!is_square_all) {
      stop("`connectivity_type='symmetric'` requires square matrices.", call. = FALSE)
    }
    rr <- min(as.integer(r_seed), as.integer(r_roi), s)
    r_seed_eff <- rr
    r_roi_eff <- rr
  } else {
    r_seed_eff <- min(as.integer(r_seed), s)
    r_roi_eff <- min(as.integer(r_roi), p)
  }
  if (r_seed_eff < 1 || r_roi_eff < 1) {
    stop("`r_seed` and `r_roi` must be >= 1 after dimension capping.", call. = FALSE)
  }
  if (mode == "seed_axis") {
    K <- min(as.integer(k_subject), n_subject, r_seed_eff * r_roi_eff)
  } else {
    K <- min(as.integer(k_subject), n_subject, r_roi_eff)
  }
  K <- max(as.integer(K), 1L)

  R_basis <- .bc_basis_col(data, r_roi_eff, laplacian_roi, lambda_roi)
  if (connectivity_type == "symmetric") {
    L_basis <- R_basis
  } else {
    L_basis <- .bc_basis_row(data, r_seed_eff, laplacian_seed, lambda_seed)
  }

  axis_fit <- NULL
  repeat_fit <- NULL

  if (mode %in% c("seed_axis", "both")) {
    core_axis <- t(vapply(data, function(X) as.vector(crossprod(L_basis, X) %*% R_basis), numeric(r_seed_eff * r_roi_eff)))
    colnames(core_axis) <- paste0("c", seq_len(ncol(core_axis)))

    if (!is.null(laplacian_repeat) && lambda_repeat > 0) {
      if (is.null(repeat_id)) {
        stop("`repeat_id` is required when laplacian_repeat is supplied.", call. = FALSE)
      }
      core_axis <- .bc_smooth_repeats(
        M = core_axis,
        subject = subject,
        repeat_id = repeat_id,
        L_repeat = laplacian_repeat,
        lambda_repeat = lambda_repeat
      )
    }

    fit_axis_core <- .bc_fit_core_model(
      M = core_axis,
      subject = subject,
      Z = Z,
      Y = Ysubj,
      k_subject = K,
      lambda_w = lambda_w,
      lambda_t = lambda_t,
      lambda_b = lambda_b,
      lambda_a = lambda_a,
      lambda_y = lambda_y,
      laplacian_design = laplacian_design,
      lambda_design_smooth = lambda_design_smooth,
      max_iter = max_iter,
      tol = tol,
      verbose = verbose
    )

    design_maps <- vector("list", Pz)
    if (Pz > 0) {
      for (j in seq_len(Pz)) {
        Bj <- matrix(fit_axis_core$B[, j], nrow = r_seed_eff, ncol = r_roi_eff)
        design_maps[[j]] <- L_basis %*% Bj %*% t(R_basis)
      }
      names(design_maps) <- colnames(Z)
    }

    trait_maps <- vector("list", K)
    seed_profiles <- matrix(NA_real_, nrow = s, ncol = K)
    roi_profiles  <- matrix(NA_real_, nrow = p, ncol = K)
    for (k in seq_len(K)) {
      Wk <- matrix(fit_axis_core$W[, k], nrow = r_seed_eff, ncol = r_roi_eff)
      trait_maps[[k]] <- L_basis %*% Wk %*% t(R_basis)
      seed_profiles[, k] <- rowMeans(trait_maps[[k]])
      roi_profiles[, k]  <- colMeans(trait_maps[[k]])
    }
    comp_names <- paste0("comp", seq_len(K))
    names(trait_maps) <- comp_names
    colnames(seed_profiles) <- comp_names
    colnames(roi_profiles)  <- comp_names

    # Propagate row/column names from data if available
    rn <- rownames(data[[1]])
    cn <- colnames(data[[1]])
    if (!is.null(rn)) rownames(seed_profiles) <- rn
    if (!is.null(cn)) rownames(roi_profiles)  <- cn

    # Name core-model outputs
    rownames(fit_axis_core$t_scores) <- subj_levels
    colnames(fit_axis_core$t_scores) <- comp_names
    colnames(fit_axis_core$W) <- comp_names
    if (!is.null(fit_axis_core$A)) {
      rownames(fit_axis_core$A) <- comp_names
    }

    axis_fit <- c(
      list(
        head = "seed_axis",
        L = L_basis,
        R = R_basis,
        design_maps = design_maps,
        trait_maps = trait_maps,
        seed_profiles = seed_profiles,
        roi_profiles = roi_profiles,
        design_names = colnames(Z)
      ),
      fit_axis_core
    )
  }

  if (mode %in% c("seed_repeat", "both")) {
    row_design_mat <- .bc_row_design(row_design, n_row = s)
    Q <- ncol(row_design_mat)

    # Build long-form projected observations: one row per (matrix, seed).
    n_long <- n_obs * s
    U_long <- matrix(0, nrow = n_long, ncol = r_roi_eff)
    subj_long <- character(n_long)
    Z_long <- if (Pz > 0) matrix(0, nrow = n_long, ncol = Pz) else matrix(numeric(0), nrow = n_long, ncol = 0)
    Q_long <- if (Q > 0) matrix(0, nrow = n_long, ncol = Q) else matrix(numeric(0), nrow = n_long, ncol = 0)
    repeat_long <- if (!is.null(repeat_id)) integer(n_long) else NULL

    idx <- 1L
    for (i in seq_len(n_obs)) {
      XiR <- data[[i]] %*% R_basis
      rng <- idx:(idx + s - 1L)
      U_long[rng, ] <- XiR
      subj_long[rng] <- as.character(subject[i])
      if (Pz > 0) Z_long[rng, ] <- Z[rep(i, s), , drop = FALSE]
      if (Q > 0) Q_long[rng, ] <- row_design_mat
      if (!is.null(repeat_long)) repeat_long[rng] <- rep.int(repeat_id[i], s)
      idx <- idx + s
    }

    D_long <- .bc_build_repeat_design(
      Z = Z_long,
      Q = Q_long,
      include_interactions = isTRUE(include_seed_interactions)
    )

    fit_repeat_core <- .bc_fit_core_model(
      M = U_long,
      subject = factor(subj_long, levels = subj_levels),
      Z = D_long,
      Y = Ysubj,
      k_subject = K,
      lambda_w = lambda_w,
      lambda_t = lambda_t,
      lambda_b = lambda_b,
      lambda_a = lambda_a,
      lambda_y = lambda_y,
      laplacian_design = NULL,
      lambda_design_smooth = 0,
      max_iter = max_iter,
      tol = tol,
      verbose = verbose
    )

    D_names <- colnames(D_long)
    roi_design_maps <- vector("list", ncol(D_long))
    for (j in seq_len(ncol(D_long))) {
      roi_design_maps[[j]] <- as.vector(R_basis %*% fit_repeat_core$B[, j, drop = FALSE])
    }
    names(roi_design_maps) <- D_names

    comp_names_rep <- paste0("comp", seq_len(K))
    roi_trait_maps <- vector("list", K)
    roi_profiles_rep <- matrix(NA_real_, nrow = p, ncol = K)
    for (k in seq_len(K)) {
      roi_trait_maps[[k]] <- as.vector(R_basis %*% fit_repeat_core$W[, k, drop = FALSE])
      roi_profiles_rep[, k] <- roi_trait_maps[[k]]
    }
    names(roi_trait_maps) <- comp_names_rep
    colnames(roi_profiles_rep) <- comp_names_rep
    cn_rep <- colnames(data[[1]])
    if (!is.null(cn_rep)) rownames(roi_profiles_rep) <- cn_rep

    # Name core-model outputs
    rownames(fit_repeat_core$t_scores) <- subj_levels
    colnames(fit_repeat_core$t_scores) <- comp_names_rep
    colnames(fit_repeat_core$W) <- comp_names_rep
    if (!is.null(fit_repeat_core$A)) {
      rownames(fit_repeat_core$A) <- comp_names_rep
    }

    repeat_fit <- c(
      list(
        head = "seed_repeat",
        R = R_basis,
        row_design = row_design_mat,
        seed_design = row_design_mat,
        design_matrix = D_long,
        design_names = D_names,
        roi_design_maps = roi_design_maps,
        roi_trait_maps = roi_trait_maps,
        roi_profiles = roi_profiles_rep
      ),
      fit_repeat_core
    )
  }

  out <- list(
    call = match.call(),
    mode = mode,
    dimensions = list(
      n_obs = n_obs,
      n_subject = n_subject,
      n_row = s,
      n_col = p,
      n_seed = s,
      n_roi = p,
      r_seed = r_seed_eff,
      r_roi = r_roi_eff,
      k_subject = K
    ),
    connectivity_type = connectivity_type,
    settings = list(
      complexity = complexity,
      lambda_w = lambda_w,
      lambda_t = lambda_t,
      lambda_b = lambda_b,
      lambda_a = lambda_a,
      lambda_y = lambda_y,
      lambda_repeat = lambda_repeat,
      lambda_seed = lambda_seed,
      lambda_roi = lambda_roi,
      lambda_design_smooth = lambda_design_smooth,
      include_seed_interactions = include_seed_interactions,
      sym_tol = sym_tol
    ),
    subject_levels = subj_levels,
    trait_names = if (!is.null(Ysubj)) colnames(Ysubj) else NULL,
    design_names = colnames(Z),
    axis = axis_fit,
    repeat_head = repeat_fit
  )
  class(out) <- "bilinear_mixed"
  out
}

#' @export
print.bilinear_mixed <- function(x, ...) {
  cat("Bilinear Repeated-Measures Mixed Model\n")
  cat("  mode:", x$mode, "\n")
  cat(
    "  observations:",
    x$dimensions$n_obs,
    "| subjects:",
    x$dimensions$n_subject,
    "| rows:",
    x$dimensions$n_row,
    "| cols:",
    x$dimensions$n_col,
    "\n"
  )
  cat(
    "  ranks (seed, roi):",
    paste0("(", x$dimensions$r_seed, ", ", x$dimensions$r_roi, ")"),
    "| subject latent K:",
    x$dimensions$k_subject,
    "\n"
  )
  cat("  connectivity type:", x$connectivity_type, "\n")
  if (!is.null(x$axis)) {
    cat("  axis head iterations:", length(x$axis$objective_trace), "\n")
  }
  if (!is.null(x$repeat_head)) {
    cat("  repeat head iterations:", length(x$repeat_head$objective_trace), "\n")
  }
  invisible(x)
}

# -----------------------------------------------------------------------------
# Internal helpers
# -----------------------------------------------------------------------------

.bc_as_design <- function(z, n_obs) {
  if (is.null(z)) {
    out <- matrix(numeric(0), nrow = n_obs, ncol = 0)
    return(out)
  }
  if (is.vector(z) && !is.list(z)) {
    z <- matrix(as.numeric(z), ncol = 1)
  }
  if (is.data.frame(z)) z <- as.matrix(z)
  if (!is.matrix(z)) {
    stop("`z` must be NULL, vector, matrix, or data.frame.", call. = FALSE)
  }
  if (nrow(z) != n_obs) {
    stop("`z` must have one row per matrix in `data`.", call. = FALSE)
  }
  storage.mode(z) <- "double"
  if (is.null(colnames(z))) colnames(z) <- paste0("z", seq_len(ncol(z)))
  z
}

.bc_row_design <- function(row_design, n_row) {
  if (is.null(row_design)) {
    return(matrix(numeric(0), nrow = n_row, ncol = 0))
  }
  if (is.vector(row_design) && !is.list(row_design)) {
    row_design <- matrix(as.numeric(row_design), ncol = 1)
  }
  if (is.data.frame(row_design)) row_design <- as.matrix(row_design)
  if (!is.matrix(row_design)) {
    stop("`row_design` must be NULL, vector, matrix, or data.frame.", call. = FALSE)
  }
  if (nrow(row_design) != n_row) {
    stop("`row_design` must have one row per matrix row.", call. = FALSE)
  }
  storage.mode(row_design) <- "double"
  if (is.null(colnames(row_design))) {
    colnames(row_design) <- paste0("row_cov", seq_len(ncol(row_design)))
  }
  row_design
}

.bc_prepare_traits <- function(y, subj_levels) {
  if (is.null(y)) return(NULL)
  n_subj <- length(subj_levels)

  if (is.vector(y) && !is.list(y)) {
    if (!is.null(names(y))) {
      y <- y[subj_levels]
    }
    if (length(y) != n_subj) {
      stop("Trait vector must have one value per subject.", call. = FALSE)
    }
    Y <- matrix(as.numeric(y), ncol = 1)
    rownames(Y) <- subj_levels
    colnames(Y) <- "trait1"
    return(Y)
  }

  if (is.data.frame(y)) y <- as.matrix(y)
  if (!is.matrix(y)) {
    stop("`y` must be NULL, vector, matrix, or data.frame.", call. = FALSE)
  }

  if (nrow(y) != n_subj) {
    stop("Trait matrix must have one row per unique subject.", call. = FALSE)
  }
  if (!is.null(rownames(y))) {
    y <- y[subj_levels, , drop = FALSE]
  } else {
    rownames(y) <- subj_levels
  }
  storage.mode(y) <- "double"
  if (is.null(colnames(y))) colnames(y) <- paste0("trait", seq_len(ncol(y)))
  y
}

.bc_basis_col <- function(data, r, L_roi = NULL, lambda_roi = 0) {
  p <- ncol(data[[1]])
  C <- matrix(0, p, p)
  for (X in data) C <- C + crossprod(X)
  C <- (C + t(C)) / (2 * length(data))
  if (!is.null(L_roi) && lambda_roi > 0) {
    if (!is.matrix(L_roi) || any(dim(L_roi) != c(p, p))) {
      stop("`laplacian_roi` must be a p x p matrix.", call. = FALSE)
    }
    C <- C - lambda_roi * ((L_roi + t(L_roi)) / 2)
  }
  ev <- eigen(C, symmetric = TRUE)
  V <- ev$vectors[, seq_len(min(r, ncol(ev$vectors))), drop = FALSE]
  qr.Q(qr(V))
}

.bc_basis_row <- function(data, r, L_seed = NULL, lambda_seed = 0) {
  s <- nrow(data[[1]])
  C <- matrix(0, s, s)
  for (X in data) C <- C + tcrossprod(X)
  C <- (C + t(C)) / (2 * length(data))
  if (!is.null(L_seed) && lambda_seed > 0) {
    if (!is.matrix(L_seed) || any(dim(L_seed) != c(s, s))) {
      stop("`laplacian_seed` must be an s x s matrix.", call. = FALSE)
    }
    C <- C - lambda_seed * ((L_seed + t(L_seed)) / 2)
  }
  ev <- eigen(C, symmetric = TRUE)
  V <- ev$vectors[, seq_len(min(r, ncol(ev$vectors))), drop = FALSE]
  qr.Q(qr(V))
}

.bc_build_repeat_design <- function(Z, Q, include_interactions = TRUE) {
  Pz <- ncol(Z)
  Qq <- ncol(Q)
  if (Pz > 0 && is.null(colnames(Z))) colnames(Z) <- paste0("z", seq_len(Pz))
  if (Qq > 0 && is.null(colnames(Q))) colnames(Q) <- paste0("row_cov", seq_len(Qq))
  blocks <- list()
  names_out <- character(0)
  if (Pz > 0) {
    blocks[[length(blocks) + 1L]] <- Z
    names_out <- c(names_out, colnames(Z))
  }
  if (Qq > 0) {
    blocks[[length(blocks) + 1L]] <- Q
    names_out <- c(names_out, colnames(Q))
  }
  if (isTRUE(include_interactions) && Pz > 0 && Qq > 0) {
    inter <- matrix(0, nrow = nrow(Z), ncol = Pz * Qq)
    nm <- character(Pz * Qq)
    idx <- 1L
    for (i in seq_len(Pz)) {
      for (j in seq_len(Qq)) {
        inter[, idx] <- Z[, i] * Q[, j]
        nm[idx] <- paste0(colnames(Z)[i], ":", colnames(Q)[j])
        idx <- idx + 1L
      }
    }
    blocks[[length(blocks) + 1L]] <- inter
    names_out <- c(names_out, nm)
  }
  if (length(blocks) == 0) {
    out <- matrix(numeric(0), nrow = nrow(Z), ncol = 0)
    return(out)
  }
  out <- do.call(cbind, blocks)
  colnames(out) <- names_out
  out
}

.bc_smooth_repeats <- function(M, subject, repeat_id, L_repeat, lambda_repeat) {
  if (!is.matrix(L_repeat)) stop("`laplacian_repeat` must be a matrix.", call. = FALSE)
  if (lambda_repeat <= 0) return(M)

  out <- M
  subj <- factor(subject)
  d <- ncol(M)
  for (sid in levels(subj)) {
    idx <- which(subj == sid)
    rid <- repeat_id[idx]
    if (anyDuplicated(rid) > 0) {
      next
    }
    if (length(rid) < 2) next
    if (max(rid) > nrow(L_repeat) || max(rid) > ncol(L_repeat)) {
      stop("repeat_id contains indices outside laplacian_repeat dimensions.", call. = FALSE)
    }
    Ls <- L_repeat[rid, rid, drop = FALSE]
    A <- diag(length(idx)) + lambda_repeat * Ls
    A <- (A + t(A)) / 2
    for (j in seq_len(d)) {
      out[idx, j] <- tryCatch(
        solve(A, M[idx, j]),
        error = function(e) as.vector(MASS::ginv(A) %*% M[idx, j])
      )
    }
  }
  out
}

.bc_fit_core_model <- function(M,
                               subject,
                               Z,
                               Y,
                               k_subject,
                               lambda_w,
                               lambda_t,
                               lambda_b,
                               lambda_a,
                               lambda_y,
                               laplacian_design = NULL,
                               lambda_design_smooth = 0,
                               max_iter = 50,
                               tol = 1e-6,
                               verbose = FALSE) {
  n <- nrow(M)
  d <- ncol(M)
  subj <- factor(subject)
  subj_idx <- as.integer(subj)
  S <- nlevels(subj)
  counts <- tabulate(subj_idx, nbins = S)
  Pz <- ncol(Z)
  K <- min(max(1L, as.integer(k_subject)), S, d)

  # Subject-level initialization from subject means in core space.
  Msub <- rowsum(M, subj_idx, reorder = FALSE)
  Msub <- Msub / counts
  Msub <- scale(Msub, center = TRUE, scale = FALSE)
  sv <- svd(Msub, nu = min(K, nrow(Msub)), nv = 0)
  tmat <- matrix(0, nrow = S, ncol = K)
  kk <- min(K, length(sv$d))
  if (kk > 0) {
    tmat[, seq_len(kk)] <- sv$u[, seq_len(kk), drop = FALSE] %*% diag(sv$d[seq_len(kk)], kk, kk)
  }

  W <- matrix(0, nrow = d, ncol = K)
  B <- matrix(0, nrow = d, ncol = Pz)
  beta0 <- colMeans(M)

  supervised <- !is.null(Y) && lambda_y > 0
  if (supervised && nrow(Y) != S) {
    stop("Trait matrix Y must have one row per subject level.", call. = FALSE)
  }
  A <- if (supervised) matrix(0, nrow = K, ncol = ncol(Y)) else NULL

  obj <- numeric(max_iter)
  prev <- Inf

  for (iter in seq_len(max_iter)) {
    Tobs <- tmat[subj_idx, , drop = FALSE]

    # Update W from current t and fixed effects.
    fixed_part <- matrix(beta0, nrow = n, ncol = d, byrow = TRUE)
    if (Pz > 0) fixed_part <- fixed_part + Z %*% t(B)
    W <- t(ls_ridge(Tobs, M - fixed_part, lambda = lambda_w))

    # Update fixed effects beta0 and B.
    resid_fixed <- M - Tobs %*% t(W)
    beta0 <- colMeans(resid_fixed)
    if (Pz > 0) {
      resid_centered <- resid_fixed - matrix(beta0, nrow = n, ncol = d, byrow = TRUE)
      G <- crossprod(Z) + (lambda_b + 1e-8) * diag(Pz)
      if (!is.null(laplacian_design) && lambda_design_smooth > 0) {
        if (!is.matrix(laplacian_design) || any(dim(laplacian_design) != c(Pz, Pz))) {
          stop("`laplacian_design` must be a P x P matrix matching design columns.", call. = FALSE)
        }
        G <- G + lambda_design_smooth * ((laplacian_design + t(laplacian_design)) / 2)
      }
      rhs <- crossprod(Z, resid_centered)
      Bt <- tryCatch(
        solve(G, rhs),
        error = function(e) MASS::ginv(G) %*% rhs
      )
      B <- t(Bt)
      beta0 <- colMeans(resid_fixed - Z %*% t(B))
    }

    if (supervised) {
      A <- ls_ridge(tmat, Y, lambda = lambda_a)
    }

    # Update subject scores t_i.
    WW <- crossprod(W)
    AA <- if (supervised) A %*% t(A) else matrix(0, K, K)
    fixed_part <- matrix(beta0, nrow = n, ncol = d, byrow = TRUE)
    if (Pz > 0) fixed_part <- fixed_part + Z %*% t(B)
    for (si in seq_len(S)) {
      idx <- which(subj_idx == si)
      Ri <- length(idx)
      e_sum <- colSums(M[idx, , drop = FALSE] - fixed_part[idx, , drop = FALSE])
      lhs <- Ri * WW + (lambda_t + 1e-8) * diag(K)
      rhs <- as.vector(crossprod(W, e_sum))
      if (supervised) {
        lhs <- lhs + lambda_y * AA
        rhs <- rhs + lambda_y * as.vector(A %*% as.numeric(Y[si, ]))
      }
      tmat[si, ] <- tryCatch(
        solve(lhs, rhs),
        error = function(e) as.vector(MASS::ginv(lhs) %*% rhs)
      )
    }

    # Objective and convergence.
    Tobs <- tmat[subj_idx, , drop = FALSE]
    Mu <- matrix(beta0, nrow = n, ncol = d, byrow = TRUE)
    if (Pz > 0) Mu <- Mu + Z %*% t(B)
    Mu <- Mu + Tobs %*% t(W)
    recon <- sum((M - Mu)^2)
    penalty <- lambda_w * sum(W^2) + lambda_t * sum(tmat^2) + lambda_b * sum(B^2)
    if (!is.null(laplacian_design) && Pz > 0 && lambda_design_smooth > 0) {
      penalty <- penalty + lambda_design_smooth * sum(diag(B %*% laplacian_design %*% t(B)))
    }
    if (supervised) {
      penalty <- penalty +
        lambda_y * sum((Y - tmat %*% A)^2) +
        lambda_a * sum(A^2)
    }
    obj[iter] <- recon + penalty

    rel <- abs(obj[iter] - prev) / (abs(prev) + 1e-8)
    if (isTRUE(verbose)) {
      message(
        sprintf(
          "bilinear core iter %d: obj=%.6g rel=%.3g",
          iter, obj[iter], rel
        )
      )
    }
    if (is.finite(prev) && rel < tol) {
      obj <- obj[seq_len(iter)]
      break
    }
    prev <- obj[iter]
  }

  list(
    beta0 = beta0,
    B = B,
    W = W,
    t_scores = tmat,
    A = A,
    objective_trace = obj
  )
}

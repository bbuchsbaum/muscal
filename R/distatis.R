#' DISTATIS for Multiple Distance Matrices
#'
#' Fits a consensus object map, an assessor map based on RV similarities, and
#' partial object maps for distances measured on the same objects.
#'
#' @param data A list of numeric square distance matrices or [stats::dist()]
#'   objects, an object-by-object-by-assessor numeric array, or a single distance
#'   matrix or `dist` object. Object names and order must agree across tables.
#' @param ncomp Positive integer; maximum number of axes retained in each map.
#'   Object and assessor ranks are determined independently.
#' @param squared Logical; are the supplied distances already squared? Defaults
#'   to `FALSE`, including for `dist` objects.
#' @param norm_method Table normalization: `"mfa"` (default) divides each centered
#'   Gram matrix by its largest eigenvalue; `"frobenius"` gives unit Frobenius
#'   norm; `"none"` preserves scale. RV similarities are used in all three cases.
#' @param labels Optional unique object labels. These override display labels,
#'   but do not bypass checks of input object names and order.
#' @param tol Relative numerical tolerance, strictly between zero and one, for
#'   input symmetry, Euclidean geometry, spectral rank, and leading-root ties.
#' @param ... Reserved for methods. Unused arguments are rejected.
#'
#' @details
#' For ordinary distance matrices \eqn{D_k}, DISTATIS computes
#' \eqn{B_k=-HD_k^2H/2}, with elementwise squaring and
#' \eqn{H=I-11^T/n}. The normalized Gram matrices \eqn{S_k} are compared using
#' RV coefficients. The leading eigenvector of the RV matrix, oriented
#' nonnegatively and rescaled to sum to one, supplies the compromise weights.
#' The weighted compromise has object coordinates \eqn{V\Lambda^{1/2}};
#' each partial map is \eqn{S_kV\Lambda^{-1/2}}. Their weighted mean equals the
#' object map. Assessor coordinates use the RV matrix's own spectrum.
#'
#' Non-Euclidean distances, zero Gram matrices, missing values, and a tied leading
#' RV eigenvalue are rejected. Small negative eigenvalues within `tol` are treated
#' as numerical zero for rank and inertia calculations. No distance correction
#' or imputation is performed. Weights describe agreement with the fitted panel,
#' rather than assessor accuracy. Inference is not currently implemented.
#'
#' @return A `distatis` object inheriting from `covstatis` and `projector`.
#'   `object_scores` (also `roi_scores`) and `subject_scores` contain the two maps;
#'   `partial_scores` contains one object map per assessor; `alpha` contains
#'   weights; `C` is the RV matrix; `compromise_full` is the normalized compromise
#'   Gram matrix. `eigenvalues` and `interstructure_eigenvalues` retain the full
#'   spectra, while `inertia` and `subject_inertia` use their complete positive
#'   inertia as denominators. `normalization_constants`, `input_labels`,
#'   `squared`, `norm_method`, `tol`, and `diagnostics` record preprocessing.
#' @references
#' Abdi, H., Valentin, D., Chollet, S., & Chrea, C. (2007). Analyzing assessors
#' and products in sorting tasks: DISTATIS, theory and applications.
#' Food Quality and Preference, 18, 627-640. \doi{10.1016/j.foodqual.2006.09.003}.
#' @seealso [covstatis()], [sorting_distances()], [project_distance()],
#'   [project_subjects()], [reconstruct.distatis()], [plot.distatis()]
#' @examples
#' partitions <- data.frame(A = c(1, 1, 2, 3), B = c(1, 1, 2, 2),
#'                          C = c(1, 2, 3, 3))
#' distances <- sorting_distances(partitions)
#' fit <- distatis(distances)
#' fit$alpha
#' multivarious::scores(fit)
#' project_distance(fit, distances[[1]])
#' @md
#' @export
distatis <- function(data, ncomp = 2, squared = FALSE,
                     norm_method = c("mfa", "frobenius", "none"),
                     labels = NULL, tol = 1e-8, ...) {
  UseMethod("distatis")
}

#' @rdname distatis
#' @md
#' @export
distatis.default <- function(data, ncomp = 2, squared = FALSE,
                             norm_method = c("mfa", "frobenius", "none"),
                             labels = NULL, tol = 1e-8, ...) {
  .distatis_unused(...)
  .distatis_count(ncomp, "ncomp")
  .distatis_flag(squared, "squared")
  if (!is.numeric(tol) || length(tol) != 1L || !is.finite(tol) || tol <= 0 || tol >= 1) {
    stop("`tol` must be a finite number strictly between zero and one.", call. = FALSE)
  }
  norm_method <- match.arg(norm_method)
  tables <- .distatis_as_tables(data)
  processed <- .distatis_preprocess(tables, squared, norm_method, tol)
  n <- nrow(processed$tables[[1L]])
  object_labels <- labels %||% processed$input_labels %||% paste0("Obs_", seq_len(n))
  object_labels <- .distatis_labels(object_labels, n, "labels")
  fit <- .fit_statis_grams(
    processed$tables, ncomp, norm_method, TRUE, names(tables), object_labels,
    interstructure = "rv", tol = tol, classes = c("distatis", "covstatis")
  )
  axes <- which(fit$eigC$values > fit$eigC$values[1L] * tol)
  axes <- head(axes, ncomp)
  fit$subject_axes <- axes
  fit$subject_scores <- sweep(fit$eigC$vectors[, axes, drop = FALSE], 2,
                              sqrt(fit$eigC$values[axes]), `*`)
  fit$interstructure_eigenvalues <- fit$eigC$values
  fit$inertia <- pmax(fit$eigenvalues, 0) / sum(pmax(fit$eigenvalues, 0))
  fit$subject_inertia <- pmax(fit$eigC$values, 0) / sum(pmax(fit$eigC$values, 0))
  fit$compromise_quality <- fit$subject_inertia[1L]
  fit$normalization_constants <- processed$constants
  fit$input_labels <- processed$input_labels
  fit$squared <- squared
  fit$tol <- tol
  fit$diagnostics <- processed$diagnostics
  fit$rv_table_matrix <- .statis_unit_table_matrix(fit$table_matrix)
  dimnames(fit$v) <- dimnames(fit$roi_scores) <- dimnames(fit$projmat) <-
    list(object_labels, paste0("Dim", seq_len(ncol(fit$v))))
  fit$object_scores <- fit$roi_scores
  dimnames(fit$subject_scores) <- list(fit$block_labels, paste0("Dim", axes))
  dimnames(fit$C) <- list(fit$block_labels, fit$block_labels)
  dimnames(fit$compromise_full) <- list(object_labels, object_labels)
  fit$partial_scores <- lapply(fit$partial_scores, function(S) {
    dimnames(S) <- dimnames(fit$object_scores)
    S
  })
  fit
}

.distatis_as_tables <- function(data) {
  if (inherits(data, "dist") || is.matrix(data)) {
    data <- list(data)
  } else if (is.array(data) && length(dim(data)) == 3L) {
    dn <- dimnames(data)
    tables <- lapply(seq_len(dim(data)[3L]), function(k) data[, , k, drop = FALSE][, , 1L])
    names(tables) <- dn[[3L]]
    data <- tables
  }
  if (!is.list(data) || is.data.frame(data) || !length(data)) {
    stop("`data` must be a distance matrix, dist object, nonempty list, or three-dimensional array.",
         call. = FALSE)
  }
  ids <- names(data) %||% paste0("Assessor_", seq_along(data))
  names(data) <- .distatis_labels(ids, length(data), "assessor names")
  data
}

.distatis_preprocess <- function(tables, squared, norm_method, tol,
                                 expected_n = NULL, expected_labels = NULL) {
  input_labels <- expected_labels
  n <- expected_n
  constants <- numeric(length(tables))
  diagnostics <- vector("list", length(tables))
  result <- lapply(seq_along(tables), function(k) {
    id <- names(tables)[k]
    fail <- function(msg) stop(sprintf("Assessor '%s': %s", id, msg), call. = FALSE)
    D <- tables[[k]]
    if (inherits(D, "dist")) D <- as.matrix(D)
    if (!is.matrix(D) || !is.numeric(D) || nrow(D) < 2L || nrow(D) != ncol(D) ||
        any(!is.finite(D))) {
      fail("distances must be a finite numeric square matrix with at least two objects.")
    }
    if (is.null(n)) n <<- nrow(D)
    if (nrow(D) != n) fail("object dimensions must agree with the other tables or fitted model.")
    rn <- rownames(D)
    cn <- colnames(D)
    if (!is.null(rn)) rn <- .distatis_labels(rn, n, "object row names")
    if (!is.null(cn)) cn <- .distatis_labels(cn, n, "object column names")
    if (!is.null(rn) && !is.null(cn) && !identical(rn, cn)) {
      fail("row and column object names/order must agree.")
    }
    ids <- rn %||% cn
    if (!is.null(ids)) {
      if (is.null(input_labels)) input_labels <<- ids
      if (!identical(ids, input_labels)) fail("object names/order must agree across tables and with the fitted model.")
    }
    scale <- max(abs(D))
    if (scale == 0) fail("zero distances give a zero Gram matrix and undefined RV similarities.")
    if (max(abs(D - t(D))) > tol * scale) fail("distance matrix must be symmetric.")
    if (max(abs(diag(D))) > tol * scale) fail("distance matrix must have a zero diagonal.")
    if (min(D) < -tol * scale) fail("distances must be nonnegative.")
    D <- (D + t(D)) / 2
    D[D < 0] <- 0
    diag(D) <- 0
    B <- -0.5 * double_center(if (squared) D else D^2)
    if (any(!is.finite(B))) fail("distance conversion overflowed; rescale the distances.")
    B <- (B + t(B)) / 2
    roots <- eigen(B, symmetric = TRUE, only.values = TRUE)$values
    spectral_scale <- max(abs(roots))
    if (spectral_scale == 0 || roots[1L] <= 0) fail("Gram matrix has no positive inertia.")
    negative_inertia <- sum(abs(roots[roots < 0])) / sum(abs(roots))
    if (min(roots) < -tol * spectral_scale) {
      fail(sprintf("non-Euclidean distances (minimum Gram eigenvalue %.6g; negative inertia %.6g).",
                   min(roots), negative_inertia))
    }
    constant <- switch(norm_method, mfa = roots[1L],
                       frobenius = .statis_frobenius(B), none = 1)
    if (!is.finite(constant) || constant <= 0) fail("normalization failed; rescale the distances.")
    constants[k] <<- constant
    diagnostics[[k]] <<- data.frame(assessor = id, minimum_eigenvalue = min(roots),
                                    negative_inertia = negative_inertia,
                                    rank = sum(roots > roots[1L] * tol))
    B / constant
  })
  names(result) <- names(constants) <- names(tables)
  list(tables = result, constants = constants, input_labels = input_labels,
       diagnostics = do.call(rbind, diagnostics))
}

.distatis_labels <- function(x, n, name) {
  if (!is.character(x) || length(x) != n || anyNA(x) || any(!nzchar(x)) || anyDuplicated(x)) {
    stop(sprintf("`%s` must contain %d unique, nonempty character labels.", name, n), call. = FALSE)
  }
  x
}

.distatis_count <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x < 1 || x != floor(x)) {
    stop(sprintf("`%s` must be a positive integer.", name), call. = FALSE)
  }
}

.distatis_flag <- function(x, name) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    stop(sprintf("`%s` must be TRUE or FALSE.", name), call. = FALSE)
  }
}

.distatis_unused <- function(...) {
  if (length(list(...))) stop("Unused arguments in `...`.", call. = FALSE)
}

#' Convert Sorting Partitions to Distance Matrices
#'
#' @param data An object-by-assessor matrix or data frame of categorical group
#'   identifiers. Identifiers may be numeric, character, logical, or factors.
#'   Every object must have exactly one nonmissing group identifier per assessor.
#' @return A named list of binary distance matrices. A pair has distance zero
#'   when the assessor grouped it together and one otherwise. Row/column names
#'   identify objects. Renaming group identifiers leaves the distances unchanged.
#' @seealso [distatis()]
#' @md
#' @export
sorting_distances <- function(data) {
  if ((!is.matrix(data) && !is.data.frame(data)) || nrow(data) < 2L || ncol(data) < 1L) {
    stop("`data` must be an object-by-assessor matrix or data frame with at least two objects.", call. = FALSE)
  }
  ids <- rownames(data) %||% paste0("Obs_", seq_len(nrow(data)))
  ids <- .distatis_labels(ids, nrow(data), "object names")
  assessors <- colnames(data) %||% paste0("Assessor_", seq_len(ncol(data)))
  assessors <- .distatis_labels(assessors, ncol(data), "assessor names")
  tables <- lapply(seq_len(ncol(data)), function(k) {
    groups <- if (is.data.frame(data)) data[[k]] else data[, k]
    if (!is.atomic(groups) || (!is.numeric(groups) && !is.character(groups) &&
        !is.logical(groups) && !is.factor(groups)) || anyNA(groups) ||
        (is.numeric(groups) && any(!is.finite(groups)))) {
      stop(sprintf("Assessor '%s': group identifiers must be nonmissing categorical values.", assessors[k]),
           call. = FALSE)
    }
    codes <- match(groups, unique(groups))
    D <- 1 * outer(codes, codes, `!=`)
    dimnames(D) <- list(ids, ids)
    D
  })
  names(tables) <- assessors
  tables
}

#' Scores from a DISTATIS Model
#'
#' @param x A fitted `distatis` model.
#' @param ... Unused arguments.
#' @return The object-by-component compromise scores.
#' @importFrom multivarious scores
#' @md
#' @export
scores.distatis <- function(x, ...) x$object_scores

#' @rdname scores.distatis
#' @return `sdev()` returns square roots of the retained object eigenvalues.
#' @importFrom multivarious sdev
#' @md
#' @export
sdev.distatis <- function(x) x$sdev

#' Project a Supplementary Assessor's Distances
#'
#' @param x A fitted `distatis` model.
#' @param new_data A distance matrix or `dist` object on the fitted objects in
#'   their original order. This projects a new assessor, not new objects.
#' @param squared Whether new distances are already squared. Defaults to the
#'   setting used when fitting `x`.
#' @param ... Unused arguments, rejected.
#' @return An object-by-component partial score matrix. The fitted compromise
#'   and assessor weights remain fixed.
#' @seealso [project_subjects()], [distatis()]
#' @md
#' @export
project_distance <- function(x, new_data, squared = x$squared, ...) {
  .distatis_unused(...)
  if (!inherits(x, "distatis")) stop("`x` must be a fitted distatis model.", call. = FALSE)
  .distatis_flag(squared, "squared")
  S <- .distatis_preprocess(list(Supplementary = new_data), squared, x$norm_method,
                            x$tol, nrow(x$v), x$input_labels)$tables[[1L]]
  out <- S %*% x$projmat
  dimnames(out) <- dimnames(x$object_scores)
  out
}

#' @rdname project_distance
#' @md
#' @export
project.distatis <- function(x, new_data, ...) project_distance(x, new_data, ...)

#' @rdname project_distance
#' @md
#' @export
project_cov.distatis <- function(x, new_data, ...) {
  stop("Use project_distance() for a DISTATIS model; project_cov() expects covariance data.", call. = FALSE)
}

#' Supplementary Assessor Maps and Diagnostics for DISTATIS
#'
#' @inheritParams project_distance
#' @param subject_ids Optional unique names for supplementary assessors.
#' @return A list with `subject_scores` in the fitted RV space, `subject_cosines`
#'   (cosines of those retained coordinates), `object_scores` (also `roi_scores`)
#'   containing partial object maps, and `scalar_summaries`. The latter reports
#'   RV against the full compromise, Frobenius `distance_to_compromise`, and
#'   `distance_to_subspace` for the retained object subspace. `new_data` accepts
#'   the same representations as [distatis()]. A training table reproduces its
#'   fitted assessor and partial object coordinates.
#' @md
#' @export
project_subjects.distatis <- function(x, new_data, subject_ids = NULL,
                                      squared = x$squared, ...) {
  .distatis_unused(...)
  .distatis_flag(squared, "squared")
  tables <- .distatis_as_tables(new_data)
  ids <- subject_ids %||% names(tables)
  ids <- .distatis_labels(ids, length(tables), "subject_ids")
  S <- .distatis_preprocess(tables, squared, x$norm_method, x$tol,
                            nrow(x$v), x$input_labels)$tables
  M <- vapply(S, c, numeric(nrow(x$v)^2))
  rv_profiles <- crossprod(.statis_unit_table_matrix(M), x$rv_table_matrix)
  axes <- x$subject_axes
  G <- sweep(rv_profiles %*% x$eigC$vectors[, axes, drop = FALSE], 2,
              sqrt(x$eigC$values[axes]), `/`)
  dimnames(G) <- list(ids, colnames(x$subject_scores))
  norms <- sqrt(rowSums(G^2))
  norms[norms == 0] <- 1
  cosines <- G / norms
  partial <- lapply(S, function(B) {
    F <- B %*% x$projmat
    dimnames(F) <- dimnames(x$object_scores)
    F
  })
  names(partial) <- ids
  P <- tcrossprod(x$v)
  unit_compromise <- .statis_unit_table_matrix(matrix(c(x$compromise_full), ncol = 1L))
  summaries <- data.frame(
    subject_id = ids,
    rv_coefficient = drop(crossprod(.statis_unit_table_matrix(M), unit_compromise)),
    distance_to_compromise = vapply(S, function(B) .statis_frobenius(B - x$compromise_full), numeric(1)),
    distance_to_subspace = vapply(S, function(B) .statis_frobenius(B - P %*% B %*% P), numeric(1)),
    row.names = ids
  )
  list(subject_scores = G, subject_cosines = cosines, scalar_summaries = summaries,
       object_scores = partial, roi_scores = partial)
}

#' Reconstruct the Normalized DISTATIS Compromise
#'
#' @param x A fitted `distatis` model.
#' @param comp Nonempty vector of distinct retained object component indices.
#' @param what Return the centered `"gram"` matrix (default) or `"distance"`.
#' @param squared For `what = "distance"`, return squared distances?
#' @param ... Unused arguments, rejected.
#' @return A square matrix in normalized compromise units, not the original
#'   assessor-specific distance scales. Distance reconstruction uses the Gram
#'   diagonal identity and optionally its elementwise square root.
#' @importFrom multivarious reconstruct
#' @md
#' @export
reconstruct.distatis <- function(x, comp = seq_len(ncol(x$v)),
                                 what = c("gram", "distance"), squared = FALSE, ...) {
  .distatis_unused(...)
  what <- match.arg(what)
  .distatis_flag(squared, "squared")
  if (!is.numeric(comp) || !length(comp) || any(!is.finite(comp)) ||
      any(comp != floor(comp) | comp < 1 | comp > ncol(x$v)) || anyDuplicated(comp)) {
    stop("`comp` must contain distinct retained object component indices.", call. = FALSE)
  }
  B <- reconstruct.covstatis(x, comp = comp)
  if (what == "gram") return(B)
  D <- outer(diag(B), diag(B), `+`) - 2 * B
  D[D < 0] <- 0
  diag(D) <- 0
  if (squared) D else sqrt(D)
}

#' Project Covariates onto DISTATIS Maps
#'
#' @param x A fitted `distatis` model.
#' @param y Finite numeric covariate, one value per assessor or object.
#' @param side `"table"` for assessors or `"feature"` for objects.
#' @param what `"dimension"` returns uncentered columnwise cosines or regression
#'   slopes in the selected map. For table covariates, `"feature"` (alias
#'   `"observation"`) returns an object-by-component pattern by projecting the
#'   partial scores across assessors onto `y`.
#' @param scale `"cosine"` or `"beta"`. For dimension projections the slope
#'   regresses `y` on each map axis. For feature patterns it regresses each
#'   partial-score coordinate on `y`. No intercept or automatic centering is used.
#' @param ... Unused arguments, rejected.
#' @return A named vector for dimension projections, otherwise an object score
#'   matrix. Assessor and object projections may have different axis counts.
#' @md
#' @export
project_covariate.distatis <- function(x, y, side = c("table", "feature"),
                                      what = "dimension", scale = c("cosine", "beta"), ...) {
  .distatis_unused(...)
  side <- match.arg(side)
  scale <- match.arg(scale)
  what <- match.arg(what, c("dimension", "feature", "observation"))
  S <- if (side == "table") x$subject_scores else x$object_scores
  if (!is.numeric(y) || !is.null(dim(y)) || length(y) != nrow(S) || any(!is.finite(y))) {
    stop("`y` must be a finite numeric vector with one value per selected entity.", call. = FALSE)
  }
  if (what == "dimension") {
    denom <- if (scale == "cosine") sqrt(colSums(S^2) * sum(y^2)) else colSums(S^2)
    denom[denom == 0] <- Inf
    return(setNames(drop(crossprod(y, S)) / denom, colnames(S)))
  }
  if (side != "table") stop("Feature-side covariates only support what = 'dimension'.", call. = FALSE)
  M <- vapply(x$partial_scores, c, numeric(length(x$object_scores)))
  denom <- if (scale == "cosine") sqrt(rowSums(M^2) * sum(y^2)) else rep(sum(y^2), nrow(M))
  denom[denom == 0] <- Inf
  matrix(drop(M %*% y) / denom, nrow(x$object_scores), dimnames = dimnames(x$object_scores))
}

#' @md
#' @export
as_multifer_adapter.distatis <- function(x, ...) {
  stop("DISTATIS inference is not implemented; its RV refits require a dedicated adapter.", call. = FALSE)
}

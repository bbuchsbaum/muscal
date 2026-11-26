#' Repeated-measures COVSTATIS
#'
#' Estimate the common covariance geometry and the covariance patterns that
#' separate conditions in a balanced, one-factor repeated-measures design.
#'
#' @param data List of covariance matrices, one per subject-condition cell.
#'   All matrices must be finite, symmetric, and describe the same variables
#'   in the same order. Supply consistent variable names throughout or use
#'   positional matching throughout. Positive semidefiniteness is assumed.
#' @param subject,condition Vectors identifying the subject and condition of
#'   each matrix, in input order. Each subject must have one matrix per
#'   condition. Fitting requires at least two subjects and two conditions.
#'   Observed factor levels specify condition order. Other labels use order
#'   of first appearance.
#' @param ncomp Maximum number of common variable axes to retain.
#' @param ncomp_condition Maximum number of condition axes to retain. The
#'   default keeps all numerically supported axes, up to one fewer than the
#'   number of conditions.
#' @param norm_method Per-cell scaling. The default, `"none"`, preserves
#'   covariance magnitude. `"frobenius"` divides each matrix by its Frobenius
#'   norm and `"mfa"` divides by its largest eigenvalue.
#' @param dcenter Apply Gower double centering to each cell before scaling.
#'   The default, `FALSE`, preserves covariance along the all-ones direction.
#' @inheritParams covstatis
#'
#' @details
#' The common map is an ordinary [covstatis()] fit on each subject's average
#' preprocessed covariance, with no further preprocessing. Variable weights
#' apply to every cell in both analyses. The estimated STATIS table weights
#' apply only to the common map.
#'
#' For the condition analysis, subtract each subject's mean matrix from their
#' cells, then average these paired deviations across subjects with equal
#' weights. Decompose the resulting condition means using the Frobenius inner
#' product. For each unit-norm loading matrix, the cell score is its inner
#' product with that cell's paired deviation. Component signs are arbitrary.
#' A tied set of components identifies a subspace up to rotation.
#'
#' For two conditions, subtract the first fitted condition effect from the
#' second to obtain the mean paired difference in the preprocessed units.
#' Eigenvectors of this signed difference give variable combinations along
#' which projected variance increased or decreased.
#'
#' @return A `covstatis_rm` list with the following fields:
#' * `common` is the fitted COVSTATIS object on subject means.
#' * `weights` contains the variable weights used for raw cells at fit and projection.
#' * `condition_effects` contains one signed mean deviation matrix per condition.
#' * `condition_maps` contains the retained Frobenius-orthonormal loading matrices.
#' * `condition_scores` gives condition centroids on the retained axes.
#' * `paired_scores` gives subject, condition, and component scores in input order.
#' * `explained` gives the fraction of total condition-effect energy on each
#'   retained axis. `sums_of_squares` partitions total variation about the
#'   equal-cell grand mean into subject, condition, and residual terms. The
#'   residual includes subject-by-condition variation and estimation error.
#'
#' A zero condition effect yields an empty list of maps and zero score columns.
#' [project_subjects()] applies the fitted maps to complete panels of new subjects.
#' @references Madssen TS et al. (2021). Repeated measures ASCA+ for analysis
#'   of longitudinal intervention studies with multivariate outcome data.
#'   \doi{10.1371/journal.pcbi.1009585}.
#' @seealso [covstatis()], [design_mfa()]
#' @examples
#' subject <- rep(c("s1", "s2", "s3"), each = 2)
#' condition <- rep(c("rest", "task"), 3)
#' cells <- lapply(seq_along(subject), function(i) {
#'   S <- diag(c(5 + match(subject[i], unique(subject)), 2, 2))
#'   S[2, 3] <- S[3, 2] <- if (condition[i] == "task") 0.5 else 0
#'   S
#' })
#' fit <- covstatis_rm(cells, subject, condition)
#' fit$common$roi_scores
#' fit$paired_scores
#' delta <- fit$condition_effects$task - fit$condition_effects$rest
#' eigen(delta, symmetric = TRUE)$values
#' @md
#' @export
covstatis_rm <- function(data, subject, condition, ncomp = 2L,
                         ncomp_condition = NULL,
                         norm_method = c("none", "frobenius", "mfa"),
                         dcenter = FALSE, weights = NULL) {
  norm_method <- match.arg(norm_method)
  .distatis_count(ncomp, "ncomp")
  if (!is.null(ncomp_condition)) .distatis_count(ncomp_condition, "ncomp_condition")
  .distatis_flag(dcenter, "dcenter")

  panel <- .covstatis_rm_panel(data, subject, condition, norm_method, dcenter, weights)
  n <- length(panel$means)
  if (n < 2L) stop("Fitting requires at least two subjects.", call. = FALSE)
  k <- length(panel$conditions)
  common <- covstatis(panel$means, ncomp = ncomp, normalize = FALSE, dcenter = FALSE)
  rownames(common$roi_scores) <- common$labels

  effect <- Reduce("+", panel$within) / n
  # Average cancellations can leave roundoff on a truly null condition effect.
  abs_mean <- Reduce("+", lapply(panel$within, abs)) / n
  effect[abs(effect) <= 8 * .Machine$double.eps * abs_mean] <- 0
  scale <- max(abs(effect))
  rank <- min(ncomp_condition %||% (k - 1L), k - 1L)
  # Rescale for a rank decision relative to the effect, independent of units.
  decomp <- .dmfa_thin_svd(if (scale > 0) effect / scale else effect,
                           weights = rep(1, k), ncomp = rank,
                           rank_tol = k * .Machine$double.eps)
  axes <- colnames(decomp$v)
  scores <- decomp$scores * scale
  dimnames(scores) <- list(panel$conditions, axes)
  maps <- lapply(seq_len(ncol(decomp$v)), function(j) panel$unpack(decomp$v[, j]))
  names(maps) <- axes
  condition_effects <- lapply(seq_len(k), function(j) panel$unpack(effect[j, ]))
  names(condition_effects) <- panel$conditions

  grand <- Reduce("+", panel$means) / n
  ss <- c(subject = k * sum(vapply(panel$means, function(m) sum((m - grand)^2), numeric(1))),
          condition = n * sum(effect^2),
          residual = sum(vapply(panel$within, function(w) sum((w - effect)^2), numeric(1))))

  structure(list(common = common, condition_effects = condition_effects,
                 condition_maps = maps, condition_scores = scores,
                 paired_scores = .covstatis_rm_scores(panel, decomp$v),
                 explained = stats::setNames(decomp$d^2 / decomp$total_energy, axes),
                 sums_of_squares = ss, conditions = panel$conditions,
                 input_labels = panel$input_labels,
                 norm_method = norm_method, dcenter = dcenter, weights = panel$weights),
            class = "covstatis_rm")
}

#' @rdname covstatis_rm
#' @param x A fitted `covstatis_rm` model.
#' @param new_data Covariance matrices for one or more new subjects, with all
#'   training conditions observed once per subject. Supply `subject` and
#'   `condition` as for the training data.
#' @param ... Reserved for methods.
#' @return For `project_subjects()`, a list with `paired_scores` in input order
#'   and `common_scores`, one partial variable map per new subject.
#' @export
project_subjects.covstatis_rm <- function(x, new_data, subject, condition, ...) {
  .distatis_unused(...)
  panel <- .covstatis_rm_panel(new_data, subject, condition, x$norm_method,
                              x$dcenter, x$weights, fit = x)
  loadings <- matrix(0, length(panel$entry_weights), length(x$condition_maps),
                     dimnames = list(NULL, names(x$condition_maps)))
  for (j in seq_along(x$condition_maps)) {
    loadings[, j] <- x$condition_maps[[j]][panel$upper] * panel$entry_weights
  }
  list(paired_scores = .covstatis_rm_scores(panel, loadings),
       common_scores = lapply(panel$means, function(m) project_cov(x$common, m)))
}

#' @export
print.covstatis_rm <- function(x, ...) {
  cat("Repeated-measures COVSTATIS\n")
  cat(length(x$common$alpha), "subjects,", length(x$conditions), "conditions\n")
  cat(ncol(x$common$roi_scores), "common axes,", length(x$condition_maps),
      "condition axes\n")
  invisible(x)
}

.covstatis_rm_scores <- function(panel, loadings) {
  out <- matrix(0, nrow(panel$design), ncol(loadings),
                dimnames = list(NULL, colnames(loadings)))
  for (i in seq_along(panel$rows)) {
    out[panel$rows[[i]], ] <- panel$within[[i]] %*% loadings
  }
  cbind(panel$design, as.data.frame(out))
}

.covstatis_rm_id <- function(x, n, name) {
  if (!is.atomic(x) || !is.null(dim(x)) || length(x) != n || anyNA(x) ||
      (is.numeric(x) && any(!is.finite(x))) || any(!nzchar(as.character(x)))) {
    stop(sprintf("Supply one nonmissing %s identifier per covariance matrix.", name),
         call. = FALSE)
  }
  as.character(x)
}

.covstatis_rm_panel <- function(data, subject, condition, norm_method, dcenter, weights, fit = NULL) {
  if (!is.list(data) || is.data.frame(data) || !length(data)) {
    stop("`data` must be a nonempty list of covariance matrices.", call. = FALSE)
  }
  n <- length(data)
  design <- data.frame(subject = .covstatis_rm_id(subject, n, "subject"),
                       condition = .covstatis_rm_id(condition, n, "condition"),
                       stringsAsFactors = FALSE)
  conditions <- if (!is.null(fit)) {
    fit$conditions
  } else if (is.factor(condition)) {
    levels(droplevels(condition))
  } else {
    unique(design$condition)
  }
  if (length(conditions) < 2L) stop("At least two conditions are required.", call. = FALSE)
  rows <- split(seq_len(n), factor(design$subject, levels = unique(design$subject)))
  for (id in names(rows)) {
    idx <- rows[[id]]
    order <- match(conditions, design$condition[idx])
    if (length(idx) != length(conditions) || anyNA(order)) {
      stop(sprintf("Subject '%s' must have exactly one matrix per condition.", id),
           call. = FALSE)
    }
    rows[[id]] <- idx[order]
  }

  first <- data[[1L]]
  if (!is.matrix(first) || nrow(first) < 1L) stop("Expected a square covariance matrix.", call. = FALSE)
  p <- if (is.null(fit)) nrow(first) else length(fit$common$labels)
  labels <- if (is.null(fit)) rownames(first) %||% colnames(first) else fit$input_labels
  if (!is.null(labels)) .distatis_labels(labels, p, "variable names")
  weights <- .covstatis_weights(weights, p, labels)
  var_names <- labels %||% paste0("Obs_", seq_len(p))
  upper <- upper.tri(matrix(0, p, p), diag = TRUE)
  entry_weights <- ifelse(diag(p)[upper] == 1, 1, sqrt(2))
  unpack <- function(v) {
    m <- matrix(0, p, p, dimnames = list(var_names, var_names))
    m[upper] <- v / entry_weights
    m[lower.tri(m)] <- t(m)[lower.tri(m)]
    m
  }
  pack <- function(i) {
    m <- data[[i]]
    if (!is.numeric(m) || !is.matrix(m) || !identical(dim(m), c(p, p)) ||
        any(!is.finite(m))) {
      stop(sprintf("Cell %d must be a finite numeric %d by %d matrix.", i, p, p), call. = FALSE)
    }
    rn <- rownames(m)
    cn <- colnames(m)
    if ((!is.null(rn) && !is.null(cn) && !identical(rn, cn)) ||
        !identical(rn %||% cn, labels)) {
      stop(sprintf("Variable names and order must match across both axes and cells (cell %d).", i),
           call. = FALSE)
    }
    if (max(abs(m - t(m))) > sqrt(.Machine$double.eps) * max(abs(m))) {
      stop(sprintf("Covariance matrix in cell %d must be symmetric.", i), call. = FALSE)
    }
    storage.mode(m) <- "double"
    m <- m / 2 + t(m) / 2
    m <- .covstatis_preprocess(m, weights, dcenter, norm_method)
    m[upper] * entry_weights
  }

  means <- within <- vector("list", length(rows))
  names(means) <- names(within) <- names(rows)
  for (i in seq_along(rows)) {
    cells <- matrix(unlist(lapply(rows[[i]], pack), use.names = FALSE),
                    nrow = length(conditions), byrow = TRUE)
    # Difference first so identical cells yield exact zero paired deviations.
    delta <- sweep(cells, 2L, cells[1L, ], "-")
    if (dcenter || norm_method != "none") {
      delta[abs(delta) <= 8 * .Machine$double.eps * abs(cells)] <- 0
    }
    offset <- colMeans(delta)
    within[[i]] <- sweep(delta, 2L, offset, "-")
    means[[i]] <- unpack(colMeans(cells))
  }
  list(design = design, conditions = conditions, rows = rows, means = means,
       within = within, unpack = unpack, upper = upper, entry_weights = entry_weights,
       weights = weights, input_labels = labels)
}

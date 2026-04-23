#' Orient Component Signs for Fitted muscal Models
#'
#' Post-process a fitted model by applying a shared per-component sign choice
#' across scores, loadings, and cached derived quantities. This does not refit
#' the model. Reconstructions remain unchanged under the sign flip; projected
#' scores follow the new orientation.
#'
#' Supported first-wave classes are `mfa`, `aligned_mfa`, `anchored_mfa` /
#' `linked_mfa` (including graph-anchored variants), `mcca`, `aligned_mcca` /
#' `anchored_mcca`, `response_aligned_mfa`, `aligned_rrr`,
#' `aligned_interbattery`, `ipca`, and `penalized_mfa`.
#'
#' The sign convention can be supplied explicitly via `signs`, inferred from a
#' basis matrix or list of matrices via `basis`, or aligned to a reference fit
#' or matrix via `reference`.
#'
#' @param x A fitted model object.
#' @param basis Basis used to choose the component signs. May be `NULL` (use a
#'   class-specific default), a character scalar naming a stored field in `x`,
#'   a numeric matrix/data frame, or a list of matrices that will be row-bound
#'   before computing the sign rule.
#' @param reference Optional reference used to align signs. May be another
#'   fitted model, a numeric matrix/data frame, or a list of matrices. When
#'   supplied, signs are chosen to maximize per-component agreement between the
#'   resolved basis for `x` and the resolved basis for `reference`.
#' @param signs Optional numeric vector of length equal to the number of
#'   components. Non-zero values are coerced to `-1` or `+1`. When supplied,
#'   `basis` and `reference` are ignored.
#' @param rule Rule used when `reference = NULL` and `signs = NULL`. `"max_abs"`
#'   makes the largest-magnitude entry in each basis column positive. `"sum"`
#'   makes the column sum positive.
#' @param tol Non-negative tolerance used when a sign decision is numerically
#'   ambiguous.
#' @param ... Unused; reserved for future extensions.
#'
#' @return An oriented copy of `x`.
#' @export
orient_components <- function(x,
                              basis = NULL,
                              reference = NULL,
                              signs = NULL,
                              rule = c("max_abs", "sum"),
                              tol = 1e-8,
                              ...) {
  UseMethod("orient_components")
}

#' @rdname orient_components
#' @export
orient_components.mfa <- function(x,
                                  basis = NULL,
                                  reference = NULL,
                                  signs = NULL,
                                  rule = c("max_abs", "sum"),
                                  tol = 1e-8,
                                  ...) {
  rule <- match.arg(rule)
  signs <- .muscal_resolve_component_signs(
    x = x,
    basis = basis,
    default_basis = "v",
    reference = reference,
    signs = signs,
    rule = rule,
    tol = tol
  )

  .muscal_apply_component_signs(
    x,
    signs = signs,
    matrix_fields = c("s", "v", "cor_loadings", "ou", "ov"),
    list_fields = c("partial_scores")
  )
}

#' @rdname orient_components
#' @export
orient_components.aligned_mfa <- function(x,
                                          basis = NULL,
                                          reference = NULL,
                                          signs = NULL,
                                          rule = c("max_abs", "sum"),
                                          tol = 1e-8,
                                          ...) {
  rule <- match.arg(rule)
  signs <- .muscal_resolve_component_signs(
    x = x,
    basis = basis,
    default_basis = "v",
    reference = reference,
    signs = signs,
    rule = rule,
    tol = tol
  )

  .muscal_apply_component_signs(
    x,
    signs = signs,
    matrix_fields = c("s", "v", "cor_loadings"),
    list_fields = c("V_list", "partial_scores")
  )
}

#' @rdname orient_components
#' @export
orient_components.linked_mfa <- function(x,
                                         basis = NULL,
                                         reference = NULL,
                                         signs = NULL,
                                         rule = c("max_abs", "sum"),
                                         tol = 1e-8,
                                         ...) {
  rule <- match.arg(rule)
  signs <- .muscal_resolve_component_signs(
    x = x,
    basis = basis,
    default_basis = "B",
    reference = reference,
    signs = signs,
    rule = rule,
    tol = tol
  )

  .muscal_apply_component_signs(
    x,
    signs = signs,
    matrix_fields = c("s", "v", "B", "cor_loadings"),
    list_fields = c("V_list", "Z_list", "partial_scores")
  )
}

#' @rdname orient_components
#' @export
orient_components.aligned_interbattery <- function(x,
                                                   basis = NULL,
                                                   reference = NULL,
                                                   signs = NULL,
                                                   rule = c("max_abs", "sum"),
                                                   tol = 1e-8,
                                                   ...) {
  rule <- match.arg(rule)
  signs <- .muscal_resolve_component_signs(
    x = x,
    basis = basis,
    default_basis = "v",
    reference = reference,
    signs = signs,
    rule = rule,
    tol = tol
  )

  out <- .muscal_apply_component_signs(
    x,
    signs = signs,
    matrix_fields = c("s", "v", "common_scores"),
    list_fields = c("T_list", "U_list", "V_list", "x_loadings", "y_loadings", "x_block_scores", "y_block_scores", "partial_scores")
  )

  D <- diag(signs, nrow = length(signs), ncol = length(signs))
  for (nm in c(
    "x_coupling",
    "y_coupling",
    "x_coupling_operator",
    "y_coupling_operator",
    "x_whitener",
    "y_whitener",
    "x_unwhitener",
    "y_unwhitener"
  )) {
    if (!is.null(out[[nm]])) {
      out[[nm]] <- D %*% as.matrix(out[[nm]]) %*% D
    }
  }

  if (!is.null(out$coupling_maps)) {
    out$coupling_maps <- lapply(out$coupling_maps, function(M) D %*% M %*% D)
  }

  out
}

#' @rdname orient_components
#' @export
orient_components.mcca <- function(x,
                                   basis = NULL,
                                   reference = NULL,
                                   signs = NULL,
                                   rule = c("max_abs", "sum"),
                                   tol = 1e-8,
                                   ...) {
  rule <- match.arg(rule)
  signs <- .muscal_resolve_component_signs(
    x = x,
    basis = basis,
    default_basis = "v",
    reference = reference,
    signs = signs,
    rule = rule,
    tol = tol
  )

  .muscal_apply_component_signs(
    x,
    signs = signs,
    matrix_fields = c("s", "v", "cor_loadings"),
    list_fields = c("partial_scores", "canonical_weights")
  )
}

#' @rdname orient_components
#' @export
orient_components.aligned_mcca <- function(x,
                                           basis = NULL,
                                           reference = NULL,
                                           signs = NULL,
                                           rule = c("max_abs", "sum"),
                                           tol = 1e-8,
                                           ...) {
  rule <- match.arg(rule)
  signs <- .muscal_resolve_component_signs(
    x = x,
    basis = basis,
    default_basis = "v",
    reference = reference,
    signs = signs,
    rule = rule,
    tol = tol
  )

  .muscal_apply_component_signs(
    x,
    signs = signs,
    matrix_fields = c("s", "v", "cor_loadings"),
    list_fields = c("partial_scores", "canonical_weights")
  )
}

#' @rdname orient_components
#' @export
orient_components.response_aligned_mfa <- function(x,
                                                   basis = NULL,
                                                   reference = NULL,
                                                   signs = NULL,
                                                   rule = c("max_abs", "sum"),
                                                   tol = 1e-8,
                                                   ...) {
  rule <- match.arg(rule)
  signs <- .muscal_resolve_component_signs(
    x = x,
    basis = basis,
    default_basis = "B",
    reference = reference,
    signs = signs,
    rule = rule,
    tol = tol
  )

  .muscal_apply_component_signs(
    x,
    signs = signs,
    matrix_fields = c("s", "v", "B", "S", "cor_loadings"),
    list_fields = c("V_list", "Z_list", "partial_scores")
  )
}

#' @rdname orient_components
#' @export
orient_components.aligned_rrr <- function(x,
                                          basis = NULL,
                                          reference = NULL,
                                          signs = NULL,
                                          rule = c("max_abs", "sum"),
                                          tol = 1e-8,
                                          ...) {
  rule <- match.arg(rule)
  signs <- .muscal_resolve_component_signs(
    x = x,
    basis = basis,
    default_basis = "B",
    reference = reference,
    signs = signs,
    rule = rule,
    tol = tol
  )

  .muscal_apply_component_signs(
    x,
    signs = signs,
    matrix_fields = c("s", "v", "B", "cor_loadings"),
    list_fields = c("W_list", "Z_list", "partial_scores")
  )
}

#' @rdname orient_components
#' @export
orient_components.ipca <- function(x,
                                   basis = NULL,
                                   reference = NULL,
                                   signs = NULL,
                                   rule = c("max_abs", "sum"),
                                   tol = 1e-8,
                                   ...) {
  rule <- match.arg(rule)
  signs <- .muscal_resolve_component_signs(
    x = x,
    basis = basis,
    default_basis = "v",
    reference = reference,
    signs = signs,
    rule = rule,
    tol = tol
  )

  out <- .muscal_apply_component_signs(
    x,
    signs = signs,
    matrix_fields = c("s", "v"),
    list_fields = c("partial_scores")
  )

  if (!is.null(out$projection_map) && !is.null(out$projection_map$coef_blocks)) {
    out$projection_map$coef_blocks <- .muscal_apply_component_signs_to_list(
      out$projection_map$coef_blocks,
      signs
    )
  }
  if (!is.null(out$Sigma_eigenvectors)) {
    out$Sigma_eigenvectors <- .muscal_apply_component_signs_to_matrix(
      out$Sigma_eigenvectors,
      signs,
      allow_extra_cols = TRUE
    )
  }
  if (!is.null(out$warm_state) && !is.null(out$warm_state$U)) {
    out$warm_state$U <- .muscal_apply_component_signs_to_matrix(
      out$warm_state$U,
      signs,
      allow_extra_cols = TRUE
    )
  }

  out
}

#' @rdname orient_components
#' @export
orient_components.penalized_mfa <- function(x,
                                            basis = NULL,
                                            reference = NULL,
                                            signs = NULL,
                                            rule = c("max_abs", "sum"),
                                            tol = 1e-8,
                                            ...) {
  rule <- match.arg(rule)
  signs <- .muscal_resolve_component_signs(
    x = x,
    basis = basis,
    default_basis = "v",
    reference = reference,
    signs = signs,
    rule = rule,
    tol = tol
  )

  out <- .muscal_apply_component_signs(
    x,
    signs = signs,
    matrix_fields = c("s", "v", "cor_loadings"),
    list_fields = c("V_list", "scores_list")
  )

  V_attr <- attr(out, "V_list", exact = TRUE)
  if (!is.null(V_attr)) {
    attr(out, "V_list") <- .muscal_apply_component_signs_to_list(V_attr, signs)
  }
  consensus_attr <- attr(out, "consensus", exact = TRUE)
  if (!is.null(consensus_attr)) {
    attr(out, "consensus") <- .muscal_apply_component_signs_to_matrix(consensus_attr, signs)
  }

  out
}

#' @rdname orient_components
#' @export
orient_components.default <- function(x,
                                      basis = NULL,
                                      reference = NULL,
                                      signs = NULL,
                                      rule = c("max_abs", "sum"),
                                      tol = 1e-8,
                                      ...) {
  stop(
    sprintf(
      "No orient_components() method is available for class %s.",
      paste(class(x), collapse = "/")
    ),
    call. = FALSE
  )
}

.muscal_resolve_component_signs <- function(x,
                                            basis = NULL,
                                            default_basis = "v",
                                            reference = NULL,
                                            signs = NULL,
                                            rule = c("max_abs", "sum"),
                                            tol = 1e-8) {
  rule <- match.arg(rule)
  ncomp <- .muscal_component_count(x)

  if (!is.null(signs)) {
    return(.muscal_normalize_component_signs(signs, ncomp))
  }

  basis_eff <- .muscal_resolve_orientation_basis(
    x = x,
    basis = basis,
    default_basis = default_basis,
    ncomp = ncomp,
    what = "x"
  )

  reference_eff <- if (is.null(reference)) {
    NULL
  } else {
    .muscal_resolve_orientation_basis(
      x = reference,
      basis = basis,
      default_basis = default_basis,
      ncomp = ncomp,
      what = "reference"
    )
  }

  .muscal_component_signs_from_basis(
    basis = basis_eff,
    rule = rule,
    reference = reference_eff,
    tol = tol
  )
}

.muscal_component_count <- function(x) {
  if (!is.null(x$s)) return(ncol(as.matrix(x$s)))
  if (!is.null(x$v)) return(ncol(as.matrix(x$v)))
  if (!is.null(x$sdev)) return(length(x$sdev))
  stop("Unable to infer the number of components from `x`.", call. = FALSE)
}

.muscal_normalize_component_signs <- function(signs, ncomp) {
  signs <- as.numeric(signs)
  if (length(signs) != ncomp) {
    stop(sprintf("`signs` must have length %d.", ncomp), call. = FALSE)
  }
  if (anyNA(signs) || any(!is.finite(signs)) || any(abs(signs) <= 0)) {
    stop("`signs` must be finite and non-zero.", call. = FALSE)
  }
  ifelse(signs > 0, 1, -1)
}

.muscal_resolve_orientation_basis <- function(x,
                                              basis = NULL,
                                              default_basis = "v",
                                              ncomp = NULL,
                                              what = "x") {
  basis_obj <- basis %||% default_basis

  if (is.character(basis_obj)) {
    if (length(basis_obj) != 1L || is.na(basis_obj) || basis_obj == "") {
      stop("`basis` must be a non-empty character scalar when supplied by name.", call. = FALSE)
    }

    basis_name <- basis_obj
    if (!is.null(x[[basis_name]])) {
      basis_obj <- x[[basis_name]]
    } else {
      basis_attr <- attr(x, basis_name, exact = TRUE)
      if (is.null(basis_attr)) {
        stop(sprintf("Could not find basis '%s' in %s.", basis_name, what), call. = FALSE)
      }
      basis_obj <- basis_attr
    }
  }

  .muscal_as_component_basis_matrix(basis_obj, ncomp = ncomp, what = what)
}

.muscal_as_component_basis_matrix <- function(x, ncomp = NULL, what = "basis") {
  if (is.data.frame(x)) {
    x <- as.matrix(x)
  }

  if (is.matrix(x)) {
    basis <- x
  } else if (is.list(x) &&
             length(x) > 0L &&
             all(vapply(x, function(el) is.matrix(el) || is.data.frame(el), logical(1)))) {
    mats <- lapply(x, as.matrix)
    ncols <- vapply(mats, ncol, integer(1))
    if (!all(ncols == ncols[1])) {
      stop(sprintf("All matrices in %s must have the same number of columns.", what), call. = FALSE)
    }
    basis <- do.call(rbind, mats)
  } else {
    stop(
      sprintf(
        "%s must resolve to a matrix/data.frame or a list of matrices/data.frames.",
        what
      ),
      call. = FALSE
    )
  }

  basis <- as.matrix(basis)
  if (nrow(basis) < 1L || ncol(basis) < 1L) {
    stop(sprintf("%s must contain at least one row and one column.", what), call. = FALSE)
  }
  if (!is.null(ncomp) && ncol(basis) != ncomp) {
    stop(sprintf("%s must have %d columns.", what, ncomp), call. = FALSE)
  }

  basis
}

.muscal_component_signs_from_basis <- function(basis,
                                               rule = c("max_abs", "sum"),
                                               reference = NULL,
                                               tol = 1e-8) {
  rule <- match.arg(rule)
  tol <- as.numeric(tol)[1]
  if (!is.finite(tol) || tol < 0) {
    stop("`tol` must be finite and non-negative.", call. = FALSE)
  }

  basis <- as.matrix(basis)
  K <- ncol(basis)
  signs <- rep.int(1, K)

  fallback_sign <- function(x) {
    x <- x[is.finite(x)]
    if (length(x) == 0L) return(1)
    if (rule == "sum") {
      s <- sum(x)
      if (!is.finite(s) || abs(s) <= tol) return(1)
      return(if (s >= 0) 1 else -1)
    }

    idx <- which.max(abs(x))
    val <- x[[idx]]
    if (!is.finite(val) || abs(val) <= tol) 1 else if (val >= 0) 1 else -1
  }

  if (!is.null(reference)) {
    reference <- as.matrix(reference)
    if (!identical(dim(reference), dim(basis))) {
      stop("`reference` must have the same dimensions as the resolved basis.", call. = FALSE)
    }

    for (j in seq_len(K)) {
      dots <- basis[, j] * reference[, j]
      dots <- dots[is.finite(dots)]
      dot_j <- if (length(dots) == 0L) 0 else sum(dots)
      signs[[j]] <- if (!is.finite(dot_j) || abs(dot_j) <= tol) {
        fallback_sign(basis[, j])
      } else if (dot_j >= 0) {
        1
      } else {
        -1
      }
    }
    return(signs)
  }

  for (j in seq_len(K)) {
    signs[[j]] <- fallback_sign(basis[, j])
  }
  signs
}

.muscal_apply_component_signs_to_matrix <- function(x,
                                                    signs,
                                                    allow_extra_cols = FALSE) {
  if (is.null(x)) return(NULL)
  x <- as.matrix(x)
  if (ncol(x) < length(signs)) {
    stop(
      sprintf("Cannot orient a matrix with %d columns using %d component signs.", ncol(x), length(signs)),
      call. = FALSE
    )
  }
  if (ncol(x) > length(signs) && !isTRUE(allow_extra_cols)) {
    stop(
      sprintf("Cannot orient a matrix with %d columns using %d component signs.", ncol(x), length(signs)),
      call. = FALSE
    )
  }

  if (ncol(x) == length(signs)) {
    return(sweep(x, 2, signs, `*`))
  }

  out <- x
  out[, seq_along(signs)] <- sweep(x[, seq_along(signs), drop = FALSE], 2, signs, `*`)
  out
}

.muscal_apply_component_signs_to_list <- function(x, signs) {
  if (is.null(x)) return(NULL)
  if (!is.list(x)) {
    stop("Expected a list of matrices when applying component signs.", call. = FALSE)
  }
  lapply(x, .muscal_apply_component_signs_to_matrix, signs = signs)
}

.muscal_apply_component_signs <- function(x,
                                          signs,
                                          matrix_fields = character(),
                                          list_fields = character()) {
  out <- x

  for (nm in matrix_fields) {
    if (!is.null(out[[nm]])) {
      out[[nm]] <- .muscal_apply_component_signs_to_matrix(out[[nm]], signs)
    }
  }

  for (nm in list_fields) {
    if (!is.null(out[[nm]])) {
      out[[nm]] <- .muscal_apply_component_signs_to_list(out[[nm]], signs)
    }
  }

  out
}

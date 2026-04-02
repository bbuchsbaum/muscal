#' Muscal Fit Object Contract
#'
#' CV-capable `muscal` methods should expose a small common contract so they can
#' participate in shared cross-validation and inference workflows.
#'
#' The current contract stores:
#' - `task`: the primary evaluation task
#' - `oos_types`: the supported `predict(type = )` outputs for out-of-sample use
#' - `fit_spec`: a small metadata record describing the fitting method and whether
#'   generic refitting is currently supported
#' - `refit`: optional metadata and callbacks for standard resampling workflows
#'
#' The contract does not guarantee refitting for every method. Methods that
#' support standard resampling attach a `refit` record describing how to rebuild
#' the fit from stored training data and how to generate bootstrap or permutation
#' replicates.
#'
#' @param method Character scalar naming the fitting method.
#' @param task Character scalar naming the primary evaluation task.
#' @param oos_types Character vector naming the supported out-of-sample
#'   prediction types.
#' @param fit_call Optional matched call captured at fit time.
#' @param refit_supported Logical; whether the object currently stores enough
#'   metadata for a generic refit path.
#' @param prediction_target Optional character scalar naming the high-level
#'   prediction target.
#' @param refit Optional refit specification created by
#'   `.muscal_make_refit_spec()`.
#'
#' @return A named list storing contract metadata.
#' @keywords internal
muscal_fit_contract <- function(method,
                                task,
                                oos_types,
                                fit_call = NULL,
                                refit_supported = FALSE,
                                prediction_target = NULL,
                                refit = NULL) {
  chk::chk_string(method)
  chk::chk_string(task)
  chk::chk_flag(refit_supported)

  oos_types <- unique(as.character(oos_types))
  if (length(oos_types) == 0L || anyNA(oos_types) || any(oos_types == "")) {
    stop("`oos_types` must contain at least one non-empty prediction type.", call. = FALSE)
  }

  if (!is.null(prediction_target)) {
    chk::chk_string(prediction_target)
  }

  if (!is.null(refit) && !is.list(refit)) {
    stop("`refit` must be NULL or a refit specification list.", call. = FALSE)
  }

  list(
    method = method,
    task = task,
    oos_types = oos_types,
    prediction_target = prediction_target,
    refit_supported = isTRUE(refit_supported) || !is.null(refit),
    refit = refit,
    call = fit_call
  )
}

.muscal_attach_fit_contract <- function(object,
                                        method,
                                        task,
                                        oos_types,
                                        fit_call = NULL,
                                        refit_supported = FALSE,
                                        prediction_target = NULL,
                                        refit = NULL) {
  contract <- muscal_fit_contract(
    method = method,
    task = task,
    oos_types = oos_types,
    fit_call = fit_call,
    refit_supported = refit_supported,
    prediction_target = prediction_target,
    refit = refit
  )

  object$task <- contract$task
  object$oos_types <- contract$oos_types
  object$fit_spec <- contract
  object
}

.muscal_make_refit_spec <- function(data,
                                    fit_fn,
                                    bootstrap_fn,
                                    permutation_fn,
                                    resample_unit = "rows") {
  chk::chk_function(fit_fn)
  chk::chk_function(bootstrap_fn)
  chk::chk_function(permutation_fn)
  chk::chk_string(resample_unit)

  list(
    data = data,
    fit_fn = fit_fn,
    bootstrap_fn = bootstrap_fn,
    permutation_fn = permutation_fn,
    resample_unit = resample_unit
  )
}

.muscal_resolve_refit_spec <- function(object, refit = NULL) {
  spec <- refit %||% object$fit_spec$refit
  if (is.null(spec) || !isTRUE(object$fit_spec$refit_supported)) {
    stop(
      "This fit does not expose refit metadata. Supply `refit =` explicitly or use a method with `fit_spec$refit_supported = TRUE`.",
      call. = FALSE
    )
  }

  required <- c("data", "fit_fn", "bootstrap_fn", "permutation_fn")
  missing <- setdiff(required, names(spec))
  if (length(missing) > 0L) {
    stop(
      sprintf("Refit spec is missing required fields: %s.", paste(missing, collapse = ", ")),
      call. = FALSE
    )
  }

  if (!is.function(spec$fit_fn) ||
      !is.function(spec$bootstrap_fn) ||
      !is.function(spec$permutation_fn)) {
    stop("Refit spec callbacks must be functions.", call. = FALSE)
  }

  spec
}

.muscal_resolve_block_id <- function(block, block_names, arg = "block") {
  if (missing(block) || is.null(block)) {
    stop(sprintf("`%s` must identify a block.", arg), call. = FALSE)
  }
  if (is.null(block_names) || length(block_names) == 0L) {
    stop("No block names are available in this fit object.", call. = FALSE)
  }

  if (is.character(block)) {
    if (length(block) != 1L || is.na(block)) {
      stop(sprintf("`%s` must be a single non-missing character block name.", arg), call. = FALSE)
    }
    bi <- match(block, block_names)
    if (is.na(bi)) {
      stop(
        sprintf("Unknown block '%s'. Available blocks: %s.", block, paste(block_names, collapse = ", ")),
        call. = FALSE
      )
    }
    return(list(index = bi, name = block_names[[bi]]))
  }

  if (!is.numeric(block) || length(block) != 1L || is.na(block) || block != as.integer(block)) {
    stop(sprintf("`%s` must be a single block name or integer index.", arg), call. = FALSE)
  }
  bi <- as.integer(block)
  if (bi < 1L || bi > length(block_names)) {
    stop(
      sprintf("`%s` index %d is outside the valid range 1:%d.", arg, bi, length(block_names)),
      call. = FALSE
    )
  }
  list(index = bi, name = block_names[[bi]])
}

.muscal_project_row_linked_scores <- function(x, new_data, block, preprocess = TRUE) {
  resolved <- .muscal_resolve_block_id(block, names(x$V_list))
  block_name <- resolved$name
  Xnew <- as.matrix(new_data)

  if (isTRUE(preprocess)) {
    Xnew <- multivarious::transform(x$block_preproc[[block_name]], Xnew)
  }

  Vk <- x$V_list[[block_name]]
  if (ncol(Xnew) != nrow(Vk)) {
    stop(
      sprintf("Block '%s' has %d columns, expected %d.", block_name, ncol(Xnew), nrow(Vk)),
      call. = FALSE
    )
  }

  G <- crossprod(Vk) + x$ridge * diag(ncol(Vk))
  Xnew %*% Vk %*% solve(G)
}

.muscal_predict_block_reconstruction <- function(object,
                                                 new_data,
                                                 block,
                                                 preprocess = TRUE) {
  resolved <- .muscal_resolve_block_id(block, names(object$V_list))
  block_name <- resolved$name
  scores_new <- .muscal_project_row_linked_scores(object, new_data, block = block, preprocess = preprocess)
  Xhat_p <- scores_new %*% t(object$V_list[[block_name]])
  multivarious::inverse_transform(object$block_preproc[[block_name]], Xhat_p)
}

.muscal_predict_anchor_response <- function(object,
                                            new_data,
                                            block,
                                            preprocess = TRUE) {
  scores_new <- .muscal_project_row_linked_scores(object, new_data, block = block, preprocess = preprocess)
  Yhat_p <- scores_new %*% t(object$B)
  multivarious::inverse_transform(object$anchor_preproc, Yhat_p)
}

#' Predict from a Multiple Factor Analysis Fit
#'
#' @param object A fitted `mfa` object.
#' @param new_data Optional matrix/data.frame of new observations in the full
#'   concatenated feature space.
#' @param type One of `"scores"` or `"reconstruction"`.
#' @param ... Additional arguments passed to the underlying projector helpers.
#'
#' @return A numeric matrix of projected scores or reconstructed observations.
#' @export
predict.mfa <- function(object, new_data = NULL,
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

#' Project New Rows into an Anchored MFA Space
#'
#' @param x A fitted `anchored_mfa` object.
#' @param new_data Numeric matrix/data.frame of new rows for a known auxiliary
#'   block.
#' @param block Character or integer identifying the auxiliary block used for
#'   projection.
#' @param preprocess Logical; if `TRUE` (default), applies the fitted block
#'   preprocessing pipeline before projection.
#' @param ... Unused.
#'
#' @return A numeric matrix of latent scores.
#' @export
project.anchored_mfa <- function(x, new_data, block, preprocess = TRUE, ...) {
  chk::chk_true(inherits(x, "anchored_mfa"))
  .muscal_project_row_linked_scores(x, new_data, block = block, preprocess = preprocess)
}

#' Predict from an Anchored MFA Fit
#'
#' @param object A fitted `anchored_mfa` object.
#' @param new_data Numeric matrix/data.frame of new rows from a known auxiliary
#'   block.
#' @param block Character or integer identifying the auxiliary block.
#' @param type One of `"response"`, `"scores"`, or `"reconstruction"`.
#' @param preprocess Logical; if `TRUE` (default), applies the fitted block
#'   preprocessing pipeline before projection.
#' @param ... Unused.
#'
#' @return A numeric matrix of predicted anchor responses, latent scores, or
#'   reconstructed block rows.
#' @export
predict.anchored_mfa <- function(object,
                                 new_data,
                                 block,
                                 type = c("response", "scores", "reconstruction"),
                                 preprocess = TRUE,
                                 ...) {
  type <- match.arg(type)

  if (type == "scores") {
    return(project.anchored_mfa(object, new_data, block = block, preprocess = preprocess, ...))
  }
  if (type == "reconstruction") {
    return(.muscal_predict_block_reconstruction(object, new_data, block = block, preprocess = preprocess))
  }

  .muscal_predict_anchor_response(object, new_data, block = block, preprocess = preprocess)
}

#' Project New Rows into an Aligned MFA Space
#'
#' @param x A fitted `aligned_mfa` object.
#' @param new_data Numeric matrix/data.frame of new rows for a known block.
#' @param block Character or integer identifying the block used for projection.
#' @param preprocess Logical; if `TRUE` (default), applies the fitted block
#'   preprocessing pipeline before projection.
#' @param ... Unused.
#'
#' @return A numeric matrix of latent scores.
#' @export
project.aligned_mfa <- function(x, new_data, block, preprocess = TRUE, ...) {
  chk::chk_true(inherits(x, "aligned_mfa"))
  .muscal_project_row_linked_scores(x, new_data, block = block, preprocess = preprocess)
}

#' Predict from an Aligned MFA Fit
#'
#' @param object A fitted `aligned_mfa` object.
#' @param new_data Numeric matrix/data.frame of new rows from a known block.
#' @param block Character or integer identifying the block.
#' @param type One of `"scores"` or `"reconstruction"`.
#' @param preprocess Logical; if `TRUE` (default), applies the fitted block
#'   preprocessing pipeline before projection.
#' @param ... Unused.
#'
#' @return A numeric matrix of latent scores or reconstructed block rows.
#' @export
predict.aligned_mfa <- function(object,
                                new_data,
                                block,
                                type = c("scores", "reconstruction"),
                                preprocess = TRUE,
                                ...) {
  type <- match.arg(type)

  if (type == "scores") {
    return(project.aligned_mfa(object, new_data, block = block, preprocess = preprocess, ...))
  }

  .muscal_predict_block_reconstruction(object, new_data, block = block, preprocess = preprocess)
}

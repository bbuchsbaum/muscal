#' Cross-Validate Muscal Models with Multidesign Folds
#'
#' A thin cross-validation wrapper that reuses `multidesign` fold objects and
#' the `muscal` fit contract. Users provide:
#' - `fit_fn(analysis)` to fit a model on the analysis split,
#' - `estimate_fn(model, assessment, ...)` to produce predictions or held-out
#'   estimates, and
#' - `truth_fn(assessment, ...)` to extract the corresponding truth object.
#'
#' The returned metrics are computed through [performance_metrics()], so fold-wise
#' outputs have a stable, task-aware tibble shape that works with
#' `multidesign::cross_validate()`.
#'
#' This function deliberately does not define a new split abstraction. Use
#' `multidesign::fold_over()` or `multidesign::cv_rows()` to construct `folds`.
#'
#' @param folds A `foldlist` object, typically created by
#'   `multidesign::fold_over()` or `multidesign::cv_rows()`.
#' @param fit_fn Function with signature `fit_fn(analysis)` returning a fitted
#'   `muscal` model.
#' @param estimate_fn Function that produces held-out estimates from
#'   `estimate_fn(model, assessment, ...)`.
#' @param truth_fn Function that extracts the truth object from
#'   `truth_fn(assessment, ...)`.
#' @param task Optional task override. If `NULL`, the task is taken from
#'   `model$task`.
#' @param metrics Optional metric vector passed to [performance_metrics()].
#' @param performance_args Optional named list of additional arguments forwarded
#'   to [performance_metrics()] on every fold.
#' @param ... Reserved for future extensions.
#'
#' @return An object inheriting from `cv_result` with additional fields `folds`,
#'   `foldframe`, `task`, `metrics`, and `call`.
#' @export
cv_muscal <- function(folds,
                      fit_fn,
                      estimate_fn,
                      truth_fn,
                      task = NULL,
                      metrics = NULL,
                      performance_args = list(),
                      ...) {
  chk::chk_s3_class(folds, "foldlist")
  chk::chk_function(fit_fn)
  chk::chk_function(estimate_fn)
  chk::chk_function(truth_fn)
  chk::chk_list(performance_args)

  call <- match.call()

  cv_result <- multidesign::cross_validate(
    folds = folds,
    fit_fn = fit_fn,
    score_fn = function(model, assessment, fold, fold_id) {
      task_eff <- if (is.null(task)) {
        if (is.null(model$task)) {
          stop("`task` was not supplied and the fitted model does not expose `model$task`.", call. = FALSE)
        }
        model$task
      } else {
        task
      }

      truth_obj <- .muscal_invoke_cv_callback(
        truth_fn,
        list(assessment = assessment, fold = fold, fold_id = fold_id)
      )
      estimate_obj <- .muscal_invoke_cv_callback(
        estimate_fn,
        list(model = model, assessment = assessment, fold = fold, fold_id = fold_id)
      )

      score <- do.call(
        performance_metrics,
        c(
          list(
            task = task_eff,
            truth = truth_obj,
            estimate = estimate_obj,
            metrics = metrics
          ),
          performance_args
        )
      )

      as.list(score[1, , drop = FALSE])
    }
  )

  structure(
    c(
      cv_result,
      list(
        folds = folds,
        foldframe = attr(folds, "foldframe"),
        task = if (is.null(task)) NULL else .muscal_normalize_metric_task(task),
        metrics = metrics
      )
    ),
    class = c("muscal_cv_result", class(cv_result))
  )
}

.muscal_invoke_cv_callback <- function(fn, available_args) {
  fn_formals <- names(formals(fn))
  has_dots <- "..." %in% fn_formals

  call_args <- list()
  for (nm in names(available_args)) {
    if (isTRUE(has_dots) || nm %in% fn_formals) {
      call_args[[nm]] <- available_args[[nm]]
    }
  }
  do.call(fn, call_args)
}

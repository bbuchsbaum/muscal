#' Metric Registry for Muscal Tasks
#'
#' Returns the task-aware performance metric registry used by
#' [performance_metrics()]. The registry records which metrics belong to which
#' evaluation task, whether larger values are better, and a short description of
#' the metric.
#'
#' Supported task families currently include:
#' - `"reconstruction"`
#' - `"response_prediction"`
#' - `"retrieval_alignment"`
#'
#' The alias `"row_alignment"` maps to `"retrieval_alignment"`.
#'
#' @param task Optional task name. If supplied, returns only the metrics for
#'   that task (after alias normalization).
#'
#' @return A tibble with one row per supported metric.
#' @export
metric_registry <- function(task = NULL) {
  registry <- tibble::tibble(
    task = c(
      rep("reconstruction", 4L),
      rep("response_prediction", 6L),
      rep("retrieval_alignment", 5L)
    ),
    metric = c(
      "mse", "rmse", "r2", "mae",
      "mse", "rmse", "r2", "mae", "mean_correlation", "mean_cosine_similarity",
      "mean_top1_similarity", "mean_topk_similarity", "oracle_gap", "recall_at_k", "mrr"
    ),
    maximize = c(
      FALSE, FALSE, TRUE, FALSE,
      FALSE, FALSE, TRUE, FALSE, TRUE, TRUE,
      TRUE, TRUE, FALSE, TRUE, TRUE
    ),
    description = c(
      "Mean squared reconstruction error",
      "Root mean squared reconstruction error",
      "Global reconstruction R-squared",
      "Mean absolute reconstruction error",
      "Mean squared response prediction error",
      "Root mean squared response prediction error",
      "Global response prediction R-squared",
      "Mean absolute response prediction error",
      "Mean row-wise Pearson correlation between truth and prediction",
      "Mean row-wise cosine similarity between truth and prediction",
      "Mean cosine similarity of the top retrieved item",
      "Mean cosine similarity averaged over the top-k retrieved items",
      "Average gap between the oracle similarity and the top retrieved similarity",
      "Fraction of queries whose true match appears in the top-k retrieved ids",
      "Mean reciprocal rank of the true match in the retrieved ids"
    )
  )

  if (is.null(task)) {
    return(registry)
  }

  task_norm <- .muscal_normalize_metric_task(task)
  registry[registry$task == task_norm, , drop = FALSE]
}

#' Default Metrics for a Muscal Task
#'
#' @param task Task name. `"row_alignment"` is treated as an alias for
#'   `"retrieval_alignment"`.
#'
#' @return A character vector of default metric names.
#' @export
default_metrics <- function(task) {
  task_norm <- .muscal_normalize_metric_task(task)

  switch(
    task_norm,
    reconstruction = c("mse", "rmse", "r2"),
    response_prediction = c("mse", "rmse", "r2", "mean_cosine_similarity"),
    retrieval_alignment = c("mean_top1_similarity", "mean_topk_similarity"),
    stop("Unsupported task: ", task, call. = FALSE)
  )
}

#' Compute Task-Aware Performance Metrics
#'
#' Computes a one-row tibble of performance metrics for a declared evaluation
#' task. This provides a stable metric interface for fold-wise cross-validation
#' and held-out evaluation.
#'
#' Supported task families currently include:
#' - `"reconstruction"`
#' - `"response_prediction"`
#' - `"retrieval_alignment"`
#'
#' The alias `"row_alignment"` maps to `"retrieval_alignment"`.
#'
#' For `"retrieval_alignment"`, `estimate` may be:
#' - a numeric matrix/data.frame interpreted as a single top-1 retrieved feature
#'   vector per query, or
#' - a list containing `retrieved_features` (a list of ranked feature matrices),
#'   plus optional `retrieved_ids`, `oracle_similarity`, or `oracle_features`.
#'
#' @param task Task name.
#' @param truth Truth object for the task. For reconstruction and response
#'   prediction this is a numeric matrix/data.frame. For retrieval/alignment this
#'   is a numeric matrix/data.frame of query-side feature vectors.
#' @param estimate Estimated or predicted object for the task.
#' @param metrics Optional character vector of metrics. If `NULL`, task-specific
#'   defaults are used.
#' @param k Integer `k` used by top-k retrieval metrics and recall@k.
#' @param truth_ids Optional vector of true ids for retrieval metrics such as
#'   `recall_at_k` and `mrr`.
#' @param by_column Logical; forwarded to reconstruction-style `r2`
#'   calculations.
#' @param ... Reserved for future task-specific options.
#'
#' @return A one-row tibble of metric values.
#' @export
performance_metrics <- function(task,
                                truth,
                                estimate,
                                metrics = NULL,
                                k = 5L,
                                truth_ids = NULL,
                                by_column = FALSE,
                                ...) {
  task_norm <- .muscal_normalize_metric_task(task)
  if (is.null(metrics)) {
    metrics <- default_metrics(task_norm)
  }
  metrics <- as.character(metrics)
  .muscal_validate_metric_names(task_norm, metrics)

  if (task_norm == "reconstruction") {
    truth_mat <- .muscal_as_numeric_matrix(truth, "truth")
    estimate_mat <- .muscal_as_numeric_matrix(estimate, "estimate")
    return(
      multivarious::measure_reconstruction_error(
        truth_mat,
        estimate_mat,
        metrics = metrics,
        by_column = by_column
      )
    )
  }

  if (task_norm == "response_prediction") {
    truth_mat <- .muscal_as_numeric_matrix(truth, "truth")
    estimate_mat <- .muscal_as_numeric_matrix(estimate, "estimate")
    .muscal_require_same_dims(truth_mat, estimate_mat, "truth", "estimate")
    return(.muscal_response_metrics(truth_mat, estimate_mat, metrics = metrics, by_column = by_column))
  }

  truth_mat <- .muscal_as_numeric_matrix(truth, "truth")
  retrieval <- .muscal_normalize_retrieval_estimate(estimate, n_query = nrow(truth_mat))
  .muscal_retrieval_metrics(
    truth = truth_mat,
    estimate = retrieval,
    metrics = metrics,
    k = as.integer(k),
    truth_ids = truth_ids
  )
}

.muscal_normalize_metric_task <- function(task) {
  chk::chk_string(task)

  task_map <- c(
    reconstruction = "reconstruction",
    response_prediction = "response_prediction",
    retrieval_alignment = "retrieval_alignment",
    row_alignment = "retrieval_alignment",
    retrieval = "retrieval_alignment",
    alignment = "retrieval_alignment"
  )

  task_norm <- unname(task_map[[task]])
  if (is.null(task_norm)) {
    stop("Unsupported task: ", task, call. = FALSE)
  }
  task_norm
}

.muscal_validate_metric_names <- function(task, metrics) {
  registry <- metric_registry(task)
  supported <- registry$metric
  bad <- setdiff(metrics, supported)
  if (length(bad) > 0L) {
    stop(
      sprintf(
        "Unsupported metrics for task '%s': %s. Supported metrics: %s.",
        task,
        paste(bad, collapse = ", "),
        paste(supported, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  invisible(metrics)
}

.muscal_as_numeric_matrix <- function(x, arg) {
  if (is.vector(x) && !is.list(x)) {
    x <- matrix(as.numeric(x), ncol = 1L)
  } else if (is.data.frame(x) || is.matrix(x)) {
    x <- as.matrix(x)
  } else {
    stop(sprintf("`%s` must be a numeric matrix, data.frame, or vector.", arg), call. = FALSE)
  }

  if (!is.numeric(x)) {
    stop(sprintf("`%s` must be numeric.", arg), call. = FALSE)
  }
  x
}

.muscal_require_same_dims <- function(x, y, x_name, y_name) {
  if (!identical(dim(x), dim(y))) {
    stop(
      sprintf(
        "`%s` has shape %s but `%s` has shape %s.",
        x_name,
        paste(dim(x), collapse = "x"),
        y_name,
        paste(dim(y), collapse = "x")
      ),
      call. = FALSE
    )
  }
}

.muscal_rowwise_cosine <- function(x, y) {
  .muscal_require_same_dims(x, y, "x", "y")
  num <- rowSums(x * y)
  den <- sqrt(rowSums(x^2)) * sqrt(rowSums(y^2))
  out <- ifelse(den > 0, num / den, NA_real_)
  as.numeric(out)
}

.muscal_rowwise_correlation <- function(x, y) {
  .muscal_require_same_dims(x, y, "x", "y")
  if (ncol(x) < 2L) {
    return(rep(NA_real_, nrow(x)))
  }

  out <- vapply(seq_len(nrow(x)), function(i) {
    xi <- x[i, ]
    yi <- y[i, ]
    if (stats::sd(xi) < 1e-12 || stats::sd(yi) < 1e-12) {
      return(NA_real_)
    }
    stats::cor(xi, yi)
  }, numeric(1))
  as.numeric(out)
}

.muscal_response_metrics <- function(truth, estimate, metrics, by_column = FALSE) {
  base_metrics <- intersect(metrics, c("mse", "rmse", "r2", "mae"))
  out <- list()

  if (length(base_metrics) > 0L) {
    base_res <- multivarious::measure_reconstruction_error(
      truth,
      estimate,
      metrics = base_metrics,
      by_column = by_column
    )
    out[names(base_res)] <- as.list(base_res[1, , drop = FALSE])
  }

  if ("mean_correlation" %in% metrics) {
    out$mean_correlation <- mean(.muscal_rowwise_correlation(truth, estimate), na.rm = TRUE)
  }
  if ("mean_cosine_similarity" %in% metrics) {
    out$mean_cosine_similarity <- mean(.muscal_rowwise_cosine(truth, estimate), na.rm = TRUE)
  }

  tibble::as_tibble(out)[metrics]
}

.muscal_normalize_retrieval_estimate <- function(estimate, n_query) {
  if (is.matrix(estimate) || is.data.frame(estimate)) {
    est_mat <- .muscal_as_numeric_matrix(estimate, "estimate")
    if (nrow(est_mat) != n_query) {
      stop(
        sprintf("Top-1 retrieval estimate has %d rows but truth has %d rows.", nrow(est_mat), n_query),
        call. = FALSE
      )
    }
    retrieved_features <- lapply(seq_len(n_query), function(i) est_mat[i, , drop = FALSE])
    return(list(retrieved_features = retrieved_features))
  }

  if (!is.list(estimate) || is.null(estimate$retrieved_features)) {
    stop(
      "`estimate` must be a numeric matrix/data.frame or a list with `retrieved_features`.",
      call. = FALSE
    )
  }

  retrieved_features <- estimate$retrieved_features
  if (!is.list(retrieved_features) || length(retrieved_features) != n_query) {
    stop(
      "`estimate$retrieved_features` must be a list with one ranked feature matrix per query.",
      call. = FALSE
    )
  }

  retrieved_features <- lapply(seq_along(retrieved_features), function(i) {
    mat <- .muscal_as_numeric_matrix(retrieved_features[[i]], sprintf("estimate$retrieved_features[[%d]]", i))
    if (nrow(mat) < 1L) {
      stop(sprintf("estimate$retrieved_features[[%d]] must contain at least one retrieved row.", i), call. = FALSE)
    }
    mat
  })

  out <- estimate
  out$retrieved_features <- retrieved_features
  out
}

.muscal_retrieval_rank <- function(retrieved_ids, truth_id, k) {
  if (is.null(retrieved_ids)) return(NA_integer_)
  pos <- match(truth_id, retrieved_ids)
  if (is.na(pos) || pos > k) return(NA_integer_)
  as.integer(pos)
}

.muscal_retrieval_metrics <- function(truth, estimate, metrics, k, truth_ids = NULL) {
  chk::chk_count(k)
  chk::chk_gte(k, 1)

  n_query <- nrow(truth)
  top1_sim <- vapply(seq_len(n_query), function(i) {
    feats <- estimate$retrieved_features[[i]]
    if (ncol(feats) != ncol(truth)) {
      stop("Retrieved feature dimensions must match `truth`.", call. = FALSE)
    }
    .muscal_rowwise_cosine(truth[i, , drop = FALSE], feats[1, , drop = FALSE])[[1]]
  }, numeric(1))

  out <- list()

  if ("mean_top1_similarity" %in% metrics) {
    out$mean_top1_similarity <- mean(top1_sim, na.rm = TRUE)
  }

  if ("mean_topk_similarity" %in% metrics) {
    out$mean_topk_similarity <- mean(vapply(seq_len(n_query), function(i) {
      feats <- estimate$retrieved_features[[i]]
      k_i <- min(k, nrow(feats))
      sims <- .muscal_rowwise_cosine(
        feats[seq_len(k_i), , drop = FALSE],
        truth[rep(i, k_i), , drop = FALSE]
      )
      mean(sims, na.rm = TRUE)
    }, numeric(1)), na.rm = TRUE)
  }

  if ("oracle_gap" %in% metrics) {
    oracle_similarity <- estimate$oracle_similarity
    if (is.null(oracle_similarity) && !is.null(estimate$oracle_features)) {
      oracle_features <- .muscal_as_numeric_matrix(estimate$oracle_features, "estimate$oracle_features")
      .muscal_require_same_dims(truth, oracle_features, "truth", "estimate$oracle_features")
      oracle_similarity <- .muscal_rowwise_cosine(truth, oracle_features)
    }
    if (is.null(oracle_similarity)) {
      stop("`oracle_gap` requires `estimate$oracle_similarity` or `estimate$oracle_features`.", call. = FALSE)
    }
    if (length(oracle_similarity) != n_query) {
      stop("`estimate$oracle_similarity` must have length equal to nrow(truth).", call. = FALSE)
    }
    out$oracle_gap <- mean(as.numeric(oracle_similarity) - top1_sim, na.rm = TRUE)
  }

  if (any(c("recall_at_k", "mrr") %in% metrics)) {
    if (is.null(truth_ids)) {
      stop("Retrieval metrics `recall_at_k` and `mrr` require `truth_ids`.", call. = FALSE)
    }
    if (length(truth_ids) != n_query) {
      stop("`truth_ids` must have length equal to nrow(truth).", call. = FALSE)
    }
    retrieved_ids <- estimate$retrieved_ids
    if (is.null(retrieved_ids) || !is.list(retrieved_ids) || length(retrieved_ids) != n_query) {
      stop("Retrieval metrics `recall_at_k` and `mrr` require `estimate$retrieved_ids` as a list per query.", call. = FALSE)
    }
    ranks <- vapply(seq_len(n_query), function(i) {
      .muscal_retrieval_rank(retrieved_ids[[i]], truth_ids[[i]], k)
    }, integer(1))

    if ("recall_at_k" %in% metrics) {
      out$recall_at_k <- mean(!is.na(ranks))
    }
    if ("mrr" %in% metrics) {
      out$mrr <- mean(ifelse(is.na(ranks), 0, 1 / ranks))
    }
  }

  tibble::as_tibble(out)[metrics]
}

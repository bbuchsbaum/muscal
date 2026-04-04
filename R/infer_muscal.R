#' Standard Resampling and Inference for Muscal Fits
#'
#' Runs bootstrap or permutation inference for `muscal` fits that expose a
#' standard refit contract through `fit_spec$refit`. This keeps resampling logic
#' generic across method families while still allowing each method to define its
#' own data reconstruction and null-generation strategy.
#'
#' The default statistics currently supported are:
#' - `"sdev"`: component-wise standard deviations / singular values
#' - `"loadings"`: the concatenated loading matrix `fit$v` (bootstrap only)
#'
#' Methods without stored refit metadata will error unless an explicit `refit`
#' specification is supplied.
#'
#' @param object A fitted `muscal` model.
#' @param method One of `"bootstrap"` or `"permutation"`.
#' @param statistic One of `"sdev"` or `"loadings"`. Ignored if
#'   `statistic_fn` is supplied.
#' @param statistic_fn Optional extractor `function(fit)`, used for custom
#'   statistics. For permutation inference this must return a numeric vector.
#' @param nrep Number of resampling replicates.
#' @param alpha Tail probability used for bootstrap intervals and permutation
#'   reference quantiles.
#' @param seed Optional RNG seed for reproducible resampling.
#' @param alternative Alternative hypothesis for permutation p-values.
#' @param refit Optional refit specification overriding `object$fit_spec$refit`.
#' @param verbose Logical; if `TRUE`, emit progress and failure counts.
#' @param ... Reserved for future extensions.
#'
#' @return A list inheriting from `muscal_inference_result` and either
#'   `muscal_bootstrap_result` or `muscal_permutation_result`.
#' @export
infer_muscal <- function(object,
                         method = c("bootstrap", "permutation"),
                         statistic = c("sdev", "loadings"),
                         statistic_fn = NULL,
                         nrep = 100,
                         alpha = 0.05,
                         seed = NULL,
                         alternative = c("greater", "less", "two.sided"),
                         refit = NULL,
                         verbose = FALSE,
                         ...) {
  method <- match.arg(method)
  alternative <- match.arg(alternative)
  chk::chk_count(nrep)
  chk::chk_number(alpha)
  chk::chk_range(alpha, c(0, 1), x_name = "alpha")
  chk::chk_flag(verbose)

  spec <- .muscal_resolve_refit_spec(object, refit = refit)
  stat_info <- .muscal_resolve_infer_statistic(statistic = statistic, statistic_fn = statistic_fn)
  observed <- stat_info$fn(object)

  if (method == "permutation" && (!is.numeric(observed) || is.matrix(observed) || is.list(observed))) {
    stop("Permutation inference currently requires a numeric vector statistic.", call. = FALSE)
  }
  if (method == "bootstrap" && stat_info$name == "loadings" && (!is.matrix(observed) || !is.numeric(observed))) {
    stop("Bootstrap loading reliability requires `fit$v`-like numeric matrix output.", call. = FALSE)
  }

  if (!is.null(seed)) {
    set.seed(seed)
  }

  replicate_stats <- vector("list", nrep)
  failures <- vector("list", 0L)

  for (i in seq_len(nrep)) {
    sampled_data <- switch(
      method,
      bootstrap = spec$bootstrap_fn(spec$data),
      permutation = spec$permutation_fn(spec$data)
    )

    fit_i <- tryCatch(
      spec$fit_fn(sampled_data),
      error = function(e) e
    )
    if (inherits(fit_i, "error")) {
      failures[[length(failures) + 1L]] <- list(index = i, stage = "fit", message = conditionMessage(fit_i))
      next
    }

    stat_i <- tryCatch(
      stat_info$fn(fit_i),
      error = function(e) e
    )
    if (inherits(stat_i, "error")) {
      failures[[length(failures) + 1L]] <- list(index = i, stage = "statistic", message = conditionMessage(stat_i))
      next
    }

    replicate_stats[[i]] <- stat_i
  }

  ok <- !vapply(replicate_stats, is.null, logical(1))
  replicate_stats <- replicate_stats[ok]
  if (length(replicate_stats) == 0L) {
    stop("All inference replicates failed.", call. = FALSE)
  }

  if (isTRUE(verbose) && length(failures) > 0L) {
    message(sprintf("infer_muscal: %d of %d replicates failed.", length(failures), nrep))
  }

  if (method == "bootstrap") {
    summary <- .muscal_summarize_bootstrap(
      fit = object,
      observed = observed,
      replicate_stats = replicate_stats,
      statistic = stat_info$name,
      alpha = alpha
    )

    return(structure(
      list(
        call = match.call(),
        method = method,
        statistic = stat_info$name,
        observed = observed,
        replicates = replicate_stats,
        summary = summary,
        alpha = alpha,
        seed = seed,
        nrep = nrep,
        n_success = length(replicate_stats),
        failures = failures
      ),
      class = c("muscal_bootstrap_result", "muscal_inference_result", "list")
    ))
  }

  perm_summary <- .muscal_summarize_permutation(
    observed = observed,
    replicate_stats = replicate_stats,
    alpha = alpha,
    alternative = alternative
  )

  structure(
    list(
      call = match.call(),
      method = method,
      statistic = stat_info$name,
      observed = observed,
      perm_values = perm_summary$perm_values,
      component_results = perm_summary$component_results,
      alpha = alpha,
      alternative = alternative,
      seed = seed,
      nrep = nrep,
      n_success = nrow(perm_summary$perm_values),
      failures = failures
    ),
    class = c("muscal_permutation_result", "muscal_inference_result", "list")
  )
}

.muscal_resolve_infer_statistic <- function(statistic, statistic_fn = NULL) {
  if (!is.null(statistic_fn)) {
    chk::chk_function(statistic_fn)
    return(list(name = "custom", fn = statistic_fn))
  }

  statistic <- match.arg(statistic, c("sdev", "loadings"))
  fn <- switch(
    statistic,
    sdev = function(fit) as.numeric(fit$sdev),
    loadings = function(fit) as.matrix(fit$v)
  )
  list(name = statistic, fn = fn)
}

.muscal_summarize_bootstrap <- function(fit,
                                        observed,
                                        replicate_stats,
                                        statistic,
                                        alpha) {
  if (statistic == "loadings") {
    return(
      loading_reliability(
        fit,
        method = "bootstrap",
        boot_loadings = replicate_stats,
        alpha = alpha,
        V_ref = observed
      )
    )
  }

  boot_mat <- do.call(rbind, lapply(replicate_stats, function(x) as.numeric(x)))
  obs <- as.numeric(observed)
  labels <- names(observed) %||% paste0("comp", seq_along(obs))

  tibble::tibble(
    component = seq_along(obs),
    label = labels,
    observed = obs,
    mean = colMeans(boot_mat),
    sd = apply(boot_mat, 2, stats::sd),
    lower = apply(boot_mat, 2, stats::quantile, probs = alpha / 2, names = FALSE),
    upper = apply(boot_mat, 2, stats::quantile, probs = 1 - alpha / 2, names = FALSE)
  )
}

.muscal_summarize_permutation <- function(observed,
                                          replicate_stats,
                                          alpha,
                                          alternative) {
  perm_mat <- do.call(rbind, lapply(replicate_stats, function(x) as.numeric(x)))
  obs <- as.numeric(observed)
  labels <- names(observed) %||% paste0("comp", seq_along(obs))

  greater <- (colSums(sweep(perm_mat, 2, obs, `>=`), na.rm = TRUE) + 1) / (nrow(perm_mat) + 1)
  less <- (colSums(sweep(perm_mat, 2, obs, `<=`), na.rm = TRUE) + 1) / (nrow(perm_mat) + 1)
  p_value <- switch(
    alternative,
    greater = greater,
    less = less,
    two.sided = pmin(1, 2 * pmin(greater, less))
  )

  list(
    perm_values = perm_mat,
    component_results = tibble::tibble(
      component = seq_along(obs),
      label = labels,
      observed = obs,
      p_value = p_value,
      lower_ci = apply(perm_mat, 2, stats::quantile, probs = alpha / 2, names = FALSE),
      upper_ci = apply(perm_mat, 2, stats::quantile, probs = 1 - alpha / 2, names = FALSE)
    )
  )
}

.muscal_bootstrap_anchored_data <- function(data) {
  Y <- as.matrix(data$Y)
  X <- data$X
  row_index <- data$row_index

  N <- nrow(Y)
  sampled <- sample.int(N, size = N, replace = TRUE)
  Y_boot <- Y[sampled, , drop = FALSE]

  X_boot <- lapply(X, function(block) block[0, , drop = FALSE])
  row_index_boot <- lapply(row_index, function(idx) integer(0))
  names(X_boot) <- names(X)
  names(row_index_boot) <- names(row_index)

  for (new_row in seq_along(sampled)) {
    old_row <- sampled[[new_row]]
    for (block_id in seq_along(X)) {
      keep <- which(row_index[[block_id]] == old_row)
      if (length(keep) == 0L) {
        next
      }
      X_boot[[block_id]] <- rbind(X_boot[[block_id]], X[[block_id]][keep, , drop = FALSE])
      row_index_boot[[block_id]] <- c(row_index_boot[[block_id]], rep.int(new_row, length(keep)))
    }
  }

  list(
    Y = Y_boot,
    X = X_boot,
    row_index = row_index_boot
  )
}

.muscal_permute_anchored_data <- function(data) {
  out <- data
  out$Y <- as.matrix(data$Y)[sample.int(nrow(data$Y)), , drop = FALSE]
  out
}

.muscal_bootstrap_aligned_data <- function(data) {
  X <- data$X
  row_index <- data$row_index
  N <- as.integer(data$N)

  sampled <- sample.int(N, size = N, replace = TRUE)

  X_boot <- lapply(X, function(block) block[0, , drop = FALSE])
  row_index_boot <- lapply(row_index, function(idx) integer(0))
  names(X_boot) <- names(X)
  names(row_index_boot) <- names(row_index)

  for (new_row in seq_along(sampled)) {
    old_row <- sampled[[new_row]]
    for (block_id in seq_along(X)) {
      keep <- which(row_index[[block_id]] == old_row)
      if (length(keep) == 0L) {
        next
      }
      X_boot[[block_id]] <- rbind(X_boot[[block_id]], X[[block_id]][keep, , drop = FALSE])
      row_index_boot[[block_id]] <- c(row_index_boot[[block_id]], rep.int(new_row, length(keep)))
    }
  }

  list(
    X = X_boot,
    row_index = row_index_boot,
    N = N
  )
}

.muscal_permute_aligned_data <- function(data) {
  out <- data
  out$row_index <- lapply(data$row_index, sample)
  out
}

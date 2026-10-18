.dmfa_bootstrap_subject_data <- function(data) {
  S <- length(data$Y)
  sampled <- sample.int(S, size = S, replace = TRUE)
  new_names <- paste0("Bootstrap_", seq_len(S))
  out <- list(
    Y = data$Y[sampled],
    design = data$design[sampled],
    trial_weights = data$trial_weights[sampled],
    block_weight = data$block_weight,
    reference_grid = data$reference_grid,
    reference_weight = data$reference_weight
  )
  if (!is.null(data$feature_map)) out$feature_map <- data$feature_map[sampled]
  names(out$Y) <- names(out$design) <- names(out$trial_weights) <- new_names
  if (!is.null(out$feature_map)) names(out$feature_map) <- new_names
  if (is.numeric(out$block_weight)) {
    out$block_weight <- out$block_weight[sampled]
    names(out$block_weight) <- new_names
  }
  attr(out, "sampled_subjects") <- names(data$Y)[sampled]
  out
}

.dmfa_refuse_implicit_permutation <- function(data) {
  stop(
    "`design_mfa` has no safe default permutation. Supply an explicit exchangeability callback to `infer_design_mfa()`.",
    call. = FALSE
  )
}

.dmfa_loading_crossproduct <- function(candidate, reference, feature_metric = NULL) {
  if (is.null(feature_metric)) return(crossprod(candidate, reference))
  crossprod(candidate, feature_metric %*% reference)
}

.dmfa_align_loading_basis <- function(candidate, reference, feature_metric = NULL) {
  candidate <- as.matrix(candidate)
  reference <- as.matrix(reference)
  K <- ncol(reference)
  if (nrow(candidate) != nrow(reference) || ncol(candidate) < K) {
    stop("Bootstrap loading basis is incompatible with the reference effect.", call. = FALSE)
  }
  candidate <- candidate[, seq_len(K), drop = FALSE]
  cross <- .dmfa_loading_crossproduct(candidate, reference, feature_metric)
  sv <- svd(cross, nu = K, nv = K)
  rotation <- sv$u %*% t(sv$v)
  candidate %*% rotation
}

.dmfa_subspace_distance <- function(candidate, reference, feature_metric = NULL) {
  K <- ncol(reference)
  if (ncol(candidate) != K || nrow(candidate) != nrow(reference)) {
    stop("Candidate and reference bases must have matching dimensions.", call. = FALSE)
  }
  overlap <- .dmfa_loading_crossproduct(candidate, reference, feature_metric)
  squared <- 1 - sum(overlap^2) / K
  sqrt(max(0, min(1, squared)))
}

.dmfa_summarize_vector_bootstrap <- function(observed, replicates, alpha, labels = NULL) {
  values <- do.call(rbind, replicates)
  if (is.null(labels)) labels <- names(observed)
  if (is.null(labels)) labels <- paste0("value_", seq_along(observed))
  data.frame(
    component = seq_along(observed),
    label = labels,
    observed = as.numeric(observed),
    mean = colMeans(values),
    sd = apply(values, 2L, stats::sd),
    lower = apply(values, 2L, stats::quantile, probs = alpha / 2, names = FALSE),
    upper = apply(values, 2L, stats::quantile, probs = 1 - alpha / 2, names = FALSE),
    row.names = NULL
  )
}

.dmfa_summarize_loading_bootstrap <- function(reference, replicates, alpha) {
  array <- abind::abind(replicates, along = 3L)
  list(
    observed = reference,
    mean = apply(array, c(1L, 2L), mean),
    sd = apply(array, c(1L, 2L), stats::sd),
    lower = apply(array, c(1L, 2L), stats::quantile, probs = alpha / 2, names = FALSE),
    upper = apply(array, c(1L, 2L), stats::quantile, probs = 1 - alpha / 2, names = FALSE),
    sign_consistency = apply(sign(array), c(1L, 2L), function(z) max(mean(z > 0), mean(z < 0)))
  )
}

.dmfa_inference_replicate <- function(i,
                                      spec,
                                      method,
                                      term,
                                      exchangeability,
                                      required_rank,
                                      needs_aligned_basis,
                                      reference_loadings,
                                      feature_metric,
                                      statistics,
                                      K) {
  sampled <- if (identical(method, "bootstrap")) {
    spec$bootstrap_fn(spec$data)
  } else {
    exchangeability(data = spec$data, term = term, replicate = i)
  }
  fit_i <- tryCatch(spec$fit_fn(sampled), error = identity)
  if (inherits(fit_i, "error")) {
    return(list(
      ok = FALSE,
      sampled_subjects = attr(sampled, "sampled_subjects"),
      failure = list(index = i, stage = "fit", message = conditionMessage(fit_i))
    ))
  }
  effect_i <- tryCatch(.dmfa_get_effect(fit_i, term), error = identity)
  if (inherits(effect_i, "error") || !isTRUE(effect_i$estimable) || effect_i$rank < required_rank) {
    message_i <- if (inherits(effect_i, "error")) {
      conditionMessage(effect_i)
    } else {
      "Effect became nonestimable or lost component rank."
    }
    return(list(
      ok = FALSE,
      sampled_subjects = attr(sampled, "sampled_subjects"),
      failure = list(index = i, stage = "effect", message = message_i)
    ))
  }
  aligned <- if (needs_aligned_basis) {
    .dmfa_align_loading_basis(
      effect_i$loadings,
      reference_loadings,
      feature_metric = feature_metric
    )
  } else {
    NULL
  }
  value_i <- list(
    partial_ss = if (is.finite(effect_i$partial_ss)) effect_i$partial_ss else effect_i$objective_delta
  )
  if ("singular_values" %in% statistics) {
    value_i$singular_values <- effect_i$singular_values[seq_len(K)]
  }
  if ("loadings" %in% statistics) value_i$loadings <- aligned
  if ("subspace" %in% statistics) {
    value_i$subspace <- .dmfa_subspace_distance(
      aligned,
      reference_loadings,
      feature_metric = feature_metric
    )
  }
  list(
    ok = TRUE,
    sampled_subjects = attr(sampled, "sampled_subjects"),
    value = value_i,
    failure = NULL
  )
}

#' Subject-Level Inference for a Design-MFA Effect
#'
#' @description
#' Resamples subjects as intact clusters for bootstrap uncertainty. Permutation
#' inference is available only when the caller supplies an explicit
#' `exchangeability` callback; there is deliberately no generic trial shuffle.
#'
#' Bootstrap loading bases are Procrustes-aligned to the fitted reference basis.
#' Nearly tied components should be interpreted through the returned subspace
#' distances rather than voxelwise loading intervals.
#'
#' @details
#' The target is conditional on the fitted reference grid, its cell weights,
#' feature maps, contrasts and chosen penalties/rank. Subject resampling keeps
#' each subject's trials, design and trial weights together and preserves the
#' fitted reference grid/weights; it does not estimate a new population grid.
#' Exchangeability of subjects or a caller-specified null transformation is a
#' scientific assumption, not established by callback schema validation. No
#' unconditional rejection-rate or interval-coverage calibration is claimed.
#'
#' Bootstrap percentile intervals on `partial_ss` (and other nonnegative
#' statistics) have a positive null expectation and must not be read as a
#' hypothesis test. Use permutation with an explicit `exchangeability`
#' callback to test a term.
#'
#' Permutation p-values use successful replicates only:
#' `(1 + b) / (n_success + 1)`. Failed refits are not treated as null draws.
#' If `n_success < nrep`, the resolution is `1 / (n_success + 1)` and a
#' warning is issued.
#'
#' When a multifer build with `effect_significance` is available, serial
#' permutation calls `multifer::infer()`. When the statistic-stability
#' sidecar is also available, bootstrap calls `multifer::bootstrap_fits()`
#' and maps percentile intervals from that sidecar. `use_future = TRUE`
#' no longer starts a local `furrr` loop on those paths.
#'
#' @param x A fitted [design_mfa()] object.
#' @param term One estimable design term.
#' @param method `"bootstrap"` or `"permutation"`.
#' @param statistics Summaries to compute. The method-dependent default is
#'   `c("partial_ss", "singular_values", "subspace")` for the bootstrap and
#'   `c("partial_ss", "singular_values")` for permutation inference. Loading
#'   intervals are optional because they can be large.
#' @param nrep Number of resampling replicates.
#' @param alpha Total error probability for central percentile intervals.
#' @param seed Optional random seed. The caller's global random-number state is
#'   restored on successful and error exits.
#' @param exchangeability For permutation, a function accepting `data`, `term`,
#'   and `replicate` and returning a complete refit-data list. The callback owns
#'   the scientific validity of the null transformation. Return aligned named
#'   `Y`, `design`, `trial_weights`, `block_weight`, `reference_grid`, and
#'   `reference_weight` fields (plus `feature_map` for soft alignment). Preserve
#'   the frozen reference grid and weights. Malformed refits are failures.
#' @param verbose Emit per-replicate progress. Failed replicates always warn.
#' @param use_future If `TRUE` and the multifer path is unavailable, run
#'   independent native refits with [furrr::future_map()]. Multifer
#'   permutation and bootstrap stay sequential.
#' @param ... Unused.
#'
#' @return A `design_mfa_inference` object containing observed statistics,
#'   replicate values, failures, and bootstrap summaries or a permutation
#'   p-value.
#'
#' @examples
#' \donttest{
#' example_data <- design_mfa_example_data()
#' example_fit <- design_mfa(example_data$Y, example_data$design, ~ saliency)
#' uncertainty <- infer_design_mfa(
#'   example_fit, "saliency", statistics = "partial_ss",
#'   nrep = 5, seed = 10
#' )
#' uncertainty$summary$partial_ss
#' }
#' @export
infer_design_mfa <- function(x,
                             term,
                             method = c("bootstrap", "permutation"),
                             statistics = NULL,
                             nrep = 200L,
                             alpha = 0.05,
                             seed = NULL,
                             exchangeability = NULL,
                             verbose = FALSE,
                             use_future = FALSE,
                             ...) {
  if (!inherits(x, "design_mfa")) stop("`x` must be a fitted `design_mfa` object.", call. = FALSE)
  effect <- .dmfa_get_effect(x, term)
  if (!effect$estimable) stop(sprintf("Term '%s' is not estimable.", term), call. = FALSE)
  method <- match.arg(method)
  if (is.null(statistics)) {
    statistics <- if (identical(method, "permutation")) {
      c("partial_ss", "singular_values")
    } else {
      c("partial_ss", "singular_values", "subspace")
    }
  }
  statistics <- unique(match.arg(
    statistics,
    c("partial_ss", "singular_values", "subspace", "loadings"),
    several.ok = TRUE
  ))
  nrep <- as.integer(nrep)
  if (length(nrep) != 1L || is.na(nrep) || nrep < 2L) stop("`nrep` must be at least two.", call. = FALSE)
  if (length(alpha) != 1L || !is.finite(alpha) || alpha <= 0 || alpha >= 1) {
    stop("`alpha` must lie strictly between zero and one.", call. = FALSE)
  }
  if (length(verbose) != 1L || is.na(verbose)) stop("`verbose` must be TRUE or FALSE.", call. = FALSE)
  if (length(use_future) != 1L || is.na(use_future)) {
    stop("`use_future` must be TRUE or FALSE.", call. = FALSE)
  }
  use_future <- isTRUE(use_future)
  if (identical(method, "permutation")) {
    if (!is.function(exchangeability)) {
      stop("Permutation inference requires an explicit `exchangeability` callback.", call. = FALSE)
    }
    if (any(statistics %in% c("subspace", "loadings"))) {
      stop("Permutation inference supports scalar/root statistics, not loading or subspace summaries.", call. = FALSE)
    }
  }

  spec <- .muscal_resolve_refit_spec(x)
  if (!identical(spec$resample_unit, "subjects")) {
    stop("The fitted refit contract is not subject-clustered.", call. = FALSE)
  }
  the_call <- match.call()
  use_multifer_permutation <- identical(method, "permutation") &&
    .dmfa_multifer_effect_contract_available()
  use_multifer_bootstrap <- identical(method, "bootstrap") &&
    .dmfa_multifer_stability_contract_available()
  if (!is.null(seed)) {
    if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed) ||
        abs(seed - round(seed)) > sqrt(.Machine$double.eps)) {
      stop("`seed` must be NULL or one integer.", call. = FALSE)
    }
    seed <- as.integer(round(seed))
    had_random_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    if (had_random_seed) {
      random_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    }
    on.exit({
      if (had_random_seed) {
        assign(".Random.seed", random_seed, envir = .GlobalEnv)
      } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        rm(".Random.seed", envir = .GlobalEnv)
      }
    }, add = TRUE)
    set.seed(seed)
  }

  if (isTRUE(use_multifer_permutation)) {
    out <- .dmfa_infer_permutation_via_multifer(
      x = x,
      term = term,
      statistics = statistics,
      nrep = nrep,
      alpha = alpha,
      seed = seed,
      exchangeability = exchangeability,
      spec = spec,
      effect = effect
    )
    out$call <- the_call
    return(out)
  }
  if (isTRUE(use_multifer_bootstrap)) {
    out <- .dmfa_infer_bootstrap_via_multifer(
      x = x,
      term = term,
      statistics = statistics,
      nrep = nrep,
      alpha = alpha,
      seed = seed,
      spec = spec,
      effect = effect,
      verbose = verbose
    )
    out$call <- the_call
    return(out)
  }

  reference_loadings <- effect$loadings
  K <- ncol(reference_loadings)
  required_rank <- if (any(statistics %in% c("singular_values", "subspace", "loadings"))) K else 0L
  needs_aligned_basis <- any(statistics %in% c("subspace", "loadings"))
  replicate_args <- list(
    spec = spec,
    method = method,
    term = term,
    exchangeability = exchangeability,
    required_rank = required_rank,
    needs_aligned_basis = needs_aligned_basis,
    reference_loadings = reference_loadings,
    feature_metric = x$feature_metric,
    statistics = statistics,
    K = K
  )
  run_one <- function(i) {
    if (isTRUE(verbose) && !isTRUE(use_future)) {
      message(sprintf("infer_design_mfa: replicate %d / %d", i, nrep))
    }
    do.call(.dmfa_inference_replicate, c(list(i = i), replicate_args))
  }
  if (isTRUE(verbose) && isTRUE(use_future)) {
    message(sprintf("infer_design_mfa: running %d replicates via future", nrep))
  }
  results <- if (isTRUE(use_future)) {
    furrr::future_map(seq_len(nrep), run_one, .options = furrr::furrr_options(seed = TRUE))
  } else {
    lapply(seq_len(nrep), run_one)
  }

  ok <- vapply(results, `[[`, logical(1), "ok")
  replicate_values <- lapply(results[ok], `[[`, "value")
  sampled_subjects <- lapply(results[ok], `[[`, "sampled_subjects")
  failures <- lapply(results[!ok], `[[`, "failure")
  if (length(replicate_values) < 2L) stop("Fewer than two inference replicates succeeded.", call. = FALSE)
  if (length(failures) > 0L) {
    warning(
      sprintf(
        "infer_design_mfa: %d of %d replicates failed; summaries use %d successful replicates.",
        length(failures), nrep, length(replicate_values)
      ),
      call. = FALSE
    )
  }

  observed <- list(
    partial_ss = if (is.finite(effect$partial_ss)) effect$partial_ss else effect$objective_delta,
    singular_values = effect$singular_values[seq_len(K)],
    loadings = reference_loadings,
    subspace = 0
  )

  if (identical(method, "permutation")) {
    summary <- list()
    if ("partial_ss" %in% statistics) {
      null <- vapply(replicate_values, `[[`, numeric(1), "partial_ss")
      summary$partial_ss <- data.frame(
        observed = observed$partial_ss,
        p_value = (1 + sum(null >= observed$partial_ss)) / (length(null) + 1),
        null_mean = mean(null),
        null_sd = stats::sd(null)
      )
    }
    if ("singular_values" %in% statistics) {
      null <- do.call(rbind, lapply(replicate_values, `[[`, "singular_values"))
      summary$singular_values <- data.frame(
        component = seq_len(K),
        observed = observed$singular_values,
        p_value = (1 + colSums(sweep(null, 2L, observed$singular_values, `>=`))) / (nrow(null) + 1)
      )
    }
  } else {
    summary <- list()
    if ("partial_ss" %in% statistics) {
      values <- lapply(replicate_values, function(z) c(partial_ss = z$partial_ss))
      summary$partial_ss <- .dmfa_summarize_vector_bootstrap(
        observed = c(partial_ss = observed$partial_ss),
        replicates = values,
        alpha = alpha,
        labels = "partial_ss"
      )
    }
    if ("singular_values" %in% statistics) {
      values <- lapply(replicate_values, `[[`, "singular_values")
      summary$singular_values <- .dmfa_summarize_vector_bootstrap(
        observed = observed$singular_values,
        replicates = values,
        alpha = alpha,
        labels = paste0("Component_", seq_len(K))
      )
    }
    if ("subspace" %in% statistics) {
      values <- vapply(replicate_values, `[[`, numeric(1), "subspace")
      summary$subspace <- data.frame(
        mean_distance = mean(values),
        sd_distance = stats::sd(values),
        lower = stats::quantile(values, alpha / 2, names = FALSE),
        upper = stats::quantile(values, 1 - alpha / 2, names = FALSE)
      )
    }
    if ("loadings" %in% statistics) {
      values <- lapply(replicate_values, `[[`, "loadings")
      summary$loadings <- .dmfa_summarize_loading_bootstrap(reference_loadings, values, alpha)
    }
  }

  structure(
    list(
      call = match.call(),
      method = method,
      term = term,
      statistics = statistics,
      observed = observed,
      replicates = replicate_values,
      sampled_subjects = sampled_subjects,
      summary = summary,
      alpha = alpha,
      seed = seed,
      nrep = nrep,
      n_success = length(replicate_values),
      failures = failures,
      validity = if (identical(method, "bootstrap")) "subject_cluster_bootstrap" else "user_supplied_exchangeability",
      engine = "native"
    ),
    class = c("design_mfa_inference", "muscal_inference_result", "list")
  )
}

#' @export
print.design_mfa_inference <- function(x, ...) {
  cat(sprintf("Design-effect %s inference for term '%s'\n", x$method, x$term))
  cat(sprintf("  Successful replicates: %d / %d\n", x$n_success, x$nrep))
  cat(sprintf("  Validity: %s\n", x$validity))
  invisible(x)
}

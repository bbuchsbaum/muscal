#' @export
as_multifer_adapter.design_mfa <- function(x, ...) {
  .design_mfa_multifer_adapter()
}

.design_mfa_data_schema <- function() {
  multifer::data_role_schema(
    Y              = multifer::role("primary", axes = c("row", "col"),
                                    container = "list"),
    design         = multifer::role("design_index", axes = "row",
                                    container = "list"),
    trial_weights  = multifer::role("metric", axes = "row",
                                    container = "list",
                                    policy = "diagonal_observation_weights"),
    feature_map    = multifer::role("constraint", axes = c("row", "col"),
                                    container = "list", optional = TRUE),
    reference_grid = multifer::role("static"),
    block_weight   = multifer::role("metric", container = "list",
                                    policy = "diagonal_observation_weights")
  )
}

.design_mfa_estimable_terms <- function(x) {
  names(x$effects)[vapply(
    x$effects,
    function(effect) isTRUE(effect$estimable),
    logical(1L)
  )]
}

.design_mfa_effect_stat <- function(x, data, effect, spec = NULL,
                                    collector = NULL,
                                    statistics = "partial_ss",
                                    required_rank = 0L,
                                    K = 1L) {
  spec <- spec %||% .muscal_resolve_refit_spec(x)
  # An identity permutation is still a null draw, even when its data match
  # the original fit. Remove the internal marker before invoking the refit.
  null_draw <- !is.null(attr(data, ".dmfa_null_replicate"))
  attr(data, ".dmfa_null_replicate") <- NULL
  original <- !null_draw && inherits(x, "design_mfa") && identical(data, spec$data)
  fit <- if (original) x else spec$fit_fn(data)
  if (!inherits(fit, "design_mfa")) {
    stop("Design-MFA refit did not return a design_mfa object.", call. = FALSE)
  }
  eff <- .dmfa_get_effect(fit, effect)
  if (!isTRUE(eff$estimable) ||
      (required_rank > 0L && (eff$rank %||% 0L) < required_rank)) {
    stop("Effect became nonestimable or lost component rank.", call. = FALSE)
  }
  stat <- if (is.finite(eff$partial_ss)) eff$partial_ss else eff$objective_delta
  if (!is.numeric(stat) || length(stat) != 1L || !is.finite(stat)) {
    stop("Effect statistic is not a finite scalar.", call. = FALSE)
  }
  if (!original && !is.null(collector)) {
    value <- list(partial_ss = as.numeric(stat))
    if ("singular_values" %in% statistics) {
      value$singular_values <- eff$singular_values[seq_len(K)]
    }
    collector$replicates[[length(collector$replicates) + 1L]] <- value
    collector$sampled_subjects[length(collector$sampled_subjects) + 1L] <-
      list(attr(data, "sampled_subjects"))
  }
  as.numeric(stat)
}

.design_mfa_named_statistics <- function(x, effect,
                                         statistics = "partial_ss",
                                         required_rank = 0L,
                                         K = 1L) {
  eff <- .dmfa_get_effect(x, effect)
  if (!isTRUE(eff$estimable) ||
      (required_rank > 0L && (eff$rank %||% 0L) < required_rank)) {
    stop("Effect became nonestimable or lost component rank.", call. = FALSE)
  }
  out <- list()
  if ("partial_ss" %in% statistics) {
    stat <- if (is.finite(eff$partial_ss)) eff$partial_ss else eff$objective_delta
    if (!is.numeric(stat) || length(stat) != 1L || !is.finite(stat)) {
      stop("Effect statistic is not a finite scalar.", call. = FALSE)
    }
    out$partial_ss <- as.numeric(stat)
  }
  if ("singular_values" %in% statistics) {
    sv <- eff$singular_values[seq_len(K)]
    if (!is.numeric(sv) || length(sv) != K || any(!is.finite(sv))) {
      stop("Singular values are not a finite vector of the required rank.",
           call. = FALSE)
    }
    out$singular_values <- as.numeric(sv)
  }
  out
}

.dmfa_ensure_unique_subject_names <- function(data) {
  if (is.null(data$Y) || !is.list(data$Y) || is.data.frame(data$Y)) {
    return(data)
  }
  nms <- names(data$Y)
  if (!is.null(nms) && !anyNA(nms) && all(nzchar(nms)) && !anyDuplicated(nms)) {
    return(data)
  }
  new_names <- paste0("Bootstrap_", seq_len(length(data$Y)))
  rename_list <- function(x) {
    if (is.null(x) || !is.list(x) || is.data.frame(x)) {
      return(x)
    }
    names(x) <- new_names
    x
  }
  data$Y <- rename_list(data$Y)
  data$design <- rename_list(data$design)
  data$trial_weights <- rename_list(data$trial_weights)
  if (!is.null(data$feature_map)) {
    data$feature_map <- rename_list(data$feature_map)
  }
  if (is.numeric(data$block_weight) && length(data$block_weight) == length(new_names)) {
    names(data$block_weight) <- new_names
  }
  data
}

.design_mfa_multifer_adapter <- function(spec = NULL,
                                         collector = NULL,
                                         statistics = "partial_ss",
                                         required_rank = 0L,
                                         K = 1L,
                                         term = NULL) {
  .dmfa_require_multifer_effect_contract()

  claimed <- "effect_significance"
  if (.dmfa_multifer_stability_contract_available()) {
    claimed <- c(claimed, "statistic_stability")
  }

  multifer::infer_adapter(
    adapter_id = "muscal_design_mfa",
    adapter_version = "0.1.0",
    shape_kinds = "adapter",
    capabilities = multifer::capability_matrix(
      list(
        geometry = "adapter",
        relation = "design_contrast",
        targets = claimed
      )
    ),
    domains = function(x, data = NULL, ...) {
      if (!is.null(term)) {
        return(term)
      }
      .design_mfa_estimable_terms(x)
    },
    roots = function(x, ...) {
      terms <- .design_mfa_estimable_terms(x)
      if (length(terms) == 0L) {
        return(numeric(0))
      }
      as.numeric(x$effects[[terms[[1L]]]]$singular_values)
    },
    loadings = function(x, domain = NULL, ...) {
      terms <- .design_mfa_estimable_terms(x)
      if (length(terms) == 0L) {
        stop("No estimable design effect has loadings.", call. = FALSE)
      }
      term <- domain %||% terms[[1L]]
      .dmfa_get_effect(x, term)$loadings
    },
    scores = function(x, domain = NULL, ...) {
      terms <- .design_mfa_estimable_terms(x)
      if (length(terms) == 0L) {
        stop("No estimable design effect has scores.", call. = FALSE)
      }
      term <- domain %||% terms[[1L]]
      .dmfa_get_effect(x, term)$scores
    },
    effects = function(x, data = NULL, ...) {
      out <- .design_mfa_estimable_terms(x)
      if (length(out) == 0L) {
        stop("No estimable design effects are available.", call. = FALSE)
      }
      out
    },
    effect_stat = function(x, data, effect, ...) {
      .design_mfa_effect_stat(
        x, data, effect,
        spec = spec,
        collector = collector,
        statistics = statistics,
        required_rank = required_rank,
        K = K
      )
    },
    statistic = function(x, data, effect, ...) {
      .design_mfa_named_statistics(
        x,
        effect = effect,
        statistics = statistics,
        required_rank = required_rank,
        K = K
      )
    },
    refit = function(x, new_data, ...) {
      spec_i <- spec %||% .muscal_resolve_refit_spec(x)
      spec_i$fit_fn(.dmfa_ensure_unique_subject_names(new_data))
    },
    validity_level = "conditional",
    declared_assumptions = c(
      "subjects_are_sampling_units",
      "caller_owned_exchangeability",
      "no_default_permutation"
    ),
    data_schema = .design_mfa_data_schema(),
    null_spec = multifer::null_spec(
      kind = "caller_owned",
      preserves = c("subject_blocks", "reference_grid", "trial_weights")
    )
  )
}

.dmfa_require_multifer_effect_contract <- function(caller = "as_multifer_adapter.design_mfa()") {
  if (!.dmfa_multifer_effect_contract_available()) {
    stop(
      paste0(
        "`", caller, "` requires a multifer build with `effect_significance`, ",
        "`clustered_units()`, and list-contained roles. Install the updated ",
        "multifer from ~/code/multifer or bbuchsbaum/multifer."
      ),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

.dmfa_multifer_effect_contract_available <- function() {
  if (!requireNamespace("multifer", quietly = TRUE)) {
    return(FALSE)
  }
  ns <- asNamespace("multifer")
  needed <- c(
    "infer", "infer_adapter", "capability_matrix", "clustered_units",
    "role", "data_role_schema", "null_spec", "valid_capability_targets"
  )
  if (!all(vapply(needed, exists, logical(1L), where = ns, inherits = FALSE))) {
    return(FALSE)
  }
  targets_ok <- tryCatch(
    "effect_significance" %in% multifer::valid_capability_targets(),
    error = function(e) FALSE
  )
  if (!isTRUE(targets_ok)) {
    return(FALSE)
  }
  role_formals <- names(formals(multifer::role))
  "container" %in% role_formals && "optional" %in% role_formals
}

.dmfa_multifer_stability_contract_available <- function() {
  if (!.dmfa_multifer_effect_contract_available()) {
    return(FALSE)
  }
  ns <- asNamespace("multifer")
  needed <- c(
    "statistic_stability_from_bootstrap", "statistic_stability",
    "bootstrap_fits", "form_effect_units"
  )
  if (!all(vapply(needed, exists, logical(1L), where = ns, inherits = FALSE))) {
    return(FALSE)
  }
  targets_ok <- tryCatch(
    "statistic_stability" %in% multifer::valid_capability_targets(),
    error = function(e) FALSE
  )
  if (!isTRUE(targets_ok)) {
    return(FALSE)
  }
  "on_error" %in% names(formals(multifer::bootstrap_fits))
}

.dmfa_translate_exchangeability <- function(exchangeability, term) {
  function(x, data, effect = term, replicate = NULL, ...) {
    out <- exchangeability(
      data = data,
      term = effect,
      replicate = if (is.null(replicate)) 0L else as.integer(replicate)
    )
    attr(out, ".dmfa_null_replicate") <- replicate %||% 0L
    out
  }
}

.dmfa_infer_permutation_via_multifer <- function(x,
                                                 term,
                                                 statistics,
                                                 nrep,
                                                 alpha,
                                                 seed,
                                                 exchangeability,
                                                 spec,
                                                 effect) {
  K <- ncol(effect$loadings)
  required_rank <- if ("singular_values" %in% statistics) K else 0L
  collector <- new.env(parent = emptyenv())
  collector$replicates <- list()
  collector$sampled_subjects <- list()
  adapter <- .design_mfa_multifer_adapter(
    spec = spec,
    collector = collector,
    statistics = statistics,
    required_rank = required_rank,
    K = K
  )
  null_action <- .dmfa_translate_exchangeability(exchangeability, term)

  bundle <- withCallingHandlers(
    multifer::infer(
      adapter = adapter,
      data = spec$data,
      geometry = "adapter",
      relation = "design_contrast",
      design = multifer::clustered_units(unit_role = "Y"),
      targets = "effect_significance",
      effects = term,
      null_action = null_action,
      model = x,
      B = nrep,
      alpha = alpha,
      seed = NULL,
      return_bundle = TRUE
    ),
    warning = function(w) {
      if (grepl("null draws failed", conditionMessage(w), fixed = TRUE)) {
        tryInvokeRestart("muffleWarning")
      }
    }
  )

  step <- bundle$artifacts$engine$ladder_result$step_results[[1L]]
  n_success <- as.integer(step$n_success %||% length(step$null_values))
  n_failed <- as.integer(step$n_failed %||% 0L)
  if (n_success < 2L) {
    stop("Fewer than two inference replicates succeeded.", call. = FALSE)
  }
  if (n_failed > 0L) {
    warning(
      sprintf(
        "infer_design_mfa: %d of %d replicates failed; summaries use %d successful replicates.",
        n_failed, nrep, n_success
      ),
      call. = FALSE
    )
  }

  replicate_values <- collector$replicates
  if (length(replicate_values) != n_success) {
    stop("Multifer returned inconsistent effect replicate statistics.", call. = FALSE)
  }
  sampled_subjects <- collector$sampled_subjects
  if (length(sampled_subjects) != length(replicate_values)) {
    sampled_subjects <- vector("list", length(replicate_values))
  }

  observed <- list(
    partial_ss = if (is.finite(effect$partial_ss)) effect$partial_ss else effect$objective_delta,
    singular_values = effect$singular_values[seq_len(K)],
    loadings = effect$loadings,
    subspace = 0
  )
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

  failures <- if (n_failed > 0L) {
    replicate(n_failed, list(index = NA_integer_, stage = "fit",
                             message = "null draw failed"), simplify = FALSE)
  } else {
    list()
  }

  structure(
    list(
      call = NULL,
      method = "permutation",
      term = term,
      statistics = statistics,
      observed = observed,
      replicates = replicate_values,
      sampled_subjects = sampled_subjects,
      summary = summary,
      alpha = alpha,
      seed = seed,
      nrep = nrep,
      n_success = n_success,
      failures = failures,
      validity = "user_supplied_exchangeability",
      engine = "multifer"
    ),
    class = c("design_mfa_inference", "muscal_inference_result", "list")
  )
}

.dmfa_bootstrap_value_from_fit <- function(fit,
                                           term,
                                           statistics,
                                           K,
                                           required_rank,
                                           reference_loadings,
                                           feature_metric) {
  effect_i <- .dmfa_get_effect(fit, term)
  if (!isTRUE(effect_i$estimable) || (effect_i$rank %||% 0L) < required_rank) {
    stop("Effect became nonestimable or lost component rank.", call. = FALSE)
  }
  needs_aligned_basis <- any(statistics %in% c("subspace", "loadings"))
  aligned <- if (needs_aligned_basis) {
    .dmfa_align_loading_basis(
      effect_i$loadings,
      reference_loadings,
      feature_metric = feature_metric
    )
  } else {
    NULL
  }
  value <- list(
    partial_ss = if (is.finite(effect_i$partial_ss)) {
      effect_i$partial_ss
    } else {
      effect_i$objective_delta
    }
  )
  if ("singular_values" %in% statistics) {
    value$singular_values <- effect_i$singular_values[seq_len(K)]
  }
  if ("loadings" %in% statistics) value$loadings <- aligned
  if ("subspace" %in% statistics) {
    value$subspace <- .dmfa_subspace_distance(
      aligned,
      reference_loadings,
      feature_metric = feature_metric
    )
  }
  value
}

.dmfa_sidecar_vector_summary <- function(stab, statistic, labels = NULL) {
  rows <- stab[stab$statistic == statistic, , drop = FALSE]
  if (nrow(rows) == 0L) {
    stop(sprintf("statistic_stability sidecar is missing '%s'.", statistic),
         call. = FALSE)
  }
  if (!is.null(labels) && length(labels) == nrow(rows)) {
    rows$label <- labels
  }
  data.frame(
    component = seq_len(nrow(rows)),
    label = rows$label,
    observed = rows$observed,
    mean = rows$mean,
    sd = rows$sd,
    lower = rows$lower,
    upper = rows$upper,
    row.names = NULL
  )
}

.dmfa_infer_bootstrap_via_multifer <- function(x,
                                               term,
                                               statistics,
                                               nrep,
                                               alpha,
                                               seed,
                                               spec,
                                               effect,
                                               verbose = FALSE) {
  K <- ncol(effect$loadings)
  required_rank <- if (any(statistics %in% c("singular_values", "subspace", "loadings"))) {
    K
  } else {
    0L
  }
  named_stats <- intersect(statistics, c("partial_ss", "singular_values"))
  adapter <- .design_mfa_multifer_adapter(
    spec = spec,
    statistics = named_stats,
    required_rank = required_rank,
    K = K,
    term = term
  )
  units <- multifer::form_effect_units(term)
  recipe <- multifer::infer_recipe(
    geometry = "adapter",
    relation = "design_contrast",
    design = multifer::clustered_units(unit_role = "Y"),
    targets = "statistic_stability",
    effects = term,
    adapter = adapter
  )
  if (isTRUE(verbose)) {
    message(sprintf("infer_design_mfa: running %d replicates via multifer", nrep))
  }
  artifact <- multifer::bootstrap_fits(
    recipe = recipe,
    adapter = adapter,
    data = spec$data,
    original_fit = x,
    units = units,
    R = nrep,
    method_align = "sign",
    seed = NULL,
    parallel = "sequential",
    store_aligned_scores = FALSE,
    on_error = "drop"
  )

  subject_ids <- names(spec$data$Y)
  if (is.null(subject_ids)) {
    subject_ids <- as.character(seq_along(spec$data$Y))
  }
  kept <- list()
  failures <- list()
  for (i in seq_along(artifact$reps)) {
    rep <- artifact$reps[[i]]
    if (isFALSE(rep$ok) || is.null(rep$fit)) {
      failures[[length(failures) + 1L]] <- list(
        index = i,
        stage = "fit",
        message = rep$error %||% "bootstrap replicate failed"
      )
      next
    }
    value <- tryCatch(
      .dmfa_bootstrap_value_from_fit(
        fit = rep$fit,
        term = term,
        statistics = statistics,
        K = K,
        required_rank = required_rank,
        reference_loadings = effect$loadings,
        feature_metric = x$feature_metric
      ),
      error = identity
    )
    if (inherits(value, "error")) {
      failures[[length(failures) + 1L]] <- list(
        index = i,
        stage = "effect",
        message = conditionMessage(value)
      )
      next
    }
    kept[[length(kept) + 1L]] <- list(
      fit = rep$fit,
      value = value,
      sampled_subjects = subject_ids[rep$resample_indices]
    )
  }

  n_success <- length(kept)
  if (n_success < 2L) {
    stop("Fewer than two inference replicates succeeded.", call. = FALSE)
  }
  if (length(failures) > 0L) {
    warning(
      sprintf(
        "infer_design_mfa: %d of %d replicates failed; summaries use %d successful replicates.",
        length(failures), nrep, n_success
      ),
      call. = FALSE
    )
  }

  observed <- list(
    partial_ss = if (is.finite(effect$partial_ss)) effect$partial_ss else effect$objective_delta,
    singular_values = effect$singular_values[seq_len(K)],
    loadings = effect$loadings,
    subspace = 0
  )
  summary <- list()
  if (length(named_stats) > 0L) {
    observed_named <- observed[named_stats]
    names(observed_named) <- named_stats
    filtered <- structure(
      list(
        reps = lapply(kept, function(z) list(fit = z$fit, ok = TRUE)),
        R = n_success
      ),
      class = "multifer_bootstrap_artifact"
    )
    stab <- multifer::statistic_stability_from_bootstrap(
      artifact = filtered,
      units = units,
      statistics = observed_named,
      nonnegative = TRUE,
      alpha = alpha,
      statistic = function(fit) {
        .design_mfa_named_statistics(
          fit,
          effect = term,
          statistics = named_stats,
          required_rank = required_rank,
          K = K
        )
      }
    )
    if ("partial_ss" %in% named_stats) {
      summary$partial_ss <- .dmfa_sidecar_vector_summary(
        stab, "partial_ss", labels = "partial_ss"
      )
    }
    if ("singular_values" %in% named_stats) {
      summary$singular_values <- .dmfa_sidecar_vector_summary(
        stab,
        "singular_values",
        labels = paste0("Component_", seq_len(K))
      )
    }
  }
  replicate_values <- lapply(kept, `[[`, "value")
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
    summary$loadings <- .dmfa_summarize_loading_bootstrap(
      effect$loadings, values, alpha
    )
  }

  structure(
    list(
      call = NULL,
      method = "bootstrap",
      term = term,
      statistics = statistics,
      observed = observed,
      replicates = replicate_values,
      sampled_subjects = lapply(kept, `[[`, "sampled_subjects"),
      summary = summary,
      alpha = alpha,
      seed = seed,
      nrep = nrep,
      n_success = n_success,
      failures = failures,
      validity = "subject_cluster_bootstrap",
      engine = "multifer"
    ),
    class = c("design_mfa_inference", "muscal_inference_result", "list")
  )
}

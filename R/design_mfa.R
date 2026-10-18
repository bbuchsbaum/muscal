#' Construct a Small Deterministic Design-MFA Example
#'
#' @description
#' Returns a balanced six-subject trialwise field dataset used by the runnable
#' examples for the Design-MFA API. The fixture is deterministic and small
#' enough for documentation and teaching.
#'
#' @return A list with named `Y` and `design` subject blocks and
#'   `feature_names`.
#'
#' @examples
#' example_data <- design_mfa_example_data()
#' vapply(example_data$Y, nrow, integer(1))
#' @export
design_mfa_example_data <- function() {
  subjects <- paste0("S", seq_len(6L))
  saliency_levels <- c(20, 40, 60, 80, 100)
  probe_levels <- c("old", "new", "lure")
  feature_names <- paste0("voxel_", seq_len(8L))
  design <- lapply(seq_along(subjects), function(s) {
    grid <- expand.grid(
      saliency = saliency_levels,
      probe = probe_levels,
      repetition = seq_len(2L),
      KEEP.OUT.ATTRS = FALSE,
      stringsAsFactors = FALSE
    )
    data.frame(
      saliency = ordered(grid$saliency, levels = saliency_levels),
      probe = factor(grid$probe, levels = probe_levels)
    )
  })
  saliency_pattern <- seq(-1, 1, length.out = length(feature_names))
  probe_pattern <- cos(seq(0, 2 * pi, length.out = length(feature_names)))
  Y <- lapply(seq_along(design), function(s) {
    d <- design[[s]]
    saliency <- as.numeric(d$saliency) - 3
    probe <- c(old = -1, new = 0, lure = 1)[as.character(d$probe)]
    deterministic_noise <- 0.03 * sin(
      outer(seq_len(nrow(d)), seq_len(length(feature_names)), "+") + s
    )
    out <- outer(saliency, saliency_pattern) +
      0.5 * outer(probe, probe_pattern) +
      matrix(0.05 * s * saliency_pattern, nrow(d), length(feature_names), byrow = TRUE) +
      deterministic_noise
    colnames(out) <- feature_names
    out
  })
  names(Y) <- names(design) <- subjects
  list(Y = Y, design = design, feature_names = feature_names)
}

#' Effect-Specific Multilevel Components for Trialwise Fields
#'
#' @description
#' `design_mfa()` fits a trial-level multivariate design model and decomposes
#' each estimable design effect into low-rank spatial components. Subjects are
#' represented by fixed intercepts, with optional ridge-pooled term deviations,
#' and trial weights are normalized within subject by default so unequal trial
#' counts do not silently change population weighting.
#'
#' @details
#' For every requested term, the method compares a containing weighted model
#' with a reduced model, evaluates their prediction difference on an explicit
#' reference grid, and performs a weighted thin SVD of that effect matrix.
#' Thus, the components maximize fitted variance associated with a named design
#' term rather than total variance in the activation fields.
#'
#' Partial sums of squares need not add across terms in an unbalanced design.
#' When `ss_type = "sequential"`, results are explicitly order-dependent.
#' Ordered factors should use planned contrasts, such as `"poly"`. Effects of
#' observed response variables such as confidence are conditional associations,
#' not causal experimental effects.
#'
#' Shared alignment supports either fixed subject intercepts or centered
#' ridge-pooled subject slopes. Soft alignment uses explicit sparse maps from a
#' reference field into each subject's feature space; it currently supports the
#' fixed-intercept subject model.
#'
#' The exported function name is retained for API continuity, but the
#' statistical method is an ASCA+/partial-RDA-style design decomposition. It
#' does not perform the block normalization that defines multiple factor
#' analysis (MFA).
#'
#' Fitted objects retain the training responses in `fit_spec$refit$data` so
#' subject-clustered refits do not need the caller to pass data again.
#' `saveRDS(fit)` therefore writes that dataset. Peak working memory during
#' the fit is approximately `O(p q)` plus the data already in memory.
#'
#' @param Y A list of subject matrices, or a `multidesign` object. Rows are
#'   trials. Columns are shared features under shared alignment and may be
#'   subject-specific under soft alignment.
#' @param design A list of data frames with one row per corresponding trial.
#'   Omit this when `Y` is a `multidesign` object.
#' @param formula A one-sided formula defining design effects.
#' @param subject Optional unique subject identifiers, one per list element.
#' @param nuisance Optional one-sided formula defining nuisance columns. Terms
#'   must not also appear in `formula`.
#' @param contrasts Optional named list of factor contrast specifications. The
#'   aliases `"poly"`, `"sum"`, `"treatment"`, and `"helmert"` are supported.
#' @param terms Optional character vector selecting formula terms to return.
#' @param reference_grid Optional data frame on which effects are evaluated.
#'   The default is the Cartesian product of declared factor levels and modest
#'   numeric grids.
#' @param reference_weight Either `"equal_cells"`, `"observed"`, or a
#'   nonnegative numeric vector with one value per reference-grid row.
#' @param trial_weights Optional list of nonnegative trial-precision vectors.
#' @param block_weight Population weighting: `"equal_subject"`, `"precision"`,
#'   or a numeric subject-weight vector.
#' @param feature_names Optional shared feature schema, required when `Y` has no
#'   column names. Under soft alignment these name the reference features.
#' @param feature_map For soft alignment, a list of sparse `p_s x p_0` maps
#'   taking reference fields into each subject's feature space.
#' @param rank Optional retained component rank, either a scalar or values named
#'   by term. By default all numerically nonzero effect components are retained.
#' @param ss_type `"partial"` for containing-versus-reduced effects or
#'   `"sequential"` for formula-order effects.
#' @param term_scope For partial effects, remove only assigned term columns or
#'   the term and all descendant interactions. The default is `"hierarchical"`
#'   for partial effects and `"assigned"` for sequential effects.
#' @param alignment `"shared"` for exact feature correspondence or `"soft"`
#'   for user-supplied sparse reference-to-subject maps.
#' @param subject_model `"fixed_intercept"` or centered `"ridge_slopes"`.
#' @param subject_penalty Positive ridge penalty for centered subject-specific
#'   deviations when `subject_model = "ridge_slopes"`.
#' @param random_terms Optional formula terms receiving subject deviations.
#'   The default is every term in `formula`.
#' @param spatial_penalty Nonnegative loading-smoothness strength.
#' @param spatial_laplacian Optional sparse feature Laplacian. Required when
#'   `spatial_penalty > 0`; named dimensions are aligned to `feature_names`.
#' @param alignment_penalty Optional reference-coefficient ridge stabilization
#'   for the soft map solver. It must be zero for shared alignment.
#' @param store `"sufficient"` stores compact component representations;
#'   `"fitted"` or `"residuals"` additionally stores reference effect matrices.
#' @param control A named list returned by [design_mfa_control()].
#' @param ... Reserved for future methods.
#'
#' @return A `design_mfa` object containing the weighted design fit, an
#'   estimability report, and one `design_mfa_effect` object per requested term.
#'   Training responses are retained in `fit_spec$refit$data` for subject-clustered
#'   refits; the stored model does not keep the dead `q x p` `XtY` product.
#'
#' @examples
#' set.seed(1)
#' subjects <- paste0("S", 1:6)
#' design <- lapply(subjects, function(id) {
#'   data.frame(
#'     saliency = ordered(rep(c(20, 40, 60, 80, 100), each = 6)),
#'     probe = factor(rep(c("old", "new", "lure"), 10))
#'   )
#' })
#' pattern <- seq(-1, 1, length.out = 8)
#' Y <- lapply(design, function(d) {
#'   outer(as.numeric(d$saliency) - 3, pattern) +
#'     matrix(rnorm(nrow(d) * 8, sd = 0.2), nrow(d), 8)
#' })
#' names(Y) <- names(design) <- subjects
#' feature_names <- paste0("voxel_", 1:8)
#' fit <- design_mfa(
#'   Y,
#'   design,
#'   ~ saliency + probe,
#'   contrasts = list(saliency = "poly", probe = "sum"),
#'   feature_names = feature_names
#' )
#' components(fit, "saliency")
#' continuous_design <- lapply(design, function(d) {
#'   d$trial_index <- seq_len(nrow(d))
#'   d
#' })
#' continuous_fit <- design_mfa(
#'   Y, continuous_design, ~ trial_index,
#'   feature_names = feature_names
#' )
#' head(effect_scores(continuous_fit, "trial_index"))
#'
#' @export
design_mfa <- function(Y, ...) {
  UseMethod("design_mfa")
}

#' @rdname design_mfa
#' @export
design_mfa.default <- function(Y, ...) {
  stop("Unsupported `Y` input. Supply a named list of subject trial-by-feature matrices.", call. = FALSE)
}

#' @rdname design_mfa
#' @export
design_mfa.multidesign <- function(Y,
                                   design = NULL,
                                   formula,
                                   subject,
                                   ...) {
  if (missing(subject)) {
    stop("`subject` must identify the subject column in `Y$design`.", call. = FALSE)
  }
  subject_quo <- rlang::enquo(subject)
  subject_name <- rlang::as_name(rlang::quo_get_expr(subject_quo))
  if (!subject_name %in% names(Y$design)) {
    stop(sprintf("Subject column '%s' is absent from `Y$design`.", subject_name), call. = FALSE)
  }
  subject_id <- as.character(Y$design[[subject_name]])
  if (anyNA(subject_id) || any(subject_id == "")) {
    stop("The subject column contains missing or empty identifiers.", call. = FALSE)
  }
  subject_levels <- unique(subject_id)
  row_index <- split(seq_len(nrow(Y$x)), factor(subject_id, levels = subject_levels))
  Y_list <- lapply(row_index, function(idx) as.matrix(Y$x[idx, , drop = FALSE]))
  design_list <- lapply(row_index, function(idx) Y$design[idx, , drop = FALSE])
  names(Y_list) <- names(design_list) <- subject_levels

  design_mfa.list(
    Y = Y_list,
    design = design_list,
    formula = formula,
    subject = subject_levels,
    ...
  )
}

#' @rdname design_mfa
#' @export
design_mfa.list <- function(Y,
                            design,
                            formula,
                            subject = NULL,
                            nuisance = NULL,
                            contrasts = NULL,
                            terms = NULL,
                            reference_grid = NULL,
                            reference_weight = c("equal_cells", "observed"),
                            trial_weights = NULL,
                            block_weight = c("equal_subject", "precision"),
                            feature_names = NULL,
                            feature_map = NULL,
                            rank = NULL,
                            ss_type = c("partial", "sequential"),
                            term_scope = NULL,
                            alignment = c("shared", "soft"),
                            subject_model = c("fixed_intercept", "ridge_slopes"),
                            subject_penalty = 1,
                            random_terms = NULL,
                            spatial_penalty = 0,
                            spatial_laplacian = NULL,
                            alignment_penalty = 0,
                            store = c("sufficient", "fitted", "residuals"),
                            control = design_mfa_control(),
                            ...) {
  fit_call <- match.call(expand.dots = FALSE)
  dots <- list(...)
  if (length(dots) > 0L) {
    stop(sprintf("Unused argument(s): %s.", paste(names(dots), collapse = ", ")), call. = FALSE)
  }
  ss_type <- match.arg(ss_type)
  if (is.null(term_scope)) {
    term_scope <- if (identical(ss_type, "partial")) "hierarchical" else "assigned"
  } else {
    term_scope <- match.arg(term_scope, c("assigned", "hierarchical"))
  }
  alignment <- match.arg(alignment)
  subject_model <- match.arg(subject_model)
  store <- match.arg(store)
  if (!is.numeric(reference_weight)) reference_weight <- match.arg(reference_weight)
  if (length(subject_penalty) != 1L || !is.finite(subject_penalty) || subject_penalty <= 0) {
    stop("`subject_penalty` must be one positive finite number.", call. = FALSE)
  }
  if (length(alignment_penalty) != 1L || !is.finite(alignment_penalty) || alignment_penalty < 0) {
    stop("`alignment_penalty` must be one nonnegative finite number.", call. = FALSE)
  }
  if (identical(alignment, "shared") && alignment_penalty != 0) {
    stop("`alignment_penalty` is only defined for soft alignment.", call. = FALSE)
  }
  if (identical(alignment, "soft") && !identical(subject_model, "fixed_intercept")) {
    stop("Soft alignment and ridge subject slopes cannot yet be combined.", call. = FALSE)
  }
  if (identical(ss_type, "sequential") && identical(term_scope, "hierarchical")) {
    stop("`term_scope = \"hierarchical\"` is only defined for partial effects.", call. = FALSE)
  }

  control <- .dmfa_merge_control(control)
  validated <- .dmfa_validate_inputs(
    Y = Y,
    design = design,
    subject = subject,
    feature_names = feature_names,
    alignment = alignment,
    feature_map = feature_map,
    feature_chunk_size = control$feature_chunk_size
  )
  Y <- validated$Y
  design <- validated$design
  subject <- validated$subject
  feature_names <- validated$feature_names
  feature_map <- validated$feature_map
  spatial <- .dmfa_validate_spatial_metric(
    spatial_laplacian = spatial_laplacian,
    spatial_penalty = spatial_penalty,
    feature_names = feature_names
  )

  weight_spec <- .dmfa_parse_weights(
    Y = Y,
    trial_weights = trial_weights,
    block_weight = block_weight
  )
  compiled <- .dmfa_compile_design(
    design = design,
    subject = subject,
    formula = formula,
    nuisance = nuisance,
    contrasts = contrasts,
    weights = weight_spec,
    terms = terms,
    reference_grid = reference_grid,
    reference_weight = reference_weight,
    control = control
  )

  if (is.null(random_terms)) {
    random_terms <- compiled$term_labels
  } else {
    random_terms <- unique(as.character(random_terms))
    unknown_random <- setdiff(random_terms, compiled$term_labels)
    if (length(unknown_random) > 0L) {
      stop(sprintf("Unknown `random_terms`: %s.", paste(unknown_random, collapse = ", ")), call. = FALSE)
    }
  }
  random_columns_global <- unlist(compiled$term_columns[random_terms], use.names = FALSE)

  Y_stacked <- if (identical(alignment, "shared") &&
      !identical(subject_model, "fixed_intercept")) {
    do.call(rbind, Y)
  } else {
    NULL
  }
  if (!is.null(Y_stacked)) colnames(Y_stacked) <- feature_names
  offsets <- cumsum(c(0L, vapply(Y, nrow, integer(1))))
  subject_rows <- lapply(seq_along(Y), function(s) seq.int(offsets[[s]] + 1L, offsets[[s + 1L]]))
  names(subject_rows) <- subject
  soft_map_cache <- if (identical(alignment, "soft")) {
    .dmfa_soft_map_cache(feature_map, control$feature_chunk_size)
  } else {
    NULL
  }
  fit_model <- function(indices) {
    if (identical(alignment, "soft")) {
      return(.dmfa_fit_soft_map(
        X = compiled$X[, indices, drop = FALSE],
        Y_list = Y,
        weight_list = weight_spec$weights,
        feature_map = feature_map,
        subject_rows = subject_rows,
        alignment_penalty = alignment_penalty,
        control = control,
        map_cache = soft_map_cache
      ))
    }
    if (identical(subject_model, "fixed_intercept")) {
      return(.dmfa_fit_shared_streamed(
        X = compiled$X[, indices, drop = FALSE],
        Y_list = Y,
        weight_list = weight_spec$weights,
        subject_rows = subject_rows,
        control = control
      ))
    }
    local_random <- match(intersect(random_columns_global, indices), indices)
    .dmfa_fit_model(
      X = compiled$X[, indices, drop = FALSE],
      Y = Y_stacked,
      weights = weight_spec$stacked,
      control = control,
      subject_model = subject_model,
      subject_rows = subject_rows,
      subject_weights = unname(weight_spec$subject_weight),
      random_columns = local_random,
      subject_penalty = subject_penalty
    )
  }
  global_full <- fit_model(seq_len(ncol(compiled$X)))
  if (!isTRUE(global_full$converged)) {
    stop(
      sprintf(
        "The model could not be fitted: %s. Increase the relevant iteration limit or relax its tolerance.",
        .dmfa_nonconvergence_message(global_full, "full", control)
      ),
      call. = FALSE
    )
  }

  effect_list <- lapply(compiled$selected_terms, function(term) {
    effect <- .dmfa_fit_effect(
      term = term,
      compiled = compiled,
      global_full = global_full,
      fit_model = fit_model,
      ss_type = ss_type,
      term_scope = term_scope,
      rank = rank,
      control = control,
      store = store,
      feature_metric = spatial$metric
    )
    effect$feature_names <- feature_names
    effect
  })
  names(effect_list) <- compiled$selected_terms

  estimability_table <- do.call(rbind, lapply(effect_list, function(effect) {
    data.frame(
      term = effect$term,
      status = effect$status,
      estimable_df = effect$estimable_df,
      grid_rows = nrow(effect$reference_grid),
      effect_nonestimable_grid_rows = sum(!effect$estimable_grid),
      prediction_nonestimable_grid_rows = sum(!effect$prediction_estimable_grid),
      estimability_scope = if (all(effect$prediction_estimable_grid)) {
        "effect_and_prediction"
      } else if (all(effect$estimable_grid)) {
        "effect_only"
      } else {
        "nonestimable_effect"
      },
      partial_ss = effect$partial_ss,
      partial_r2 = effect$partial_r2,
      objective_delta = effect$objective_delta,
      objective_fraction = effect$objective_fraction,
      comparison = effect$comparison,
      component_rank = effect$rank,
      stringsAsFactors = FALSE
    )
  }))
  rownames(estimability_table) <- NULL

  fit <- list(
    call = fit_call,
    formula = formula,
    nuisance = nuisance,
    subjects = subject,
    subject_model = subject_model,
    subject_penalty = subject_penalty,
    random_terms = random_terms,
    subject_deviations = global_full$subject_deviations,
    alignment = alignment,
    computation = global_full$computation,
    feature_map = feature_map,
    subject_feature_counts = validated$subject_feature_counts,
    alignment_coverage = validated$alignment_coverage,
    alignment_penalty = alignment_penalty,
    spatial_penalty = spatial_penalty,
    spatial_laplacian = spatial$laplacian,
    feature_metric = spatial$metric,
    feature_names = feature_names,
    n_features = length(feature_names),
    n_trials = vapply(Y, nrow, integer(1)),
    weights = weight_spec$weights,
    raw_trial_weights = weight_spec$raw,
    subject_weights = weight_spec$subject_weight,
    weight_scheme = weight_spec$scheme,
    design = compiled,
    full_model = .dmfa_compact_fitted_model(global_full),
    effects = effect_list,
    estimability = estimability_table,
    ss_type = ss_type,
    term_scope = term_scope,
    store = store,
    control = control
  )
  class(fit) <- c("design_mfa", "list")
  block_weight_refit <- if (is.numeric(block_weight)) weight_spec$subject_weight else weight_spec$scheme
  frozen_reference_grid <- compiled$reference_grid
  frozen_reference_weight <- compiled$reference_weights
  refit_data <- list(
    Y = Y,
    design = design,
    trial_weights = weight_spec$raw,
    block_weight = block_weight_refit,
    reference_grid = frozen_reference_grid,
    reference_weight = frozen_reference_weight
  )
  if (!is.null(feature_map)) refit_data$feature_map <- feature_map
  .muscal_attach_fit_contract(
    fit,
    method = "design_mfa",
    task = "design_decomposition",
    oos_types = "response",
    fit_call = fit_call,
    refit_supported = TRUE,
    prediction_target = "activation_field",
    refit = .muscal_make_refit_spec(
      data = refit_data,
      fit_fn = .dmfa_make_refit_fit_fn(
        formula = formula,
        nuisance = nuisance,
        contrasts = contrasts,
        terms = terms,
        reference_grid = frozen_reference_grid,
        reference_weight = frozen_reference_weight,
        feature_names = feature_names,
        rank = rank,
        ss_type = ss_type,
        term_scope = term_scope,
        alignment = alignment,
        subject_model = subject_model,
        subject_penalty = subject_penalty,
        random_terms = random_terms,
        spatial_penalty = spatial_penalty,
        spatial_laplacian = spatial_laplacian,
        alignment_penalty = alignment_penalty,
        store = store,
        control = control
      ),
      bootstrap_fn = .dmfa_bootstrap_subject_data,
      permutation_fn = .dmfa_refuse_implicit_permutation,
      resample_unit = "subjects"
    )
  )
}

#' Control Numerical Design-MFA Calculations
#'
#' @param rank_tol Relative eigenvalue threshold for numerical design rank and
#'   effect-component rank.
#' @param estimability_tol Relative tolerance for reference-grid estimability.
#' @param max_numeric_grid Maximum number of distinct numeric values included
#'   automatically in a reference grid. Variables with more distinct values
#'   use this many bounded sample quantiles.
#' @param max_grid_rows Maximum allowed reference-grid size.
#' @param feature_chunk_size Maximum number of response features processed at
#'   once while accumulating shared-grid sufficient statistics.
#' @param multilevel_max_iter Maximum block-coordinate iterations for ridge
#'   subject deviations. The default leaves headroom for weak-pooling fits,
#'   which can converge substantially more slowly than strongly pooled fits.
#' @param multilevel_tol Relative penalized-objective convergence tolerance.
#' @param soft_max_iter Maximum conjugate-gradient iterations for soft alignment.
#' @param soft_tol Relative normal-equation residual tolerance for soft alignment.
#'
#' @return A named control list for [design_mfa()].
#'
#' @examples
#' control <- design_mfa_control(feature_chunk_size = 1000L, soft_tol = 1e-7)
#' control[c("feature_chunk_size", "soft_tol")]
#' @export
design_mfa_control <- function(rank_tol = sqrt(.Machine$double.eps),
                               estimability_tol = 1e-8,
                               max_numeric_grid = 20L,
                               max_grid_rows = 5000L,
                               feature_chunk_size = 25000L,
                               multilevel_max_iter = 5000L,
                               multilevel_tol = 1e-8,
                               soft_max_iter = 500L,
                               soft_tol = 1e-8) {
  .dmfa_control(
    rank_tol,
    estimability_tol,
    max_numeric_grid,
    max_grid_rows,
    feature_chunk_size,
    multilevel_max_iter,
    multilevel_tol,
    soft_max_iter,
    soft_tol
  )
}

.dmfa_make_refit_fit_fn <- function(formula,
                                    nuisance,
                                    contrasts,
                                    terms,
                                    reference_grid,
                                    reference_weight,
                                    feature_names,
                                    rank,
                                    ss_type,
                                    term_scope,
                                    alignment,
                                    subject_model,
                                    subject_penalty,
                                    random_terms,
                                    spatial_penalty,
                                    spatial_laplacian,
                                    alignment_penalty,
                                    store,
                                    control) {
  function(data) {
    required <- c("Y", "design", "trial_weights", "block_weight", "reference_grid", "reference_weight")
    if (!is.list(data) || !all(required %in% names(data))) {
      stop("Callbacks must return a complete refit-data list including weights and frozen reference fields.", call. = FALSE)
    }
    if (!identical(data$reference_grid, reference_grid) ||
        !identical(data$reference_weight, reference_weight)) {
      stop("Callbacks must preserve the frozen reference grid and weights.", call. = FALSE)
    }
    ids <- names(data$Y)
    if (!is.list(data$Y) || !length(data$Y) || is.null(ids) || anyNA(ids) ||
        any(!nzchar(ids)) || anyDuplicated(ids) || !is.list(data$design) ||
        !identical(names(data$design), ids) || !is.list(data$trial_weights) ||
        !identical(names(data$trial_weights), ids) ||
        (!is.null(data$feature_map) && !identical(names(data$feature_map), ids)) ||
        (is.numeric(data$block_weight) && !identical(names(data$block_weight), ids))) {
      stop("Refit-data subject identities must align across responses, design, maps and weights.", call. = FALSE)
    }
    design_mfa(
      Y = data$Y,
      design = data$design,
      formula = formula,
      subject = names(data$Y),
      nuisance = nuisance,
      contrasts = contrasts,
      terms = terms,
      reference_grid = reference_grid,
      reference_weight = reference_weight,
      trial_weights = data$trial_weights,
      block_weight = data$block_weight,
      feature_names = feature_names,
      feature_map = data$feature_map,
      rank = rank,
      ss_type = ss_type,
      term_scope = term_scope,
      alignment = alignment,
      subject_model = subject_model,
      subject_penalty = subject_penalty,
      random_terms = random_terms,
      spatial_penalty = spatial_penalty,
      spatial_laplacian = spatial_laplacian,
      alignment_penalty = alignment_penalty,
      store = store,
      control = control
    )
  }
}

.dmfa_get_effect <- function(x, term) {
  if (!inherits(x, "design_mfa")) stop("`x` must be a fitted `design_mfa` object.", call. = FALSE)
  if (missing(term) || length(term) != 1L || is.na(term)) {
    stop("`term` must identify one fitted design term.", call. = FALSE)
  }
  term <- as.character(term)
  if (!term %in% names(x$effects)) {
    stop(sprintf("Unknown term '%s'. Available terms: %s.", term, paste(names(x$effects), collapse = ", ")), call. = FALSE)
  }
  x$effects[[term]]
}

#' Extract Design Effects
#'
#' @param x A fitted object.
#' @param term Optional design term.
#' @param ... Unused.
#'
#' @return A list of effects, or one `design_mfa_effect` when `term` is supplied.
#'
#' @examples
#' example_data <- design_mfa_example_data()
#' example_fit <- design_mfa(example_data$Y, example_data$design, ~ saliency + probe)
#' effects(example_fit, "saliency")
#' @importFrom stats effects
#' @export
stats::effects

#' Extract Design Effects from a Design-MFA Fit
#'
#' @param object A fitted `design_mfa` object.
#' @param term Optional design term.
#' @param ... Unused.
#'
#' @return A list of effects, or one `design_mfa_effect` when `term` is supplied.
#'
#' @examples
#' example_data <- design_mfa_example_data()
#' example_fit <- design_mfa(example_data$Y, example_data$design, ~ saliency + probe)
#' effects(example_fit, "probe")
#' @export
effects.design_mfa <- function(object, term = NULL, ...) {
  if (is.null(term)) object$effects else .dmfa_get_effect(object, term)
}

#' Summarize Components for a Design Effect
#'
#' @param x A fitted `design_mfa` object.
#' @param term A fitted design term.
#' @param ... Unused.
#'
#' @return A data frame of singular values and fractions of total weighted
#'   reference-grid effect energy. Fractions need not sum to one when `rank`
#'   truncates the effect decomposition.
#' @importFrom multivarious components
#' @export
multivarious::components

#' Summarize Components for a Design-MFA Effect
#'
#' @param x A fitted `design_mfa` object.
#' @param term A fitted design term.
#' @param ... Unused.
#'
#' @return A data frame of singular values and fractions of total weighted
#'   reference-grid effect energy. Fractions need not sum to one when `rank`
#'   truncates the effect decomposition.
#'
#' @examples
#' example_data <- design_mfa_example_data()
#' example_fit <- design_mfa(example_data$Y, example_data$design, ~ saliency)
#' components(example_fit, "saliency")
#' @export
components.design_mfa <- function(x, term, ...) {
  effect <- .dmfa_get_effect(x, term)
  data.frame(
    component = seq_along(effect$singular_values),
    singular_value = effect$singular_values,
    explained = effect$explained
  )
}

#' Extract Effect Scores or Spatial Loadings
#'
#' @param x A fitted `design_mfa` object.
#' @param term A fitted design term.
#' @param components Optional component indices.
#'
#' @return A numeric matrix.
#'
#' @examples
#' example_data <- design_mfa_example_data()
#' example_fit <- design_mfa(example_data$Y, example_data$design, ~ saliency)
#' head(effect_scores(example_fit, "saliency"))
#' head(effect_loadings(example_fit, "saliency"))
#' @export
effect_scores <- function(x, term, components = NULL) {
  effect <- .dmfa_get_effect(x, term)
  if (is.null(effect$scores)) stop(sprintf("Term '%s' has no estimable component scores.", term), call. = FALSE)
  if (is.null(components)) return(effect$scores)
  effect$scores[, components, drop = FALSE]
}

#' @rdname effect_scores
#' @export
effect_loadings <- function(x, term, components = NULL) {
  effect <- .dmfa_get_effect(x, term)
  if (is.null(effect$loadings)) stop(sprintf("Term '%s' has no estimable spatial loadings.", term), call. = FALSE)
  if (is.null(components)) return(effect$loadings)
  effect$loadings[, components, drop = FALSE]
}

#' Extract One Effect Loading Map
#'
#' @param x A fitted `design_mfa` object.
#' @param term A fitted design term.
#' @param component A single component index.
#'
#' @return A named numeric loading vector.
#'
#' @examples
#' example_data <- design_mfa_example_data()
#' example_fit <- design_mfa(example_data$Y, example_data$design, ~ saliency)
#' effect_map(example_fit, "saliency", component = 1L)
#' @export
effect_map <- function(x, term, component = 1L) {
  loading <- effect_loadings(x, term, components = component)
  out <- as.numeric(loading[, 1L])
  names(out) <- rownames(loading)
  out
}

#' Report Design-MFA Estimability
#'
#' @param x A fitted object.
#' @param ... Unused.
#'
#' @return A data frame with term-level rank and reference-grid diagnostics.
#'   `effect_nonestimable_grid_rows` concerns the tested term contrast, whereas
#'   `prediction_nonestimable_grid_rows` concerns the complete fitted response,
#'   always against the global model and independently of `ss_type`.
#'   These can differ when an unrelated between-subject covariate is aliased
#'   with subject fixed effects.
#'
#' @examples
#' example_data <- design_mfa_example_data()
#' ages <- seq(20, 70, length.out = length(example_data$design))
#' adjusted_design <- Map(function(d, age) transform(d, age = age),
#'                        example_data$design, ages)
#' names(adjusted_design) <- names(example_data$design)
#' adjusted <- design_mfa(
#'   example_data$Y, adjusted_design, ~ saliency + age,
#'   contrasts = list(saliency = "poly"), terms = "saliency"
#' )
#' estimability(adjusted)
#' @export
estimability <- function(x, ...) {
  UseMethod("estimability")
}

#' @rdname estimability
#' @export
estimability.design_mfa <- function(x, ...) {
  x$estimability
}

#' Report Reference-Feature Alignment Coverage
#'
#' @param x A fitted `design_mfa` object.
#'
#' @return A data frame with the number and fraction of subjects supporting
#'   each reference feature and the aggregate sparse-map weight. Shared
#'   alignment reports complete support.
#'
#' @examples
#' example_data <- design_mfa_example_data()
#' identity_maps <- lapply(example_data$Y, function(y) {
#'   P <- Matrix::Diagonal(ncol(y))
#'   dimnames(P) <- list(colnames(y), colnames(y))
#'   P
#' })
#' names(identity_maps) <- names(example_data$Y)
#' soft_fit <- design_mfa(
#'   example_data$Y, example_data$design, ~ saliency,
#'   alignment = "soft", feature_map = identity_maps
#' )
#' head(alignment_coverage(soft_fit))
#' @export
alignment_coverage <- function(x) {
  if (!inherits(x, "design_mfa")) {
    stop("`x` must be a fitted `design_mfa` object.", call. = FALSE)
  }
  x$alignment_coverage
}

#' Extract Subject-Specific Design Effects
#'
#' @param x A fitted `design_mfa` object.
#' @param term Optional term.
#'
#' @return Subject deviations when fitted by a partial-pooling model.
#'
#' @examples
#' example_data <- design_mfa_example_data()
#' pooled <- design_mfa(
#'   example_data$Y, example_data$design, ~ saliency,
#'   subject_model = "ridge_slopes", random_terms = "saliency"
#' )
#' lapply(subject_effects(pooled, "saliency"), dim)
#' @export
subject_effects <- function(x, term = NULL) {
  if (!inherits(x, "design_mfa")) stop("`x` must be a fitted `design_mfa` object.", call. = FALSE)
  if (is.null(x$subject_deviations)) {
    stop(
      "Subject-specific deviations are unavailable for `subject_model = \"fixed_intercept\"`.",
      call. = FALSE
    )
  }
  if (is.null(term)) return(x$subject_deviations)
  if (!term %in% x$design$term_labels) {
    stop(sprintf("Unknown term '%s'.", term), call. = FALSE)
  }
  wanted <- colnames(x$design$X)[x$design$term_columns[[term]]]
  lapply(x$subject_deviations, function(U) {
    keep <- intersect(wanted, rownames(U))
    U[keep, , drop = FALSE]
  })
}

.dmfa_prediction_matrix <- function(object, new_design) {
  E <- .dmfa_build_effect_matrix(object$design, new_design)
  m <- nrow(E)
  subject_part <- matrix(
    rep(object$subject_weights, each = m),
    nrow = m,
    byrow = FALSE
  )
  colnames(subject_part) <- colnames(object$design$X)[object$design$subject_columns]
  if (length(object$design$nuisance_columns) > 0L) {
    nuisance_part <- matrix(
      rep(object$design$nuisance_mean, each = m),
      nrow = m,
      byrow = FALSE
    )
    colnames(nuisance_part) <- colnames(object$design$X)[object$design$nuisance_columns]
  } else {
    nuisance_part <- matrix(0, m, 0L)
  }
  L <- cbind(subject_part, nuisance_part, E)
  colnames(L) <- colnames(object$design$X)
  L
}

#' Predict Group-Average Activation Fields
#'
#' @param object A fitted `design_mfa` object.
#' @param new_design Optional design rows. Defaults to the fitted reference grid.
#' @param type Currently `"response"`.
#' @param ... Unused.
#'
#' @details
#' Full-response prediction is stricter than extracting a fitted design effect.
#' Every coefficient functional needed for each requested prediction row must
#' be estimable under the complete fitted model. A within-subject term effect
#' and its partial sum of squares can therefore remain estimable when an
#' unrelated between-subject covariate is aliased with the subject intercepts,
#' even though `predict()` refuses the corresponding full response. Reported
#' term effects use only the estimable contribution of the columns removed by
#' that term's hypothesis. Use [estimability()] to compare the
#' `effect_nonestimable_grid_rows` and
#' `prediction_nonestimable_grid_rows` diagnostics before predicting.
#'
#' @return A matrix of group-average fitted activation fields.
#'
#' @examples
#' example_data <- design_mfa_example_data()
#' example_fit <- design_mfa(example_data$Y, example_data$design, ~ saliency + probe)
#' predict(example_fit, new_design = data.frame(
#'   saliency = ordered(c(20, 100), levels = c(20, 40, 60, 80, 100)),
#'   probe = factor(c("old", "lure"), levels = c("old", "new", "lure"))
#' ))
#' @export
predict.design_mfa <- function(object,
                               new_design = NULL,
                               type = c("response"),
                               ...) {
  type <- match.arg(type)
  if (is.null(new_design)) new_design <- object$design$reference_grid
  L <- .dmfa_prediction_matrix(object, new_design)
  check <- .dmfa_estimable_rows(L, object$full_model, object$control$estimability_tol)
  if (!all(check$estimable)) {
    stop(
      sprintf(
        paste0(
          "%d requested design rows are not estimable for full-response prediction. ",
          "A term contrast may still be estimable when the failure comes from unrelated aliased columns; ",
          "see `estimability(object)$prediction_nonestimable_grid_rows`."
        ),
        sum(!check$estimable)
      ),
      call. = FALSE
    )
  }
  out <- L %*% object$full_model$coefficient
  colnames(out) <- object$feature_names
  out
}

#' Project New Fields onto a Design-Effect Basis
#'
#' @param x A `design_mfa_effect` object.
#' @param new_data A numeric row-by-feature matrix.
#' @param ... Unused.
#'
#' @return Component scores for the new fields.
#'
#' @examples
#' example_data <- design_mfa_example_data()
#' example_fit <- design_mfa(example_data$Y, example_data$design, ~ saliency)
#' project(effects(example_fit, "saliency"), example_data$Y[[1L]][1:3, ])
#' @export
project.design_mfa_effect <- function(x, new_data, ...) {
  if (!inherits(x, "design_mfa_effect")) stop("`x` must be a `design_mfa_effect`.", call. = FALSE)
  if (is.null(x$loadings)) stop("This effect has no estimable loading basis.", call. = FALSE)
  new_data <- as.matrix(new_data)
  if (!is.numeric(new_data) || any(!is.finite(new_data))) {
    stop("`new_data` must be a finite numeric matrix.", call. = FALSE)
  }
  if (!is.null(colnames(new_data))) {
    if (!setequal(colnames(new_data), x$feature_names)) {
      stop("`new_data` does not match the fitted feature schema.", call. = FALSE)
    }
    new_data <- new_data[, x$feature_names, drop = FALSE]
  } else if (ncol(new_data) != length(x$feature_names)) {
    stop(sprintf("`new_data` has %d features; expected %d.", ncol(new_data), length(x$feature_names)), call. = FALSE)
  }
  new_data %*% x$loadings
}

#' @export
print.design_mfa <- function(x, ...) {
  cat("Design effect decomposition (ASCA+/partial RDA)\n")
  cat(sprintf("  Subjects: %d\n", length(x$subjects)))
  cat(sprintf("  Trials: %d total (%s)\n", sum(x$n_trials), paste(x$n_trials, collapse = ", ")))
  cat(sprintf("  Features: %d\n", x$n_features))
  cat(sprintf("  Design rank: %d / %d\n", x$full_model$rank, ncol(x$design$X)))
  cat(sprintf("  Weighting: %s\n", x$weight_scheme))
  cat(sprintf("  Alignment: %s\n", x$alignment))
  cat(sprintf("  Computation: %s\n", x$computation))
  cat(sprintf("  Subject model: %s\n", x$subject_model))
  if (x$full_model$iterations > 0L) {
    cat(sprintf(
      "  Solver: %s after %d iterations\n",
      if (isTRUE(x$full_model$converged)) "converged" else "not converged",
      x$full_model$iterations
    ))
  }
  cat(sprintf("  Effects: %d estimable of %d\n", sum(x$estimability$status == "estimable"), nrow(x$estimability)))
  invisible(x)
}

#' @export
print.design_mfa_effect <- function(x, ...) {
  cat(sprintf("Design effect: %s\n", x$term))
  cat(sprintf("  Status: %s\n", x$status))
  cat(sprintf("  Estimable df: %d\n", x$estimable_df))
  cat(sprintf("  Components: %d\n", x$rank))
  if (identical(x$comparison, "penalized_objective")) {
    cat(sprintf("  Objective fraction: %s\n", format(x$objective_fraction, digits = 4)))
  } else {
    cat(sprintf("  Partial R2: %s\n", format(x$partial_r2, digits = 4)))
  }
  if (!is.null(x$reason)) cat(sprintf("  Note: %s\n", x$reason))
  invisible(x)
}

#' @export
summary.design_mfa <- function(object, ...) {
  out <- object$estimability
  attr(out, "design_rank") <- object$full_model$rank
  attr(out, "n_subjects") <- length(object$subjects)
  class(out) <- c("summary_design_mfa", class(out))
  out
}

#' @export
print.summary_design_mfa <- function(x, ...) {
  cat(sprintf(
    "Design-effect summary: %d subjects, design rank %d\n",
    attr(x, "n_subjects"), attr(x, "design_rank")
  ))
  print.data.frame(x, row.names = FALSE, ...)
  invisible(x)
}

#' Plot a Design-Effect Component Trajectory
#'
#' @param x A fitted `design_mfa` object.
#' @param term Optional fitted design term. By default, the first estimable
#'   effect with a retained component is plotted.
#' @param component Component index.
#' @param xvar Optional reference-grid variable for the horizontal axis.
#' @param ... Additional arguments passed to [graphics::plot()].
#'
#' @return The plotted data, invisibly.
#'
#' @examples
#' example_data <- design_mfa_example_data()
#' example_fit <- design_mfa(example_data$Y, example_data$design, ~ saliency)
#' plot(example_fit)
#' @export
plot.design_mfa <- function(x, term = NULL, component = 1L, xvar = NULL, ...) {
  if (is.null(term)) {
    eligible <- names(x$effects)[vapply(
      x$effects,
      function(effect) isTRUE(effect$estimable) && effect$rank > 0L,
      logical(1)
    )]
    if (length(eligible) == 0L) {
      stop("No estimable fitted effect has a retained component to plot.", call. = FALSE)
    }
    term <- eligible[[1L]]
  }
  effect <- .dmfa_get_effect(x, term)
  score <- effect_scores(x, term, component)[, 1L]
  grid <- effect$reference_grid
  varying <- names(grid)[vapply(grid, function(z) length(unique(z)) > 1L, logical(1))]
  if (is.null(xvar)) {
    xvar <- if (length(varying) == 0L) NULL else varying[[1L]]
  }
  if (!is.null(xvar) && !xvar %in% names(grid)) {
    stop(sprintf("Unknown reference-grid variable '%s'.", xvar), call. = FALSE)
  }
  xvalue <- if (is.null(xvar)) seq_len(nrow(grid)) else grid[[xvar]]
  xplot <- if (is.factor(xvalue)) as.numeric(xvalue) else xvalue
  graphics::plot(
    xplot,
    score,
    xlab = if (is.null(xvar)) "reference row" else xvar,
    ylab = sprintf("%s component %d score", term, component),
    xaxt = if (is.factor(xvalue)) "n" else "s",
    ...
  )
  if (!is.null(xvar) && is.factor(xvalue)) {
    graphics::axis(1, at = seq_along(levels(xvalue)), labels = levels(xvalue))
  }
  invisible(cbind(grid, score = score))
}

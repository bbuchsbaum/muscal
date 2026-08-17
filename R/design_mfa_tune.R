.dmfa_reconstruct_effect <- function(effect, feature_metric = NULL) {
  if (is.null(effect$scores) || is.null(effect$loadings)) {
    stop("The effect has no retained component representation.", call. = FALSE)
  }
  dual_loadings <- if (is.null(feature_metric)) {
    effect$loadings
  } else {
    as.matrix(feature_metric %*% effect$loadings)
  }
  effect$scores %*% t(dual_loadings)
}

.dmfa_tuning_grid <- function(grid, spatial_laplacian, subject_model, alignment) {
  if (is.null(grid)) {
    grid <- expand.grid(
      rank = 1:3,
      spatial_penalty = if (is.null(spatial_laplacian)) 0 else c(0, 1),
      subject_penalty = if (identical(subject_model, "ridge_slopes")) c(0.1, 1, 10) else 1,
      alignment_penalty = if (identical(alignment, "soft")) c(0, 0.01, 0.1) else 0,
      KEEP.OUT.ATTRS = FALSE
    )
  } else if (is.list(grid) && !is.data.frame(grid)) {
    if (is.null(names(grid)) || any(names(grid) == "")) {
      stop("A list `grid` must be named.", call. = FALSE)
    }
    grid <- expand.grid(grid, KEEP.OUT.ATTRS = FALSE)
  } else if (!is.data.frame(grid)) {
    stop("`grid` must be NULL, a data frame, or a named list.", call. = FALSE)
  }

  allowed <- c("rank", "spatial_penalty", "subject_penalty", "alignment_penalty")
  unknown <- setdiff(names(grid), allowed)
  if (length(unknown) > 0L) {
    stop(sprintf("Unknown tuning-grid field(s): %s.", paste(unknown, collapse = ", ")), call. = FALSE)
  }
  defaults <- list(rank = 1L, spatial_penalty = 0, subject_penalty = 1, alignment_penalty = 0)
  for (name in allowed) if (!name %in% names(grid)) grid[[name]] <- defaults[[name]]
  grid <- grid[, allowed, drop = FALSE]
  if (nrow(grid) < 1L) stop("`grid` must contain at least one candidate.", call. = FALSE)

  raw_rank <- as.numeric(grid$rank)
  if (any(!is.finite(raw_rank)) || any(raw_rank < 1) || any(abs(raw_rank - round(raw_rank)) > 1e-8)) {
    stop("Tuning `rank` values must be positive integers.", call. = FALSE)
  }
  grid$rank <- as.integer(round(raw_rank))
  grid$spatial_penalty <- as.numeric(grid$spatial_penalty)
  grid$subject_penalty <- as.numeric(grid$subject_penalty)
  grid$alignment_penalty <- as.numeric(grid$alignment_penalty)
  if (any(!is.finite(grid$spatial_penalty)) || any(grid$spatial_penalty < 0)) {
    stop("Tuning `spatial_penalty` values must be nonnegative and finite.", call. = FALSE)
  }
  if (any(!is.finite(grid$subject_penalty)) || any(grid$subject_penalty <= 0)) {
    stop("Tuning `subject_penalty` values must be positive and finite.", call. = FALSE)
  }
  if (any(!is.finite(grid$alignment_penalty)) || any(grid$alignment_penalty < 0)) {
    stop("Tuning `alignment_penalty` values must be nonnegative and finite.", call. = FALSE)
  }
  if (any(grid$spatial_penalty > 0) && is.null(spatial_laplacian)) {
    stop("A positive tuning value for `spatial_penalty` requires `spatial_laplacian`.", call. = FALSE)
  }
  if (!identical(subject_model, "ridge_slopes") && length(unique(grid$subject_penalty)) > 1L) {
    stop("Varying `subject_penalty` requires `subject_model = \"ridge_slopes\"`.", call. = FALSE)
  }
  if (!identical(alignment, "soft") && any(grid$alignment_penalty > 0)) {
    stop("A positive tuning value for `alignment_penalty` requires `alignment = \"soft\"`.", call. = FALSE)
  }
  unique(grid)
}

.dmfa_subject_folds <- function(subject, n_folds, seed) {
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
  shuffled <- sample(sort(subject), length(subject), replace = FALSE)
  assessment <- split(shuffled, rep(seq_len(n_folds), length.out = length(shuffled)))
  names(assessment) <- paste0("Fold_", seq_along(assessment))
  lapply(assessment, function(test) {
    list(analysis = setdiff(sort(subject), test), assessment = sort(test))
  })
}

.dmfa_tuning_subject_names <- function(Y) {
  subject_names <- names(Y)
  if (is.null(subject_names)) subject_names <- paste0("Subject_", seq_along(Y))
  subject_names
}

.dmfa_canonicalize_subject_arg <- function(value,
                                           argument,
                                           subject_names,
                                           type = c("list", "numeric")) {
  type <- match.arg(type)
  if (is.null(value)) return(value)
  if (identical(type, "list") && !is.list(value)) return(value)
  if (identical(type, "numeric") && !is.numeric(value)) return(value)

  n_subjects <- length(subject_names)
  value_names <- names(value)
  if (is.null(value_names)) {
    if (length(value) != n_subjects) {
      stop(
        sprintf(
          "Unnamed `%s` must have exactly one entry per Y block (%d); got %d.",
          argument, n_subjects, length(value)
        ),
        call. = FALSE
      )
    }
    names(value) <- subject_names
    return(value)
  }

  invalid_names <- is.na(value_names) | value_names == ""
  duplicated_names <- unique(value_names[duplicated(value_names) & !invalid_names])
  if (any(invalid_names) || length(duplicated_names) > 0L) {
    details <- character(0)
    if (any(invalid_names)) details <- c(details, "missing or empty keys")
    if (length(duplicated_names) > 0L) {
      details <- c(details, sprintf("duplicated keys {%s}", paste(duplicated_names, collapse = ", ")))
    }
    stop(
      sprintf("Named `%s` has %s.", argument, paste(details, collapse = " and ")),
      call. = FALSE
    )
  }

  missing_keys <- setdiff(subject_names, value_names)
  extra_keys <- setdiff(value_names, subject_names)
  if (length(missing_keys) > 0L || length(extra_keys) > 0L) {
    format_keys <- function(x) if (length(x) == 0L) "none" else paste(x, collapse = ", ")
    stop(
      sprintf(
        "Named `%s` must match the Y block keys; missing {%s}; extra {%s}.",
        argument,
        format_keys(missing_keys),
        format_keys(extra_keys)
      ),
      call. = FALSE
    )
  }
  value[subject_names]
}

.dmfa_canonicalize_tuning_args <- function(dots, Y) {
  subject_names <- .dmfa_tuning_subject_names(Y)
  dots$trial_weights <- .dmfa_canonicalize_subject_arg(
    dots$trial_weights, "trial_weights", subject_names, "list"
  )
  dots$feature_map <- .dmfa_canonicalize_subject_arg(
    dots$feature_map, "feature_map", subject_names, "list"
  )
  if (is.numeric(dots$block_weight)) {
    dots$block_weight <- .dmfa_canonicalize_subject_arg(
      dots$block_weight, "block_weight", subject_names, "numeric"
    )
  }
  dots
}

.dmfa_subset_tuning_args <- function(base_args,
                                     indices,
                                     Y,
                                     design,
                                     subject,
                                     total_subjects) {
  out <- base_args
  out$Y <- Y[indices]
  out$design <- design[indices]
  out$subject <- subject[indices]
  subset_subject_list <- function(value) {
    if (!is.list(value)) return(value)
    if (length(value) != total_subjects || is.null(names(value)) ||
        !identical(names(value), names(Y))) {
      stop("Internal tuning arguments were not canonicalized before fold subsetting.", call. = FALSE)
    }
    value <- value[indices]
    names(value) <- names(out$Y)
    value
  }
  out$trial_weights <- subset_subject_list(out$trial_weights)
  out$feature_map <- subset_subject_list(out$feature_map)
  if (is.numeric(out$block_weight) && length(out$block_weight) == total_subjects) {
    if (is.null(names(out$block_weight)) || !identical(names(out$block_weight), names(Y))) {
      stop("Internal tuning block weights were not canonicalized before fold subsetting.", call. = FALSE)
    }
    out$block_weight <- out$block_weight[indices]
    names(out$block_weight) <- names(out$Y)
  }
  out
}

.dmfa_effect_holdout_score <- function(training, target, terms) {
  term_scores <- stats::setNames(rep(NA_real_, length(terms)), terms)
  for (term in terms) {
    train_effect <- training$effects[[term]]
    target_effect <- target$effects[[term]]
    if (is.null(train_effect) || is.null(target_effect) ||
        !isTRUE(train_effect$estimable) || !isTRUE(target_effect$estimable) ||
        is.null(target_effect$effect_matrix)) {
      next
    }
    estimate <- .dmfa_reconstruct_effect(train_effect, training$feature_metric)
    truth <- target_effect$effect_matrix
    if (!identical(dim(estimate), dim(truth))) next
    weights <- target_effect$reference_weights
    error <- sum(rowSums((truth - estimate)^2) * weights)
    energy <- sum(rowSums(truth^2) * weights)
    if (is.finite(error) && is.finite(energy) && energy > .Machine$double.eps) {
      term_scores[[term]] <- error / energy
    }
  }
  finite <- is.finite(term_scores)
  list(
    score = if (any(finite)) mean(term_scores[finite]) else NA_real_,
    n_terms = sum(finite),
    term_scores = term_scores
  )
}

.dmfa_assert_tuning_subjects <- function(actual, expected, stage, fold_id) {
  actual <- sort(as.character(actual))
  expected <- sort(as.character(expected))
  if (!identical(actual, expected)) {
    stop(
      sprintf(
        "Subject leakage detected in fold %d %s fit: expected {%s}, got {%s}.",
        fold_id,
        stage,
        paste(expected, collapse = ", "),
        paste(actual, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  actual
}

.dmfa_tuning_failure_reason <- function(message) {
  if (grepl("did not converge", message, fixed = TRUE)) return("nonconvergence")
  if (grepl("not estimable", message, ignore.case = TRUE) ||
      grepl("No requested term was estimable", message, fixed = TRUE)) {
    return("nonestimable_effect")
  }
  if (grepl("nonzero held-out effect energy", message, fixed = TRUE) ||
      grepl("zero target energy", message, ignore.case = TRUE)) {
    return("zero_target_energy")
  }
  if (grepl("Invalid|must |requires |Unknown ", message)) return("invalid_candidate")
  "other_fit_error"
}

.dmfa_numeric_fingerprint <- function(x) {
  values <- as.numeric(x)
  if (length(values) == 0L) return(c(sum = 0, sumsq = 0, indexed = 0))
  c(
    sum = sum(values),
    sumsq = sum(values^2),
    indexed = sum(values * ((seq_along(values) %% 104729L) + 1L))
  )
}

#' Tune Design-MFA Components with Held-Out Subjects
#'
#' @description
#' Selects component rank, spatial smoothness, fixed partial-pooling penalties,
#' and soft-alignment stabilization using subject-disjoint folds. For each
#' fold, a candidate effect is learned from analysis subjects and compared with
#' the same design-effect matrix estimated independently from assessment
#' subjects.
#'
#' @details
#' The score is the mean, over estimable requested terms, of reference-weighted
#' squared disagreement divided by the held-out effect energy. Spatial
#' generalized components are reconstructed with their dual metric loadings.
#' Every candidate must succeed on every fold to be eligible.
#' Failed candidate-fold fits remain visible in `fold_scores` and are aggregated
#' in `failure_summary`; an otherwise successful tuning run warns when any
#' candidate was excluded by a fit failure.
#'
#' This procedure selects hyperparameters; its minimum score is not an unbiased
#' estimate of final performance. For an unbiased performance claim, run
#' `tune_design_mfa()` inside an outer subject-level resampling loop.
#'
#' @param Y A named list of subject trial-by-feature matrices.
#' @param design A matching list of trial-design data frames.
#' @param formula A one-sided design formula passed to [design_mfa()].
#' @param subject Optional unique subject identifiers.
#' @param terms Optional terms scored during tuning. The default scores every
#'   requested formula term.
#' @param grid Candidate data frame or named list with `rank`,
#'   `spatial_penalty`, `subject_penalty`, and `alignment_penalty` fields.
#' @param n_folds Number of subject-disjoint folds. Every assessment fold must
#'   contain at least two subjects.
#' @param seed Seed for deterministic named-subject fold construction. The
#'   caller's global random-number state is restored on successful and error
#'   exits.
#' @param refit Refit [design_mfa()] on all subjects with the selected values.
#' @param verbose Report failed candidate-fold fits as they occur.
#' @param ... Fixed arguments passed to [design_mfa()]. Tuned fields and
#'   `store` are controlled by this function.
#'
#' @return A `design_mfa_tuning` object containing candidate summaries,
#'   fold-level scores, retained failure diagnostics, subject membership,
#'   selected parameters, and the full refit when `refit = TRUE`.
#'
#' @examples
#' \donttest{
#' example_data <- design_mfa_example_data()
#' tuned <- tune_design_mfa(
#'   example_data$Y, example_data$design, ~ saliency,
#'   terms = "saliency", grid = data.frame(rank = c(1L, 2L)),
#'   n_folds = 3, seed = 11, refit = FALSE
#' )
#' tuned$best_params
#' }
#' @export
tune_design_mfa <- function(Y,
                            design,
                            formula,
                            subject = NULL,
                            terms = NULL,
                            grid = NULL,
                            n_folds = 3L,
                            seed = 1L,
                            refit = TRUE,
                            verbose = FALSE,
                            ...) {
  n_folds_defaulted <- missing(n_folds)
  dots <- list(...)
  forbidden <- intersect(
    names(dots),
    c("rank", "spatial_penalty", "subject_penalty", "alignment_penalty", "store", "terms")
  )
  if (length(forbidden) > 0L) {
    stop(
      sprintf("Argument(s) %s are controlled by `grid` or the tuning workflow.", paste(forbidden, collapse = ", ")),
      call. = FALSE
    )
  }
  if (!is.numeric(n_folds) || length(n_folds) != 1L || !is.finite(n_folds) ||
      n_folds < 2L || abs(n_folds - round(n_folds)) > sqrt(.Machine$double.eps)) {
    stop("`n_folds` must be an integer of at least two.", call. = FALSE)
  }
  if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed) ||
      abs(seed - round(seed)) > sqrt(.Machine$double.eps)) {
    stop("`seed` must be one integer.", call. = FALSE)
  }
  n_folds <- as.integer(round(n_folds))
  seed <- as.integer(round(seed))
  if (!is.logical(refit) || length(refit) != 1L || is.na(refit) ||
      !is.logical(verbose) || length(verbose) != 1L || is.na(verbose)) {
    stop("`refit` and `verbose` must each be TRUE or FALSE.", call. = FALSE)
  }

  dots <- .dmfa_canonicalize_tuning_args(dots, Y)
  alignment <- if (is.null(dots$alignment)) "shared" else match.arg(dots$alignment, c("shared", "soft"))
  validated <- .dmfa_validate_inputs(
    Y = Y,
    design = design,
    subject = subject,
    feature_names = dots$feature_names,
    alignment = alignment,
    feature_map = dots$feature_map
  )
  Y <- validated$Y
  design <- validated$design
  subject <- validated$subject
  n_subjects <- length(Y)
  if (n_subjects < 4L) stop("Held-out tuning requires at least four subjects.", call. = FALSE)
  maximum_folds <- floor(n_subjects / 2L)
  if (n_folds > maximum_folds) {
    if (isTRUE(n_folds_defaulted)) {
      n_folds <- maximum_folds
    } else {
      stop("Every assessment fold must contain at least two subjects; reduce `n_folds`.", call. = FALSE)
    }
  }

  dots$feature_names <- validated$feature_names
  dots$feature_map <- validated$feature_map
  subject_model <- if (is.null(dots$subject_model)) "fixed_intercept" else {
    match.arg(dots$subject_model, c("fixed_intercept", "ridge_slopes"))
  }
  grid <- .dmfa_tuning_grid(grid, dots$spatial_laplacian, subject_model, alignment)
  folds <- .dmfa_subject_folds(subject, n_folds, seed)
  index_by_subject <- stats::setNames(seq_along(subject), subject)
  fold_reference <- vector("list", n_folds)
  fold_target <- vector("list", n_folds)
  fold_target_error <- vector("list", n_folds)
  fit_subject_audit <- vector("list", n_folds)

  base_args <- c(list(formula = formula, terms = terms), dots)
  results <- cbind(
    data.frame(candidate = seq_len(nrow(grid)), grid, stringsAsFactors = FALSE),
    score = Inf,
    score_sd = NA_real_,
    score_se = NA_real_,
    n_success = 0L
  )
  fold_rows <- vector("list", nrow(grid) * n_folds)
  row_cursor <- 0L

  for (candidate in seq_len(nrow(grid))) {
    candidate_scores <- numeric(0)
    for (fold_id in seq_len(n_folds)) {
      row_cursor <- row_cursor + 1L
      fold <- folds[[fold_id]]
      analysis_index <- unname(index_by_subject[fold$analysis])
      assessment_index <- unname(index_by_subject[fold$assessment])
      train_args <- .dmfa_subset_tuning_args(
        base_args,
        analysis_index,
        Y,
        design,
        subject,
        n_subjects
      )
      train_args$rank <- grid$rank[[candidate]]
      train_args$spatial_penalty <- grid$spatial_penalty[[candidate]]
      train_args$subject_penalty <- grid$subject_penalty[[candidate]]
      train_args$alignment_penalty <- grid$alignment_penalty[[candidate]]
      train_args$store <- "fitted"
      if (!is.null(fold_reference[[fold_id]])) {
        train_args$reference_grid <- fold_reference[[fold_id]]$grid
        train_args$reference_weight <- fold_reference[[fold_id]]$weight
      }
      training <- tryCatch(do.call(design_mfa, train_args), error = identity)
      failure_message <- NULL
      failure_stage <- ""
      failure_reason <- ""
      failure_class <- ""
      analysis_fingerprint <- rep(NA_real_, 3L)
      analysis_design_fingerprint <- rep(NA_real_, 3L)
      analysis_rank <- NA_integer_
      score_info <- NULL
      if (inherits(training, "error")) {
        failure_message <- conditionMessage(training)
        failure_stage <- "analysis_fit"
        failure_reason <- .dmfa_tuning_failure_reason(failure_message)
        failure_class <- class(training)[[1L]]
      } else {
        analysis_fingerprint <- .dmfa_numeric_fingerprint(training$full_model$coefficient)
        analysis_design_fingerprint <- .dmfa_numeric_fingerprint(training$design$X)
        analysis_rank <- training$full_model$rank
        analysis_subjects <- .dmfa_assert_tuning_subjects(
          training$subjects, fold$analysis, "analysis", fold_id
        )
        if (is.null(fold_reference[[fold_id]])) {
          first_effect <- training$effects[[1L]]
          fold_reference[[fold_id]] <- list(
            grid = first_effect$reference_grid,
            weight = first_effect$reference_weights
          )
        }
        if (is.null(fold_target[[fold_id]]) && is.null(fold_target_error[[fold_id]])) {
          target_args <- .dmfa_subset_tuning_args(
            base_args,
            assessment_index,
            Y,
            design,
            subject,
            n_subjects
          )
          target_args$rank <- 1L
          target_args$subject_model <- "fixed_intercept"
          target_args$subject_penalty <- 1
          target_args$random_terms <- NULL
          target_args$spatial_penalty <- 0
          target_args$alignment_penalty <- 0
          target_args$store <- "fitted"
          target_args$reference_grid <- fold_reference[[fold_id]]$grid
          target_args$reference_weight <- fold_reference[[fold_id]]$weight
          target <- tryCatch(do.call(design_mfa, target_args), error = identity)
          if (inherits(target, "error")) {
            fold_target_error[[fold_id]] <- list(
              message = conditionMessage(target),
              condition_class = class(target)[[1L]]
            )
          } else {
            assessment_subjects <- .dmfa_assert_tuning_subjects(
              target$subjects, fold$assessment, "assessment", fold_id
            )
            fold_target[[fold_id]] <- target
            fit_subject_audit[[fold_id]] <- list(
              analysis = analysis_subjects,
              assessment = assessment_subjects,
              reference_grid = fold_reference[[fold_id]]$grid,
              reference_weights = fold_reference[[fold_id]]$weight,
              target_coefficient_fingerprint = .dmfa_numeric_fingerprint(
                target$full_model$coefficient
              )
            )
          }
        }
        if (!is.null(fold_target_error[[fold_id]])) {
          failure_message <- fold_target_error[[fold_id]]$message
          failure_stage <- "target_fit"
          failure_reason <- .dmfa_tuning_failure_reason(failure_message)
          failure_class <- fold_target_error[[fold_id]]$condition_class
        } else {
          scored_terms <- if (is.null(terms)) names(training$effects) else as.character(terms)
          score_info <- .dmfa_effect_holdout_score(
            training,
            fold_target[[fold_id]],
            scored_terms
          )
          if (!is.finite(score_info$score)) {
            failure_message <- "No requested term was estimable with nonzero held-out effect energy."
            failure_stage <- "scoring"
            failure_reason <- "zero_target_energy"
            failure_class <- "design_mfa_scoring_failure"
          } else {
            candidate_scores <- c(candidate_scores, score_info$score)
          }
        }
      }
      if (isTRUE(verbose) && !is.null(failure_message)) {
        message(sprintf("candidate %d fold %d failed: %s", candidate, fold_id, failure_message))
      }
      fold_rows[[row_cursor]] <- data.frame(
        candidate = candidate,
        fold = fold_id,
        score = if (is.null(score_info)) NA_real_ else score_info$score,
        n_terms = if (is.null(score_info)) 0L else score_info$n_terms,
        status = if (is.null(failure_message)) "ok" else "failed",
        stage = failure_stage,
        reason = failure_reason,
        condition_class = failure_class,
        message = if (is.null(failure_message)) "" else failure_message,
        analysis_rank = analysis_rank,
        analysis_coefficient_sum = analysis_fingerprint[[1L]],
        analysis_coefficient_sumsq = analysis_fingerprint[[2L]],
        analysis_coefficient_indexed = analysis_fingerprint[[3L]],
        analysis_design_sum = analysis_design_fingerprint[[1L]],
        analysis_design_sumsq = analysis_design_fingerprint[[2L]],
        analysis_design_indexed = analysis_design_fingerprint[[3L]],
        stringsAsFactors = FALSE
      )
    }
    results$n_success[[candidate]] <- length(candidate_scores)
    if (length(candidate_scores) == n_folds) {
      results$score[[candidate]] <- mean(candidate_scores)
      results$score_sd[[candidate]] <- if (n_folds > 1L) stats::sd(candidate_scores) else 0
      results$score_se[[candidate]] <- results$score_sd[[candidate]] / sqrt(n_folds)
    }
  }

  fold_scores <- do.call(rbind, fold_rows)
  results$n_failures <- n_folds - results$n_success
  failed_scores <- fold_scores[fold_scores$status == "failed", , drop = FALSE]
  if (nrow(failed_scores) > 0L) {
    category <- interaction(
      failed_scores$stage,
      failed_scores$reason,
      failed_scores$condition_class,
      drop = TRUE,
      lex.order = TRUE
    )
    groups <- split(seq_len(nrow(failed_scores)), category)
    failure_summary <- do.call(rbind, lapply(groups, function(rows) {
      data.frame(
        stage = failed_scores$stage[rows[[1L]]],
        reason = failed_scores$reason[rows[[1L]]],
        condition_class = failed_scores$condition_class[rows[[1L]]],
        message = failed_scores$message[rows[[1L]]],
        n_failures = length(rows),
        stringsAsFactors = FALSE
      )
    }))
    rownames(failure_summary) <- NULL
  } else {
    failure_summary <- data.frame(
      stage = character(0),
      reason = character(0),
      condition_class = character(0),
      message = character(0),
      n_failures = integer(0),
      stringsAsFactors = FALSE
    )
  }

  eligible <- results$n_success == n_folds & is.finite(results$score)
  if (!any(eligible)) {
    examples <- head(
      sprintf(
        "candidate %d, fold %d: %s",
        failed_scores$candidate,
        failed_scores$fold,
        gsub("[\r\n]+", " ", failed_scores$message)
      ),
      5L
    )
    detail <- if (length(examples) == 0L) {
      "No fold-level error was recorded."
    } else {
      paste(examples, collapse = "; ")
    }
    stop(
      sprintf(
        "No tuning candidate succeeded on every subject fold. Fold failures: %s",
        detail
      ),
      call. = FALSE
    )
  }
  ordering <- order(
    results$score,
    results$rank,
    -results$spatial_penalty,
    -results$subject_penalty,
    -results$alignment_penalty
  )
  best_index <- ordering[eligible[ordering]][[1L]]
  best_params <- list(
    rank = results$rank[[best_index]],
    spatial_penalty = results$spatial_penalty[[best_index]],
    subject_penalty = results$subject_penalty[[best_index]],
    alignment_penalty = results$alignment_penalty[[best_index]]
  )

  final_fit <- NULL
  if (isTRUE(refit)) {
    final_args <- c(
      list(Y = Y, design = design, formula = formula, subject = subject, terms = terms),
      dots
    )
    final_args$rank <- best_params$rank
    final_args$spatial_penalty <- best_params$spatial_penalty
    final_args$subject_penalty <- best_params$subject_penalty
    final_args$alignment_penalty <- best_params$alignment_penalty
    final_fit <- do.call(design_mfa, final_args)
  }

  out <- list(
    call = match.call(),
    metric = "heldout_relative_effect_mse",
    results = results,
    fold_scores = fold_scores,
    failure_summary = failure_summary,
    folds = folds,
    fit_subject_audit = fit_subject_audit,
    best_params = best_params,
    best_score = results$score[[best_index]],
    selected_candidate = best_index,
    fit = final_fit,
    seed = seed,
    n_folds = n_folds,
    terms = terms
  )
  class(out) <- c("design_mfa_tuning", "list")
  if (nrow(failed_scores) > 0L) {
    warning(
      sprintf(
        "%d of %d candidate-fold fits failed; inspect `$fold_scores` and `$failure_summary` for the retained diagnostics.",
        nrow(failed_scores),
        nrow(fold_scores)
      ),
      call. = FALSE
    )
  }
  out
}

#' @export
print.design_mfa_tuning <- function(x, ...) {
  cat("Design-effect Held-out-subject tuning\n")
  cat(sprintf("  Folds: %d\n", x$n_folds))
  cat(sprintf("  Candidates: %d\n", nrow(x$results)))
  cat(sprintf("  Best score: %s\n", format(x$best_score, digits = 4)))
  cat(sprintf(
    "  Selected: rank=%d, spatial_penalty=%s, subject_penalty=%s, alignment_penalty=%s\n",
    x$best_params$rank,
    format(x$best_params$spatial_penalty, digits = 4),
    format(x$best_params$subject_penalty, digits = 4),
    format(x$best_params$alignment_penalty, digits = 4)
  ))
  invisible(x)
}

#' Nested Held-Out-Subject Evaluation for Design-MFA Tuning
#'
#' @description
#' Runs [tune_design_mfa()] entirely within each outer analysis set and scores
#' the selected model against an independently estimated design effect in the
#' outer assessment subjects. The resulting outer scores evaluate the complete
#' tuning procedure without reusing an outer subject for candidate selection.
#'
#' @param Y A named list of subject trial-by-feature matrices.
#' @param design A matching list of trial-design data frames.
#' @param formula A one-sided design formula.
#' @param subject Optional unique subject identifiers.
#' @param terms Optional terms scored in inner and outer folds.
#' @param grid Candidate grid passed to [tune_design_mfa()].
#' @param outer_folds Number of outer subject folds.
#' @param inner_folds Number of inner subject folds. This is reduced when
#'   necessary to retain at least two assessment subjects per fold.
#' @param seed Seed controlling both named-subject fold levels. The caller's
#'   global random-number state is restored on return, including error exits.
#' @param refit Run a final tuning pass on all subjects for deployment.
#' @param verbose Report inner fitting failures.
#' @param ... Fixed [design_mfa()] arguments passed through both levels.
#'
#' @return A `nested_design_mfa_tuning` object with outer scores and folds,
#'   every inner tuning result, and an optional final full-data tuning/refit.
#'
#' @examples
#' \donttest{
#' example_data <- design_mfa_example_data()
#' nested <- nested_tune_design_mfa(
#'   example_data$Y, example_data$design, ~ saliency,
#'   terms = "saliency", grid = data.frame(rank = 1L),
#'   outer_folds = 3, inner_folds = 2, seed = 12, refit = FALSE
#' )
#' nested$outer_results
#' }
#' @export
nested_tune_design_mfa <- function(Y,
                                   design,
                                   formula,
                                   subject = NULL,
                                   terms = NULL,
                                   grid = NULL,
                                   outer_folds = 3L,
                                   inner_folds = 3L,
                                   seed = 1L,
                                   refit = TRUE,
                                   verbose = FALSE,
                                   ...) {
  outer_defaulted <- missing(outer_folds)
  dots <- list(...)
  valid_fold_count <- function(x) {
    is.numeric(x) && length(x) == 1L && is.finite(x) && x >= 2L &&
      abs(x - round(x)) <= sqrt(.Machine$double.eps)
  }
  if (!valid_fold_count(outer_folds) || !valid_fold_count(inner_folds)) {
    stop("`outer_folds` and `inner_folds` must be integers of at least two.", call. = FALSE)
  }
  if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed) ||
      abs(seed - round(seed)) > sqrt(.Machine$double.eps)) {
    stop("`seed` must be one integer.", call. = FALSE)
  }
  outer_folds <- as.integer(round(outer_folds))
  inner_folds <- as.integer(round(inner_folds))
  seed <- as.integer(round(seed))
  if (!is.logical(refit) || length(refit) != 1L || is.na(refit) ||
      !is.logical(verbose) || length(verbose) != 1L || is.na(verbose)) {
    stop("`refit` and `verbose` must each be TRUE or FALSE.", call. = FALSE)
  }

  dots <- .dmfa_canonicalize_tuning_args(dots, Y)
  alignment <- if (is.null(dots$alignment)) "shared" else match.arg(dots$alignment, c("shared", "soft"))
  validated <- .dmfa_validate_inputs(
    Y = Y,
    design = design,
    subject = subject,
    feature_names = dots$feature_names,
    alignment = alignment,
    feature_map = dots$feature_map
  )
  Y <- validated$Y
  design <- validated$design
  subject <- validated$subject
  dots$feature_names <- validated$feature_names
  dots$feature_map <- validated$feature_map
  n_subjects <- length(Y)
  if (n_subjects < 6L) {
    stop("Nested tuning requires at least six subjects.", call. = FALSE)
  }
  maximum_outer <- floor(n_subjects / 2L)
  if (outer_folds > maximum_outer) {
    if (isTRUE(outer_defaulted)) {
      outer_folds <- maximum_outer
    } else {
      stop("Every outer assessment fold must contain at least two subjects.", call. = FALSE)
    }
  }

  folds <- .dmfa_subject_folds(subject, outer_folds, seed)
  index_by_subject <- stats::setNames(seq_along(subject), subject)
  base_args <- c(list(formula = formula, terms = terms), dots)
  outer_tuning <- vector("list", outer_folds)
  outer_rows <- vector("list", outer_folds)
  outer_fit_subject_audit <- vector("list", outer_folds)

  for (outer_id in seq_len(outer_folds)) {
    fold <- folds[[outer_id]]
    analysis_index <- unname(index_by_subject[fold$analysis])
    assessment_index <- unname(index_by_subject[fold$assessment])
    if (length(analysis_index) < 4L) {
      stop("Each outer analysis set must retain at least four subjects for inner tuning.", call. = FALSE)
    }
    inner_folds_effective <- min(inner_folds, floor(length(analysis_index) / 2L))
    training_args <- .dmfa_subset_tuning_args(
      base_args,
      analysis_index,
      Y,
      design,
      subject,
      n_subjects
    )
    training_args$grid <- grid
    training_args$n_folds <- inner_folds_effective
    training_args$seed <- seed + outer_id
    training_args$refit <- TRUE
    training_args$verbose <- verbose
    tuned <- do.call(tune_design_mfa, training_args)
    outer_tuning[[outer_id]] <- tuned

    selected <- tuned$fit
    first_effect <- selected$effects[[1L]]
    target_args <- .dmfa_subset_tuning_args(
      base_args,
      assessment_index,
      Y,
      design,
      subject,
      n_subjects
    )
    target_args$rank <- 1L
    target_args$subject_model <- "fixed_intercept"
    target_args$subject_penalty <- 1
    target_args$random_terms <- NULL
    target_args$spatial_penalty <- 0
    target_args$alignment_penalty <- 0
    target_args$store <- "fitted"
    target_args$reference_grid <- first_effect$reference_grid
    target_args$reference_weight <- first_effect$reference_weights
    target <- do.call(design_mfa, target_args)
    analysis_subjects <- .dmfa_assert_tuning_subjects(
      selected$subjects, fold$analysis, "outer analysis", outer_id
    )
    assessment_subjects <- .dmfa_assert_tuning_subjects(
      target$subjects, fold$assessment, "outer assessment", outer_id
    )
    outer_fit_subject_audit[[outer_id]] <- list(
      analysis = analysis_subjects,
      assessment = assessment_subjects,
      reference_grid = first_effect$reference_grid,
      reference_weights = first_effect$reference_weights,
      analysis_coefficient_fingerprint = .dmfa_numeric_fingerprint(
        selected$full_model$coefficient
      ),
      analysis_design_fingerprint = .dmfa_numeric_fingerprint(selected$design$X),
      target_coefficient_fingerprint = .dmfa_numeric_fingerprint(
        target$full_model$coefficient
      )
    )
    scored_terms <- if (is.null(terms)) names(selected$effects) else as.character(terms)
    score <- .dmfa_effect_holdout_score(selected, target, scored_terms)
    if (!is.finite(score$score)) {
      stop(sprintf("Outer fold %d had no estimable nonzero target effect.", outer_id), call. = FALSE)
    }
    outer_rows[[outer_id]] <- data.frame(
      fold = outer_id,
      score = score$score,
      n_terms = score$n_terms,
      rank = tuned$best_params$rank,
      spatial_penalty = tuned$best_params$spatial_penalty,
      subject_penalty = tuned$best_params$subject_penalty,
      alignment_penalty = tuned$best_params$alignment_penalty,
      stringsAsFactors = FALSE
    )
  }
  outer_results <- do.call(rbind, outer_rows)

  final_tuning <- NULL
  if (isTRUE(refit)) {
    final_inner_folds <- min(inner_folds, floor(n_subjects / 2L))
    final_args <- c(
      list(
        Y = Y,
        design = design,
        formula = formula,
        subject = subject,
        terms = terms,
        grid = grid,
        n_folds = final_inner_folds,
        seed = seed + outer_folds + 1L,
        refit = TRUE,
        verbose = verbose
      ),
      dots
    )
    final_tuning <- do.call(tune_design_mfa, final_args)
  }

  out <- list(
    call = match.call(),
    metric = "outer_heldout_relative_effect_mse",
    outer_results = outer_results,
    outer_folds = folds,
    outer_tuning = outer_tuning,
    outer_fit_subject_audit = outer_fit_subject_audit,
    mean_score = mean(outer_results$score),
    score_sd = stats::sd(outer_results$score),
    score_se = stats::sd(outer_results$score) / sqrt(outer_folds),
    final_tuning = final_tuning,
    fit = if (is.null(final_tuning)) NULL else final_tuning$fit,
    seed = seed,
    n_outer_folds = outer_folds,
    inner_folds = inner_folds
  )
  class(out) <- c("nested_design_mfa_tuning", "list")
  out
}

#' @export
print.nested_design_mfa_tuning <- function(x, ...) {
  cat("Nested Design MFA tuning\n")
  cat(sprintf("  Outer folds: %d\n", x$n_outer_folds))
  cat(sprintf("  Mean outer score: %s\n", format(x$mean_score, digits = 4)))
  cat(sprintf("  Outer score SE: %s\n", format(x$score_se, digits = 4)))
  if (!is.null(x$final_tuning)) {
    cat(sprintf(
      "  Final selection: rank=%d, spatial_penalty=%s, subject_penalty=%s, alignment_penalty=%s\n",
      x$final_tuning$best_params$rank,
      format(x$final_tuning$best_params$spatial_penalty, digits = 4),
      format(x$final_tuning$best_params$subject_penalty, digits = 4),
      format(x$final_tuning$best_params$alignment_penalty, digits = 4)
    ))
  }
  invisible(x)
}

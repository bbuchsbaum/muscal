.dmfa_control <- function(rank_tol = sqrt(.Machine$double.eps),
                          estimability_tol = 1e-8,
                          max_numeric_grid = 20L,
                          max_grid_rows = 5000L,
                          feature_chunk_size = 25000L,
                          multilevel_max_iter = 5000L,
                          multilevel_tol = 1e-8,
                          soft_max_iter = 500L,
                          soft_tol = 1e-8) {
  rank_tol <- as.numeric(rank_tol)
  estimability_tol <- as.numeric(estimability_tol)
  max_numeric_grid <- as.integer(max_numeric_grid)
  max_grid_rows <- as.integer(max_grid_rows)
  feature_chunk_size <- as.integer(feature_chunk_size)
  multilevel_max_iter <- as.integer(multilevel_max_iter)
  multilevel_tol <- as.numeric(multilevel_tol)
  soft_max_iter <- as.integer(soft_max_iter)
  soft_tol <- as.numeric(soft_tol)

  if (length(rank_tol) != 1L || !is.finite(rank_tol) || rank_tol <= 0) {
    stop("`rank_tol` must be one positive finite number.", call. = FALSE)
  }
  if (length(estimability_tol) != 1L ||
      !is.finite(estimability_tol) || estimability_tol <= 0) {
    stop("`estimability_tol` must be one positive finite number.", call. = FALSE)
  }
  if (length(max_numeric_grid) != 1L || is.na(max_numeric_grid) || max_numeric_grid < 1L) {
    stop("`max_numeric_grid` must be a positive integer.", call. = FALSE)
  }
  if (length(max_grid_rows) != 1L || is.na(max_grid_rows) || max_grid_rows < 1L) {
    stop("`max_grid_rows` must be a positive integer.", call. = FALSE)
  }
  if (length(feature_chunk_size) != 1L || is.na(feature_chunk_size) || feature_chunk_size < 1L) {
    stop("`feature_chunk_size` must be a positive integer.", call. = FALSE)
  }
  if (length(multilevel_max_iter) != 1L || is.na(multilevel_max_iter) || multilevel_max_iter < 1L) {
    stop("`multilevel_max_iter` must be a positive integer.", call. = FALSE)
  }
  if (length(multilevel_tol) != 1L || !is.finite(multilevel_tol) || multilevel_tol <= 0) {
    stop("`multilevel_tol` must be one positive finite number.", call. = FALSE)
  }
  if (length(soft_max_iter) != 1L || is.na(soft_max_iter) || soft_max_iter < 1L) {
    stop("`soft_max_iter` must be a positive integer.", call. = FALSE)
  }
  if (length(soft_tol) != 1L || !is.finite(soft_tol) || soft_tol <= 0) {
    stop("`soft_tol` must be one positive finite number.", call. = FALSE)
  }

  list(
    rank_tol = rank_tol,
    estimability_tol = estimability_tol,
    max_numeric_grid = max_numeric_grid,
    max_grid_rows = max_grid_rows,
    feature_chunk_size = feature_chunk_size,
    multilevel_max_iter = multilevel_max_iter,
    multilevel_tol = multilevel_tol,
    soft_max_iter = soft_max_iter,
    soft_tol = soft_tol
  )
}

.dmfa_merge_control <- function(control) {
  defaults <- .dmfa_control()
  if (is.null(control)) return(defaults)
  if (!is.list(control) || is.null(names(control))) {
    stop("`control` must be a named list returned by `design_mfa_control()`.", call. = FALSE)
  }
  unknown <- setdiff(names(control), names(defaults))
  if (length(unknown) > 0L) {
    stop(sprintf("Unknown control field(s): %s.", paste(unknown, collapse = ", ")), call. = FALSE)
  }
  do.call(.dmfa_control, utils::modifyList(defaults, control))
}

.dmfa_validate_formula <- function(x, what) {
  if (!inherits(x, "formula") || length(x) != 2L) {
    stop(sprintf("`%s` must be a one-sided formula.", what), call. = FALSE)
  }
  x
}

.dmfa_validate_inputs <- function(Y,
                                  design,
                                  subject = NULL,
                                  feature_names = NULL,
                                  alignment = "shared",
                                  feature_map = NULL,
                                  feature_chunk_size = 25000L) {
  if (!is.list(Y) || length(Y) < 2L) {
    stop("`Y` must be a list containing at least two subject matrices.", call. = FALSE)
  }
  if (!is.list(design) || length(design) != length(Y)) {
    stop("`design` must be a list with one table per element of `Y`.", call. = FALSE)
  }

  if (is.null(names(Y))) names(Y) <- paste0("Subject_", seq_along(Y))
  if (anyNA(names(Y)) || any(names(Y) == "") || anyDuplicated(names(Y))) {
    stop("`Y` must have unique, non-empty subject names.", call. = FALSE)
  }
  if (is.null(names(design))) {
    names(design) <- names(Y)
  } else {
    if (anyNA(names(design)) || any(names(design) == "") || anyDuplicated(names(design))) {
      stop("`design` must have unique, non-empty subject names.", call. = FALSE)
    }
    missing_design <- setdiff(names(Y), names(design))
    extra_design <- setdiff(names(design), names(Y))
    if (length(missing_design) > 0L || length(extra_design) > 0L) {
      stop("Named `Y` and `design` lists must contain the same subjects.", call. = FALSE)
    }
    design <- design[names(Y)]
  }

  Y <- lapply(seq_along(Y), function(s) {
    ys <- Y[[s]]
    if (!is.matrix(ys) && !is.data.frame(ys)) {
      stop(sprintf("`Y[[%d]]` must be a numeric matrix or data frame.", s), call. = FALSE)
    }
    ys <- as.matrix(ys)
    if (!is.numeric(ys) || nrow(ys) < 1L || ncol(ys) < 1L) {
      stop(sprintf("`Y[[%d]]` must be a nonempty numeric matrix.", s), call. = FALSE)
    }
    has_nonfinite <- FALSE
    for (chunk_start in seq.int(1L, ncol(ys), by = feature_chunk_size)) {
      columns <- seq.int(chunk_start, min(ncol(ys), chunk_start + feature_chunk_size - 1L))
      if (any(!is.finite(ys[, columns, drop = FALSE]))) {
        has_nonfinite <- TRUE
        break
      }
    }
    if (has_nonfinite) {
      stop(sprintf("`Y[[%d]]` contains non-finite values.", s), call. = FALSE)
    }
    ys
  })
  names(Y) <- names(design)

  p_subject <- vapply(Y, ncol, integer(1))
  supplied_names <- if (is.null(feature_names)) NULL else as.character(feature_names)
  if (identical(alignment, "shared")) {
    p <- p_subject[[1L]]
    if (!all(p_subject == p)) {
      stop("`alignment = \"shared\"` requires the same number of features for every subject.", call. = FALSE)
    }
    if (!is.null(supplied_names) &&
        (length(supplied_names) != p || anyNA(supplied_names) ||
         any(supplied_names == "") || anyDuplicated(supplied_names))) {
      stop("`feature_names` must contain one unique, non-empty name per feature.", call. = FALSE)
    }

    block_colnames <- lapply(Y, colnames)
    has_names <- !vapply(block_colnames, is.null, logical(1))
    if (all(has_names)) {
      reference_names <- block_colnames[[1L]]
      if (anyNA(reference_names) || any(reference_names == "") || anyDuplicated(reference_names)) {
        stop("Feature column names must be unique, non-missing, and non-empty.", call. = FALSE)
      }
      for (s in seq_along(Y)) {
        if (!setequal(colnames(Y[[s]]), reference_names)) {
          stop("Shared-grid subjects must have the same named feature schema.", call. = FALSE)
        }
        if (!identical(colnames(Y[[s]]), reference_names)) {
          Y[[s]] <- Y[[s]][, reference_names, drop = FALSE]
        }
      }
      if (!is.null(supplied_names) && !identical(supplied_names, reference_names)) {
        stop("`feature_names` does not match the feature names stored in `Y`.", call. = FALSE)
      }
      feature_names <- reference_names
    } else if (any(has_names)) {
      stop("Either every `Y` block must have feature names or none may have them.", call. = FALSE)
    } else {
      if (is.null(supplied_names)) {
        stop(
          "Shared-grid feature identity is ambiguous. Supply identical column names or `feature_names`.",
          call. = FALSE
        )
      }
      feature_names <- supplied_names
      Y <- lapply(Y, function(ys) {
        colnames(ys) <- feature_names
        ys
      })
    }
    feature_map <- NULL
    alignment_coverage <- data.frame(
      feature = feature_names,
      subject_support = length(Y),
      coverage_fraction = 1,
      map_weight = length(Y),
      stringsAsFactors = FALSE
    )
  } else {
    if (!is.list(feature_map) || length(feature_map) != length(Y)) {
      stop("`alignment = \"soft\"` requires one sparse `feature_map` per subject.", call. = FALSE)
    }
    if (!is.null(names(feature_map))) {
      if (!setequal(names(feature_map), names(Y))) {
        stop("Named `feature_map` entries must match the subjects in `Y`.", call. = FALSE)
      }
      feature_map <- feature_map[names(Y)]
    }
    map_colnames <- lapply(feature_map, colnames)
    if (is.null(supplied_names)) {
      if (any(vapply(map_colnames, is.null, logical(1)))) {
        stop("Soft alignment requires `feature_names` or named reference columns on every map.", call. = FALSE)
      }
      supplied_names <- map_colnames[[1L]]
    }
    if (length(supplied_names) < 1L || anyNA(supplied_names) ||
        any(supplied_names == "") || anyDuplicated(supplied_names)) {
      stop("Soft-alignment `feature_names` must be unique, non-missing, and non-empty.", call. = FALSE)
    }
    p_reference <- length(supplied_names)
    row_coverage <- numeric(p_reference)
    subject_coverage <- numeric(p_reference)
    map_weight <- numeric(p_reference)
    for (s in seq_along(feature_map)) {
      P <- Matrix::Matrix(feature_map[[s]], sparse = TRUE)
      if (nrow(P) != p_subject[[s]] || ncol(P) != p_reference) {
        stop(sprintf("`feature_map[[%d]]` must be %d x %d.", s, p_subject[[s]], p_reference), call. = FALSE)
      }
      if (length(P@x) > 0L && (any(!is.finite(P@x)) || any(P@x < 0))) {
        stop(sprintf("`feature_map[[%d]]` must contain finite nonnegative weights.", s), call. = FALSE)
      }
      if (!is.null(colnames(P))) {
        if (!setequal(colnames(P), supplied_names)) {
          stop("Reference columns of every `feature_map` must match `feature_names`.", call. = FALSE)
        }
        P <- P[, supplied_names, drop = FALSE]
      } else {
        colnames(P) <- supplied_names
      }
      if (!is.null(rownames(P)) && !is.null(colnames(Y[[s]]))) {
        if (!setequal(rownames(P), colnames(Y[[s]]))) {
          stop(sprintf("Rows of `feature_map[[%d]]` do not match subject feature names.", s), call. = FALSE)
        }
        P <- P[colnames(Y[[s]]), , drop = FALSE]
      } else if (!is.null(rownames(P)) && is.null(colnames(Y[[s]]))) {
        colnames(Y[[s]]) <- rownames(P)
      } else if (is.null(colnames(Y[[s]]))) {
        colnames(Y[[s]]) <- paste0("subject_feature_", seq_len(p_subject[[s]]))
        rownames(P) <- colnames(Y[[s]])
      }
      row_total <- as.numeric(Matrix::rowSums(P))
      if (any(row_total <= 0)) {
        stop(sprintf("Every row of `feature_map[[%d]]` must have positive support.", s), call. = FALSE)
      }
      subject_rows <- as.numeric(Matrix::colSums(P != 0))
      row_coverage <- row_coverage + subject_rows
      subject_coverage <- subject_coverage + as.numeric(subject_rows > 0)
      map_weight <- map_weight + as.numeric(Matrix::colSums(P))
      feature_map[[s]] <- P
    }
    if (any(row_coverage == 0)) {
      stop("One or more reference features have no support in any `feature_map`.", call. = FALSE)
    }
    names(feature_map) <- names(Y)
    feature_names <- supplied_names
    alignment_coverage <- data.frame(
      feature = feature_names,
      subject_support = subject_coverage,
      coverage_fraction = subject_coverage / length(Y),
      map_weight = map_weight,
      stringsAsFactors = FALSE
    )
  }

  design <- lapply(seq_along(design), function(s) {
    ds <- design[[s]]
    if (!is.data.frame(ds)) ds <- as.data.frame(ds)
    if (nrow(ds) != nrow(Y[[s]])) {
      stop(sprintf("`design[[%d]]` must have one row per trial in `Y[[%d]]`.", s, s), call. = FALSE)
    }
    ds
  })
  names(design) <- names(Y)

  if (is.null(subject)) subject <- names(Y)
  subject <- as.character(subject)
  if (length(subject) != length(Y) || anyNA(subject) ||
      any(subject == "") || anyDuplicated(subject)) {
    stop("`subject` must contain one unique, non-empty identifier per `Y` block.", call. = FALSE)
  }

  list(
    Y = Y,
    design = design,
    subject = subject,
    feature_names = feature_names,
    feature_map = feature_map,
    subject_feature_counts = p_subject,
    alignment_coverage = alignment_coverage
  )
}

.dmfa_parse_weights <- function(Y,
                                trial_weights = NULL,
                                block_weight = c("equal_subject", "precision")) {
  S <- length(Y)
  n <- vapply(Y, nrow, integer(1))
  block_names <- names(Y)

  if (is.null(trial_weights)) {
    raw <- lapply(n, function(ni) rep(1, ni))
  } else if (is.list(trial_weights)) {
    if (length(trial_weights) != S) {
      stop("`trial_weights` must have one vector per subject.", call. = FALSE)
    }
    if (!is.null(names(trial_weights))) {
      if (!setequal(names(trial_weights), block_names)) {
        stop("Named `trial_weights` must match the subjects in `Y`.", call. = FALSE)
      }
      trial_weights <- trial_weights[block_names]
    }
    raw <- trial_weights
  } else {
    stop("`trial_weights` must be NULL or a list of subject-specific vectors.", call. = FALSE)
  }

  raw <- lapply(seq_len(S), function(s) {
    ws <- as.numeric(raw[[s]])
    if (length(ws) == 1L) ws <- rep(ws, n[[s]])
    if (length(ws) != n[[s]] || any(!is.finite(ws)) || any(ws < 0)) {
      stop(sprintf("Invalid trial weights for subject '%s'.", block_names[[s]]), call. = FALSE)
    }
    if (sum(ws) <= 0) {
      stop(sprintf("Subject '%s' must have positive total trial weight.", block_names[[s]]), call. = FALSE)
    }
    ws
  })
  names(raw) <- block_names

  if (is.character(block_weight)) {
    if (length(block_weight) > 1L) block_weight <- block_weight[[1L]]
    block_weight <- match.arg(block_weight, c("equal_subject", "precision"))
    if (identical(block_weight, "equal_subject")) {
      subject_weight <- rep(1 / S, S)
      final <- Map(function(w, sw) sw * w / sum(w), raw, subject_weight)
    } else {
      total <- sum(vapply(raw, sum, numeric(1)))
      final <- lapply(raw, function(w) w / total)
      subject_weight <- vapply(final, sum, numeric(1))
    }
  } else {
    if (!is.null(names(block_weight))) {
      if (anyNA(names(block_weight)) || any(names(block_weight) == "") ||
          anyDuplicated(names(block_weight)) ||
          !setequal(names(block_weight), block_names)) {
        stop("Named `block_weight` values must match the subjects in `Y`.", call. = FALSE)
      }
      block_weight <- block_weight[block_names]
    }
    subject_weight <- as.numeric(block_weight)
    if (length(subject_weight) != S || any(!is.finite(subject_weight)) ||
        any(subject_weight < 0) || sum(subject_weight) <= 0) {
      stop("Numeric `block_weight` must contain one nonnegative finite value per subject.", call. = FALSE)
    }
    subject_weight <- subject_weight / sum(subject_weight)
    final <- Map(function(w, sw) sw * w / sum(w), raw, subject_weight)
    block_weight <- "custom"
  }
  names(final) <- names(raw)
  names(subject_weight) <- names(raw)

  list(
    weights = final,
    stacked = unlist(final, use.names = FALSE),
    raw = raw,
    subject_weight = subject_weight,
    scheme = block_weight
  )
}

.dmfa_harmonize_design <- function(design, variables, subject) {
  for (var in variables) {
    missing_var <- !vapply(design, function(ds) var %in% names(ds), logical(1))
    if (any(missing_var)) {
      stop(sprintf("Design variable '%s' is missing for one or more subjects.", var), call. = FALSE)
    }
    vals <- lapply(design, `[[`, var)
    if (any(vapply(vals, is.character, logical(1)))) {
      stop(
        sprintf("Design variable '%s' is character. Convert it to a factor with explicit level order.", var),
        call. = FALSE
      )
    }
    is_factor <- vapply(vals, is.factor, logical(1))
    if (any(is_factor) && !all(is_factor)) {
      stop(sprintf("Design variable '%s' must have the same type for every subject.", var), call. = FALSE)
    }
    if (all(is_factor)) {
      reference_levels <- levels(vals[[1L]])
      reference_ordered <- is.ordered(vals[[1L]])
      compatible <- vapply(vals, function(x) {
        identical(levels(x), reference_levels) && identical(is.ordered(x), reference_ordered)
      }, logical(1))
      if (!all(compatible)) {
        stop(
          sprintf("Factor '%s' must declare identical levels in identical order for every subject.", var),
          call. = FALSE
        )
      }
    }
    valid_atomic <- vapply(vals, function(x) {
      is.factor(x) || is.numeric(x) || is.logical(x)
    }, logical(1))
    if (!all(valid_atomic)) {
      stop(sprintf("Unsupported design variable type for '%s'.", var), call. = FALSE)
    }
    if (any(vapply(vals, anyNA, logical(1)))) {
      stop(sprintf("Design variable '%s' contains missing values.", var), call. = FALSE)
    }
  }

  stacked <- do.call(rbind, Map(function(ds, sid) {
    ds$.subject <- sid
    ds
  }, design, subject))
  rownames(stacked) <- NULL
  stacked$.subject <- factor(stacked$.subject, levels = subject)
  stacked
}

.dmfa_resolve_contrasts <- function(contrasts, data) {
  if (is.null(contrasts)) return(NULL)
  if (!is.list(contrasts) || is.null(names(contrasts)) || any(names(contrasts) == "")) {
    stop("`contrasts` must be a named list.", call. = FALSE)
  }
  unknown <- setdiff(names(contrasts), names(data))
  if (length(unknown) > 0L) {
    stop(sprintf("Unknown contrast variable(s): %s.", paste(unknown, collapse = ", ")), call. = FALSE)
  }
  aliases <- c(
    poly = "contr.poly",
    sum = "contr.sum",
    treatment = "contr.treatment",
    helmert = "contr.helmert"
  )
  out <- lapply(names(contrasts), function(var) {
    value <- contrasts[[var]]
    if (!is.factor(data[[var]])) {
      stop(sprintf("Contrasts can only be assigned to factors; '%s' is not a factor.", var), call. = FALSE)
    }
    if (is.character(value) && length(value) == 1L && value %in% names(aliases)) {
      value <- unname(aliases[[value]])
    }
    if (!(is.character(value) || is.function(value) || is.matrix(value))) {
      stop(sprintf("Invalid contrast specification for '%s'.", var), call. = FALSE)
    }
    value
  })
  stats::setNames(out, names(contrasts))
}

.dmfa_drop_intercept <- function(X) {
  assignment <- attr(X, "assign")
  keep <- assignment != 0L
  list(X = X[, keep, drop = FALSE], assign = assignment[keep], keep = keep)
}

.dmfa_make_reference_grid <- function(stacked,
                                      formula,
                                      xlev,
                                      reference_grid,
                                      control) {
  variables <- all.vars(formula)
  if (!is.null(reference_grid)) {
    ref <- as.data.frame(reference_grid)
    missing_var <- setdiff(variables, names(ref))
    if (length(missing_var) > 0L) {
      stop(
        sprintf("`reference_grid` is missing variable(s): %s.", paste(missing_var, collapse = ", ")),
        call. = FALSE
      )
    }
    ref <- ref[, variables, drop = FALSE]
  } else {
    values <- lapply(variables, function(var) {
      x <- stacked[[var]]
      if (is.factor(x)) {
        declared <- if (!is.null(xlev) && !is.null(xlev[[var]])) xlev[[var]] else levels(x)
        return(declared)
      }
      ux <- sort(unique(x))
      if (is.numeric(x) && length(ux) > control$max_numeric_grid) {
        probs <- seq(0, 1, length.out = control$max_numeric_grid)
        return(unique(as.numeric(stats::quantile(x, probs, names = FALSE, type = 8))))
      }
      ux
    })
    names(values) <- variables
    ref <- expand.grid(values, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  }

  if (nrow(ref) < 1L || nrow(ref) > control$max_grid_rows) {
    stop(
      sprintf("Reference grid has %d rows; the allowed range is 1:%d.", nrow(ref), control$max_grid_rows),
      call. = FALSE
    )
  }

  for (var in variables) {
    train <- stacked[[var]]
    if (is.factor(train)) {
      declared <- if (!is.null(xlev) && !is.null(xlev[[var]])) xlev[[var]] else levels(train)
      ref[[var]] <- factor(ref[[var]], levels = declared, ordered = is.ordered(train))
      if (anyNA(ref[[var]])) {
        stop(sprintf("`reference_grid` contains an unknown level for '%s'.", var), call. = FALSE)
      }
    } else if (is.numeric(train)) {
      ref[[var]] <- as.numeric(ref[[var]])
      if (any(!is.finite(ref[[var]]))) {
        stop(sprintf("`reference_grid` contains non-finite values for '%s'.", var), call. = FALSE)
      }
    } else if (is.logical(train)) {
      ref[[var]] <- as.logical(ref[[var]])
      if (anyNA(ref[[var]])) {
        stop(sprintf("`reference_grid` contains invalid values for '%s'.", var), call. = FALSE)
      }
    }
  }
  rownames(ref) <- NULL
  ref
}

.dmfa_row_key <- function(data, variables) {
  pieces <- lapply(variables, function(var) {
    x <- data[[var]]
    if (is.numeric(x)) format(x, digits = 17L, scientific = FALSE, trim = TRUE) else as.character(x)
  })
  do.call(paste, c(pieces, sep = "\r"))
}

.dmfa_reference_weights <- function(reference_weight,
                                    reference_grid,
                                    stacked,
                                    variables,
                                    weights) {
  m <- nrow(reference_grid)
  if (is.numeric(reference_weight)) {
    out <- as.numeric(reference_weight)
    if (length(out) != m || any(!is.finite(out)) || any(out < 0) || sum(out) <= 0) {
      stop("Numeric `reference_weight` must contain one nonnegative value per grid row.", call. = FALSE)
    }
    return(out / sum(out))
  }
  reference_weight <- match.arg(reference_weight, c("equal_cells", "observed"))
  if (identical(reference_weight, "equal_cells")) return(rep(1 / m, m))

  ref_key <- .dmfa_row_key(reference_grid, variables)
  obs_key <- .dmfa_row_key(stacked, variables)
  matched <- match(obs_key, ref_key)
  if (anyNA(matched)) {
    stop("Observed design rows do not map exactly to `reference_grid`; use explicit numeric weights.", call. = FALSE)
  }
  out <- numeric(m)
  for (i in seq_along(matched)) out[[matched[[i]]]] <- out[[matched[[i]]]] + weights[[i]]
  if (sum(out) <= 0) stop("Observed reference-grid weight is zero.", call. = FALSE)
  out / sum(out)
}

.dmfa_compile_design <- function(design,
                                 subject,
                                 formula,
                                 nuisance,
                                 contrasts,
                                 weights,
                                 terms,
                                 reference_grid,
                                 reference_weight,
                                 control) {
  formula <- .dmfa_validate_formula(formula, "formula")
  if (!is.null(nuisance)) nuisance <- .dmfa_validate_formula(nuisance, "nuisance")
  effect_variables <- all.vars(formula)
  nuisance_variables <- if (is.null(nuisance)) character(0) else all.vars(nuisance)
  variables <- unique(c(effect_variables, nuisance_variables))
  if (length(effect_variables) == 0L) {
    stop("`formula` must contain at least one design variable.", call. = FALSE)
  }

  stacked <- .dmfa_harmonize_design(design, variables, subject)
  contrast_arg <- .dmfa_resolve_contrasts(contrasts, stacked)
  effect_contrast_arg <- contrast_arg[intersect(names(contrast_arg), effect_variables)]
  nuisance_contrast_arg <- contrast_arg[intersect(names(contrast_arg), nuisance_variables)]
  if (length(effect_contrast_arg) == 0L) effect_contrast_arg <- NULL
  if (length(nuisance_contrast_arg) == 0L) nuisance_contrast_arg <- NULL

  effect_frame <- stats::model.frame(
    formula,
    data = stacked,
    na.action = stats::na.fail,
    drop.unused.levels = FALSE
  )
  # Use the terms object attached to the training model frame.  It carries
  # predvars (for poly(), scale(), splines, and similar stateful transforms)
  # into reference-grid construction and prediction.
  effect_terms <- attr(effect_frame, "terms")
  effect_matrix_raw <- stats::model.matrix(effect_terms, effect_frame, contrasts.arg = effect_contrast_arg)
  effect_contrasts <- attr(effect_matrix_raw, "contrasts")
  effect_xlev <- stats::.getXlevels(effect_terms, effect_frame)
  effect_drop <- .dmfa_drop_intercept(effect_matrix_raw)
  effect_matrix <- effect_drop$X
  effect_assign <- effect_drop$assign
  term_labels <- attr(effect_terms, "term.labels")
  if (length(term_labels) == 0L || ncol(effect_matrix) == 0L) {
    stop("`formula` must contain at least one non-intercept term.", call. = FALSE)
  }

  if (!is.null(nuisance)) {
    nuisance_frame <- stats::model.frame(
      nuisance,
      data = stacked,
      na.action = stats::na.fail,
      drop.unused.levels = FALSE
    )
    nuisance_terms <- attr(nuisance_frame, "terms")
    duplicate_terms <- intersect(term_labels, attr(nuisance_terms, "term.labels"))
    if (length(duplicate_terms) > 0L) {
      stop(
        sprintf("Term(s) cannot appear in both `formula` and `nuisance`: %s.", paste(duplicate_terms, collapse = ", ")),
        call. = FALSE
      )
    }
    nuisance_raw <- stats::model.matrix(nuisance_terms, nuisance_frame, contrasts.arg = nuisance_contrast_arg)
    nuisance_drop <- .dmfa_drop_intercept(nuisance_raw)
    nuisance_matrix <- nuisance_drop$X
    nuisance_contrasts <- attr(nuisance_raw, "contrasts")
    nuisance_xlev <- stats::.getXlevels(nuisance_terms, nuisance_frame)
  } else {
    nuisance_terms <- NULL
    nuisance_matrix <- matrix(0, nrow(stacked), 0L)
    nuisance_contrasts <- NULL
    nuisance_xlev <- NULL
  }

  subject_matrix <- stats::model.matrix(~ 0 + .subject, data = stacked)
  colnames(subject_matrix) <- paste0("subject::", subject)
  if (ncol(nuisance_matrix) > 0L) {
    colnames(nuisance_matrix) <- paste0("nuisance::", colnames(nuisance_matrix))
  }
  colnames(effect_matrix) <- paste0("effect::", colnames(effect_matrix))
  X <- cbind(subject_matrix, nuisance_matrix, effect_matrix)
  if (anyDuplicated(colnames(X))) colnames(X) <- make.unique(colnames(X))

  base_count <- ncol(subject_matrix) + ncol(nuisance_matrix)
  term_columns <- lapply(seq_along(term_labels), function(j) {
    base_count + which(effect_assign == j)
  })
  names(term_columns) <- term_labels

  if (is.null(terms)) {
    selected_terms <- term_labels
  } else {
    selected_terms <- as.character(terms)
    unknown <- setdiff(selected_terms, term_labels)
    if (length(unknown) > 0L) {
      stop(sprintf("Unknown requested term(s): %s.", paste(unknown, collapse = ", ")), call. = FALSE)
    }
    selected_terms <- unique(selected_terms)
  }

  factors <- attr(effect_terms, "factors")
  term_order <- attr(effect_terms, "order")
  descendants <- lapply(seq_along(term_labels), function(j) {
    focal <- factors[, j] > 0
    which(vapply(seq_along(term_labels), function(k) {
      all(!focal | factors[, k] > 0) && term_order[[k]] >= term_order[[j]]
    }, logical(1)))
  })
  names(descendants) <- term_labels

  reference_grid <- .dmfa_make_reference_grid(
    stacked = stacked,
    formula = formula,
    xlev = effect_xlev,
    reference_grid = reference_grid,
    control = control
  )
  reference_frame <- stats::model.frame(
    effect_terms,
    data = reference_grid,
    xlev = effect_xlev,
    na.action = stats::na.fail
  )
  effect_reference_raw <- stats::model.matrix(
    effect_terms,
    reference_frame,
    contrasts.arg = effect_contrasts
  )
  effect_reference <- effect_reference_raw[, effect_drop$keep, drop = FALSE]
  colnames(effect_reference) <- colnames(effect_matrix)

  subject_reference <- matrix(
    rep(weights$subject_weight, each = nrow(reference_grid)),
    nrow = nrow(reference_grid),
    byrow = FALSE
  )
  colnames(subject_reference) <- colnames(subject_matrix)
  if (ncol(nuisance_matrix) > 0L) {
    nuisance_mean <- colSums(nuisance_matrix * weights$stacked) / sum(weights$stacked)
    nuisance_reference <- matrix(
      rep(nuisance_mean, each = nrow(reference_grid)),
      nrow = nrow(reference_grid),
      byrow = FALSE
    )
    colnames(nuisance_reference) <- colnames(nuisance_matrix)
  } else {
    nuisance_mean <- numeric(0)
    nuisance_reference <- matrix(0, nrow(reference_grid), 0L)
  }
  X_reference <- cbind(subject_reference, nuisance_reference, effect_reference)
  colnames(X_reference) <- colnames(X)

  reference_weights <- .dmfa_reference_weights(
    reference_weight = reference_weight,
    reference_grid = reference_grid,
    stacked = stacked,
    variables = effect_variables,
    weights = weights$stacked
  )

  list(
    X = X,
    X_reference = X_reference,
    stacked_design = stacked,
    reference_grid = reference_grid,
    reference_weights = reference_weights,
    effect_terms = effect_terms,
    effect_contrasts = effect_contrasts,
    effect_xlev = effect_xlev,
    effect_keep = effect_drop$keep,
    nuisance_terms = nuisance_terms,
    nuisance_contrasts = nuisance_contrasts,
    nuisance_xlev = nuisance_xlev,
    nuisance_mean = nuisance_mean,
    subject_columns = seq_len(ncol(subject_matrix)),
    nuisance_columns = if (ncol(nuisance_matrix) == 0L) integer(0) else ncol(subject_matrix) + seq_len(ncol(nuisance_matrix)),
    term_labels = term_labels,
    term_columns = term_columns,
    descendants = descendants,
    selected_terms = selected_terms,
    effect_variables = effect_variables,
    formula = formula,
    nuisance = nuisance
  )
}

.dmfa_term_indices <- function(compiled,
                               term,
                               ss_type = c("partial", "sequential"),
                               term_scope = c("assigned", "hierarchical")) {
  ss_type <- match.arg(ss_type)
  term_scope <- match.arg(term_scope)
  term_id <- match(term, compiled$term_labels)
  base <- c(compiled$subject_columns, compiled$nuisance_columns)

  if (identical(ss_type, "sequential")) {
    previous <- if (term_id <= 1L) integer(0) else seq_len(term_id - 1L)
    reduced_terms <- previous
    full_terms <- seq_len(term_id)
    return(list(
      full = c(base, unlist(compiled$term_columns[full_terms], use.names = FALSE)),
      reduced = c(base, unlist(compiled$term_columns[reduced_terms], use.names = FALSE)),
      removed_terms = term
    ))
  }

  remove_ids <- if (identical(term_scope, "hierarchical")) {
    compiled$descendants[[term]]
  } else {
    term_id
  }
  remove_cols <- unlist(compiled$term_columns[remove_ids], use.names = FALSE)
  list(
    full = seq_len(ncol(compiled$X)),
    reduced = setdiff(seq_len(ncol(compiled$X)), remove_cols),
    removed_terms = compiled$term_labels[remove_ids]
  )
}

.dmfa_build_effect_matrix <- function(compiled, new_design) {
  new_design <- as.data.frame(new_design)
  missing_var <- setdiff(compiled$effect_variables, names(new_design))
  if (length(missing_var) > 0L) {
    stop(sprintf("`new_design` is missing variable(s): %s.", paste(missing_var, collapse = ", ")), call. = FALSE)
  }
  frame <- stats::model.frame(
    compiled$effect_terms,
    data = new_design,
    xlev = compiled$effect_xlev,
    na.action = stats::na.fail
  )
  raw <- stats::model.matrix(
    compiled$effect_terms,
    frame,
    contrasts.arg = compiled$effect_contrasts
  )
  out <- raw[, compiled$effect_keep, drop = FALSE]
  colnames(out) <- colnames(compiled$X)[unlist(compiled$term_columns, use.names = FALSE)]
  out
}

.dmfa_validate_spatial_metric <- function(spatial_laplacian,
                                          spatial_penalty,
                                          feature_names,
                                          tol = 1e-8) {
  if (length(spatial_penalty) != 1L || !is.finite(spatial_penalty) || spatial_penalty < 0) {
    stop("`spatial_penalty` must be one nonnegative finite number.", call. = FALSE)
  }
  if (spatial_penalty == 0) return(list(laplacian = NULL, metric = NULL))
  if (is.null(spatial_laplacian)) {
    stop("A positive `spatial_penalty` requires `spatial_laplacian`.", call. = FALSE)
  }
  L <- Matrix::Matrix(spatial_laplacian, sparse = TRUE)
  p <- length(feature_names)
  if (nrow(L) != p || ncol(L) != p) {
    stop(sprintf("`spatial_laplacian` must be a %d x %d matrix.", p, p), call. = FALSE)
  }
  if (length(L@x) > 0L && any(!is.finite(L@x))) {
    stop("`spatial_laplacian` contains non-finite values.", call. = FALSE)
  }
  if (!is.null(rownames(L)) || !is.null(colnames(L))) {
    if (is.null(rownames(L)) || is.null(colnames(L)) ||
        !setequal(rownames(L), feature_names) || !setequal(colnames(L), feature_names)) {
      stop("Named `spatial_laplacian` dimensions must match `feature_names`.", call. = FALSE)
    }
    L <- L[feature_names, feature_names, drop = FALSE]
  }
  asymmetry <- L - Matrix::t(L)
  if (length(asymmetry@x) > 0L && max(abs(asymmetry@x)) > tol) {
    stop("`spatial_laplacian` must be symmetric.", call. = FALSE)
  }
  L <- Matrix::forceSymmetric((L + Matrix::t(L)) / 2)
  probes <- cbind(
    rep(1, p),
    seq(-1, 1, length.out = p),
    sin(seq(0, 2 * pi, length.out = p))
  )
  quadratic <- vapply(seq_len(ncol(probes)), function(j) {
    sum(probes[, j] * as.numeric(L %*% probes[, j]))
  }, numeric(1))
  scale <- max(1, max(abs(quadratic)))
  if (any(quadratic < -tol * scale)) {
    stop("`spatial_laplacian` is not positive semidefinite on validation probes.", call. = FALSE)
  }
  spectral_scale <- max(1, max(abs(Matrix::diag(L))))
  diagonal <- as.numeric(Matrix::diag(L))
  off_diagonal <- L - Matrix::Diagonal(x = diagonal)
  off_values <- off_diagonal@x
  graph_laplacian_certificate <-
    all(diagonal >= -tol * spectral_scale) &&
    (length(off_values) == 0L || all(off_values <= tol * spectral_scale)) &&
    all(diagonal + tol * spectral_scale >= as.numeric(Matrix::rowSums(abs(off_diagonal))))
  minimum_eigenvalue <- if (graph_laplacian_certificate) {
    0
  } else if (p <= 300L) {
    min(eigen(as.matrix(L), symmetric = TRUE, only.values = TRUE)$values)
  } else {
    tryCatch(
      as.numeric(RSpectra::eigs_sym(
        L,
        k = 1L,
        which = "SA",
        opts = list(retvec = FALSE, maxitr = max(1000L, 5L * p))
      )$values),
      error = function(e) NA_real_
    )
  }
  if (!is.finite(minimum_eigenvalue)) {
    stop("Could not verify that `spatial_laplacian` is positive semidefinite.", call. = FALSE)
  }
  if (minimum_eigenvalue < -tol * spectral_scale) {
    stop("`spatial_laplacian` must be positive semidefinite.", call. = FALSE)
  }
  Q <- Matrix::Diagonal(p) + spatial_penalty * L
  dimnames(Q) <- list(feature_names, feature_names)
  list(laplacian = L, metric = Q)
}

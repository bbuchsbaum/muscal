.dmfa_wls_from_stats <- function(XtX,
                                 XtY,
                                 response_ss,
                                 coefficient_names,
                                 feature_names,
                                 rank_tol,
                                 sse = NA_real_) {
  column_scales <- sqrt(pmax(0, diag(XtX)))
  column_scales[!is.finite(column_scales) | column_scales == 0] <- 1
  inverse_scales <- 1 / column_scales
  XtX_equilibrated <- XtX * tcrossprod(inverse_scales)
  XtY_equilibrated <- XtY * inverse_scales

  eig <- eigen(XtX_equilibrated, symmetric = TRUE)
  values <- pmax(eig$values, 0)
  threshold <- rank_tol * max(1, values[[1L]])
  keep <- which(values > threshold)
  if (length(keep) == 0L) {
    coefficient <- matrix(0, nrow(XtX), ncol(XtY))
    row_projector <- matrix(0, nrow(XtX), nrow(XtX))
    condition <- Inf
  } else {
    basis <- eig$vectors[, keep, drop = FALSE]
    rotated <- crossprod(basis, XtY_equilibrated)
    rotated <- rotated / values[keep]
    coefficient <- (basis %*% rotated) * inverse_scales
    row_projector <- tcrossprod(basis)
    condition <- max(values[keep]) / min(values[keep])
  }
  dimnames(coefficient) <- list(coefficient_names, feature_names)

  if (length(sse) == 1L && is.na(sse)) {
    fitted_ss <- NA_real_
  } else {
    sse <- as.numeric(sse)
    if (length(sse) != 1L || !is.finite(sse) || sse < 0) {
      stop(
        "Internal residual SSE must be supplied directly as one nonnegative finite number, or as NA while residuals are pending.",
        call. = FALSE
      )
    }
    fitted_ss <- response_ss - sse
  }

  list(
    coefficient = coefficient,
    rank = length(keep),
    condition = condition,
    row_projector = row_projector,
    column_scales = column_scales,
    XtX = XtX,
    XtX_equilibrated = XtX_equilibrated,
    XtY = XtY,
    response_ss = response_ss,
    fitted_ss = fitted_ss,
    sse = sse,
    rank_threshold = threshold
  )
}

.dmfa_compact_fitted_model <- function(model) {
  model$XtY <- NULL
  model$XtX_equilibrated <- NULL
  model
}

.dmfa_centered_response_ss <- function(Y, weights, groups = NULL) {
  if (is.null(groups)) return(sum(Y^2 * weights))
  if (length(groups) != nrow(Y)) {
    stop("Internal response-centering dimension mismatch.", call. = FALSE)
  }
  total <- 0
  for (rows in split(seq_len(nrow(Y)), groups)) {
    ws <- weights[rows]
    Ys <- Y[rows, , drop = FALSE]
    center <- colSums(Ys * ws) / sum(ws)
    total <- total + sum(sweep(Ys, 2L, center, "-")^2 * ws)
  }
  total
}

.dmfa_wls_sufficient <- function(X,
                                 Y,
                                 weights,
                                 rank_tol,
                                 subject_rows = NULL) {
  if (nrow(X) != nrow(Y) || length(weights) != nrow(X)) {
    stop("Internal weighted-model dimension mismatch.", call. = FALSE)
  }
  Xw <- X * sqrt(weights)
  groups <- NULL
  if (!is.null(subject_rows)) {
    groups <- integer(nrow(X))
    for (s in seq_along(subject_rows)) groups[subject_rows[[s]]] <- s
    if (any(groups == 0L)) stop("Internal subject rows do not cover the stacked response.", call. = FALSE)
  }
  out <- .dmfa_wls_from_stats(
    XtX = crossprod(Xw),
    XtY = crossprod(X, Y * weights),
    response_ss = .dmfa_centered_response_ss(Y, weights, groups),
    coefficient_names = colnames(X),
    feature_names = colnames(Y),
    rank_tol = rank_tol,
    sse = NA_real_
  )
  residual <- Y - X %*% out$coefficient
  out$sse <- sum(residual^2 * weights)
  out$fitted_ss <- out$response_ss - out$sse
  out
}

.dmfa_wls_list_sufficient <- function(X,
                                      Y_list,
                                      weight_list,
                                      subject_rows,
                                      rank_tol,
                                      feature_chunk_size = 25000L) {
  if (length(Y_list) == 0L || length(Y_list) != length(weight_list) ||
      length(Y_list) != length(subject_rows)) {
    stop("Internal subject-wise weighted-model dimension mismatch.", call. = FALSE)
  }
  p <- ncol(Y_list[[1L]])
  feature_names <- colnames(Y_list[[1L]])
  XtX <- matrix(0, ncol(X), ncol(X))
  XtY <- matrix(0, ncol(X), p)
  response_ss <- 0

  for (s in seq_along(Y_list)) {
    rows <- subject_rows[[s]]
    Ys <- Y_list[[s]]
    ws <- weight_list[[s]]
    if (length(rows) != nrow(Ys) || length(ws) != nrow(Ys) ||
        ncol(Ys) != p) {
      stop("Internal subject-wise weighted-model dimension mismatch.", call. = FALSE)
    }
    Xs <- X[rows, , drop = FALSE]
    XtX <- XtX + crossprod(Xs * sqrt(ws))
    chunk_starts <- seq.int(1L, p, by = feature_chunk_size)
    for (chunk_start in chunk_starts) {
      columns <- seq.int(chunk_start, min(p, chunk_start + feature_chunk_size - 1L))
      Ys_chunk <- Ys[, columns, drop = FALSE]
      XtY[, columns] <- XtY[, columns, drop = FALSE] + crossprod(Xs, Ys_chunk * ws)
      center <- colSums(Ys_chunk * ws) / sum(ws)
      response_ss <- response_ss + sum(sweep(Ys_chunk, 2L, center, "-")^2 * ws)
    }
  }
  dimnames(XtX) <- list(colnames(X), colnames(X))
  dimnames(XtY) <- list(colnames(X), feature_names)

  out <- .dmfa_wls_from_stats(
    XtX = XtX,
    XtY = XtY,
    response_ss = response_ss,
    coefficient_names = colnames(X),
    feature_names = feature_names,
    rank_tol = rank_tol,
    sse = NA_real_
  )
  sse <- 0
  for (s in seq_along(Y_list)) {
    rows <- subject_rows[[s]]
    Xs <- X[rows, , drop = FALSE]
    Ys <- Y_list[[s]]
    ws <- weight_list[[s]]
    for (chunk_start in seq.int(1L, p, by = feature_chunk_size)) {
      columns <- seq.int(chunk_start, min(p, chunk_start + feature_chunk_size - 1L))
      residual <- Ys[, columns, drop = FALSE] -
        Xs %*% out$coefficient[, columns, drop = FALSE]
      sse <- sse + sum(residual^2 * ws)
    }
  }
  out$sse <- sse
  out$fitted_ss <- out$response_ss - out$sse
  out
}

.dmfa_fit_shared_streamed <- function(X,
                                      Y_list,
                                      weight_list,
                                      subject_rows,
                                      control) {
  out <- .dmfa_wls_list_sufficient(
    X = X,
    Y_list = Y_list,
    weight_list = weight_list,
    subject_rows = subject_rows,
    rank_tol = control$rank_tol,
    feature_chunk_size = control$feature_chunk_size
  )
  out$penalized <- FALSE
  out$objective <- out$sse
  out$penalty_value <- 0
  out$subject_deviations <- NULL
  out$converged <- TRUE
  out$iterations <- 0L
  out$computation <- "subject_streamed_sufficient"
  out
}

.dmfa_estimable_rows <- function(L, model, tol) {
  if (ncol(L) != nrow(model$row_projector)) {
    stop("Internal estimability dimension mismatch.", call. = FALSE)
  }
  column_scales <- model$column_scales
  if (is.null(column_scales)) column_scales <- rep(1, ncol(L))
  if (length(column_scales) != ncol(L)) {
    stop("Internal estimability scaling mismatch.", call. = FALSE)
  }
  L_equilibrated <- sweep(L, 2L, column_scales, "/")
  residual <- L_equilibrated - L_equilibrated %*% model$row_projector
  residual_norm <- sqrt(rowSums(residual^2))
  scale <- pmax(1, sqrt(rowSums(L_equilibrated^2)))
  list(
    estimable = residual_norm <= tol * scale,
    residual_norm = residual_norm,
    relative_residual = residual_norm / scale
  )
}

.dmfa_ridge_solve <- function(G, C, ridge, rank_tol) {
  A <- G + diag(ridge, nrow(G))
  eig <- eigen(A, symmetric = TRUE)
  values <- pmax(eig$values, 0)
  threshold <- rank_tol * max(1, values[[1L]])
  keep <- which(values > threshold)
  if (length(keep) == 0L) return(matrix(0, nrow(G), ncol(C)))
  basis <- eig$vectors[, keep, drop = FALSE]
  basis %*% (crossprod(basis, C) / values[keep])
}

.dmfa_fit_partial_pool <- function(X,
                                   Y,
                                   weights,
                                   subject_rows,
                                   subject_weights,
                                   random_columns,
                                   subject_penalty,
                                   control) {
  group <- .dmfa_wls_sufficient(
    X, Y, weights, control$rank_tol, subject_rows = subject_rows
  )
  q_random <- length(random_columns)
  subject_names <- names(subject_rows)
  deviations <- lapply(subject_rows, function(rows) matrix(0, q_random, ncol(Y)))
  names(deviations) <- subject_names
  if (q_random > 0L) {
    for (s in seq_along(deviations)) {
      rownames(deviations[[s]]) <- colnames(X)[random_columns]
      colnames(deviations[[s]]) <- colnames(Y)
    }
  }

  objective_trace <- group$sse
  if (q_random == 0L) {
    group$penalized <- TRUE
    group$penalty_value <- 0
    group$objective <- group$sse
    group$objective_trace <- objective_trace
    group$subject_deviations <- deviations
    group$converged <- TRUE
    group$iterations <- 0L
    group$final_relative_change <- 0
    return(group)
  }

  previous <- objective_trace[[1L]]
  converged <- FALSE
  for (iteration in seq_len(control$multilevel_max_iter)) {
    residual_group <- Y - X %*% group$coefficient
    uncentered <- lapply(seq_along(subject_rows), function(s) {
      rows <- subject_rows[[s]]
      E <- X[rows, random_columns, drop = FALSE]
      ws <- weights[rows]
      G <- crossprod(E * sqrt(ws))
      C <- crossprod(E, residual_group[rows, , drop = FALSE] * ws)
      .dmfa_ridge_solve(
        G,
        C,
        ridge = subject_penalty * subject_weights[[s]],
        rank_tol = control$rank_tol
      )
    })

    mean_deviation <- Reduce(`+`, Map(function(U, sw) sw * U, uncentered, subject_weights))
    deviations <- lapply(uncentered, function(U) U - mean_deviation)
    names(deviations) <- subject_names
    for (s in seq_along(deviations)) {
      rownames(deviations[[s]]) <- colnames(X)[random_columns]
      colnames(deviations[[s]]) <- colnames(Y)
    }

    deviation_fit <- matrix(0, nrow(Y), ncol(Y))
    for (s in seq_along(subject_rows)) {
      rows <- subject_rows[[s]]
      deviation_fit[rows, ] <- X[rows, random_columns, drop = FALSE] %*% deviations[[s]]
    }
    group <- .dmfa_wls_sufficient(
      X,
      Y - deviation_fit,
      weights,
      control$rank_tol,
      subject_rows = subject_rows
    )
    residual <- Y - X %*% group$coefficient - deviation_fit
    sse <- sum(residual^2 * weights)
    penalty_value <- subject_penalty * sum(vapply(seq_along(deviations), function(s) {
      subject_weights[[s]] * sum(deviations[[s]]^2)
    }, numeric(1)))
    objective <- sse + penalty_value
    increase_tol <- 1000 * .Machine$double.eps * max(1, abs(previous), abs(objective))
    if (objective > previous + increase_tol) {
      stop("Partial-pooling objective increased materially during block-coordinate updates.", call. = FALSE)
    }
    objective_trace <- c(objective_trace, objective)
    relative_change <- abs(previous - objective) / max(1, abs(previous))
    previous <- objective
    if (relative_change < control$multilevel_tol) {
      converged <- TRUE
      break
    }
  }

  group$sse <- sse
  group$penalized <- TRUE
  group$penalty_value <- penalty_value
  group$objective <- objective
  group$objective_trace <- objective_trace
  group$subject_deviations <- deviations
  group$converged <- converged
  group$iterations <- length(objective_trace) - 1L
  group$final_relative_change <- relative_change
  group$random_columns <- colnames(X)[random_columns]
  group
}

.dmfa_fit_model <- function(X,
                            Y,
                            weights,
                            control,
                            subject_model,
                            subject_rows,
                            subject_weights,
                            random_columns,
                            subject_penalty) {
  if (identical(subject_model, "fixed_intercept")) {
    out <- .dmfa_wls_sufficient(
      X, Y, weights, control$rank_tol, subject_rows = subject_rows
    )
    out$penalized <- FALSE
    out$objective <- out$sse
    out$penalty_value <- 0
    out$subject_deviations <- NULL
    out$converged <- TRUE
    out$iterations <- 0L
    out$computation <- "stacked_sufficient"
    return(out)
  }
  out <- .dmfa_fit_partial_pool(
    X = X,
    Y = Y,
    weights = weights,
    subject_rows = subject_rows,
    subject_weights = subject_weights,
    random_columns = random_columns,
    subject_penalty = subject_penalty,
    control = control
  )
  out$computation <- "stacked_partial_pooling"
  out
}

.dmfa_soft_map_cache <- function(feature_map, feature_chunk_size) {
  lapply(feature_map, function(P) {
    list(
      diagonal = as.numeric(Matrix::colSums(P^2)),
      chunks = lapply(
        seq.int(1L, nrow(P), by = feature_chunk_size),
        function(first) seq.int(first, min(nrow(P), first + feature_chunk_size - 1L))
      )
    )
  })
}

.dmfa_fit_soft_map <- function(X,
                               Y_list,
                               weight_list,
                               feature_map,
                               subject_rows,
                               alignment_penalty,
                               control,
                               map_cache = NULL) {
  q <- ncol(X)
  p_reference <- ncol(feature_map[[1L]])
  feature_names <- colnames(feature_map[[1L]])
  G_list <- vector("list", length(Y_list))
  C <- matrix(0, q, p_reference)
  raw_response_ss <- 0
  preconditioner <- matrix(alignment_penalty, q, p_reference)
  if (is.null(map_cache)) {
    map_cache <- .dmfa_soft_map_cache(feature_map, control$feature_chunk_size)
  }
  if (length(map_cache) != length(feature_map)) {
    stop("Internal soft-map cache dimension mismatch.", call. = FALSE)
  }

  for (s in seq_along(Y_list)) {
    rows <- subject_rows[[s]]
    Xs <- X[rows, , drop = FALSE]
    Ys <- Y_list[[s]]
    ws <- weight_list[[s]]
    Ps <- feature_map[[s]]
    G_list[[s]] <- crossprod(Xs * sqrt(ws))
    for (columns in map_cache[[s]]$chunks) {
      Ys_chunk <- Ys[, columns, drop = FALSE]
      Ps_chunk <- Ps[columns, , drop = FALSE]
      C <- C + as.matrix(crossprod(Xs, Ys_chunk * ws) %*% Ps_chunk)
      raw_response_ss <- raw_response_ss + sum(Ys_chunk^2 * ws)
    }
    preconditioner <- preconditioner + outer(
      as.numeric(Matrix::diag(G_list[[s]])),
      map_cache[[s]]$diagonal
    )
  }
  preconditioner_floor <- sqrt(.Machine$double.eps) * max(1, max(preconditioner))
  preconditioner[preconditioner < preconditioner_floor] <- preconditioner_floor

  apply_operator <- function(B, include_ridge = TRUE) {
    out <- matrix(0, q, p_reference)
    for (s in seq_along(G_list)) {
      GB <- G_list[[s]] %*% B
      Ps <- feature_map[[s]]
      for (columns in map_cache[[s]]$chunks) {
        Ps_chunk <- Ps[columns, , drop = FALSE]
        mapped <- as.matrix(GB %*% Matrix::t(Ps_chunk))
        out <- out + as.matrix(mapped %*% Ps_chunk)
      }
    }
    if (isTRUE(include_ridge) && alignment_penalty > 0) {
      out <- out + alignment_penalty * B
    }
    out
  }

  coefficient <- matrix(0, q, p_reference)
  residual <- C
  preconditioned_residual <- residual / preconditioner
  direction <- preconditioned_residual
  residual_sq <- sum(residual^2)
  residual_inner <- sum(residual * preconditioned_residual)
  rhs_scale <- max(1, sqrt(residual_sq))
  # The additive ||Y||^2 constant is omitted while iterating.  This preserves
  # the exact CG objective ordering without subtracting nearly equal O(Y^2)
  # quantities when responses contain a large subject offset.
  objective_trace <- 0
  converged <- residual_sq == 0
  iterations <- 0L
  operator_applications <- 0L

  if (!converged) {
    for (iteration in seq_len(control$soft_max_iter)) {
      Ad <- apply_operator(direction)
      operator_applications <- operator_applications + 1L
      denominator <- sum(direction * Ad)
      if (!is.finite(denominator) || denominator <= 0) {
        stop("Soft-alignment normal operator is not positive definite on the preconditioned CG search direction.", call. = FALSE)
      }
      step <- residual_inner / denominator
      coefficient <- coefficient + step * direction
      residual_new <- residual - step * Ad
      residual_sq_new <- sum(residual_new^2)
      objective_relative <- -sum(C * coefficient) - sum(residual_new * coefficient)
      objective_trace <- c(objective_trace, objective_relative)
      iterations <- iteration
      if (sqrt(residual_sq_new) <= control$soft_tol * rhs_scale) {
        residual <- residual_new
        residual_sq <- residual_sq_new
        converged <- TRUE
        break
      }
      preconditioned_residual_new <- residual_new / preconditioner
      residual_inner_new <- sum(residual_new * preconditioned_residual_new)
      beta <- residual_inner_new / residual_inner
      direction <- preconditioned_residual_new + beta * direction
      residual <- residual_new
      residual_sq <- residual_sq_new
      residual_inner <- residual_inner_new
    }
  }

  response_ss <- 0
  sse <- 0
  for (s in seq_along(Y_list)) {
    rows <- subject_rows[[s]]
    Xs <- X[rows, , drop = FALSE]
    Ys <- Y_list[[s]]
    ws <- weight_list[[s]]
    Ps <- feature_map[[s]]
    for (chunk_start in seq.int(1L, ncol(Ys), by = control$feature_chunk_size)) {
      columns <- seq.int(
        chunk_start,
        min(ncol(Ys), chunk_start + control$feature_chunk_size - 1L)
      )
      Ys_chunk <- Ys[, columns, drop = FALSE]
      center <- colSums(Ys_chunk * ws) / sum(ws)
      response_ss <- response_ss + sum(sweep(Ys_chunk, 2L, center, "-")^2 * ws)
      predicted <- Xs %*% coefficient %*% Matrix::t(Ps[columns, , drop = FALSE])
      sse <- sse + sum((Ys_chunk - as.matrix(predicted))^2 * ws)
    }
  }
  penalty_value <- alignment_penalty * sum(coefficient^2)
  objective <- sse + penalty_value
  objective_trace <- objective_trace + (objective - tail(objective_trace, 1L))

  design_model <- .dmfa_wls_sufficient(
    X,
    matrix(0, nrow(X), 1L),
    unlist(weight_list, use.names = FALSE),
    control$rank_tol,
    subject_rows = subject_rows
  )
  dimnames(coefficient) <- list(colnames(X), feature_names)
  design_model$coefficient <- coefficient
  design_model$response_ss <- response_ss
  design_model$fitted_ss <- response_ss - sse
  design_model$sse <- sse
  design_model$penalized <- alignment_penalty > 0
  design_model$penalty_value <- penalty_value
  design_model$objective <- objective
  design_model$objective_trace <- objective_trace
  design_model$subject_deviations <- NULL
  design_model$converged <- converged
  design_model$iterations <- iterations
  design_model$solver_residual <- sqrt(residual_sq)
  design_model$solver_relative_residual <- sqrt(residual_sq) / rhs_scale
  design_model$solver_rhs_norm <- rhs_scale
  design_model$raw_response_ss <- raw_response_ss
  design_model$preconditioner_range <- range(preconditioner)
  design_model$operator_applications <- operator_applications
  design_model$map_gram_materialized <- FALSE
  design_model$operator_chunk_size <- control$feature_chunk_size
  design_model$soft_alignment <- TRUE
  design_model$computation <- "matrix_free_soft_alignment"
  design_model
}

.dmfa_comparison_tolerance <- function(full, reduced, scale, control) {
  numeric_tol <- 100 * .Machine$double.eps * max(1, scale)
  iterative <- isTRUE(full$soft_alignment) || isTRUE(reduced$soft_alignment)
  if (!iterative) return(numeric_tol)

  residual_bound <- function(model) {
    if (is.null(model$solver_residual) || is.null(model$coefficient)) return(0)
    model$solver_residual * sqrt(sum(model$coefficient^2))
  }
  max(
    numeric_tol,
    10 * (residual_bound(full) + residual_bound(reduced)),
    10 * control$soft_tol * max(1, scale)
  )
}

.dmfa_nonconvergence_message <- function(model, stage, control) {
  if (isTRUE(model$soft_alignment)) {
    return(sprintf(
      paste0(
        "%s soft-alignment solver did not converge after %d iterations ",
        "(relative residual %.6g; tolerance %.6g; RHS norm %.6g)"
      ),
      stage,
      model$iterations,
      model$solver_relative_residual,
      control$soft_tol,
      model$solver_rhs_norm
    ))
  }
  if (!is.null(model$final_relative_change)) {
    return(sprintf(
      paste0(
        "%s partial-pooling solver did not converge after %d iterations ",
        "(final relative objective change %.6g; tolerance %.6g)"
      ),
      stage,
      model$iterations,
      model$final_relative_change,
      control$multilevel_tol
    ))
  }
  sprintf("%s solver did not converge after %d iterations", stage, model$iterations)
}

.dmfa_resolve_effect_rank <- function(rank, term, maximum) {
  if (maximum <= 0L) return(0L)
  if (is.null(rank)) return(maximum)
  value <- rank
  if (is.list(rank)) {
    if (is.null(names(rank)) || !term %in% names(rank)) return(maximum)
    value <- rank[[term]]
  } else if (length(rank) > 1L) {
    if (is.null(names(rank)) || !term %in% names(rank)) {
      stop("Multi-valued `rank` must be named by design term.", call. = FALSE)
    }
    value <- rank[[term]]
  }
  value <- as.integer(value)
  if (length(value) != 1L || is.na(value) || value < 1L) {
    stop(sprintf("Invalid component rank for term '%s'.", term), call. = FALSE)
  }
  min(value, maximum)
}

.dmfa_thin_svd <- function(H,
                           weights,
                           ncomp = NULL,
                           rank_tol = sqrt(.Machine$double.eps),
                           feature_metric = NULL) {
  if (!is.matrix(H)) H <- as.matrix(H)
  if (length(weights) != nrow(H) || any(weights < 0) || sum(weights) <= 0) {
    stop("Invalid effect-grid weights.", call. = FALSE)
  }
  Hw <- H * sqrt(weights)
  nr <- nrow(Hw)
  nc <- ncol(Hw)
  total_energy <- sum(Hw^2)
  if (!is.null(ncomp) && ncomp <= 0L) {
    return(list(
      d = numeric(0),
      u = matrix(0, nr, 0),
      v = matrix(0, nc, 0),
      scores = matrix(0, nr, 0),
      total_energy = total_energy
    ))
  }

  if (!is.null(feature_metric)) {
    if (nrow(feature_metric) != nc || ncol(feature_metric) != nc) {
      stop("Feature metric dimension does not match the effect matrix.", call. = FALSE)
    }
    solved <- Matrix::solve(feature_metric, Matrix::Matrix(t(Hw), sparse = FALSE))
    gram <- as.matrix(Hw %*% solved)
    gram <- (gram + t(gram)) / 2
    total_energy <- max(0, sum(diag(gram)))
    eig <- eigen(gram, symmetric = TRUE)
    values <- pmax(eig$values, 0)
    threshold <- rank_tol * max(1, values[[1L]])
    keep <- which(values > threshold)
    if (length(keep) == 0L) {
      return(list(
        d = numeric(0),
        u = matrix(0, nr, 0),
        v = matrix(0, nc, 0),
        scores = matrix(0, nr, 0),
        total_energy = total_energy
      ))
    }
    if (!is.null(ncomp)) keep <- keep[seq_len(min(length(keep), ncomp))]
    d <- sqrt(values[keep])
    u <- eig$vectors[, keep, drop = FALSE]
    v <- as.matrix(solved %*% u)
    v <- sweep(v, 2L, d, "/")
    raw_scores <- H %*% v
    for (j in seq_len(ncol(v))) {
      anchor <- which.max(abs(v[, j]))
      if (v[anchor, j] < 0) {
        v[, j] <- -v[, j]
        u[, j] <- -u[, j]
        raw_scores[, j] <- -raw_scores[, j]
      }
    }
    colnames(v) <- colnames(u) <- colnames(raw_scores) <- paste0("Component_", seq_len(ncol(v)))
    rownames(v) <- colnames(H)
    return(list(
      d = d,
      u = u,
      v = v,
      scores = raw_scores,
      total_energy = total_energy
    ))
  }

  if (nr <= nc) {
    eig <- eigen(tcrossprod(Hw), symmetric = TRUE)
    values <- pmax(eig$values, 0)
    threshold <- rank_tol * max(1, values[[1L]])
    keep <- which(values > threshold)
    if (length(keep) == 0L) return(list(d = numeric(0), u = matrix(0, nr, 0), v = matrix(0, nc, 0), scores = matrix(0, nr, 0), total_energy = total_energy))
    if (!is.null(ncomp)) keep <- keep[seq_len(min(length(keep), ncomp))]
    d <- sqrt(values[keep])
    u <- eig$vectors[, keep, drop = FALSE]
    v <- crossprod(Hw, u)
    v <- sweep(v, 2L, d, "/")
  } else {
    eig <- eigen(crossprod(Hw), symmetric = TRUE)
    values <- pmax(eig$values, 0)
    threshold <- rank_tol * max(1, values[[1L]])
    keep <- which(values > threshold)
    if (length(keep) == 0L) return(list(d = numeric(0), u = matrix(0, nr, 0), v = matrix(0, nc, 0), scores = matrix(0, nr, 0), total_energy = total_energy))
    if (!is.null(ncomp)) keep <- keep[seq_len(min(length(keep), ncomp))]
    d <- sqrt(values[keep])
    v <- eig$vectors[, keep, drop = FALSE]
    u <- H %*% v
    u <- u * sqrt(weights)
    u <- sweep(u, 2L, d, "/")
  }

  raw_scores <- H %*% v
  for (j in seq_len(ncol(v))) {
    anchor <- which.max(abs(v[, j]))
    if (v[anchor, j] < 0) {
      v[, j] <- -v[, j]
      u[, j] <- -u[, j]
      raw_scores[, j] <- -raw_scores[, j]
    }
  }
  colnames(v) <- colnames(u) <- colnames(raw_scores) <- paste0("Component_", seq_len(ncol(v)))
  rownames(v) <- colnames(H)

  list(
    d = d,
    u = u,
    v = v,
    scores = raw_scores,
    total_energy = total_energy,
    path = "dense",
    gram_dim = c(min(nr, nc), min(nr, nc))
  )
}

.dmfa_empty_svd <- function(nr, nc, total_energy = 0, path = "factored_small_gram") {
  list(
    d = numeric(0),
    u = matrix(0, nr, 0),
    v = matrix(0, nc, 0),
    scores = matrix(0, nr, 0),
    total_energy = total_energy,
    path = path,
    gram_dim = c(0L, 0L)
  )
}

.dmfa_unpivoted_thin_qr <- function(X) {
  qrd <- qr(X)
  Q <- qr.Q(qrd, complete = FALSE)
  R_pivoted <- qr.R(qrd)
  R <- matrix(0, nrow(R_pivoted), ncol(X))
  R[, qrd$pivot] <- R_pivoted
  list(Q = Q, R = R)
}

.dmfa_finalize_svd <- function(d, u, v, scores, feature_names) {
  for (j in seq_len(ncol(v))) {
    anchor <- which.max(abs(v[, j]))
    if (v[anchor, j] < 0) {
      v[, j] <- -v[, j]
      u[, j] <- -u[, j]
      scores[, j] <- -scores[, j]
    }
  }
  colnames(v) <- colnames(u) <- colnames(scores) <- paste0("Component_", seq_len(ncol(v)))
  rownames(v) <- feature_names
  list(d = d, u = u, v = v, scores = scores)
}

.dmfa_thin_svd_factored <- function(left,
                                    right,
                                    weights,
                                    ncomp = NULL,
                                    rank_tol = sqrt(.Machine$double.eps),
                                    feature_metric = NULL) {
  if (!is.matrix(left)) left <- as.matrix(left)
  if (!is.matrix(right)) right <- as.matrix(right)
  if (nrow(right) != ncol(left)) {
    stop("Factored SVD factors have incompatible inner dimensions.", call. = FALSE)
  }
  if (length(weights) != nrow(left) || any(weights < 0) || sum(weights) <= 0) {
    stop("Invalid effect-grid weights.", call. = FALSE)
  }
  nr <- nrow(left)
  nc <- ncol(right)
  if (!is.null(ncomp) && ncomp <= 0L) {
    return(.dmfa_empty_svd(nr, nc, path = "factored_small_gram"))
  }
  if (ncol(left) == 0L || nr == 0L || nc == 0L) {
    return(.dmfa_empty_svd(nr, nc, path = "factored_small_gram"))
  }

  Lw <- left * sqrt(weights)
  qr_lw <- .dmfa_unpivoted_thin_qr(Lw)
  Q <- qr_lw$Q
  R <- qr_lw$R

  if (!is.null(feature_metric)) {
    if (nrow(feature_metric) != nc || ncol(feature_metric) != nc) {
      stop("Feature metric dimension does not match the effect matrix.", call. = FALSE)
    }
    Z <- as.matrix(Matrix::solve(feature_metric, Matrix::Matrix(t(right), sparse = FALSE)))
    mid <- right %*% Z
    mid <- (mid + t(mid)) / 2
    gram <- R %*% mid %*% t(R)
  } else {
    gram <- R %*% tcrossprod(right) %*% t(R)
  }
  gram <- (gram + t(gram)) / 2
  total_energy <- max(0, sum(diag(gram)))
  eig <- eigen(gram, symmetric = TRUE)
  values <- pmax(eig$values, 0)
  if (length(values) == 0L) {
    return(.dmfa_empty_svd(nr, nc, total_energy, path = "factored_small_gram"))
  }
  threshold <- rank_tol * max(1, values[[1L]])
  keep <- which(values > threshold)
  if (length(keep) == 0L) {
    return(.dmfa_empty_svd(nr, nc, total_energy, path = "factored_small_gram"))
  }
  if (!is.null(ncomp)) keep <- keep[seq_len(min(length(keep), ncomp))]
  d <- sqrt(values[keep])
  u_small <- eig$vectors[, keep, drop = FALSE]
  u <- Q %*% u_small
  if (!is.null(feature_metric)) {
    v <- Z %*% crossprod(R, u_small)
  } else {
    v <- crossprod(right, crossprod(R, u_small))
  }
  v <- sweep(v, 2L, d, "/")
  raw_scores <- left %*% (right %*% v)
  aligned <- .dmfa_finalize_svd(d, u, v, raw_scores, colnames(right))
  list(
    d = aligned$d,
    u = aligned$u,
    v = aligned$v,
    scores = aligned$scores,
    total_energy = total_energy,
    path = "factored_small_gram",
    gram_dim = dim(gram)
  )
}

.dmfa_fit_effect <- function(term,
                             compiled,
                             global_full,
                             fit_model,
                             ss_type,
                             term_scope,
                             rank,
                             control,
                             store,
                             feature_metric = NULL) {
  indices <- .dmfa_term_indices(compiled, term, ss_type, term_scope)
  if (identical(ss_type, "partial")) {
    full <- global_full
  } else {
    full <- fit_model(indices$full)
  }
  reduced <- fit_model(indices$reduced)
  if (!isTRUE(full$converged) || !isTRUE(reduced$converged)) {
    failures <- character(0)
    if (!isTRUE(full$converged)) {
      failures <- c(failures, .dmfa_nonconvergence_message(full, "containing", control))
    }
    if (!isTRUE(reduced$converged)) {
      failures <- c(failures, .dmfa_nonconvergence_message(reduced, "reduced", control))
    }
    stop(
      sprintf(
        "Term '%s' could not be decomposed: %s. Increase the relevant iteration limit or relax its tolerance.",
        term,
        paste(failures, collapse = "; ")
      ),
      call. = FALSE
    )
  }

  L_full <- compiled$X_reference[, indices$full, drop = FALSE]
  L_reduced <- compiled$X_reference[, indices$reduced, drop = FALSE]
  estimable_full <- .dmfa_estimable_rows(L_full, full, control$estimability_tol)
  estimable_reduced <- .dmfa_estimable_rows(L_reduced, reduced, control$estimability_tol)
  effect_columns <- setdiff(indices$full, indices$reduced)
  local_effect_columns <- match(effect_columns, indices$full)
  L_effect <- matrix(0, nrow(L_full), ncol(L_full))
  L_effect[, local_effect_columns] <- L_full[, local_effect_columns, drop = FALSE]
  effect_estimable <- .dmfa_estimable_rows(
    L_effect, full, control$estimability_tol
  )
  estimable_grid <- effect_estimable$estimable
  raw_estimable_grid <- estimable_full$estimable & estimable_reduced$estimable
  prediction_check <- .dmfa_estimable_rows(
    compiled$X_reference,
    global_full,
    control$estimability_tol
  )
  prediction_estimable_grid <- prediction_check$estimable
  df <- full$rank - reduced$rank

  if (isTRUE(full$penalized)) {
    delta_objective_raw <- reduced$objective - full$objective
    comparison_tol <- .dmfa_comparison_tolerance(
      full, reduced, max(reduced$objective, full$objective), control
    )
    if (delta_objective_raw < -comparison_tol) {
      stop(sprintf("Reduced penalized model for term '%s' has a lower objective than its containing model.", term), call. = FALSE)
    }
    partial_ss <- NA_real_
    partial_r2 <- NA_real_
    objective_delta <- max(0, delta_objective_raw)
    objective_fraction <- if (reduced$objective > 0) objective_delta / reduced$objective else NA_real_
    comparison <- "penalized_objective"
  } else {
    delta_ss_raw <- reduced$sse - full$sse
    comparison_tol <- .dmfa_comparison_tolerance(
      full, reduced, max(reduced$sse, full$sse), control
    )
    if (delta_ss_raw < -comparison_tol) {
      stop(sprintf("Reduced model for term '%s' fits better than its containing model.", term), call. = FALSE)
    }
    partial_ss <- max(0, delta_ss_raw)
    partial_r2 <- if (reduced$sse > 0) partial_ss / reduced$sse else NA_real_
    objective_delta <- partial_ss
    objective_fraction <- partial_r2
    comparison <- "partial_sse"
  }

  status <- "estimable"
  reason <- NULL
  decomposition <- NULL
  H <- NULL
  if (df <= 0L) {
    status <- "nonestimable"
    reason <- "Removing the term does not reduce model rank."
  } else if (!all(estimable_grid)) {
    status <- "nonestimable_reference"
    reason <- sprintf(
      "%d of %d reference-grid rows are not estimable under the full and reduced models.",
      sum(!estimable_grid), length(estimable_grid)
    )
  } else {
    A <- L_full[, local_effect_columns, drop = FALSE]
    rownames(A) <- NULL
    B <- full$coefficient[local_effect_columns, , drop = FALSE]
    center_A <- as.numeric(crossprod(compiled$reference_weights, A))
    A <- sweep(A, 2L, center_A, "-")
    if (store %in% c("fitted", "residuals")) {
      H <- A %*% B
    }
    maximum <- min(df, nrow(A) - 1L, ncol(B))
    requested_rank <- .dmfa_resolve_effect_rank(rank, term, maximum)
    decomposition <- .dmfa_thin_svd_factored(
      left = A,
      right = B,
      weights = compiled$reference_weights,
      ncomp = requested_rank,
      rank_tol = control$rank_tol,
      feature_metric = feature_metric
    )
    if (length(decomposition$d) == 0L) {
      status <- "degenerate"
      reason <- "The estimable reference-grid effect is numerically zero."
    }
  }

  effect <- list(
    term = term,
    removed_terms = indices$removed_terms,
    status = status,
    reason = reason,
    estimable = identical(status, "estimable"),
    estimable_df = max(0L, df),
    estimable_grid = estimable_grid,
    raw_estimable_grid = raw_estimable_grid,
    prediction_estimable_grid = prediction_estimable_grid,
    estimability_residual = effect_estimable$relative_residual,
    raw_estimability_residual = pmax(
      estimable_full$relative_residual,
      estimable_reduced$relative_residual
    ),
    reference_grid = compiled$reference_grid,
    reference_weights = compiled$reference_weights,
    partial_ss = partial_ss,
    partial_r2 = partial_r2,
    objective_delta = objective_delta,
    objective_fraction = objective_fraction,
    comparison = comparison,
    full_sse = full$sse,
    reduced_sse = reduced$sse,
    rank = if (is.null(decomposition)) 0L else length(decomposition$d),
    singular_values = if (is.null(decomposition)) numeric(0) else decomposition$d,
    explained = if (is.null(decomposition) || decomposition$total_energy == 0) {
      numeric(0)
    } else {
      decomposition$d^2 / decomposition$total_energy
    },
    scores = if (is.null(decomposition)) NULL else decomposition$scores,
    loadings = if (is.null(decomposition)) NULL else decomposition$v,
    ss_type = ss_type,
    term_scope = term_scope,
    loading_metric = if (is.null(feature_metric)) "euclidean" else "spatial",
    svd_path = if (is.null(decomposition)) NULL else decomposition$path,
    model_columns = list(full = colnames(compiled$X)[indices$full], reduced = colnames(compiled$X)[indices$reduced])
  )
  if (store %in% c("fitted", "residuals") && !is.null(H)) effect$effect_matrix <- H
  class(effect) <- c("design_mfa_effect", "list")
  effect
}

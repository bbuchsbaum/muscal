#!/usr/bin/env Rscript

source("dev/benchmark_qualification.R")
suppressPackageStartupMessages(bench_load_local())

simulate_design_mfa <- function(n_subjects, trials, p, seed) {
  set.seed(seed)
  ids <- paste0("S", seq_len(n_subjects))
  feature_names <- paste0("v", seq_len(p))
  saliency_levels <- c(20, 40, 60, 80, 100)
  probe_levels <- c("old", "new", "lure")
  map1 <- scale(seq(-1, 1, length.out = p), scale = FALSE)[, 1]
  map2 <- scale(sin(seq(0, 2 * pi, length.out = p)), scale = FALSE)[, 1]
  design <- Y <- vector("list", n_subjects)
  for (s in seq_len(n_subjects)) {
    d <- data.frame(
      saliency = ordered(sample(saliency_levels, trials[s], TRUE),
                         levels = saliency_levels),
      probe = factor(sample(probe_levels, trials[s], TRUE), levels = probe_levels)
    )
    sal <- as.numeric(d$saliency) - 3
    probe <- c(old = -1, new = 0, lure = 1)[as.character(d$probe)]
    y <- outer(sal, map1) + 0.6 * outer(probe, map2) +
      matrix(rnorm(trials[s] * p, sd = 0.5), trials[s], p)
    colnames(y) <- feature_names
    design[[s]] <- d
    Y[[s]] <- y
  }
  names(design) <- names(Y) <- ids
  list(Y = Y, design = design, feature_names = feature_names)
}

make_laplacian <- function(p, feature_names) {
  adjacency <- Matrix::bandSparse(
    p,
    k = c(-1, 1),
    diagonals = list(rep(1, p - 1L), rep(1, p - 1L))
  )
  out <- Matrix::Diagonal(x = Matrix::rowSums(adjacency)) - adjacency
  dimnames(out) <- list(feature_names, feature_names)
  out
}

make_soft_data <- function(sim, seed) {
  set.seed(seed)
  p <- length(sim$feature_names)
  angle <- seq(0, 6 * pi, length.out = p)
  reference_coords <- cbind(
    x = seq(-3, 3, length.out = p),
    y = sin(angle),
    z = cos(angle)
  )
  rownames(reference_coords) <- sim$feature_names
  subject_coords <- vector("list", length(sim$Y))
  for (s in seq_along(sim$Y)) {
    coords <- reference_coords + matrix(rnorm(p * 3L, sd = 0.002), p, 3L)
    coords <- coords[sample.int(p), , drop = FALSE]
    rownames(coords) <- paste0("S", s, "_v", seq_len(p))
    subject_coords[[s]] <- coords
  }
  names(subject_coords) <- names(sim$Y)
  maps <- coordinate_feature_maps(
    reference_coords,
    subject_coords,
    k = 4,
    bandwidth = 0.01,
    kernel = "gaussian"
  )
  Y <- Map(function(y, P) {
    out <- y %*% Matrix::t(P)
    colnames(out) <- rownames(P)
    as.matrix(out)
  }, sim$Y, maps)
  names(Y) <- names(sim$Y)
  list(Y = Y, maps = maps)
}

max_profiled_allocation <- function(path) {
  if (!file.exists(path)) return(NA_real_)
  lines <- readLines(path, warn = FALSE)
  bytes <- suppressWarnings(as.numeric(sub(" .*", "", lines)))
  if (!any(is.finite(bytes))) NA_real_ else max(bytes, na.rm = TRUE)
}

profile_fit <- function(label, expression, n_trials, p_reference) {
  profile <- tempfile("design-mfa-memory-")
  profiling <- capabilities("profmem")
  if (profiling) utils::Rprofmem(profile)
  gc(reset = TRUE)
  started <- proc.time()[["elapsed"]]
  fit <- tryCatch(force(expression), error = function(e) e)
  elapsed <- proc.time()[["elapsed"]] - started
  heap <- gc()
  if (profiling) utils::Rprofmem(NULL)
  max_used <- if ("max used" %in% colnames(heap)) {
    sum(heap[, "max used"] * c(56, 8))
  } else {
    NA_real_
  }
  if (inherits(fit, "error")) {
    row <- data.frame(engine = label, ok = FALSE, error = conditionMessage(fit),
      computation = NA_character_, subjects = NA_integer_, trials = n_trials,
      reference_features = p_reference, elapsed_sec = elapsed,
      trial_features_per_sec = NA_real_, fit_bytes = NA_real_, approx_max_r_heap_bytes = max_used,
      max_single_allocation_bytes = NA_real_, design_rank = NA_integer_, solver_iterations = NA_integer_,
      solver_converged = NA, solver_residual = NA_real_, solver_relative_residual = NA_real_)
    unlink(profile)
    return(list(fit = NULL, row = row))
  }
  row <- data.frame(
    engine = label,
    ok = !inherits(fit, "error"),
    error = if (inherits(fit, "error")) conditionMessage(fit) else NA_character_,
    computation = if (inherits(fit, "error")) NA_character_ else fit$computation,
    subjects = length(fit$subjects),
    trials = n_trials,
    reference_features = p_reference,
    elapsed_sec = elapsed,
    trial_features_per_sec = n_trials * p_reference / max(elapsed, .Machine$double.eps),
    fit_bytes = as.numeric(object.size(fit)),
    approx_max_r_heap_bytes = max_used,
    max_single_allocation_bytes = if (profiling) max_profiled_allocation(profile) else NA_real_,
    design_rank = fit$full_model$rank,
    solver_iterations = fit$full_model$iterations,
    solver_converged = if (is.null(fit$full_model$converged)) NA else fit$full_model$converged,
    solver_residual = if (is.null(fit$full_model$solver_residual)) {
      NA_real_
    } else {
      fit$full_model$solver_residual
    },
    solver_relative_residual = if (is.null(fit$full_model$solver_relative_residual)) {
      NA_real_
    } else {
      fit$full_model$solver_relative_residual
    },
    stringsAsFactors = FALSE
  )
  unlink(profile)
  list(fit = fit, row = row)
}

compare_sufficient_solvers <- function(fit, Y, repetitions = 5L) {
  X <- fit$design$X
  offsets <- cumsum(c(0L, vapply(Y, nrow, integer(1))))
  rows <- lapply(seq_along(Y), function(s) seq.int(offsets[[s]] + 1L, offsets[[s + 1L]]))
  names(rows) <- names(Y)
  streamed_time <- system.time({
    for (iteration in seq_len(repetitions)) {
      streamed <- muscal:::.dmfa_wls_list_sufficient(
        X,
        Y,
        fit$weights,
        rows,
        fit$control$rank_tol
      )
    }
  })[["elapsed"]] / repetitions
  dense_time <- system.time({
    for (iteration in seq_len(repetitions)) {
      dense <- muscal:::.dmfa_wls_sufficient(
        X,
        do.call(rbind, Y),
        unlist(fit$weights, use.names = FALSE),
        fit$control$rank_tol
      )
    }
  })[["elapsed"]] / repetitions
  data.frame(
    streamed_stats_sec = streamed_time,
    dense_stats_sec = dense_time,
    streamed_dense_runtime_ratio = streamed_time / max(dense_time, .Machine$double.eps),
    streamed_dense_coefficient_max_abs = max(abs(streamed$coefficient - dense$coefficient)),
    streamed_dense_sse_abs = abs(streamed$sse - dense$sse),
    stringsAsFactors = FALSE
  )
}

benchmark_coordinate_maps <- function(p, n_subjects, seed, subjects_to_map = n_subjects) {
  set.seed(seed)
  reference <- cbind(
    x = runif(p, -90, 90),
    y = runif(p, -126, 90),
    z = runif(p, -72, 108)
  )
  rownames(reference) <- paste0("v", seq_len(p))
  subject_coords <- lapply(seq_len(subjects_to_map), function(s) {
    out <- reference + matrix(rnorm(p * 3L, sd = 0.35), p, 3L)
    rownames(out) <- paste0("S", s, "_v", seq_len(p))
    out
  })
  names(subject_coords) <- paste0("S", seq_len(subjects_to_map))
  elapsed <- system.time({
    maps <- coordinate_feature_maps(
      reference,
      subject_coords,
      k = 4,
      bandwidth = 1,
      kernel = "gaussian"
    )
  })[["elapsed"]]
  data.frame(
    coordinate_map_subjects = subjects_to_map,
    coordinate_map_elapsed_sec = elapsed,
    coordinate_map_nonzeros = sum(vapply(maps, function(P) length(P@x), integer(1))),
    coordinate_map_bytes = sum(vapply(maps, function(P) as.numeric(object.size(P)), numeric(1))),
    stringsAsFactors = FALSE
  )
}

args <- commandArgs(trailingOnly = TRUE)
quick <- any(args %in% c("--quick", "-q"))
imaging <- any(args == "--imaging")
imaging_soft <- any(args == "--imaging-soft")
replicate_arg <- grep("^--replicate=", args, value = TRUE)
benchmark_replicate <- if (length(replicate_arg) == 0L) {
  1L
} else {
  as.integer(sub("^--replicate=", "", replicate_arg[[1L]]))
}
output_arg <- grep("^--output=", args, value = TRUE)
output_path <- if (length(output_arg) == 0L) NULL else sub("^--output=", "", output_arg[[1L]])
if (length(benchmark_replicate) != 1L || is.na(benchmark_replicate) || benchmark_replicate < 1L) {
  stop("`--replicate` must be a positive integer.", call. = FALSE)
}

if (imaging_soft) {
  n_subjects <- 4L
  trials <- 24:27
  p <- 100000L
} else if (imaging) {
  n_subjects <- 20L
  trials <- 32:51
  p <- 100000L
} else if (quick) {
  n_subjects <- 4L
  trials <- c(60L, 70L, 80L, 90L)
  p <- 256L
} else {
  n_subjects <- 20L
  trials <- seq.int(60L, 98L, length.out = n_subjects)
  p <- 4000L
}

sim <- simulate_design_mfa(n_subjects, trials, p, seed = 20260815 + benchmark_replicate)
formula <- ~ saliency + probe
contrasts <- list(saliency = "poly", probe = "sum")

shared <- profile_fit(
  "shared",
  design_mfa(sim$Y, sim$design, formula, contrasts = contrasts, rank = 2),
  sum(trials),
  p
)

if (imaging && !imaging_soft) {
  rows <- shared$row
  solver_comparison <- data.frame(
    streamed_stats_sec = NA_real_,
    dense_stats_sec = NA_real_,
    streamed_dense_runtime_ratio = NA_real_,
    streamed_dense_coefficient_max_abs = NA_real_,
    streamed_dense_sse_abs = NA_real_
  )
  map_benchmark <- benchmark_coordinate_maps(p, n_subjects, 20260817 + benchmark_replicate)
  prediction_difference <- NA_real_
} else {
  laplacian <- make_laplacian(p, sim$feature_names)
  spatial <- profile_fit(
    "shared_spatial",
    design_mfa(
      sim$Y,
      sim$design,
      formula,
      contrasts = contrasts,
      spatial_penalty = 2,
      spatial_laplacian = laplacian,
      rank = 2
    ),
    sum(trials),
    p
  )

  soft_data <- make_soft_data(sim, seed = 20260816 + benchmark_replicate)
  soft <- profile_fit(
    "soft_sparse_map",
    design_mfa(
      soft_data$Y,
      sim$design,
      formula,
      contrasts = contrasts,
      alignment = "soft",
      feature_map = soft_data$maps,
      alignment_penalty = 0.05,
      rank = 2,
      control = design_mfa_control(soft_max_iter = 1000, soft_tol = 1e-7)
    ),
    sum(trials),
    p
  )
  rows <- rbind(shared$row, spatial$row, soft$row)
  solver_comparison <- if (imaging_soft) {
    data.frame(
      streamed_stats_sec = NA_real_,
      dense_stats_sec = NA_real_,
      streamed_dense_runtime_ratio = NA_real_,
      streamed_dense_coefficient_max_abs = NA_real_,
      streamed_dense_sse_abs = NA_real_
    )
  } else if (!is.null(shared$fit)) {
    compare_sufficient_solvers(shared$fit, sim$Y)
  } else data.frame()
  map_benchmark <- benchmark_coordinate_maps(
    p, n_subjects, 20260817 + benchmark_replicate
  )
  prediction_difference <- if (isTRUE(shared$row$solver_converged) && isTRUE(spatial$row$solver_converged)) max(abs(predict(shared$fit) - predict(spatial$fit))) else NA_real_
}
mode <- if (imaging_soft) "imaging_soft" else if (imaging) "imaging" else if (quick) "quick" else "standard"
saveRDS(bench_provenance(20260815 + benchmark_replicate), paste0(if (is.null(output_path)) sprintf("dev/design_mfa_benchmark_results_%s.csv", mode) else output_path, ".provenance.rds"))
rows$replicate <- benchmark_replicate
rows$benchmark_mode <- mode
rows$input_response_bytes <- sum(vapply(sim$Y, function(y) as.numeric(object.size(y)), numeric(1)))
rows$estimated_dense_stack_bytes <- 8 * sum(trials) * p
rows$r_version <- R.version.string
rows$platform <- R.version$platform
rows$os <- unname(Sys.info()[["sysname"]])
rows$os_release <- unname(Sys.info()[["release"]])
rows$matrix_version <- as.character(utils::packageVersion("Matrix"))
rows$rspectra_version <- as.character(utils::packageVersion("RSpectra"))
rows$blas <- as.character(extSoftVersion()[["BLAS"]])
rows$compiler <- if (is.null(R.version$compiler)) NA_character_ else R.version$compiler
rows$shared_spatial_prediction_max_abs <- NA_real_
if (!imaging) rows$shared_spatial_prediction_max_abs[[2L]] <- prediction_difference
for (field in names(solver_comparison)) rows[[field]] <- solver_comparison[[field]][[1L]]
for (field in names(map_benchmark)) rows[[field]] <- map_benchmark[[field]][[1L]]
print(rows, row.names = FALSE)

if (is.null(output_path)) {
  output_path <- sprintf("dev/design_mfa_benchmark_results_%s.csv", mode)
}
utils::write.csv(rows, output_path, row.names = FALSE)
cat(sprintf("\nWrote benchmark table to %s\n", output_path))
cat(sprintf("Mode: %s\n", mode))
print(sessionInfo())

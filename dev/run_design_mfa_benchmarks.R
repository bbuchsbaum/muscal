#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
mode_arg <- grep("^--mode=", args, value = TRUE)
mode <- if (length(mode_arg) == 0L) "quick" else sub("^--mode=", "", mode_arg[[1L]])
repetition_arg <- grep("^--repetitions=", args, value = TRUE)
repetitions <- if (length(repetition_arg) == 0L) {
  5L
} else {
  as.integer(sub("^--repetitions=", "", repetition_arg[[1L]]))
}
if (!mode %in% c("quick", "standard", "imaging", "imaging_soft")) {
  stop("`--mode` must be quick, standard, imaging, or imaging_soft.", call. = FALSE)
}
if (length(repetitions) != 1L || is.na(repetitions) || repetitions < 5L) {
  stop("`--repetitions` must be an integer of at least five.", call. = FALSE)
}

mode_flag <- switch(
  mode,
  quick = "--quick",
  standard = character(0),
  imaging = "--imaging",
  imaging_soft = "--imaging-soft"
)
temporary <- tempfile("design-mfa-benchmark-")
dir.create(temporary)
on.exit(unlink(temporary, recursive = TRUE), add = TRUE)

paths <- file.path(temporary, sprintf("replicate-%02d.csv", seq_len(repetitions)))
for (replicate in seq_len(repetitions)) {
  command <- c(
    "dev/benchmark_design_mfa.R",
    mode_flag,
    sprintf("--replicate=%d", replicate),
    sprintf("--output=%s", paths[[replicate]])
  )
  log_path <- file.path(temporary, sprintf("replicate-%02d.log", replicate))
  status <- system2(
    file.path(R.home("bin"), "Rscript"), command,
    stdout = log_path, stderr = log_path
  )
  if (!identical(status, 0L) || !file.exists(paths[[replicate]])) {
    if (file.exists(log_path)) cat(readLines(log_path, warn = FALSE), sep = "\n")
    stop(sprintf("Fresh-process benchmark replicate %d failed.", replicate), call. = FALSE)
  }
  message(sprintf("completed fresh-process replicate %d of %d", replicate, repetitions))
}

raw <- do.call(rbind, lapply(paths, utils::read.csv, check.names = FALSE))
rownames(raw) <- NULL
numeric_summary <- c(
  "elapsed_sec",
  "trial_features_per_sec",
  "fit_bytes",
  "approx_max_r_heap_bytes",
  "max_single_allocation_bytes",
  "solver_iterations",
  "solver_residual",
  "solver_relative_residual",
  "coordinate_map_elapsed_sec",
  "coordinate_map_bytes"
)
numeric_summary <- intersect(numeric_summary, names(raw))
groups <- split(seq_len(nrow(raw)), raw$engine)
summary <- do.call(rbind, lapply(groups, function(rows) {
  values <- raw[rows, numeric_summary, drop = FALSE]
  out <- data.frame(
    engine = raw$engine[rows[[1L]]],
    repetitions = length(rows),
    stringsAsFactors = FALSE
  )
  for (field in numeric_summary) {
    x <- values[[field]]
    finite <- x[is.finite(x)]
    statistics <- if (length(finite) == 0L) {
      rep(NA_real_, 5L)
    } else {
      c(
        stats::median(finite),
        stats::quantile(finite, 0.25, names = FALSE),
        stats::quantile(finite, 0.75, names = FALSE),
        min(finite),
        max(finite)
      )
    }
    names(statistics) <- c("median", "q1", "q3", "min", "max")
    for (statistic in names(statistics)) {
      out[[paste0(field, "_", statistic)]] <- statistics[[statistic]]
    }
  }
  out
}))
rownames(summary) <- NULL

raw_path <- sprintf("dev/design_mfa_benchmark_%s_repeated_raw.csv", mode)
summary_path <- sprintf("dev/design_mfa_benchmark_%s_repeated_summary.csv", mode)
utils::write.csv(raw, raw_path, row.names = FALSE)
utils::write.csv(summary, summary_path, row.names = FALSE)
print(summary, row.names = FALSE)
cat(sprintf("Wrote %s and %s\n", raw_path, summary_path))

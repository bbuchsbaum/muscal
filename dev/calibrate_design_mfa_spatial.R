#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  if (requireNamespace("devtools", quietly = TRUE)) {
    devtools::load_all(".", quiet = TRUE)
  }
  library(muscal)
})
source("tests/testthat/helper-design-mfa-sim.R")

args <- commandArgs(trailingOnly = TRUE)
quick <- any(args %in% c("--quick", "-q"))
rep_arg <- grep("^--replicates=", args, value = TRUE)
n_replicates <- if (length(rep_arg) > 0L) {
  as.integer(sub("^--replicates=", "", rep_arg[[1L]]))
} else if (quick) {
  8L
} else {
  60L
}
if (length(n_replicates) != 1L || is.na(n_replicates) || n_replicates < 4L) {
  stop("`--replicates` must specify an integer of at least four.", call. = FALSE)
}

reference_effect_truth <- function(effect, map) {
  score <- as.numeric(effect$reference_grid$saliency) - 3
  H <- outer(score, map)
  sweep(H, 2L, colSums(H * effect$reference_weights), "-")
}

relative_effect_error <- function(effect, truth) {
  sum(rowSums((effect$effect_matrix - truth)^2) * effect$reference_weights) /
    sum(rowSums(truth^2) * effect$reference_weights)
}

effect_energy <- function(effect) {
  sum(rowSums(effect$effect_matrix^2) * effect$reference_weights)
}

fit_pair <- function(sim) {
  common <- list(
    design = sim$design,
    formula = ~ saliency,
    contrasts = list(saliency = "poly"),
    rank = 1,
    store = "fitted"
  )
  soft <- do.call(design_mfa, c(
    list(Y = sim$Y),
    common,
    list(
      alignment = "soft",
      feature_map = sim$maps,
      alignment_penalty = 0.1,
      control = design_mfa_control(soft_max_iter = 1000, soft_tol = 1e-8)
    )
  ))
  naive <- do.call(design_mfa, c(list(Y = sim$Y_naive), common))
  list(soft = soft, naive = naive)
}

run_replicate <- function(replicate, null) {
  seed <- 202608150L + replicate + if (null) 10000L else 0L
  sim <- .sim_design_mfa_misregistration(
    n_subjects = if (quick) 6L else 10L,
    trials = if (quick) 36L else seq.int(48L, 75L, length.out = 10L),
    p = if (quick) 20L else 40L,
    effect_scale = if (null) 0 else 1,
    noise_sd = 0.35,
    seed = seed
  )
  fitted <- fit_pair(sim)
  soft_effect <- effects(fitted$soft, "saliency")
  naive_effect <- effects(fitted$naive, "saliency")
  if (null) {
    data.frame(
      scenario = "null",
      replicate = replicate,
      seed = seed,
      soft_error = NA_real_,
      naive_error = NA_real_,
      soft_loading_correlation = NA_real_,
      naive_loading_correlation = NA_real_,
      soft_energy = effect_energy(soft_effect),
      naive_energy = effect_energy(naive_effect)
    )
  } else {
    truth <- reference_effect_truth(soft_effect, sim$truth$saliency)
    data.frame(
      scenario = "misregistered_effect",
      replicate = replicate,
      seed = seed,
      soft_error = relative_effect_error(soft_effect, truth),
      naive_error = relative_effect_error(naive_effect, truth),
      soft_loading_correlation = abs(cor(effect_map(fitted$soft, "saliency"), sim$truth$saliency)),
      naive_loading_correlation = abs(cor(effect_map(fitted$naive, "saliency"), sim$truth$saliency)),
      soft_energy = effect_energy(soft_effect),
      naive_energy = effect_energy(naive_effect)
    )
  }
}

started <- proc.time()[["elapsed"]]
effect_rows <- lapply(seq_len(n_replicates), run_replicate, null = FALSE)
null_rows <- lapply(seq_len(n_replicates), run_replicate, null = TRUE)
results <- do.call(rbind, c(effect_rows, null_rows))
elapsed <- proc.time()[["elapsed"]] - started

effect <- results[results$scenario == "misregistered_effect", , drop = FALSE]
null <- results[results$scenario == "null", , drop = FALSE]
summary <- data.frame(
  gate = c(
    "mean_same_index_error_below_0.03",
    "mean_same_index_error_below_mean_soft_error",
    "mean_soft_error_below_0.20",
    "median_soft_loading_correlation_above_0.95"
  ),
  value = c(
    mean(effect$naive_error),
    mean(effect$naive_error) / mean(effect$soft_error),
    mean(effect$soft_error),
    stats::median(effect$soft_loading_correlation)
  ),
  threshold = c(0.03, 1, 0.20, 0.95),
  direction = c("below", "below", "below", "above"),
  stringsAsFactors = FALSE
)
summary$passed <- ifelse(
  summary$direction == "below",
  summary$value < summary$threshold,
  summary$value > summary$threshold
)
summary$replicates_per_scenario <- n_replicates
summary$mode <- if (quick) "quick" else "standard"
summary$protocol <- "direct_coordinate_fields_preserved_index"
summary$elapsed_sec <- elapsed
summary$r_version <- R.version.string
summary$platform <- R.version$platform

utils::write.csv(results, "dev/design_mfa_spatial_calibration_fair.csv", row.names = FALSE)
utils::write.csv(summary, "dev/design_mfa_spatial_calibration_fair_summary.csv", row.names = FALSE)
print(summary, row.names = FALSE)
cat(sprintf("\nCompleted %d replicates per scenario in %.2f seconds.\n", n_replicates, elapsed))
if (!all(summary$passed)) {
  stop("One or more corrected spatial calibration gates failed.", call. = FALSE)
}

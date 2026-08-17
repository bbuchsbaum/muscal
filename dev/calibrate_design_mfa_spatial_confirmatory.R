#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  if (requireNamespace("devtools", quietly = TRUE)) {
    devtools::load_all(".", quiet = TRUE)
  }
  library(muscal)
})

protocol_path <- "dev/design_mfa_spatial_protocol.md"
script_path <- "dev/calibrate_design_mfa_spatial_confirmatory.R"
helper_path <- "tests/testthat/helper-design-mfa-sim.R"
args <- commandArgs(trailingOnly = TRUE)
quick <- any(args %in% c("--quick", "-q"))
recovery_replicates <- if (quick) 6L else 60L
inference_replicates <- if (quick) 8L else 100L
n_permutations <- if (quick) 19L else 199L

if (!quick) {
  tracked <- system2(
    "git", c("ls-files", "--error-unmatch", protocol_path),
    stdout = FALSE, stderr = FALSE
  )
  changed <- system2(
    "git", c("diff", "--quiet", "HEAD", "--", protocol_path),
    stdout = FALSE, stderr = FALSE
  )
  if (!identical(tracked, 0L) || !identical(changed, 0L)) {
    stop(
      "Standard calibration requires the protocol to be tracked, committed, and unchanged.",
      call. = FALSE
    )
  }
}

sha256 <- function(path) {
  output <- system2(
    "shasum", c("-a", "256", path), stdout = TRUE,
    env = c("LC_ALL=C", "LANG=C")
  )
  if (length(output) != 1L) stop(sprintf("Could not hash %s.", path), call. = FALSE)
  sub("[[:space:]].*$", "", output)
}

protocol_git_sha <- system2(
  "git", c("log", "-1", "--format=%H", "--", protocol_path), stdout = TRUE
)
if (length(protocol_git_sha) == 0L) protocol_git_sha <- NA_character_
implementation_files <- sort(Sys.glob("R/design_mfa*.R"))
implementation_hashes <- paste(
  paste0(basename(implementation_files), "=", vapply(implementation_files, sha256, character(1))),
  collapse = ";"
)
provenance <- list(
  evidence_status = if (quick) "developmental_smoke" else "confirmatory",
  protocol_git_sha = protocol_git_sha[[1L]],
  protocol_sha256 = sha256(protocol_path),
  script_sha256 = sha256(script_path),
  helper_sha256 = sha256(helper_path),
  implementation_sha256 = implementation_hashes,
  r_version = R.version.string,
  platform = R.version$platform,
  os = paste(unname(Sys.info()[c("sysname", "release")]), collapse = " "),
  matrix_version = as.character(utils::packageVersion("Matrix")),
  blas = as.character(extSoftVersion()[["BLAS"]])
)

make_reference <- function(p = 30L) {
  t <- seq(-3, 3, length.out = p)
  out <- cbind(x = t, y = sin(1.4 * t), z = cos(1.4 * t))
  rownames(out) <- paste0("reference_", seq_len(p))
  out
}

saliency_field <- function(coords) {
  exp(-rowSums(sweep(coords, 2L, c(-1.1, -0.2, 0.4))^2) / 1.6) -
    0.75 * exp(-rowSums(sweep(coords, 2L, c(1.3, 0.3, -0.3))^2) / 1.2)
}

probe_field <- function(coords) {
  sin(0.8 * coords[, "x"]) + 0.25 * coords[, "y"] - 0.15 * coords[, "z"]
}

simulate_spatial_data <- function(scenario,
                                  seed,
                                  effect_scale = 1,
                                  n_subjects = 8L,
                                  p = 30L) {
  set.seed(seed)
  reference <- make_reference(p)
  reference_saliency <- saliency_field(reference)
  reference_saliency <- reference_saliency - mean(reference_saliency)
  ids <- paste0("S", seq_len(n_subjects))
  trials <- seq.int(36L, 50L, length.out = n_subjects)
  saliency_levels <- c(20, 40, 60, 80, 100)
  probe_levels <- c("old", "new", "lure")
  subject_coords <- design <- Y <- vector("list", n_subjects)

  for (s in seq_len(n_subjects)) {
    coords <- switch(
      scenario,
      identity = reference,
      preserved_index = reference + matrix(rnorm(p * 3L, sd = 0.14), p, 3L),
      independent_coordinates = {
        t <- sort(seq(-3, 3, length.out = p) + rnorm(p, sd = 0.07))
        cbind(
          x = t,
          y = sin(1.4 * t) + rnorm(p, sd = 0.08),
          z = cos(1.4 * t) + rnorm(p, sd = 0.08)
        )
      },
      incomplete_coverage = {
        t <- sort(runif(p, -3, -0.2))
        cbind(
          x = t,
          y = sin(1.4 * t) + rnorm(p, sd = 0.05),
          z = cos(1.4 * t) + rnorm(p, sd = 0.05)
        )
      },
      stop(sprintf("Unknown scenario '%s'.", scenario), call. = FALSE)
    )
    rownames(coords) <- paste0(ids[[s]], "_feature_", seq_len(p))
    subject_coords[[s]] <- coords

    cells <- expand.grid(
      saliency = saliency_levels,
      probe = probe_levels,
      KEEP.OUT.ATTRS = FALSE,
      stringsAsFactors = FALSE
    )
    d <- cells[sample.int(nrow(cells), trials[[s]], replace = TRUE), , drop = FALSE]
    d$saliency <- ordered(d$saliency, levels = saliency_levels)
    d$probe <- factor(d$probe, levels = probe_levels)
    saliency <- as.numeric(d$saliency) - 3
    probe <- c(old = -1, new = 0, lure = 1)[as.character(d$probe)]
    subject_saliency <- saliency_field(coords) - mean(saliency_field(reference))
    subject_probe <- probe_field(coords) - mean(probe_field(reference))
    intercept <- rnorm(p, sd = 0.2)
    response <- effect_scale * outer(saliency, subject_saliency) +
      0.45 * outer(probe, subject_probe) +
      matrix(intercept, nrow(d), p, byrow = TRUE) +
      matrix(rnorm(nrow(d) * p, sd = 0.35), nrow(d), p)
    colnames(response) <- rownames(coords)
    design[[s]] <- d
    Y[[s]] <- response
  }
  names(subject_coords) <- names(design) <- names(Y) <- ids
  list(
    Y = Y,
    design = design,
    subject_coords = subject_coords,
    reference_coords = reference,
    reference_saliency = reference_saliency
  )
}

make_maps <- function(sim, kernel, k, bandwidth = NULL, require_full_reference = TRUE) {
  coordinate_feature_maps(
    sim$reference_coords,
    sim$subject_coords,
    kernel = kernel,
    k = k,
    bandwidth = bandwidth,
    require_full_reference = require_full_reference
  )
}

fit_shared <- function(sim, rank = 1L) {
  Y <- lapply(sim$Y, function(y) {
    colnames(y) <- rownames(sim$reference_coords)
    y
  })
  design_mfa(
    Y, sim$design, ~ saliency + probe,
    contrasts = list(saliency = "poly", probe = "sum"),
    terms = "saliency", rank = rank, store = "fitted"
  )
}

fit_mapped <- function(sim, maps, alignment_penalty = 0.05) {
  design_mfa(
    sim$Y, sim$design, ~ saliency + probe,
    contrasts = list(saliency = "poly", probe = "sum"),
    terms = "saliency", rank = 1L, store = "fitted",
    alignment = "soft", feature_map = maps,
    alignment_penalty = alignment_penalty,
    control = design_mfa_control(
      soft_max_iter = 1500L, soft_tol = 1e-8, feature_chunk_size = 10L
    )
  )
}

truth_effect <- function(effect, reference_saliency) {
  score <- as.numeric(effect$reference_grid$saliency) - 3
  H <- outer(score, reference_saliency)
  sweep(H, 2L, colSums(H * effect$reference_weights), "-")
}

relative_error <- function(effect, truth) {
  numerator <- sum(rowSums((effect$effect_matrix - truth)^2) * effect$reference_weights)
  denominator <- sum(rowSums(truth^2) * effect$reference_weights)
  numerator / denominator
}

run_recovery <- function(replicate) {
  common <- function(scenario_id, scenario) {
    seed <- 202609010L + replicate + 10000L * scenario_id
    list(seed = seed, sim = simulate_spatial_data(scenario, seed))
  }

  identity_case <- common(1L, "identity")
  identity_maps <- lapply(identity_case$sim$Y, function(y) {
    P <- Matrix::Diagonal(ncol(y))
    dimnames(P) <- list(colnames(y), rownames(identity_case$sim$reference_coords))
    P
  })
  names(identity_maps) <- names(identity_case$sim$Y)
  identity_shared <- fit_shared(identity_case$sim)
  identity_soft <- fit_mapped(identity_case$sim, identity_maps, alignment_penalty = 0)
  identity_row <- data.frame(
    scenario = "identity",
    replicate = replicate,
    seed = identity_case$seed,
    primary_error = max(abs(predict(identity_shared) - predict(identity_soft))),
    comparator_error = max(abs(
      effects(identity_shared, "saliency")$effect_matrix -
        effects(identity_soft, "saliency")$effect_matrix
    )),
    loading_correlation = max(abs(
      identity_shared$full_model$coefficient - identity_soft$full_model$coefficient
    )),
    rejected = FALSE,
    message = "",
    stringsAsFactors = FALSE
  )

  preserved <- common(2L, "preserved_index")
  preserved_maps <- make_maps(preserved$sim, "gaussian", 3L, 0.35)
  preserved_soft <- fit_mapped(preserved$sim, preserved_maps)
  preserved_shared <- fit_shared(preserved$sim)
  preserved_effect <- effects(preserved_soft, "saliency")
  preserved_truth <- truth_effect(preserved_effect, preserved$sim$reference_saliency)
  preserved_row <- data.frame(
    scenario = "preserved_index",
    replicate = replicate,
    seed = preserved$seed,
    primary_error = relative_error(preserved_effect, preserved_truth),
    comparator_error = relative_error(
      effects(preserved_shared, "saliency"), preserved_truth
    ),
    loading_correlation = abs(cor(
      effect_map(preserved_soft, "saliency"), preserved$sim$reference_saliency
    )),
    rejected = FALSE,
    message = "",
    stringsAsFactors = FALSE
  )

  independent <- common(3L, "independent_coordinates")
  gaussian_maps <- make_maps(independent$sim, "gaussian", 3L, 0.35)
  nearest_maps <- make_maps(independent$sim, "nearest", 1L)
  misspecified_maps <- lapply(gaussian_maps, function(P) {
    out <- P[, rev(seq_len(ncol(P))), drop = FALSE]
    colnames(out) <- colnames(P)
    out
  })
  names(misspecified_maps) <- names(gaussian_maps)
  gaussian <- fit_mapped(independent$sim, gaussian_maps)
  nearest <- fit_mapped(independent$sim, nearest_maps)
  misspecified <- fit_mapped(independent$sim, misspecified_maps)
  gaussian_effect <- effects(gaussian, "saliency")
  independent_truth <- truth_effect(gaussian_effect, independent$sim$reference_saliency)
  independent_row <- data.frame(
    scenario = "independent_coordinates",
    replicate = replicate,
    seed = independent$seed,
    primary_error = relative_error(gaussian_effect, independent_truth),
    comparator_error = relative_error(effects(nearest, "saliency"), independent_truth),
    loading_correlation = abs(cor(
      effect_map(gaussian, "saliency"), independent$sim$reference_saliency
    )),
    rejected = FALSE,
    message = "",
    stringsAsFactors = FALSE
  )
  misspecified_row <- data.frame(
    scenario = "misspecified_map",
    replicate = replicate,
    seed = independent$seed,
    primary_error = relative_error(effects(misspecified, "saliency"), independent_truth),
    comparator_error = relative_error(gaussian_effect, independent_truth),
    loading_correlation = abs(cor(
      effect_map(misspecified, "saliency"), independent$sim$reference_saliency
    )),
    rejected = FALSE,
    message = "",
    stringsAsFactors = FALSE
  )

  incomplete <- common(4L, "incomplete_coverage")
  incomplete_maps <- make_maps(
    incomplete$sim, "gaussian", 3L, 0.35, require_full_reference = FALSE
  )
  incomplete_attempt <- tryCatch(
    fit_mapped(incomplete$sim, incomplete_maps),
    error = identity
  )
  incomplete_message <- if (inherits(incomplete_attempt, "error")) {
    conditionMessage(incomplete_attempt)
  } else {
    ""
  }
  incomplete_row <- data.frame(
    scenario = "incomplete_coverage",
    replicate = replicate,
    seed = incomplete$seed,
    primary_error = NA_real_,
    comparator_error = NA_real_,
    loading_correlation = NA_real_,
    rejected = inherits(incomplete_attempt, "error"),
    message = incomplete_message,
    stringsAsFactors = FALSE
  )

  rbind(identity_row, preserved_row, independent_row, misspecified_row, incomplete_row)
}

wilson_interval <- function(successes, trials, confidence = 0.95) {
  z <- stats::qnorm(1 - (1 - confidence) / 2)
  estimate <- successes / trials
  denominator <- 1 + z^2 / trials
  center <- (estimate + z^2 / (2 * trials)) / denominator
  half <- z * sqrt(estimate * (1 - estimate) / trials + z^2 / (4 * trials^2)) / denominator
  c(lower = max(0, center - half), upper = min(1, center + half))
}

permute_saliency_within_probe <- function(data, term, replicate) {
  out <- data
  out$design <- lapply(out$design, function(d) {
    for (probe in levels(d$probe)) {
      rows <- which(d$probe == probe)
      d$saliency[rows] <- sample(d$saliency[rows], length(rows), replace = FALSE)
    }
    d
  })
  out
}

run_inference <- function(replicate, null) {
  seed <- 202619010L + replicate + if (null) 10000L else 0L
  sim <- simulate_spatial_data(
    "identity", seed, effect_scale = if (null) 0 else 0.7
  )
  fit <- fit_shared(sim, rank = 2L)
  inferred <- infer_design_mfa(
    fit,
    "saliency",
    method = "permutation",
    statistics = c("partial_ss", "singular_values"),
    nrep = n_permutations,
    seed = seed + 500000L,
    exchangeability = permute_saliency_within_probe
  )
  null_roots <- do.call(rbind, lapply(inferred$replicates, `[[`, "singular_values"))
  max_null_root <- apply(null_roots, 1L, max)
  observed_roots <- inferred$observed$singular_values
  max_root_p <- vapply(observed_roots, function(root) {
    (1 + sum(max_null_root >= root)) / (length(max_null_root) + 1)
  }, numeric(1))
  data.frame(
    scenario = if (null) "null" else "planted_effect",
    replicate = replicate,
    seed = seed,
    permutations = n_permutations,
    term_p_value = inferred$summary$partial_ss$p_value,
    minimum_component_p_value = min(inferred$summary$singular_values$p_value),
    minimum_max_root_p_value = min(max_root_p),
    term_reject = inferred$summary$partial_ss$p_value <= 0.05,
    max_root_family_reject = any(max_root_p <= 0.05),
    stringsAsFactors = FALSE
  )
}

append_provenance <- function(data) {
  for (field in names(provenance)) data[[field]] <- provenance[[field]]
  data
}

run_replicates <- function(n, fun, label, ...) {
  out <- vector("list", n)
  report_every <- max(1L, min(10L, n))
  for (i in seq_len(n)) {
    out[[i]] <- fun(i, ...)
    if (i %% report_every == 0L || i == n) {
      message(sprintf("%s: %d/%d replicates complete", label, i, n))
    }
  }
  do.call(rbind, out)
}

started <- proc.time()[["elapsed"]]
recovery <- run_replicates(recovery_replicates, run_recovery, "Recovery")
inference <- rbind(
  run_replicates(
    inference_replicates, run_inference, "Null inference", null = TRUE
  ),
  run_replicates(
    inference_replicates, run_inference, "Power inference", null = FALSE
  )
)
elapsed <- proc.time()[["elapsed"]] - started
recovery <- append_provenance(recovery)
inference <- append_provenance(inference)

identity <- recovery[recovery$scenario == "identity", , drop = FALSE]
preserved <- recovery[recovery$scenario == "preserved_index", , drop = FALSE]
independent <- recovery[recovery$scenario == "independent_coordinates", , drop = FALSE]
misspecified <- recovery[recovery$scenario == "misspecified_map", , drop = FALSE]
incomplete <- recovery[recovery$scenario == "incomplete_coverage", , drop = FALSE]
null <- inference[inference$scenario == "null", , drop = FALSE]
planted <- inference[inference$scenario == "planted_effect", , drop = FALSE]
term_null_interval <- wilson_interval(sum(null$term_reject), nrow(null))
root_null_interval <- wilson_interval(sum(null$max_root_family_reject), nrow(null))

gates <- data.frame(
  gate = c(
    "identity_prediction",
    "identity_effect",
    "identity_coefficient",
    "preserved_index_negative_control",
    "independent_gaussian_noninferiority",
    "independent_loading_recovery",
    "misspecification_degrades",
    "incomplete_coverage_rejected",
    "null_term_wilson",
    "null_max_root_wilson",
    "term_power",
    "max_root_power"
  ),
  value = c(
    max(identity$primary_error),
    max(identity$comparator_error),
    max(identity$loading_correlation),
    mean(preserved$comparator_error) / mean(preserved$primary_error),
    mean(independent$primary_error) / mean(independent$comparator_error),
    stats::median(independent$loading_correlation),
    mean(misspecified$primary_error) / mean(misspecified$comparator_error),
    mean(incomplete$rejected),
    mean(null$term_reject),
    mean(null$max_root_family_reject),
    mean(planted$term_reject),
    mean(planted$max_root_family_reject)
  ),
  passed = c(
    max(identity$primary_error) <= 1e-7,
    max(identity$comparator_error) <= 1e-7,
    max(identity$loading_correlation) <= 1e-7,
    mean(preserved$comparator_error) < mean(preserved$primary_error),
    mean(independent$primary_error) <= 1.25 * mean(independent$comparator_error),
    stats::median(independent$loading_correlation) >= 0.90,
    mean(misspecified$primary_error) >= mean(misspecified$comparator_error),
    all(incomplete$rejected) && all(grepl("support", incomplete$message, ignore.case = TRUE)),
    term_null_interval[["lower"]] <= 0.05 &&
      term_null_interval[["upper"]] >= 0.05 && term_null_interval[["upper"]] <= 0.12,
    root_null_interval[["lower"]] <= 0.05 &&
      root_null_interval[["upper"]] >= 0.05 && root_null_interval[["upper"]] <= 0.12,
    mean(planted$term_reject) >= 0.80,
    mean(planted$max_root_family_reject) >= 0.70
  ),
  stringsAsFactors = FALSE
)
gates$recovery_replicates <- recovery_replicates
gates$inference_replicates <- inference_replicates
gates$permutations <- n_permutations
gates$term_null_wilson_lower <- term_null_interval[["lower"]]
gates$term_null_wilson_upper <- term_null_interval[["upper"]]
gates$max_root_null_wilson_lower <- root_null_interval[["lower"]]
gates$max_root_null_wilson_upper <- root_null_interval[["upper"]]
gates$elapsed_sec <- elapsed
gates <- append_provenance(gates)

prefix <- if (quick) {
  "dev/design_mfa_spatial_confirmatory_smoke"
} else {
  "dev/design_mfa_spatial_confirmatory"
}
utils::write.csv(recovery, paste0(prefix, "_recovery.csv"), row.names = FALSE)
utils::write.csv(inference, paste0(prefix, "_inference.csv"), row.names = FALSE)
utils::write.csv(gates, paste0(prefix, "_gates.csv"), row.names = FALSE)
print(gates[, c("gate", "value", "passed")], row.names = FALSE)
cat(sprintf("Completed in %.2f seconds.\n", elapsed))
if (!quick && !all(gates$passed)) {
  stop("One or more preregistered confirmatory gates failed.", call. = FALSE)
}

.sim_design_mfa <- function(n_subjects = 4L,
                             trials = c(45L, 52L, 39L, 48L),
                             p = 12L,
                             noise_sd = 0.08,
                             seed = 1L) {
  set.seed(seed)
  if (length(trials) == 1L) trials <- rep(trials, n_subjects)
  stopifnot(length(trials) == n_subjects)

  subject_ids <- paste0("S", seq_len(n_subjects))
  feature_names <- paste0("voxel_", seq_len(p))
  saliency_levels <- c(20, 40, 60, 80, 100)
  probe_levels <- c("old", "new", "lure")

  saliency_map <- scale(seq(-1, 1, length.out = p), center = TRUE, scale = FALSE)[, 1]
  probe_map <- scale(cos(seq(0, 2 * pi, length.out = p)), center = TRUE, scale = FALSE)[, 1]
  nuisance_map <- scale(sin(seq(0, 2 * pi, length.out = p)), center = TRUE, scale = FALSE)[, 1]

  design <- vector("list", n_subjects)
  Y <- vector("list", n_subjects)
  for (s in seq_len(n_subjects)) {
    grid <- expand.grid(
      saliency = saliency_levels,
      probe = probe_levels,
      KEEP.OUT.ATTRS = FALSE,
      stringsAsFactors = FALSE
    )
    idx <- sample.int(nrow(grid), trials[[s]], replace = TRUE)
    ds <- grid[idx, , drop = FALSE]
    ds$saliency <- ordered(ds$saliency, levels = saliency_levels)
    ds$probe <- factor(ds$probe, levels = probe_levels)
    ds$run <- factor(sample(c("run1", "run2"), nrow(ds), replace = TRUE), levels = c("run1", "run2"))
    sal <- as.numeric(ds$saliency) - 3
    probe <- c(old = -1, new = 0, lure = 1)[as.character(ds$probe)]
    run <- ifelse(ds$run == "run2", 1, -1)
    subject_intercept <- rnorm(p, sd = 0.3)
    mu <- outer(sal, saliency_map) +
      0.65 * outer(probe, probe_map) +
      0.15 * outer(run, nuisance_map) +
      matrix(subject_intercept, nrow(ds), p, byrow = TRUE)
    ys <- mu + matrix(rnorm(nrow(ds) * p, sd = noise_sd), nrow(ds), p)
    colnames(ys) <- feature_names
    design[[s]] <- ds
    Y[[s]] <- ys
  }
  names(design) <- names(Y) <- subject_ids

  list(
    Y = Y,
    design = design,
    subjects = subject_ids,
    feature_names = feature_names,
    truth = list(saliency = saliency_map, probe = probe_map, nuisance = nuisance_map)
  )
}

.dmfa_projector <- function(V) {
  if (is.null(V) || ncol(V) == 0L) return(matrix(0, nrow(V), nrow(V)))
  V %*% t(V)
}

.sim_design_mfa_misregistration <- function(n_subjects = 8L,
                                             trials = seq.int(42L, 63L, length.out = n_subjects),
                                             p = 30L,
                                             effect_scale = 1,
                                             noise_sd = 0.35,
                                             coordinate_sd = 0.16,
                                             seed = 1L) {
  set.seed(seed)
  if (length(trials) == 1L) trials <- rep(trials, n_subjects)
  stopifnot(length(trials) == n_subjects)
  ids <- paste0("S", seq_len(n_subjects))
  reference_names <- paste0("reference_", seq_len(p))
  angle <- seq(0, 4 * pi, length.out = p)
  reference_coords <- cbind(
    x = seq(-3, 3, length.out = p),
    y = sin(angle),
    z = cos(angle)
  )
  rownames(reference_coords) <- reference_names
  saliency_field <- function(coords) {
    exp(-rowSums(sweep(coords, 2, c(-1.2, -0.3, 0.5))^2) / 1.8) -
      0.8 * exp(-rowSums(sweep(coords, 2, c(1.4, 0.4, -0.4))^2) / 1.4)
  }
  probe_field <- function(coords) {
    sin(0.7 * coords[, "x"]) + 0.35 * coords[, "y"] - 0.2 * coords[, "z"]
  }
  saliency_center <- mean(saliency_field(reference_coords))
  probe_center <- mean(probe_field(reference_coords))
  saliency_map <- saliency_field(reference_coords) - saliency_center
  probe_map <- probe_field(reference_coords) - probe_center

  subject_coords <- design <- Y <- vector("list", n_subjects)
  saliency_levels <- c(20, 40, 60, 80, 100)
  probe_levels <- c("old", "new", "lure")
  for (s in seq_len(n_subjects)) {
    observed <- reference_coords + matrix(rnorm(p * 3L, sd = coordinate_sd), p, 3L)
    # Preserve the known same-index correspondence.  Reordering features here
    # would destroy the baseline by construction rather than simulate spatial
    # uncertainty.
    rownames(observed) <- paste0(ids[[s]], "_feature_", seq_len(p))
    subject_coords[[s]] <- observed

    base_grid <- expand.grid(
      saliency = saliency_levels,
      probe = probe_levels,
      KEEP.OUT.ATTRS = FALSE,
      stringsAsFactors = FALSE
    )
    d <- base_grid[sample.int(nrow(base_grid), trials[[s]], replace = TRUE), , drop = FALSE]
    d$saliency <- ordered(d$saliency, levels = saliency_levels)
    d$probe <- factor(d$probe, levels = probe_levels)
    design[[s]] <- d
  }
  names(subject_coords) <- names(design) <- ids
  maps <- coordinate_feature_maps(
    reference_coords,
    subject_coords,
    k = 3,
    bandwidth = 0.32,
    kernel = "gaussian"
  )

  for (s in seq_len(n_subjects)) {
    d <- design[[s]]
    saliency <- as.numeric(d$saliency) - 3
    probe <- c(old = -1, new = 0, lure = 1)[as.character(d$probe)]
    observed_saliency <- saliency_field(subject_coords[[s]]) - saliency_center
    observed_probe <- probe_field(subject_coords[[s]]) - probe_center
    subject_intercept <- rnorm(p, sd = 0.2)
    # Generate fields directly from coordinates.  The fitted interpolation map
    # is not used to manufacture the data it will later be judged on.
    observed_mean <- effect_scale * outer(saliency, observed_saliency) +
      0.45 * effect_scale * outer(probe, observed_probe)
    ys <- observed_mean + matrix(subject_intercept, nrow(d), p, byrow = TRUE) +
      matrix(rnorm(nrow(d) * p, sd = noise_sd), nrow(d), p)
    colnames(ys) <- rownames(maps[[s]])
    Y[[s]] <- ys
  }
  names(Y) <- ids
  Y_naive <- lapply(Y, function(ys) {
    colnames(ys) <- reference_names
    ys
  })

  list(
    Y = Y,
    Y_naive = Y_naive,
    design = design,
    maps = maps,
    reference_coords = reference_coords,
    truth = list(saliency = saliency_map, probe = probe_map)
  )
}

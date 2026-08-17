test_that("design_mfa bootstrap resamples intact subjects and is reproducible", {
  sim <- .sim_design_mfa(n_subjects = 4, trials = c(24, 30, 27, 33), p = 7, seed = 40)
  fit <- design_mfa(
    sim$Y,
    sim$design,
    ~ saliency + probe,
    nuisance = ~ run,
    contrasts = list(saliency = "poly", probe = "sum", run = "sum"),
    rank = c(saliency = 2, probe = 2)
  )

  set.seed(9041)
  caller_seed <- .Random.seed
  boot1 <- infer_design_mfa(fit, "saliency", nrep = 6, seed = 123)
  expect_identical(.Random.seed, caller_seed)
  boot2 <- infer_design_mfa(fit, "saliency", nrep = 6, seed = 123)

  expect_s3_class(boot1, "design_mfa_inference")
  expect_identical(boot1$validity, "subject_cluster_bootstrap")
  expect_equal(boot1$n_success, 6)
  expect_equal(boot1$summary$partial_ss, boot2$summary$partial_ss)
  expect_equal(boot1$summary$singular_values, boot2$summary$singular_values)
  expect_equal(boot1$summary$subspace, boot2$summary$subspace)
  expect_true(all(vapply(boot1$sampled_subjects, length, integer(1)) == 4L))
  expect_true(all(vapply(boot1$replicates, function(z) is.finite(z$subspace), logical(1))))
})

test_that("optional bootstrap loading summaries respect the reference shape", {
  sim <- .sim_design_mfa(n_subjects = 4, trials = 30, p = 6, seed = 41)
  fit <- design_mfa(
    sim$Y,
    sim$design,
    ~ saliency,
    contrasts = list(saliency = "poly"),
    rank = 1
  )
  boot <- infer_design_mfa(
    fit,
    "saliency",
    statistics = c("loadings", "subspace"),
    nrep = 5,
    seed = 44
  )

  expect_equal(dim(boot$summary$loadings$mean), c(6L, 1L))
  expect_equal(dim(boot$summary$loadings$sd), c(6L, 1L))
  expect_true(all(boot$summary$loadings$sign_consistency >= 0.5))
})

test_that("permutation inference requires and records explicit exchangeability", {
  sim <- .sim_design_mfa(n_subjects = 4, trials = 28, p = 6, seed = 42)
  fit <- design_mfa(
    sim$Y,
    sim$design,
    ~ saliency,
    contrasts = list(saliency = "poly"),
    rank = 1
  )
  expect_error(
    infer_design_mfa(fit, "saliency", method = "permutation", nrep = 4),
    "explicit `exchangeability`"
  )

  within_subject_shuffle <- function(data, term, replicate) {
    out <- data
    out$design <- lapply(out$design, function(d) {
      d$saliency <- sample(d$saliency)
      d
    })
    out
  }
  perm <- infer_design_mfa(
    fit,
    "saliency",
    method = "permutation",
    exchangeability = within_subject_shuffle,
    nrep = 5,
    seed = 45
  )

  expect_identical(perm$validity, "user_supplied_exchangeability")
  expect_identical(perm$statistics, c("partial_ss", "singular_values"))
  expect_equal(perm$n_success, 5)
  null <- vapply(perm$replicates, `[[`, numeric(1), "partial_ss")
  expected_p <- (1 + sum(null >= perm$observed$partial_ss)) / (length(null) + 1)
  expect_equal(perm$summary$partial_ss$p_value, expected_p)
  expect_true(
    perm$summary$partial_ss$p_value %in%
      (seq_len(length(null) + 1) / (length(null) + 1))
  )
})

test_that("bootstrap refits freeze discriminating observed-cell weights", {
  set.seed(4201)
  ids <- paste0("S", 1:4)
  counts <- list(c(2, 10), c(4, 8), c(7, 5), c(11, 1))
  design <- lapply(seq_along(ids), function(s) {
    data.frame(condition = factor(
      rep(c("A", "B"), counts[[s]]),
      levels = c("A", "B")
    ))
  })
  Y <- lapply(seq_along(ids), function(s) {
    condition <- ifelse(design[[s]]$condition == "B", 1, -1)
    out <- outer(condition, c(-1, 0.5, 1.5, 2)) +
      matrix(rnorm(length(condition) * 4, sd = 0.05), length(condition), 4)
    colnames(out) <- paste0("v", 1:4)
    out
  })
  names(Y) <- names(design) <- ids
  fit <- design_mfa(
    Y, design, ~ condition, rank = 1, reference_weight = "observed"
  )
  spec <- muscal:::.muscal_resolve_refit_spec(fit)

  expect_identical(spec$data$reference_grid, fit$design$reference_grid)
  expect_identical(spec$data$reference_weight, fit$design$reference_weights)

  sampled <- NULL
  for (seed in seq_len(100L)) {
    set.seed(seed)
    candidate <- spec$bootstrap_fn(spec$data)
    if (!identical(sort(attr(candidate, "sampled_subjects")), sort(ids))) {
      sampled <- candidate
      break
    }
  }
  expect_false(is.null(sampled))
  refit <- spec$fit_fn(sampled)
  expect_identical(refit$design$reference_grid, fit$design$reference_grid)
  expect_equal(refit$design$reference_weights, fit$design$reference_weights, tolerance = 1e-15)

  rederived <- design_mfa(
    sampled$Y,
    sampled$design,
    ~ condition,
    subject = names(sampled$Y),
    trial_weights = sampled$trial_weights,
    block_weight = sampled$block_weight,
    reference_grid = fit$design$reference_grid,
    reference_weight = "observed",
    rank = 1
  )
  expect_false(isTRUE(all.equal(
    rederived$design$reference_weights,
    fit$design$reference_weights,
    tolerance = 1e-15
  )))
})

test_that("generic infer_muscal can use the subject refit contract for scalar effects", {
  sim <- .sim_design_mfa(n_subjects = 4, trials = 24, p = 6, seed = 43)
  fit <- design_mfa(
    sim$Y,
    sim$design,
    ~ saliency,
    contrasts = list(saliency = "poly"),
    rank = 1
  )
  out <- infer_muscal(
    fit,
    method = "bootstrap",
    statistic_fn = function(z) c(partial_ss = effects(z, "saliency")$partial_ss),
    nrep = 4,
    seed = 46
  )

  expect_s3_class(out, "muscal_bootstrap_result")
  expect_equal(out$n_success, 4)
  expect_equal(out$summary$label, "partial_ss")
})

test_that("subspace distance uses small metric-aware cross-products", {
  set.seed(44)
  q1 <- qr.Q(qr(matrix(rnorm(30 * 3), 30, 3)))
  q2 <- qr.Q(qr(matrix(rnorm(30 * 3), 30, 3)))
  dense_distance <- norm(tcrossprod(q1) - tcrossprod(q2), type = "F") / sqrt(6)
  expect_equal(.dmfa_subspace_distance(q1, q2), dense_distance, tolerance = 1e-12)

  metric <- crossprod(matrix(rnorm(30 * 30), 30, 30)) + diag(30)
  metric_chol <- chol(metric)
  v1 <- backsolve(metric_chol, q1)
  v2 <- backsolve(metric_chol, q2)
  expect_equal(crossprod(v1, metric %*% v1), diag(3), tolerance = 1e-10)
  expect_equal(
    .dmfa_subspace_distance(v1, v2, feature_metric = metric),
    dense_distance,
    tolerance = 1e-10
  )
})

test_that("seeded inference restores present and absent RNG state on every exit", {
  sim <- .sim_design_mfa(n_subjects = 4, trials = 24, p = 5, seed = 441)
  fit <- design_mfa(sim$Y, sim$design, ~ saliency, rank = 1)
  original_exists <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (original_exists) original_seed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit({
    if (original_exists) {
      assign(".Random.seed", original_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)

  check_state <- function(seed_present, callback, errors = FALSE) {
    if (seed_present) {
      set.seed(9941)
      before <- .Random.seed
    } else {
      if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        rm(".Random.seed", envir = .GlobalEnv)
      }
      before <- NULL
    }
    if (errors) expect_error(callback()) else expect_no_error(callback())
    if (seed_present) {
      expect_identical(.Random.seed, before)
    } else {
      expect_false(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
    }
  }

  for (present in c(FALSE, TRUE)) {
    check_state(present, function() {
      out <- infer_design_mfa(
        fit, "saliency", statistics = "partial_ss", nrep = 2, seed = 12.0
      )
      expect_identical(out$seed, 12L)
    })
    check_state(present, function() {
      infer_design_mfa(
        fit,
        "saliency",
        method = "permutation",
        statistics = "partial_ss",
        nrep = 2,
        seed = 13.0,
        exchangeability = function(data, term, replicate) stop("forced exchangeability failure")
      )
    }, errors = TRUE)
  }
})

test_that("failed inference replicates warn and enter the permutation denominator", {
  sim <- .sim_design_mfa(n_subjects = 4, trials = 24, p = 5, seed = 47)
  fit <- design_mfa(sim$Y, sim$design, ~ saliency, rank = 1)
  orig <- fit$fit_spec$refit$fit_fn
  patch_second_failure <- function(object) {
    counter <- new.env(parent = emptyenv())
    counter$i <- 0L
    object$fit_spec$refit$fit_fn <- function(data) {
      counter$i <- counter$i + 1L
      if (identical(counter$i, 2L)) stop("injected failure")
      orig(data)
    }
    object
  }

  boot_fit <- patch_second_failure(fit)
  boot <- expect_warning(
    infer_design_mfa(
      boot_fit, "saliency", statistics = "partial_ss", nrep = 4, seed = 48
    ),
    "1 of 4 replicates failed"
  )
  expect_equal(boot$n_success, 3)
  expect_equal(length(boot$failures), 1L)
  expect_identical(boot$failures[[1L]]$stage, "fit")
  expect_match(boot$failures[[1L]]$message, "injected failure")

  within_subject_shuffle <- function(data, term, replicate) {
    out <- data
    out$design <- lapply(out$design, function(d) {
      d$saliency <- sample(d$saliency)
      d
    })
    out
  }
  perm_fit <- patch_second_failure(fit)
  perm <- expect_warning(
    infer_design_mfa(
      perm_fit,
      "saliency",
      method = "permutation",
      statistics = "partial_ss",
      exchangeability = within_subject_shuffle,
      nrep = 4,
      seed = 49
    ),
    "successful replicates"
  )
  expect_equal(perm$n_success, 3)
  null <- vapply(perm$replicates, `[[`, numeric(1), "partial_ss")
  expect_equal(length(null), 3L)
  expect_equal(
    perm$summary$partial_ss$p_value,
    (1 + sum(null >= perm$observed$partial_ss)) / (perm$n_success + 1)
  )
})

test_that("serial inference progress and optional future refits keep the contract", {
  sim <- .sim_design_mfa(n_subjects = 4, trials = 22, p = 5, seed = 50)
  fit <- design_mfa(sim$Y, sim$design, ~ saliency, rank = 1)

  expect_message(
    infer_design_mfa(
      fit, "saliency", statistics = "partial_ss", nrep = 2, seed = 51, verbose = TRUE
    ),
    "replicate 1 / 2"
  )

  skip_if_not_installed("future")
  old_plan <- future::plan(future::sequential)
  on.exit(future::plan(old_plan), add = TRUE)
  fut <- infer_design_mfa(
    fit,
    "saliency",
    statistics = "partial_ss",
    nrep = 3,
    seed = 52,
    use_future = TRUE
  )
  expect_equal(fut$n_success, 3)
  expect_s3_class(fut, "design_mfa_inference")
  expect_error(
    infer_design_mfa(fit, "saliency", nrep = 3, use_future = NA),
    "TRUE or FALSE"
  )
})

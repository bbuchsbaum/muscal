test_that("variable weights reproduce rescaled data at fit and projection", {
  set.seed(391)
  weights <- c(a = 4, b = 0.25, c = 1)
  X <- replicate(8, matrix(rnorm(60), 20, 3, dimnames = list(NULL, names(weights))),
                 simplify = FALSE)
  raw <- lapply(X, cov)
  scaled <- lapply(X, function(x) cov(sweep(x, 2, sqrt(weights), "*")))
  train <- 1:6
  held <- 7:8
  subject <- rep(1:4, each = 2)
  condition <- rep(c("rest", "task"), 4)

  for (center in c(FALSE, TRUE)) {
    for (method in c("none", "frobenius", "mfa")) {
      fit <- covstatis(raw[train], weights = rev(weights), dcenter = center,
                       norm_method = method)
      reference <- covstatis(scaled[train], dcenter = center, norm_method = method)
      expect_identical(fit$weights, weights)
      expect_equal(fit$table_matrix, reference$table_matrix)
      expect_equal(fit$compromise_full, reference$compromise_full)
      expect_equal(project_cov(fit, raw[[1L]]), fit$partial_scores[[1L]])
      expect_equal(project_cov(fit, raw[[7L]]), project_cov(reference, scaled[[7L]]))
      expect_equal(project_subjects(fit, raw[held]), project_subjects(reference, scaled[held]))

      paired <- covstatis_rm(raw[train], subject[train], condition[train], weights = weights,
                             dcenter = center, norm_method = method)
      paired_ref <- covstatis_rm(scaled[train], subject[train], condition[train],
                                 dcenter = center, norm_method = method)
      expect_equal(paired$common$compromise_full, paired_ref$common$compromise_full)
      expect_equal(paired$condition_effects, paired_ref$condition_effects)
      expect_equal(project_subjects(paired, raw[held], subject[held], condition[held]),
                   project_subjects(paired_ref, scaled[held], subject[held], condition[held]))
    }
  }

  # The normalization shortcut and explicit override keep their usual precedence.
  for (normalize in c(FALSE, TRUE)) {
    fit <- covstatis(raw[train], normalize = normalize, weights = weights)
    ref <- covstatis(scaled[train], normalize = normalize)
    expect_equal(project_cov(fit, raw[[7L]]), project_cov(ref, scaled[[7L]]))
  }
  fit <- covstatis(raw[train], normalize = FALSE, norm_method = "mfa", weights = weights)
  ref <- covstatis(scaled[train], norm_method = "mfa")
  expect_equal(project_cov(fit, raw[[7L]]), project_cov(ref, scaled[[7L]]))
})

test_that("weights enforce variable identities and preserve positional inputs", {
  S <- list(diag(c(2, 3, 4)), diag(c(7, 2, 1)))
  for (bad in list(1, c(1, -1, 0), c(1, Inf, 2), c(1, NA, 2), c(0, 0, 0),
                   matrix(1, 3, 1), letters[1:3])) {
    expect_error(covstatis(S, weights = bad), "weight")
  }
  w <- c(a = 4, b = 1, c = 0)
  expect_error(covstatis(S, weights = w), "variable names")
  positional <- covstatis(S, weights = unname(w), normalize = FALSE, dcenter = FALSE)
  expect_equal(positional$compromise_full[3, ], c(0, 0, 0))
  named <- covstatis(S, labels = names(w), weights = rev(w), normalize = FALSE, dcenter = FALSE)
  expect_equal(named$compromise_full, positional$compromise_full)
  expect_equal(project_cov(named, S[[1L]]), named$partial_scores[[1L]])

  S <- lapply(S, function(s) { dimnames(s) <- list(names(w), names(w)); s })
  expect_error(covstatis(S, weights = c(a = 1, a = 1, c = 1)), "variable names")
  expect_error(covstatis(S, weights = c(a = 1, b = 1, d = 1)), "variable names")
  wrong_order <- S[[2L]][3:1, 3:1]
  expect_error(covstatis(list(S[[1L]], wrong_order), weights = w), "names and order")
  fit <- covstatis(S, labels = LETTERS[1:3], weights = w)
  expect_identical(fit$labels, LETTERS[1:3])
  expect_error(project_cov(fit, wrong_order), "names and order")
  expect_equal(project_cov(fit, S[[1L]]), fit$partial_scores[[1L]])

  ordinary <- covstatis(S)
  unit <- covstatis(S, weights = rep(1, 3))
  expect_equal(ordinary$compromise_full, unit$compromise_full)
  expect_equal(project_subjects(ordinary, S), project_subjects(unit, S))
})

test_that("normalization cancels a common weight multiplier and rejects zero tables", {
  S <- list(diag(c(2, 3, 4)), diag(c(7, 2, 1)))
  held <- matrix(c(3, 1, 0, 1, 2, 0, 0, 0, 1), 3)
  w <- c(4, 1, 0.25)
  for (method in c("frobenius", "mfa")) {
    reference <- covstatis(S, weights = w, norm_method = method)
    for (multiplier in c(2^-600, 2^600)) {
      fit <- covstatis(S, weights = w * multiplier, norm_method = method)
      expect_equal(fit$table_matrix, reference$table_matrix)
      expect_equal(project_cov(fit, held), project_cov(reference, held))
    }
    zero_after_weighting <- diag(c(0, 0, 2))
    expect_error(covstatis(list(zero_after_weighting), weights = c(1, 1, 0),
                           norm_method = method), "positive")
    expect_error(project_cov(reference, matrix(0, 3, 3)), "positive")
  }
})

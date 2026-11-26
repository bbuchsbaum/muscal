test_that("a paired contrast recovers a covariance change outside the common axis", {
  subject <- rep(c("a", "b", "c"), each = 2)
  condition <- rep(c("rest", "task"), 3)
  changes <- c(0.25, 0.5, 0.75)
  S <- lapply(seq_along(subject), function(i) {
    s <- match(subject[i], unique(subject))
    m <- diag(c(10 + s, 2, 2))
    dimnames(m) <- list(letters[1:3], letters[1:3])
    if (condition[i] == "task") m[2, 3] <- m[3, 2] <- changes[s]
    m
  })
  fit <- covstatis_rm(S, subject, condition, ncomp = 1)
  expect_s3_class(fit$common, "covstatis")
  expect_equal(rownames(fit$common$roi_scores), letters[1:3])
  expect_equal(abs(as.numeric(coef(fit$common))), c(1, 0, 0))
  means <- lapply(seq(1, 5, 2), function(i) (S[[i]] + S[[i+1]]) / 2)
  ordinary <- covstatis(means, normalize = FALSE, dcenter = FALSE)
  expect_equal(fit$common$compromise_full, ordinary$compromise_full)

  delta <- fit$condition_effects$task - fit$condition_effects$rest
  target <- matrix(0, 3, 3, dimnames = dimnames(S[[1]]))
  target[2, 3] <- target[3, 2] <- mean(changes)
  expect_equal(delta, target)
  expect_equal(eigen(delta, symmetric = TRUE)$values, c(0.5, 0, -0.5))
  expect_length(fit$condition_maps, 1L)
  expect_equal(sum(fit$condition_maps[[1]]^2), 1)
  expect_equal(unname(abs(diff(fit$condition_scores[, 1]))), sqrt(sum(target^2)))
  expect_equal(abs(fit$paired_scores[condition == "task", 3]), changes / sqrt(2))
  expect_equal(fit$sums_of_squares[["condition"]], 0.75)
  expect_equal(fit$sums_of_squares[["residual"]], 0.125)

  train <- subject != "c"
  trained <- covstatis_rm(S[train], subject[train], condition[train])
  new <- project_subjects(trained, rev(S[!train]), subject[!train], rev(condition[!train]))
  expected <- sum((S[[6]] - S[[5]]) * trained$condition_maps[[1]]) / 2
  expect_equal(new$paired_scores[, 3], c(expected, -expected))
  expect_equal(new$paired_scores$condition, c("task", "rest"))
  expect_equal(new$common_scores$c, project_cov(trained$common, means[[3]]))
})

test_that("balanced condition axes agree with a full-matrix design model", {
  set.seed(371)
  p <- 4
  subject <- rep(LETTERS[1:5], each = 3)
  condition <- factor(rep(c("lo", "mid", "hi"), 5), levels = c("hi", "lo", "mid"))
  S <- lapply(subject, function(id) crossprod(matrix(rnorm(8 * p), 8, p)) + diag(p))
  shuffle <- sample(seq_along(S))
  S <- S[shuffle]
  subject <- subject[shuffle]
  condition <- condition[shuffle]
  fit <- covstatis_rm(S, subject, condition)

  # Ordinary multivariate least squares on full (duplicated off-diagonal) entries.
  Y <- t(vapply(S, as.numeric, numeric(p^2)))
  reduced <- lm.fit(model.matrix(~ subject), Y)
  full <- lm.fit(model.matrix(~ subject + condition), Y)
  effect <- rowsum(full$fitted.values - reduced$fitted.values, condition)[levels(condition), ] / 5
  oracle <- svd(effect)
  maps <- vapply(fit$condition_maps, as.numeric, numeric(p^2))
  expect_equal(unname(abs(crossprod(oracle$v[, 1:2], maps))), diag(2))
  expect_equal(fit$condition_scores %*% t(maps), effect)
  expect_equal(unname(as.matrix(fit$paired_scores[, 3:4])), unname(reduced$residuals %*% maps))
  expect_equal(unname(fit$explained), (oracle$d^2 / sum(oracle$d^2))[1:2])
  expect_equal(sum(fit$sums_of_squares), sum(scale(Y, scale = FALSE)^2))
  expect_equal(fit$sums_of_squares[["residual"]], sum(full$residuals^2))
  limited <- covstatis_rm(S, subject, condition, ncomp_condition = 1)
  expect_equal(limited$explained, fit$explained[1])
  back <- project_subjects(fit, S, subject, condition)
  expect_equal(back$paired_scores, fit$paired_scores)
  expect_equal(back$common_scores, fit$common$partial_scores)
})

test_that("common baselines and cancellation leave no spurious condition axis", {
  m <- matrix(c(2, 0.3, 0.3, 1), 2)
  s <- rep(c("one", "two"), each = 2)
  cond <- rep(c("A", "B"), 2)
  opposite <- covstatis_rm(list(m, 2*m, 2*m, m), s, cond)
  expect_length(opposite$condition_maps, 0)
  expect_equal(dim(opposite$condition_scores), c(2, 0))
  expect_named(opposite$paired_scores, c("subject", "condition"))
  expect_gt(opposite$sums_of_squares[["residual"]], 0)
  expect_equal(opposite$sums_of_squares[["condition"]], 0)
  expect_equal(project_subjects(opposite, list(m, m), c(3, 3), cond[1:2])$paired_scores,
               data.frame(subject = c("3","3"), condition = cond[1:2]))

  cells <- list(m, 3*m, 2*m, 5*m)
  raw <- covstatis_rm(cells, s, cond)
  for (method in c("frobenius", "mfa")) {
    shape <- covstatis_rm(cells, s, cond, norm_method = method)
    expect_length(shape$condition_maps, 0)
  }
  tiny <- covstatis_rm(lapply(cells, function(x) 1e-10 * x), s, cond)
  expect_equal(tiny$condition_maps, raw$condition_maps)
  expect_equal(tiny$paired_scores[, 3] / 1e-10, raw$paired_scores[, 3])
  projected <- project_subjects(tiny, cells[1:2], c(3, 3), cond[1:2])
  expect_equal(projected$paired_scores[, 3], raw$paired_scores[1:2, 3])
  centered <- covstatis_rm(list(m, m + 2, m, m + 2), s, cond, dcenter = TRUE)
  expect_length(centered$condition_maps, 0)
})

test_that("scalar covariances preserve the subject weights and pairing", {
  cells <- lapply(c(2, 3, 3, 5, 4, 7), matrix)
  fit <- covstatis_rm(cells, rep(1:3, each = 2), rep(c("a", "b"), 3))
  expect_equal(dim(fit$common$C), c(3, 3))
  expect_equal(fit$common$alpha, c(2.5, 4, 5.5)/12)
  expect_equal(unname(fit$condition_effects$b - fit$condition_effects$a),
               matrix(2, 1, 1))
  expect_equal(unname(fit$sums_of_squares), c(9, 6, 1))
})

test_that("a weak orthogonal condition contrast is retained above numerical rank", {
  cells <- lapply(seq_len(3), function(i) {
    diag(c(5 + c(-1, 0, 1)[i], 2 + 1e-5 * c(1, -2, 1)[i]))
  })
  fit <- covstatis_rm(rep(cells, 2), rep(c("a", "b"), each = 3), rep(1:3, 2))
  expect_length(fit$condition_maps, 2)
  expect_equal(unname(fit$explained[2] / fit$explained[1]) / 1e-10, 3,
               tolerance = 1e-6)
})

test_that("the repeated-measures design and variable identities are enforced", {
  m <- diag(3)
  dimnames(m) <- list(letters[1:3], letters[1:3])
  S <- rep(list(m), 4)
  s <- c(1, 1, 2, 2)
  cond <- c("a", "b", "a", "b")
  expect_error(covstatis_rm(S[-1], s[-1], cond[-1]), "exactly one")
  expect_error(covstatis_rm(S, s, c("a", "b", "a", "a")), "exactly one")
  expect_error(covstatis_rm(S[1:2], s[1:2], cond[1:2]), "two subjects")
  expect_error(covstatis_rm(S, s, rep("a", 4)), "two conditions")
  expect_error(covstatis_rm(S, s, c(NA, cond[-1])), "nonmissing")
  expect_error(covstatis_rm(S, s, cond, ncomp_condition = 1.5), "positive integer")
  expect_error(covstatis_rm(list(m,m,m,diag(2)), s, cond), "3 by 3")
  asymmetric <- m
  asymmetric[1, 2] <- 1
  expect_error(covstatis_rm(list(m,m,m,asymmetric), s, cond), "symmetric")
  fit <- covstatis_rm(S, s, cond)
  expect_error(project_subjects(fit, list(m,m), c(3,3), c("a","c")), "exactly one")
  expect_error(project_subjects(fit, list(m,m[3:1,3:1]), c(3,3), cond[1:2]), "names and order")
  expect_error(project_subjects(fit, list(unname(m),unname(m)), c(3,3), cond[1:2]), "names and order")
})

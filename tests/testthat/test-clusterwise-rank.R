.clusterwise_rank_fixture <- function() {
  list(
    data = list(small = sweep(diag(10)[, 1:2], 2, c(3, 1), "*"),
                wide = sweep(diag(10)[, 1:5], 2, 5:1, "*")),
    coords = list(cbind(1:2, 0, 0), cbind(1:5, 1, 0))
  )
}

test_that("clusterwise common rank is checked before preprocessing or gradients", {
  f <- .clusterwise_rank_fixture()
  touched <- FALSE
  expect_error(penalized_mfa_clusterwise(f$data, f$coords, ncomp = 3, lambda = 1,
    preproc = {touched <- TRUE; multivarious::center()}), "common rank")
  expect_false(touched)
  reduced <- penalized_mfa_clusterwise(f$data, f$coords, ncomp = 3, lambda = 1,
    rank_action = "reduce", max_iter = 2)
  expect_equal(reduced$requested_ncomp, 3)
  expect_equal(reduced$ncomp_block, rep(2, 2))
  expect_equal(ncol(reduced$v), 2)
  expect_equal(unname(reduced$effective_rank), c(2, 5))
})

test_that("common-rank updates are dimensionally coherent for both optimizers", {
  f <- .clusterwise_rank_fixture()
  for (optimizer in c("gradient", "adam")) {
    for (lambda in c(0, 1)) {
      fit <- penalized_mfa_clusterwise(f$data, f$coords, ncomp = 2,
        lambda = lambda, optimizer = optimizer, max_iter = 2)
      expect_equal(unname(vapply(fit$V_list, ncol, integer(1))), c(2L, 2L))
      expect_true(all(is.finite(fit$obj_values)))
      for (s in seq_along(fit$V_list)) {
        V <- fit$V_list[[s]]
        expect_equal(crossprod(V), diag(2), tolerance = 1e-10)
        if (lambda == 0) {
          reference <- stats::prcomp(f$data[[s]], center = TRUE, scale. = FALSE,
                                     rank. = 2)$rotation
          expect_equal(tcrossprod(V), tcrossprod(reference), tolerance = 1e-9)
        }
      }
    }
  }
})

test_that("effective rank prevents unsupported and random components", {
  f <- .clusterwise_rank_fixture()
  f$data[[1]] <- cbind(1:10, 2 * (1:10))
  expect_error(penalized_mfa_clusterwise(f$data, f$coords, ncomp = 2), "effective rank")
  reduced <- penalized_mfa_clusterwise(f$data, f$coords, ncomp = 2,
    rank_action = "reduce", lambda = 0, max_iter = 1)
  expect_equal(reduced$ncomp_block, c(1, 1))
  f$data[[1]][, ] <- 1
  expect_error(penalized_mfa_clusterwise(f$data, f$coords, ncomp = 1,
    rank_action = "reduce"), "effective rank")
})

test_that("rank and finite-data boundary errors are explicit", {
  f <- .clusterwise_rank_fixture()
  for (rank in list(1.5, NA_real_, Inf, c(1, 2))) {
    expect_error(penalized_mfa_clusterwise(f$data, f$coords, ncomp = rank), "ncomp")
  }
  f$data[[1]][1, 1] <- Inf
  expect_error(penalized_mfa_clusterwise(f$data, f$coords, ncomp = 1), "finite")
})

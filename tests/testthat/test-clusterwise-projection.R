.clusterwise_projection_fixture <- function() {
  Z <- cbind(rep(c(-1, 1), 4) * 3, rep(c(-1, -1, 1, 1), 2) * 2,
              rep(c(-1, 1), each = 4))
  A <- sweep(Z, 2, c(10, 20, 30), "+")
  B <- sweep(Z[, 1:2], 2, c(100, -40), "+")
  colnames(A) <- paste0("a", 1:3)
  colnames(B) <- paste0("b", 1:2)
  list(data = list(A = A, B = B),
       coords = list(cbind(1:3, 0, 0), cbind(1:2, 1, 0)))
}

test_that("clusterwise projections apply fitted training preprocessing exactly once", {
  f <- .clusterwise_projection_fixture()
  for (kind in c("default", "pass", "standardize", "mixed")) {
    preproc <- switch(kind, default = NULL, pass = multivarious::pass(),
      standardize = multivarious::standardize(),
      mixed = list(multivarious::center(), multivarious::standardize()))
    fit <- penalized_mfa_clusterwise(f$data, f$coords, ncomp = 2, lambda = 0,
                                     max_iter = 1, preproc = preproc)
    new <- lapply(f$data, function(X) X[1:3, , drop = FALSE] + 7)
    expected <- lapply(names(f$data), function(block) {
      train <- f$data[[block]]
      X <- new[[block]]
      if (kind != "pass") X <- sweep(X, 2, colMeans(train), "-")
      if (kind == "standardize" || (kind == "mixed" && block == "B")) {
        X <- sweep(X, 2, apply(train, 2, stats::sd), "/")
      }
      X %*% fit$V_list[[block]]
    })
    names(expected) <- names(f$data)
    actual <- multivarious::project(fit, new[2:1])
    expect_equal(actual, expected, tolerance = 1e-10)
    expect_equal(multivarious::project(fit, new$A, block = "A"), expected$A)
    expect_equal(predict(fit, new, type = "scores"), expected)
    # The stored concatenated processor must also retain its fitted transforms.
    joined <- do.call(cbind, new)
    transformed <- multivarious::transform(fit$preproc, joined)
    expect_equal(transformed[, fit$block_indices$A, drop = FALSE] %*% fit$V_list$A,
                 expected$A)
  }
})

test_that("clusterwise reconstruction uses fitted inverse transforms and survives serialization", {
  f <- .clusterwise_projection_fixture()
  fit <- penalized_mfa_clusterwise(f$data, f$coords, ncomp = 2, lambda = 0,
    preproc = multivarious::standardize(), max_iter = 1)
  new <- lapply(f$data, function(X) X[1:2, , drop = FALSE] + 5)
  scores <- multivarious::project(fit, new)
  expected <- lapply(names(new), function(block) {
    processed <- scores[[block]] %*% t(fit$V_list[[block]])
    X <- sweep(processed, 2, apply(f$data[[block]], 2, stats::sd), "*")
    X <- sweep(X, 2, colMeans(f$data[[block]]), "+")
    colnames(X) <- colnames(new[[block]])
    X
  })
  names(expected) <- names(new)
  expect_equal(predict(fit, new, type = "reconstruction"), expected, tolerance = 1e-10)
  path <- tempfile(fileext = ".rds")
  saveRDS(fit, path)
  restored <- readRDS(path)
  unlink(path)
  expect_equal(multivarious::project(restored, new), scores)
  expect_equal(predict(restored, new$A, block = "A", type = "reconstruction"), expected$A)
})

test_that("clusterwise projection refuses ambiguous blocks and feature schemas", {
  f <- .clusterwise_projection_fixture()
  fit <- penalized_mfa_clusterwise(f$data, f$coords, ncomp = 2, lambda = 0)
  expect_error(multivarious::project(fit, f$data$A), "block")
  expect_error(multivarious::project(fit, f$data$A, block = "unknown"), "Unknown block")
  expect_error(multivarious::project(fit, f$data[1]), "complete")
  wrong <- f$data$A[, c(2, 1, 3)]
  expect_error(multivarious::project(fit, wrong, block = "A"), "feature")
  bad <- f$data$A
  bad[1, 1] <- Inf
  expect_error(multivarious::project(fit, bad, block = "A"), "finite")
})

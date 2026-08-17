test_that("spatial generalized SVD agrees with a dense whitening oracle", {
  set.seed(60)
  H <- matrix(rnorm(6 * 10), 6, 10)
  H <- sweep(H, 2, colMeans(H), "-")
  w <- seq_len(6)
  w <- w / sum(w)
  A <- Matrix::bandSparse(
    10,
    k = c(-1, 1),
    diagonals = list(rep(1, 9), rep(1, 9))
  )
  L <- Matrix::Diagonal(x = Matrix::rowSums(A)) - A
  Q <- Matrix::Diagonal(10) + 3 * L

  out <- muscal:::.dmfa_thin_svd(
    H,
    weights = w,
    ncomp = 4,
    rank_tol = 1e-12,
    feature_metric = Q
  )

  qe <- eigen(as.matrix(Q), symmetric = TRUE)
  Q_inv_sqrt <- qe$vectors %*% diag(1 / sqrt(qe$values)) %*% t(qe$vectors)
  Q_sqrt <- qe$vectors %*% diag(sqrt(qe$values)) %*% t(qe$vectors)
  whitened <- H * sqrt(w)
  ref <- svd(whitened %*% Q_inv_sqrt, nu = 4, nv = 4)

  expect_equal(out$d, ref$d[1:4], tolerance = 1e-9)
  expect_equal(
    .dmfa_projector(Q_sqrt %*% out$v),
    .dmfa_projector(ref$v[, 1:4, drop = FALSE]),
    tolerance = 1e-8
  )
  expect_equal(unname(as.matrix(crossprod(out$v, Q %*% out$v))), diag(4), tolerance = 1e-8)
})

test_that("factored spatial SVD matches the dense metric path on a p-by-r solve", {
  set.seed(6011)
  m <- 25L
  p <- 40L
  r <- 4L
  A <- matrix(rnorm(m * r), m, r)
  B <- matrix(rnorm(r * p), r, p)
  colnames(B) <- paste0("f", seq_len(p))
  w <- seq_len(m)
  w <- w / sum(w)
  A <- sweep(A, 2L, as.numeric(crossprod(w, A)), "-")
  H <- A %*% B
  adjacency <- Matrix::bandSparse(
    p,
    k = c(-1, 1),
    diagonals = list(rep(1, p - 1L), rep(1, p - 1L))
  )
  laplacian <- Matrix::Diagonal(x = Matrix::rowSums(adjacency)) - adjacency
  Q <- Matrix::Diagonal(p) + 2 * laplacian

  dense <- muscal:::.dmfa_thin_svd(H, w, ncomp = 3, rank_tol = 1e-12, feature_metric = Q)
  factored <- muscal:::.dmfa_thin_svd_factored(A, B, w, ncomp = 3, rank_tol = 1e-12, feature_metric = Q)

  expect_identical(factored$path, "factored_small_gram")
  expect_identical(factored$gram_dim, c(r, r))
  expect_equal(factored$d, dense$d, tolerance = 1e-9)
  expect_equal(factored$total_energy, dense$total_energy, tolerance = 1e-9)
  expect_equal(
    .dmfa_projector(as.matrix(Q %*% factored$v)),
    .dmfa_projector(as.matrix(Q %*% dense$v)),
    tolerance = 1e-8
  )
  expect_equal(unname(as.matrix(crossprod(factored$v, Q %*% factored$v))), diag(3), tolerance = 1e-8)
})

test_that("shared-grid spatial penalty produces metric-orthonormal smoother maps", {
  sim <- .sim_design_mfa(n_subjects = 4, trials = 35, p = 20, noise_sd = 0.25, seed = 61)
  A <- Matrix::bandSparse(
    20,
    k = c(-1, 1),
    diagonals = list(rep(1, 19), rep(1, 19))
  )
  L <- Matrix::Diagonal(x = Matrix::rowSums(A)) - A
  dimnames(L) <- list(sim$feature_names, sim$feature_names)

  plain <- design_mfa(
    sim$Y,
    sim$design,
    ~ saliency,
    contrasts = list(saliency = "poly"),
    rank = 1
  )
  smooth <- design_mfa(
    sim$Y,
    sim$design,
    ~ saliency,
    contrasts = list(saliency = "poly"),
    spatial_penalty = 15,
    spatial_laplacian = L,
    rank = 1
  )

  v_plain <- effect_loadings(plain, "saliency")
  v_smooth <- effect_loadings(smooth, "saliency")
  roughness <- function(v) as.numeric(crossprod(v, L %*% v) / crossprod(v))

  expect_lt(roughness(v_smooth), roughness(v_plain))
  expect_equal(
    unname(as.matrix(crossprod(v_smooth, smooth$feature_metric %*% v_smooth))),
    matrix(1, 1, 1),
    tolerance = 1e-8
  )
  expect_identical(effects(smooth, "saliency")$loading_metric, "spatial")
})

test_that("spatial metric validation rejects invalid Laplacians", {
  sim <- .sim_design_mfa(n_subjects = 3, trials = 25, p = 6, seed = 62)
  expect_error(
    design_mfa(sim$Y, sim$design, ~ saliency, spatial_penalty = 1),
    "requires `spatial_laplacian`"
  )
  expect_error(
    design_mfa(sim$Y, sim$design, ~ saliency, spatial_penalty = -1),
    "nonnegative finite"
  )
  asymmetric <- matrix(0, 6, 6)
  asymmetric[1, 2] <- 1
  expect_error(
    design_mfa(
      sim$Y,
      sim$design,
      ~ saliency,
      spatial_penalty = 1,
      spatial_laplacian = asymmetric
    ),
    "symmetric"
  )
  negative <- -diag(6)
  expect_error(
    design_mfa(
      sim$Y,
      sim$design,
      ~ saliency,
      spatial_penalty = 1,
      spatial_laplacian = negative
    ),
    "positive semidefinite"
  )
})

test_that("large graph Laplacians use a deterministic sparse PSD certificate", {
  p <- 600L
  adjacency <- Matrix::bandSparse(
    p,
    k = c(-1, 1),
    diagonals = list(rep(1, p - 1L), rep(1, p - 1L))
  )
  laplacian <- Matrix::Diagonal(x = Matrix::rowSums(adjacency)) - adjacency
  feature_names <- paste0("v", seq_len(p))
  dimnames(laplacian) <- list(feature_names, feature_names)

  validated <- muscal:::.dmfa_validate_spatial_metric(
    laplacian,
    spatial_penalty = 2,
    feature_names = feature_names
  )

  expect_s4_class(validated$laplacian, "symmetricMatrix")
  expect_equal(nrow(validated$metric), p)
  expect_lte(length(validated$metric@x), 3L * p)
})

test_that("named spatial Laplacians are restored to fitted feature order", {
  sim <- .sim_design_mfa(n_subjects = 3, trials = 28, p = 7, seed = 63)
  A <- Matrix::bandSparse(
    7,
    k = c(-1, 1),
    diagonals = list(rep(1, 6), rep(1, 6))
  )
  L <- Matrix::Diagonal(x = Matrix::rowSums(A)) - A
  dimnames(L) <- list(sim$feature_names, sim$feature_names)
  order <- c(4, 1, 7, 2, 6, 3, 5)

  fit <- design_mfa(
    sim$Y,
    sim$design,
    ~ saliency,
    spatial_penalty = 2,
    spatial_laplacian = L,
    rank = 1
  )
  fit_permuted_L <- design_mfa(
    sim$Y,
    sim$design,
    ~ saliency,
    spatial_penalty = 2,
    spatial_laplacian = L[order, order],
    rank = 1
  )

  expect_equal(
    effect_loadings(fit, "saliency"),
    effect_loadings(fit_permuted_L, "saliency"),
    tolerance = 1e-9
  )
})

test_that("identity soft maps reduce to the shared-grid engine", {
  sim <- .sim_design_mfa(n_subjects = 3, trials = 32, p = 8, seed = 64)
  maps <- lapply(sim$Y, function(y) {
    P <- Matrix::Diagonal(ncol(y))
    dimnames(P) <- list(colnames(y), colnames(y))
    P
  })
  names(maps) <- names(sim$Y)
  args <- list(
    formula = ~ saliency + probe,
    contrasts = list(saliency = "poly", probe = "sum"),
    rank = c(saliency = 2, probe = 2)
  )
  shared <- do.call(design_mfa, c(list(Y = sim$Y, design = sim$design), args))
  soft <- do.call(design_mfa, c(list(
    Y = sim$Y,
    design = sim$design,
    alignment = "soft",
    feature_map = maps,
    control = design_mfa_control(soft_tol = 1e-10)
  ), args))

  expect_equal(predict(soft), predict(shared), tolerance = 1e-8)
  expect_equal(soft$full_model$sse, shared$full_model$sse, tolerance = 1e-8)
  expect_equal(
    .dmfa_projector(effect_loadings(soft, "saliency")),
    .dmfa_projector(effect_loadings(shared, "saliency")),
    tolerance = 1e-7
  )
  expect_true(soft$full_model$converged)
  expect_true(all(diff(soft$full_model$objective_trace) <= 1e-9))
  expect_equal(alignment_coverage(soft)$subject_support, rep(3, 8))
  expect_equal(alignment_coverage(soft)$coverage_fraction, rep(1, 8))
})

test_that("soft residual SSE is stable under large subject offsets", {
  sim <- .sim_design_mfa(n_subjects = 3, trials = c(31, 35, 39), p = 7, seed = 641)
  maps <- lapply(sim$Y, function(y) {
    P <- Matrix::Diagonal(ncol(y))
    dimnames(P) <- list(colnames(y), colnames(y))
    P
  })
  names(maps) <- names(sim$Y)
  args <- list(
    design = sim$design,
    formula = ~ saliency + probe,
    contrasts = list(saliency = "poly", probe = "sum"),
    alignment = "soft",
    feature_map = maps,
    control = design_mfa_control(soft_tol = 1e-10)
  )
  baseline <- do.call(design_mfa, c(list(Y = sim$Y), args))
  shifted_Y <- Map(function(y, offset) y + offset, sim$Y, 1e5 + seq_along(sim$Y))
  names(shifted_Y) <- names(sim$Y)
  shifted <- do.call(design_mfa, c(list(Y = shifted_Y), args))

  expect_equal(shifted$full_model$sse, baseline$full_model$sse, tolerance = 2e-5)
  expect_equal(
    effects(shifted, "saliency")$partial_ss,
    effects(baseline, "saliency")$partial_ss,
    tolerance = 2e-5
  )
})

test_that("soft model comparisons use solver-scale tolerance", {
  model <- list(
    soft_alignment = TRUE,
    solver_residual = 1e-6,
    coefficient = matrix(1, 2, 2)
  )
  tol <- muscal:::.dmfa_comparison_tolerance(
    model,
    model,
    scale = 1,
    control = design_mfa_control(soft_tol = 1e-6)
  )

  expect_gte(tol, 1e-5)
  expect_gt(tol, 100 * .Machine$double.eps)
})

test_that("soft-map conjugate gradients agree with a dense Kronecker oracle", {
  sim <- .sim_design_mfa(n_subjects = 3, trials = 24, p = 5, seed = 65)
  set.seed(66)
  maps <- lapply(seq_along(sim$Y), function(s) {
    ps <- 4L + s
    first <- sample.int(5, ps, replace = TRUE)
    second <- (first %% 5L) + 1L
    P <- Matrix::sparseMatrix(
      i = rep(seq_len(ps), 2L),
      j = c(first, second),
      x = rep(c(0.7, 0.3), each = ps),
      dims = c(ps, 5)
    )
    colnames(P) <- sim$feature_names
    rownames(P) <- paste0("subject_", s, "_feature_", seq_len(ps))
    P
  })
  names(maps) <- names(sim$Y)
  Y_soft <- Map(function(y, P) {
    out <- y %*% t(as.matrix(P))
    colnames(out) <- rownames(P)
    out
  }, sim$Y, maps)
  names(Y_soft) <- names(sim$Y)

  fit <- design_mfa(
    Y_soft,
    sim$design,
    ~ saliency + probe,
    contrasts = list(saliency = "poly", probe = "sum"),
    alignment = "soft",
    feature_map = maps,
    alignment_penalty = 1e-4,
    control = design_mfa_control(soft_tol = 1e-11, soft_max_iter = 1000)
  )

  X <- fit$design$X
  q <- ncol(X)
  p <- length(sim$feature_names)
  normal <- matrix(0, q * p, q * p)
  C <- matrix(0, q, p)
  offsets <- cumsum(c(0L, vapply(Y_soft, nrow, integer(1))))
  for (s in seq_along(Y_soft)) {
    rows <- seq.int(offsets[[s]] + 1L, offsets[[s + 1L]])
    Xs <- X[rows, , drop = FALSE]
    ws <- fit$weights[[s]]
    G <- crossprod(Xs * sqrt(ws))
    H <- as.matrix(crossprod(maps[[s]]))
    normal <- normal + kronecker(H, G)
    C <- C + as.matrix(crossprod(Xs, Y_soft[[s]] * ws) %*% maps[[s]])
  }
  normal <- normal + 1e-4 * diag(q * p)
  oracle <- matrix(solve(normal, as.vector(C)), q, p)

  expect_equal(unname(fit$full_model$coefficient), oracle, tolerance = 2e-7)
  expect_lt(fit$full_model$solver_residual, 1e-8)
  expect_lt(fit$full_model$solver_relative_residual, 1e-10)
  expect_identical(
    fit$full_model$operator_applications,
    fit$full_model$iterations
  )
  expect_false(fit$full_model$map_gram_materialized)

  chunked <- design_mfa(
    Y_soft,
    sim$design,
    ~ saliency + probe,
    contrasts = list(saliency = "poly", probe = "sum"),
    alignment = "soft",
    feature_map = maps,
    alignment_penalty = 1e-4,
    control = design_mfa_control(
      soft_tol = 1e-11,
      soft_max_iter = 1000,
      feature_chunk_size = 2L
    )
  )
  expect_equal(
    chunked$full_model$coefficient,
    fit$full_model$coefficient,
    tolerance = 2e-9
  )
  expect_equal(chunked$full_model$sse, fit$full_model$sse, tolerance = 2e-9)
  expect_identical(chunked$full_model$operator_chunk_size, 2L)
})

test_that("subject-feature permutations are absorbed by matching soft maps", {
  sim <- .sim_design_mfa(n_subjects = 4, trials = 30, p = 9, seed = 67)
  maps <- vector("list", length(sim$Y))
  Y_permuted <- vector("list", length(sim$Y))
  set.seed(68)
  for (s in seq_along(sim$Y)) {
    permutation <- sample.int(9)
    Y_permuted[[s]] <- sim$Y[[s]][, permutation, drop = FALSE]
    subject_names <- paste0("observed_", seq_len(9))
    colnames(Y_permuted[[s]]) <- subject_names
    P <- Matrix::sparseMatrix(
      i = seq_len(9),
      j = permutation,
      x = 1,
      dims = c(9, 9),
      dimnames = list(subject_names, sim$feature_names)
    )
    maps[[s]] <- P
  }
  names(Y_permuted) <- names(maps) <- names(sim$Y)

  reference <- design_mfa(
    sim$Y,
    sim$design,
    ~ saliency,
    contrasts = list(saliency = "poly"),
    rank = 1
  )
  soft <- design_mfa(
    Y_permuted,
    sim$design,
    ~ saliency,
    contrasts = list(saliency = "poly"),
    alignment = "soft",
    feature_map = maps,
    rank = 1,
    control = design_mfa_control(soft_tol = 1e-10)
  )

  expect_equal(predict(soft), predict(reference), tolerance = 1e-8)
  expect_equal(effect_loadings(soft, "saliency"), effect_loadings(reference, "saliency"), tolerance = 1e-7)
})

test_that("soft map validation rejects unsupported or malformed references", {
  sim <- .sim_design_mfa(n_subjects = 3, trials = 20, p = 6, seed = 69)
  maps <- lapply(sim$Y, function(y) {
    P <- Matrix::Diagonal(6)
    dimnames(P) <- list(colnames(y), colnames(y))
    P
  })
  names(maps) <- names(sim$Y)

  bad <- maps
  bad[[1]] <- bad[[1]][, -1]
  expect_error(
    design_mfa(
      sim$Y,
      sim$design,
      ~ saliency,
      alignment = "soft",
      feature_map = bad,
      feature_names = sim$feature_names
    ),
    "must be 6 x 6"
  )
  unsupported <- maps
  unsupported <- lapply(unsupported, function(P) {
    P[, 6] <- 0
    P[6, 5] <- 1
    P
  })
  expect_error(
    design_mfa(sim$Y, sim$design, ~ saliency, alignment = "soft", feature_map = unsupported),
    "no support"
  )
  expect_error(
    design_mfa(
      sim$Y,
      sim$design,
      ~ saliency,
      alignment = "soft",
      feature_map = maps,
      subject_model = "ridge_slopes"
    ),
    "cannot yet be combined"
  )

  expect_error(
    design_mfa(
      sim$Y,
      sim$design,
      ~ saliency + probe,
      alignment = "soft",
      feature_map = maps,
      control = design_mfa_control(soft_max_iter = 1, soft_tol = 1e-15)
    ),
    "relative residual.*tolerance.*RHS norm"
  )
})

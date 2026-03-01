library(muscal)

.safe_center_scale <- function(x) {
  x <- as.numeric(x)
  xc <- x - mean(x)
  sdx <- stats::sd(xc)
  if (!is.finite(sdx) || sdx <= sqrt(.Machine$double.eps)) {
    return(matrix(xc, ncol = 1))
  }
  matrix(xc / sdx, ncol = 1)
}

.sim_bilinear_mixed_data <- function(
    n_subject = 8,
    n_repeat = 3,
    n_seed = 4,
    n_roi = 12,
    r_seed = 2,
    r_roi = 3,
    k_subject = 2) {
  subj_levels <- paste0("S", seq_len(n_subject))
  n_obs <- n_subject * n_repeat

  subject <- rep(subj_levels, each = n_repeat)
  repeat_id <- rep(seq_len(n_repeat), times = n_subject)
  z <- cbind(
    rep_linear = .safe_center_scale(repeat_id),
    rep_quad = .safe_center_scale((repeat_id - mean(repeat_id))^2)
  )

  L_true <- qr.Q(qr(matrix(stats::rnorm(n_seed * r_seed), n_seed, r_seed)))
  R_true <- qr.Q(qr(matrix(stats::rnorm(n_roi * r_roi), n_roi, r_roi)))
  W_true <- matrix(stats::rnorm(r_seed * r_roi * k_subject, sd = 0.4), r_seed * r_roi, k_subject)
  B_true <- matrix(stats::rnorm(r_seed * r_roi * ncol(z), sd = 0.3), r_seed * r_roi, ncol(z))
  t_true <- matrix(stats::rnorm(n_subject * k_subject), n_subject, k_subject)

  X <- vector("list", n_obs)
  for (ii in seq_len(n_obs)) {
    si <- ((ii - 1L) %/% n_repeat) + 1L
    m_vec <- B_true %*% z[ii, ] + W_true %*% t_true[si, ]
    M <- matrix(m_vec, nrow = r_seed, ncol = r_roi)
    X[[ii]] <- L_true %*% M %*% t(R_true) +
      matrix(stats::rnorm(n_seed * n_roi, sd = 0.12), n_seed, n_roi)
  }

  A_true <- matrix(c(0.7, -0.2, 0.4, 0.6), nrow = k_subject, byrow = TRUE)
  y <- t_true %*% A_true + matrix(stats::rnorm(n_subject * ncol(A_true), sd = 0.05), n_subject, ncol(A_true))
  rownames(y) <- subj_levels
  colnames(y) <- paste0("trait", seq_len(ncol(y)))

  row_design <- cbind(
    ap = seq(-1, 1, length.out = n_seed),
    hemi = rep(c(-1, 1), length.out = n_seed)
  )
  colnames(row_design) <- c("ap", "hemi")

  list(
    data = X,
    subject = subject,
    repeat_id = repeat_id,
    z = z,
    y = y,
    row_design = row_design
  )
}

.sim_bilinear_symmetric_data <- function(
    n_subject = 6,
    n_repeat = 3,
    n_node = 10,
    r = 3,
    k_subject = 2) {
  subj_levels <- paste0("S", seq_len(n_subject))
  n_obs <- n_subject * n_repeat

  subject <- rep(subj_levels, each = n_repeat)
  repeat_id <- rep(seq_len(n_repeat), times = n_subject)
  z <- cbind(
    rep_linear = .safe_center_scale(repeat_id),
    rep_quad = .safe_center_scale((repeat_id - mean(repeat_id))^2)
  )

  U_true <- qr.Q(qr(matrix(stats::rnorm(n_node * r), n_node, r)))
  W_true <- matrix(stats::rnorm(r * r * k_subject, sd = 0.35), r * r, k_subject)
  B_true <- matrix(stats::rnorm(r * r * ncol(z), sd = 0.25), r * r, ncol(z))
  t_true <- matrix(stats::rnorm(n_subject * k_subject), n_subject, k_subject)

  X <- vector("list", n_obs)
  for (ii in seq_len(n_obs)) {
    si <- ((ii - 1L) %/% n_repeat) + 1L
    m_vec <- B_true %*% z[ii, ] + W_true %*% t_true[si, ]
    S <- matrix(m_vec, nrow = r, ncol = r)
    S <- (S + t(S)) / 2
    E <- matrix(stats::rnorm(n_node * n_node, sd = 0.08), n_node, n_node)
    E <- (E + t(E)) / 2
    X[[ii]] <- U_true %*% S %*% t(U_true) + E
    X[[ii]] <- (X[[ii]] + t(X[[ii]])) / 2
  }

  A_true <- matrix(c(0.6, -0.3, 0.2, 0.7), nrow = k_subject, byrow = TRUE)
  y <- t_true %*% A_true + matrix(stats::rnorm(n_subject * ncol(A_true), sd = 0.05), n_subject, ncol(A_true))
  rownames(y) <- subj_levels
  colnames(y) <- paste0("trait", seq_len(ncol(y)))

  row_design <- cbind(
    axis1 = seq(-1, 1, length.out = n_node),
    axis2 = sin(seq(0, pi, length.out = n_node))
  )
  colnames(row_design) <- c("axis1", "axis2")

  list(
    data = X,
    subject = subject,
    repeat_id = repeat_id,
    z = z,
    y = y,
    row_design = row_design
  )
}

.chain_laplacian <- function(k) {
  L <- matrix(0, k, k)
  for (i in seq_len(k - 1L)) {
    L[i, i] <- L[i, i] + 1
    L[i + 1L, i + 1L] <- L[i + 1L, i + 1L] + 1
    L[i, i + 1L] <- -1
    L[i + 1L, i] <- -1
  }
  L
}

test_that("bilinear_mixed seed_axis fit returns expected structure", {
  set.seed(10)
  sim <- .sim_bilinear_mixed_data()

  fit <- bilinear_mixed(
    data = sim$data,
    subject = sim$subject,
    z = sim$z,
    y = sim$y,
    mode = "seed_axis",
    r_seed = 2,
    r_roi = 3,
    k_subject = 2,
    lambda_y = 0.5,
    max_iter = 30,
    tol = 1e-5
  )

  expect_s3_class(fit, "bilinear_mixed")
  expect_identical(fit$mode, "seed_axis")
  expect_true(!is.null(fit$axis))
  expect_true(is.null(fit$repeat_head))

  expect_equal(dim(fit$axis$L), c(4, 2))
  expect_equal(dim(fit$axis$R), c(12, 3))
  expect_equal(length(fit$axis$design_maps), ncol(sim$z))
  expect_equal(length(fit$axis$trait_maps), 2)
  expect_equal(dim(fit$axis$t_scores), c(length(unique(sim$subject)), 2))
  expect_true(all(is.finite(fit$axis$objective_trace)))
})

test_that("bilinear_mixed seed_repeat fit includes seed effects and interactions", {
  set.seed(20)
  sim <- .sim_bilinear_mixed_data()

  fit <- bilinear_mixed(
    data = sim$data,
    subject = sim$subject,
    z = sim$z,
    y = sim$y,
    row_design = sim$row_design,
    mode = "seed_repeat",
    include_seed_interactions = TRUE,
    r_seed = 2,
    r_roi = 3,
    k_subject = 2,
    lambda_y = 0.5,
    max_iter = 30,
    tol = 1e-5
  )

  expect_s3_class(fit, "bilinear_mixed")
  expect_identical(fit$mode, "seed_repeat")
  expect_true(is.null(fit$axis))
  expect_true(!is.null(fit$repeat_head))

  expected_cols <- ncol(sim$z) + ncol(sim$row_design) + ncol(sim$z) * ncol(sim$row_design)
  expect_equal(ncol(fit$repeat_head$design_matrix), expected_cols)
  expect_equal(length(fit$repeat_head$roi_design_maps), expected_cols)
  expect_equal(length(fit$repeat_head$roi_design_maps[[1]]), ncol(sim$data[[1]]))
  expect_equal(length(fit$repeat_head$roi_trait_maps), 2)
  expect_true(all(is.finite(fit$repeat_head$objective_trace)))
})

test_that("bilinear_mixed both mode fits both heads with shared ROI basis", {
  set.seed(30)
  sim <- .sim_bilinear_mixed_data()
  L_rep <- .chain_laplacian(3)
  L_des <- diag(ncol(sim$z))

  fit <- bilinear_mixed(
    data = sim$data,
    subject = sim$subject,
    z = sim$z,
    y = sim$y,
    row_design = sim$row_design,
    mode = "both",
    repeat_id = sim$repeat_id,
    laplacian_repeat = L_rep,
    laplacian_design = L_des,
    lambda_repeat = 0.1,
    lambda_design_smooth = 0.1,
    r_seed = 2,
    r_roi = 3,
    k_subject = 2,
    lambda_y = 0.5,
    max_iter = 25,
    tol = 1e-5
  )

  expect_s3_class(fit, "bilinear_mixed")
  expect_true(!is.null(fit$axis))
  expect_true(!is.null(fit$repeat_head))
  expect_equal(fit$axis$R, fit$repeat_head$R, tolerance = 1e-8)
})

test_that("bilinear_mixed supports symmetric n x n connectivity with row covariates", {
  set.seed(40)
  sim <- .sim_bilinear_symmetric_data()

  fit <- bilinear_mixed(
    data = sim$data,
    subject = sim$subject,
    z = sim$z,
    y = sim$y,
    row_design = sim$row_design,
    mode = "both",
    connectivity_type = "symmetric",
    r_seed = 3,
    r_roi = 3,
    k_subject = 2,
    lambda_y = 0.5,
    max_iter = 25,
    tol = 1e-5
  )

  fit_auto <- bilinear_mixed(
    data = sim$data,
    subject = sim$subject,
    z = sim$z,
    y = sim$y,
    mode = "seed_axis",
    connectivity_type = "auto",
    r_seed = 3,
    r_roi = 3,
    k_subject = 2,
    lambda_y = 0.5,
    max_iter = 15,
    tol = 1e-5
  )

  expect_s3_class(fit, "bilinear_mixed")
  expect_identical(fit$connectivity_type, "symmetric")
  expect_true(!is.null(fit$axis))
  expect_true(!is.null(fit$repeat_head))
  expect_equal(fit$axis$L, fit$axis$R, tolerance = 1e-8)
  expect_equal(nrow(fit$repeat_head$row_design), nrow(sim$data[[1]]))
  expect_equal(nrow(fit$axis$design_maps[[1]]), nrow(sim$data[[1]]))
  expect_equal(ncol(fit$axis$design_maps[[1]]), ncol(sim$data[[1]]))

  expect_identical(fit_auto$connectivity_type, "symmetric")
})

test_that("bilinear_mixed_recommend returns usable defaults", {
  set.seed(50)
  sim <- .sim_bilinear_mixed_data(
    n_subject = 6,
    n_repeat = 2,
    n_seed = 4,
    n_roi = 10
  )

  rec <- bilinear_mixed_recommend(
    data = sim$data,
    subject = sim$subject,
    z = sim$z,
    y = sim$y,
    row_design = sim$row_design,
    mode = "auto",
    connectivity_type = "auto",
    profile = "fast"
  )

  expect_true(is.list(rec))
  expect_true(rec$mode %in% c("seed_axis", "seed_repeat", "both"))
  expect_true(rec$connectivity_type %in% c("cross", "symmetric"))
  expect_true(rec$r_seed >= 1)
  expect_true(rec$r_roi >= 1)
  expect_true(rec$k_subject >= 1)
  expect_true(rec$lambda_t > 0)
})

test_that("bilinear_mixed_easy fits with compact API", {
  set.seed(51)
  sim <- .sim_bilinear_mixed_data(
    n_subject = 6,
    n_repeat = 2,
    n_seed = 4,
    n_roi = 10
  )

  fit <- bilinear_mixed_easy(
    data = sim$data,
    subject = sim$subject,
    z = sim$z,
    y = sim$y,
    row_design = sim$row_design,
    profile = "fast"
  )

  expect_s3_class(fit, "bilinear_mixed")
  expect_true(is.list(fit$settings))
  expect_true(fit$settings$lambda_t > 0)
})

test_that("bilinear_mixed_tune returns selection table and refit", {
  set.seed(52)
  sim <- .sim_bilinear_mixed_data(
    n_subject = 6,
    n_repeat = 2,
    n_seed = 4,
    n_roi = 10
  )

  tuned <- bilinear_mixed_tune(
    data = sim$data,
    subject = sim$subject,
    z = sim$z,
    y = sim$y,
    row_design = sim$row_design,
    mode = "both",
    connectivity_type = "cross",
    metric = "reconstruction",
    n_folds = 2,
    grid = data.frame(
      r_seed = c(2, 3),
      r_roi = c(3, 4),
      k_subject = c(1, 2),
      lambda_y = c(0.5, 0.5)
    )
  )

  expect_s3_class(tuned, "bilinear_mixed_tuning")
  expect_true(is.data.frame(tuned$results))
  expect_true(nrow(tuned$results) == 2)
  expect_s3_class(tuned$fit, "bilinear_mixed")
  expect_true(all(is.finite(tuned$results$score[tuned$results$n_success > 0])))
})

test_that("bilinear_mixed complexity shortcut sets adaptive defaults", {
  set.seed(53)
  sim <- .sim_bilinear_mixed_data(
    n_subject = 6,
    n_repeat = 2,
    n_seed = 4,
    n_roi = 10
  )

  fit <- bilinear_mixed(
    data = sim$data,
    subject = sim$subject,
    z = sim$z,
    y = sim$y,
    row_design = sim$row_design,
    mode = "both",
    connectivity_type = "cross",
    complexity = "fast",
    max_iter = 20,
    tol = 1e-5
  )

  expect_s3_class(fit, "bilinear_mixed")
  expect_identical(fit$settings$complexity, "fast")
  expect_true(fit$settings$lambda_t > 0)
  expect_true(fit$dimensions$r_seed >= 1)
  expect_true(fit$dimensions$r_roi >= 1)
})

test_that("print.bilinear_mixed_tuning shows compact summary", {
  set.seed(54)
  sim <- .sim_bilinear_mixed_data(
    n_subject = 6,
    n_repeat = 2,
    n_seed = 4,
    n_roi = 10
  )

  tuned <- bilinear_mixed_tune(
    data = sim$data,
    subject = sim$subject,
    z = sim$z,
    y = sim$y,
    row_design = sim$row_design,
    metric = "reconstruction",
    n_folds = 2,
    grid = data.frame(
      r_seed = c(2, 3),
      r_roi = c(3, 4),
      k_subject = c(1, 2),
      lambda_y = c(0.5, 0.5)
    )
  )

  txt <- capture.output(print(tuned, n = 1))
  expect_true(any(grepl("Bilinear Mixed Tuning", txt)))
  expect_true(any(grepl("Top candidates", txt)))
})

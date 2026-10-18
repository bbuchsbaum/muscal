#!/usr/bin/env Rscript

# Read-only regression reproducer for findings in muscal-review-design.md.
pkgload::load_all(".", compile = FALSE, quiet = TRUE)

# 1. Heterogeneous feature domains are recycled as if coordinatewise comparable.
heterogeneous <- list(
  V_list = list(
    matrix(c(1, 0, 0), 3L, 1L),
    matrix(c(1, 0, 0, 0, 0, 0), 6L, 1L)
  ),
  sdev = 1
)
icc <- suppressWarnings(significant_components(
  heterogeneous, n = 20L, k_vec = c(3L, 6L), check_rmt = FALSE
))
stopifnot(
  identical(unname(icc$icc), 1),
  icc$icc_pvalue < 0.05,
  identical(icc$keep, 1L)
)

# 2. True feature dimensions and estimate_components() take opposite RMT
# branches. The latter substitutes ncol(V), i.e. the component count.
set.seed(1)
fit <- list(
  V_list = list(matrix(rnorm(60), 20L, 3L), matrix(rnorm(60), 20L, 3L)),
  sdev = sqrt(c(9, 1, 0.5)),
  s = matrix(rnorm(30), 10L, 3L)
)
true_dimensions <- suppressWarnings(significant_components(
  fit, n = 10L, k_vec = c(20L, 20L), check_rmt = TRUE
))
wrong_dimensions <- estimate_components(fit, method = "rmt")
stopifnot(
  identical(true_dimensions$rmt_pass, rep(TRUE, 3L)),
  !identical(wrong_dimensions$criterion, true_dimensions$rmt_pass)
)

# 3. In the non-skipped branch, changing gamma does not change the reported
# MP edge because the implementation's scale factors cancel algebraically.
edge_a <- suppressWarnings(significant_components(
  fit, n = 100L, k_vec = c(20L, 20L), check_rmt = TRUE
))
edge_b <- suppressWarnings(significant_components(
  fit, n = 50L, k_vec = c(20L, 20L), check_rmt = TRUE
))
stopifnot(
  !is.na(edge_a$mp_edge),
  identical(edge_a$mp_edge, edge_b$mp_edge)
)

# 4. alpha is a public formal but never read by estimate_components().
stopifnot("alpha" %in% names(formals(estimate_components)))
stopifnot(!("alpha" %in% all.names(body(estimate_components), functions = FALSE)))

# 5. Generic inference accepts an empty observed sdev and counts refits as
# successful, returning a 0-row bootstrap summary.
fake <- list(fit_spec = list(
  refit_supported = TRUE,
  refit = list(
    data = list(),
    fit_fn = function(data) list(sdev = 1),
    bootstrap_fn = function(data) data,
    permutation_fn = function(data) data
  )
))
empty <- infer_muscal(fake, method = "bootstrap", nrep = 2L, seed = 4L)
stopifnot(empty$n_success == 2L, length(empty$observed) == 0L, nrow(empty$summary) == 0L)

# 6. The current live multifer effect-significance route handles rank one with
# both scalar partial sums of squares and a length-one singular-value vector.
source("tests/testthat/helper-design-mfa-sim.R")
sim <- .sim_design_mfa(n_subjects = 4L, trials = 24L, p = 5L, seed = 91L)
rank_one <- design_mfa(sim$Y, sim$design, ~ saliency, rank = 1L)
identity_null <- function(data, term, replicate) data
rank_one_inference <- infer_design_mfa(
  rank_one, "saliency", method = "permutation",
  statistics = c("partial_ss", "singular_values"), nrep = 3L, seed = 92L,
  exchangeability = identity_null
)
stopifnot(
  identical(rank_one_inference$engine, "multifer"),
  rank_one_inference$n_success == 3L,
  identical(names(rank_one_inference$summary), c("partial_ss", "singular_values")),
  nrow(rank_one_inference$summary$singular_values) == 1L
)

cat("All review reproductions passed.\n")

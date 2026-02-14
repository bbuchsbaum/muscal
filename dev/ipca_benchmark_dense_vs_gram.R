#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  if (requireNamespace("devtools", quietly = TRUE)) {
    devtools::load_all(".", quiet = TRUE)
  }
  library(muscal)
})

sim_ipca_blocks <- function(n, r, p_vec, noise_sd = 0.05, seed = 1) {
  set.seed(seed)
  Z <- matrix(rnorm(n * r), nrow = n, ncol = r)
  blocks <- lapply(p_vec, function(p) {
    A <- matrix(rnorm(p * r), nrow = p, ncol = r)
    Z %*% t(A) + matrix(rnorm(n * p, sd = noise_sd), nrow = n, ncol = p)
  })
  names(blocks) <- paste0("B", seq_along(blocks))
  blocks
}

proj_diff <- function(S1, S2) {
  P1 <- S1 %*% solve(crossprod(S1), t(S1))
  P2 <- S2 %*% solve(crossprod(S2), t(S2))
  norm(P1 - P2, type = "F") / (norm(P2, type = "F") + 1e-12)
}

run_once <- function(blocks, method, ncomp = 3, lambda = 1, max_iter = 30, tol = 1e-4) {
  gc()
  t0 <- proc.time()[["elapsed"]]
  fit <- ipca(
    data = blocks,
    ncomp = ncomp,
    lambda = lambda,
    method = method,
    max_iter = max_iter,
    tol = tol
  )
  elapsed <- proc.time()[["elapsed"]] - t0
  list(
    fit = fit,
    elapsed_sec = elapsed,
    converged = isTRUE(fit$converged),
    iter = fit$iter,
    method_used = fit$method_used
  )
}

args <- commandArgs(trailingOnly = TRUE)
quick <- any(args %in% c("--quick", "-q"))

scenarios <- list(
  list(name = "small_dense", n = 120, K = 3, p_vec = c(40, 50, 60), r = 3, noise = 0.05),
  list(name = "medium_dense", n = 200, K = 4, p_vec = c(80, 100, 120, 90), r = 4, noise = 0.05),
  list(name = "highdim_mixed", n = 150, K = 4, p_vec = c(300, 500, 700, 450), r = 4, noise = 0.08),
  list(name = "highdim_pggn", n = 120, K = 3, p_vec = c(800, 1000, 1200), r = 3, noise = 0.08)
)

if (quick) {
  scenarios <- scenarios[seq_len(2)]
}

rows <- list()
for (i in seq_along(scenarios)) {
  sc <- scenarios[[i]]
  cat(sprintf("\nScenario %d/%d: %s (n=%d, K=%d, p=[%s])\n",
              i, length(scenarios), sc$name, sc$n, sc$K, paste(sc$p_vec, collapse = ",")))

  blocks <- sim_ipca_blocks(
    n = sc$n,
    r = sc$r,
    p_vec = sc$p_vec,
    noise_sd = sc$noise,
    seed = 100 + i
  )

  max_iter <- if (quick) 12 else 30
  dense_res <- tryCatch(run_once(blocks, "dense", ncomp = min(3, sc$n), lambda = 1, max_iter = max_iter), error = function(e) e)
  gram_res <- tryCatch(run_once(blocks, "gram", ncomp = min(3, sc$n), lambda = 1, max_iter = max_iter), error = function(e) e)

  dense_time <- if (inherits(dense_res, "error")) NA_real_ else dense_res$elapsed_sec
  gram_time <- if (inherits(gram_res, "error")) NA_real_ else gram_res$elapsed_sec
  speedup <- if (is.finite(dense_time) && is.finite(gram_time) && gram_time > 0) dense_time / gram_time else NA_real_
  subspace_rel <- if (!inherits(dense_res, "error") && !inherits(gram_res, "error")) {
    proj_diff(multivarious::scores(dense_res$fit), multivarious::scores(gram_res$fit))
  } else {
    NA_real_
  }

  rows[[i]] <- data.frame(
    scenario = sc$name,
    n = sc$n,
    K = sc$K,
    p_min = min(sc$p_vec),
    p_max = max(sc$p_vec),
    dense_sec = dense_time,
    gram_sec = gram_time,
    speedup_dense_over_gram = speedup,
    dense_converged = if (inherits(dense_res, "error")) FALSE else dense_res$converged,
    gram_converged = if (inherits(gram_res, "error")) FALSE else gram_res$converged,
    dense_iter = if (inherits(dense_res, "error")) NA_integer_ else dense_res$iter,
    gram_iter = if (inherits(gram_res, "error")) NA_integer_ else gram_res$iter,
    score_subspace_rel_diff = subspace_rel,
    stringsAsFactors = FALSE
  )
}

res <- do.call(rbind, rows)
print(res, row.names = FALSE)

out_path <- "dev/ipca_benchmark_results.csv"
utils::write.csv(res, out_path, row.names = FALSE)
cat(sprintf("\nWrote benchmark table to %s\n", out_path))
if (quick) cat("Quick mode was enabled (--quick).\n")

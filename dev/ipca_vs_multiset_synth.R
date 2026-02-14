#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  if (requireNamespace("muscal", quietly = TRUE)) {
    library(muscal)
  }
})

# Development fallback for environments where load_all/compilation is unavailable.
if (!exists("ipca", mode = "function")) {
  source("R/all_generic.R")
  source("R/utils.R")
  source("R/mfa.R")
  source("R/mcca.R")
  source("R/ipca.R")
  `%>%` <- dplyr::`%>%`
}

parse_args <- function(args) {
  out <- list(quick = FALSE, reps = NULL, seed = 1L, delta = 0.04)
  if (length(args) == 0) return(out)

  i <- 1L
  while (i <= length(args)) {
    a <- args[[i]]
    if (a %in% c("--quick", "-q")) {
      out$quick <- TRUE
      i <- i + 1L
    } else if (a == "--reps" && i < length(args)) {
      out$reps <- as.integer(args[[i + 1L]])
      i <- i + 2L
    } else if (a == "--seed" && i < length(args)) {
      out$seed <- as.integer(args[[i + 1L]])
      i <- i + 2L
    } else if (a == "--delta" && i < length(args)) {
      out$delta <- as.numeric(args[[i + 1L]])
      i <- i + 2L
    } else {
      i <- i + 1L
    }
  }
  out
}

orthonormal_basis <- function(S, r = ncol(S)) {
  Q <- qr.Q(qr(S))
  Q[, seq_len(min(r, ncol(Q))), drop = FALSE]
}

subspace_similarity <- function(S_est, S_true, r) {
  if (is.null(S_est) || is.null(S_true)) return(NA_real_)
  if (!is.matrix(S_est) || !is.matrix(S_true)) return(NA_real_)
  if (nrow(S_est) != nrow(S_true)) return(NA_real_)
  r_eff <- min(r, ncol(S_est), ncol(S_true))
  if (r_eff < 1) return(NA_real_)

  Qe <- orthonormal_basis(S_est, r_eff)
  Qt <- orthonormal_basis(S_true, r_eff)
  sv <- svd(crossprod(Qt, Qe), nu = 0, nv = 0)$d
  mean(pmin(pmax(sv[seq_len(r_eff)], 0), 1))
}

make_spd <- function(p, decay = 0.9) {
  Z <- matrix(rnorm(p * p), p, p)
  Q <- qr.Q(qr(Z))
  vals <- decay^(0:(p - 1))
  Q %*% (vals * t(Q))
}

sim_matrix_normal <- function(n, p_vec, r_signal = 3, decay_sigma = 0.85, decay_delta = 0.95) {
  Sigma <- make_spd(n, decay = decay_sigma)
  eig_sigma <- eigen(Sigma, symmetric = TRUE)
  U_true <- eig_sigma$vectors[, seq_len(r_signal), drop = FALSE]
  Ls <- eig_sigma$vectors %*% diag(sqrt(pmax(eig_sigma$values, 1e-12)), n, n)

  blocks <- lapply(p_vec, function(pk) {
    Delta_k <- make_spd(pk, decay = decay_delta)
    eig_d <- eigen(Delta_k, symmetric = TRUE)
    Ld <- eig_d$vectors %*% diag(sqrt(pmax(eig_d$values, 1e-12)), pk, pk)
    Z <- matrix(rnorm(n * pk), n, pk)
    Ls %*% Z %*% t(Ld)
  })
  names(blocks) <- paste0("B", seq_along(blocks))
  list(blocks = blocks, S_true = U_true)
}

sim_latent <- function(n, p_vec, r = 3, noise_sd = 0.1, heavy_tail = FALSE, imbalance = FALSE) {
  Z <- matrix(rnorm(n * r), n, r)
  Z <- scale(Z, center = TRUE, scale = FALSE)
  Z_true <- orthonormal_basis(Z, r)

  blocks <- vector("list", length(p_vec))
  for (k in seq_along(p_vec)) {
    pk <- p_vec[[k]]
    Bk <- matrix(rnorm(pk * r), pk, r)
    signal <- Z %*% t(Bk)
    sd_k <- noise_sd
    if (imbalance && k == length(p_vec)) {
      sd_k <- noise_sd * 3
    }
    E <- if (heavy_tail) {
      matrix(rt(n * pk, df = 3), n, pk) * sd_k
    } else {
      matrix(rnorm(n * pk, sd = sd_k), n, pk)
    }
    blocks[[k]] <- signal + E
  }
  names(blocks) <- paste0("B", seq_along(blocks))
  list(blocks = blocks, S_true = Z_true)
}

safe_fit <- function(method, blocks, r) {
  t0 <- proc.time()[["elapsed"]]
  fit <- NULL
  err <- NULL

  if (method == "ipca") {
    fit <- tryCatch(
      ipca(blocks, ncomp = r, lambda = 1, method = "auto", max_iter = 80, tol = 1e-5),
      error = function(e) e
    )
  } else if (method == "mcca") {
    fit <- tryCatch(
      mcca(blocks, preproc = multivarious::center(), ncomp = r, ridge = 1e-6),
      error = function(e) e
    )
    if (inherits(fit, "error")) {
      fit <- tryCatch(
        mcca(blocks, preproc = multivarious::center(), ncomp = r, ridge = 1e-3),
        error = function(e) e
      )
    }
  } else if (method == "mfa") {
    fit <- tryCatch(
      mfa(
        blocks,
        preproc = multivarious::center(),
        ncomp = r,
        normalization = "MFA"
      ),
      error = function(e) e
    )
  } else {
    fit <- simpleError("Unknown method.")
  }

  elapsed <- proc.time()[["elapsed"]] - t0

  if (inherits(fit, "error")) {
    return(list(ok = FALSE, fit = NULL, elapsed = elapsed, error = conditionMessage(fit)))
  }
  list(ok = TRUE, fit = fit, elapsed = elapsed, error = NA_character_)
}

bootstrap_ci <- function(x, B = 1000, conf = 0.95) {
  x <- x[is.finite(x)]
  n <- length(x)
  if (n < 2) return(c(mean = if (n == 1) x else NA_real_, lwr = NA_real_, upr = NA_real_))
  bmeans <- replicate(B, mean(sample(x, n, replace = TRUE)))
  alpha <- (1 - conf) / 2
  c(mean = mean(x), lwr = unname(stats::quantile(bmeans, alpha)), upr = unname(stats::quantile(bmeans, 1 - alpha)))
}

main <- function() {
  opts <- parse_args(commandArgs(trailingOnly = TRUE))
  set.seed(opts$seed)

  reps <- if (!is.null(opts$reps)) opts$reps else if (opts$quick) 5L else 20L
  methods <- c("ipca", "mcca", "mfa")
  scenarios <- list(
    list(name = "matrix_normal", type = "matrix_normal", n = 80, p = c(120, 140, 160), r = 2, decay_sigma = 0.98, decay_delta = 0.99),
    list(name = "latent_balanced", type = "latent", n = 120, p = c(70, 80, 90), r = 3, noise = 0.10),
    list(name = "highdim_pggn", type = "latent", n = 90, p = c(350, 450, 550), r = 3, noise = 0.15),
    list(name = "block_imbalance", type = "latent_imbalance", n = 120, p = c(60, 80, 500), r = 3, noise = 0.10),
    list(name = "heavy_tail", type = "latent_heavytail", n = 120, p = c(70, 80, 90), r = 3, noise = 0.15)
  )
  if (opts$quick) scenarios <- scenarios[seq_len(3)]

  if ("mfa" %in% methods && !requireNamespace("genpca", quietly = TRUE)) {
    message("genpca not installed; mfa will be skipped.")
  }

  raw_rows <- list()
  row_idx <- 1L

  for (sc in scenarios) {
    cat(sprintf("\nScenario: %s | reps=%d\n", sc$name, reps))
    for (rep in seq_len(reps)) {
      sim <- if (sc$type == "matrix_normal") {
        sim_matrix_normal(
          n = sc$n,
          p_vec = sc$p,
          r_signal = sc$r,
          decay_sigma = if (!is.null(sc$decay_sigma)) sc$decay_sigma else 0.85,
          decay_delta = if (!is.null(sc$decay_delta)) sc$decay_delta else 0.95
        )
      } else if (sc$type == "latent") {
        sim_latent(n = sc$n, p_vec = sc$p, r = sc$r, noise_sd = sc$noise)
      } else if (sc$type == "latent_imbalance") {
        sim_latent(n = sc$n, p_vec = sc$p, r = sc$r, noise_sd = sc$noise, imbalance = TRUE)
      } else {
        sim_latent(n = sc$n, p_vec = sc$p, r = sc$r, noise_sd = sc$noise, heavy_tail = TRUE)
      }

      blocks <- sim$blocks
      S_true <- sim$S_true
      y <- as.numeric(S_true[, 1] + rnorm(nrow(S_true), sd = 0.5))

      for (m in methods) {
        if (m == "mfa" && !requireNamespace("genpca", quietly = TRUE)) {
          raw_rows[[row_idx]] <- data.frame(
            scenario = sc$name,
            rep = rep,
            method = m,
            ok = FALSE,
            similarity = NA_real_,
            r2_y = NA_real_,
            elapsed_sec = NA_real_,
            converged = NA,
            error = "genpca not installed",
            stringsAsFactors = FALSE
          )
          row_idx <- row_idx + 1L
          next
        }

        fitres <- safe_fit(m, blocks, sc$r)
        if (!fitres$ok) {
          raw_rows[[row_idx]] <- data.frame(
            scenario = sc$name,
            rep = rep,
            method = m,
            ok = FALSE,
            similarity = NA_real_,
            r2_y = NA_real_,
            elapsed_sec = fitres$elapsed,
            converged = NA,
            error = fitres$error,
            stringsAsFactors = FALSE
          )
          row_idx <- row_idx + 1L
          next
        }

        S_est <- multivarious::scores(fitres$fit)
        sim_score <- subspace_similarity(S_est, S_true, sc$r)
        df <- data.frame(y = y, S_est[, seq_len(min(sc$r, ncol(S_est))), drop = FALSE])
        names(df) <- c("y", paste0("s", seq_len(ncol(df) - 1)))
        mod <- stats::lm(y ~ ., data = df)
        r2 <- summary(mod)$adj.r.squared

        raw_rows[[row_idx]] <- data.frame(
          scenario = sc$name,
          rep = rep,
          method = m,
          ok = TRUE,
          similarity = sim_score,
          r2_y = r2,
          elapsed_sec = fitres$elapsed,
          converged = if (!is.null(fitres$fit$converged)) isTRUE(fitres$fit$converged) else NA,
          error = NA_character_,
          stringsAsFactors = FALSE
        )
        row_idx <- row_idx + 1L
      }
    }
  }

  raw <- do.call(rbind, raw_rows)
  raw_path <- "dev/ipca_vs_multiset_raw.csv"
  utils::write.csv(raw, raw_path, row.names = FALSE)

  ok_raw <- raw[raw$ok, , drop = FALSE]
  if (nrow(ok_raw) == 0) {
    cat("\nNo successful fits in synthetic battery.\n")
    cat(sprintf("Raw results (with errors): %s\n", raw_path))
    quit(save = "no", status = 1)
  }
  groups <- unique(ok_raw[c("scenario", "method")])
  summary_rows <- vector("list", nrow(groups))
  for (i in seq_len(nrow(groups))) {
    sc <- groups$scenario[[i]]
    md <- groups$method[[i]]
    sub <- ok_raw[ok_raw$scenario == sc & ok_raw$method == md, , drop = FALSE]
    summary_rows[[i]] <- data.frame(
      scenario = sc,
      method = md,
      similarity_mean = mean(sub$similarity, na.rm = TRUE),
      similarity_sd = stats::sd(sub$similarity, na.rm = TRUE),
      similarity_n = sum(is.finite(sub$similarity)),
      r2_mean = mean(sub$r2_y, na.rm = TRUE),
      r2_sd = stats::sd(sub$r2_y, na.rm = TRUE),
      r2_n = sum(is.finite(sub$r2_y)),
      elapsed_mean = mean(sub$elapsed_sec, na.rm = TRUE),
      elapsed_sd = stats::sd(sub$elapsed_sec, na.rm = TRUE),
      elapsed_n = sum(is.finite(sub$elapsed_sec)),
      stringsAsFactors = FALSE
    )
  }
  summary_df <- do.call(rbind, summary_rows)

  summary_path <- "dev/ipca_vs_multiset_summary.csv"
  utils::write.csv(summary_df, summary_path, row.names = FALSE)

  compare_metric <- function(metric_col, comparator) {
    out <- list()
    idx <- 1L
    for (sc in unique(raw$scenario)) {
      ip <- raw[raw$scenario == sc & raw$method == "ipca", c("rep", metric_col)]
      cp <- raw[raw$scenario == sc & raw$method == comparator, c("rep", metric_col)]
      merged <- merge(ip, cp, by = "rep", suffixes = c("_ipca", "_cmp"))
      d <- merged[[paste0(metric_col, "_ipca")]] - merged[[paste0(metric_col, "_cmp")]]
      d <- d[is.finite(d)]
      ci <- bootstrap_ci(d, B = if (opts$quick) 300 else 1200)
      out[[idx]] <- data.frame(
        scenario = sc,
        metric = metric_col,
        comparator = comparator,
        mean_delta = ci["mean"],
        lwr = ci["lwr"],
        upr = ci["upr"],
        n = length(d),
        noninferior = if (is.finite(ci["lwr"])) ci["lwr"] > (-opts$delta) else NA,
        stringsAsFactors = FALSE
      )
      idx <- idx + 1L
    }
    do.call(rbind, out)
  }

  ni_rows <- list(
    compare_metric("similarity", "mcca"),
    compare_metric("similarity", "mfa"),
    compare_metric("r2_y", "mcca"),
    compare_metric("r2_y", "mfa")
  )
  ni_df <- do.call(rbind, ni_rows)

  # vs best comparator (mcca/mfa) per replicate
  best_rows <- list()
  bidx <- 1L
  for (sc in unique(raw$scenario)) {
    sub <- raw[raw$scenario == sc & raw$method %in% c("ipca", "mcca", "mfa"), c("rep", "method", "similarity", "r2_y")]
    reps_sc <- sort(unique(sub$rep))
    for (metric in c("similarity", "r2_y")) {
      deltas <- c()
      for (rp in reps_sc) {
        ss <- sub[sub$rep == rp, , drop = FALSE]
        ip <- ss[ss$method == "ipca", metric]
        cm <- ss[ss$method %in% c("mcca", "mfa"), metric]
        if (length(ip) == 1 && is.finite(ip) && any(is.finite(cm))) {
          deltas <- c(deltas, ip - max(cm[is.finite(cm)]))
        }
      }
      ci <- bootstrap_ci(deltas, B = if (opts$quick) 300 else 1200)
      best_rows[[bidx]] <- data.frame(
        scenario = sc,
        metric = metric,
        comparator = "best(mcca,mfa)",
        mean_delta = ci["mean"],
        lwr = ci["lwr"],
        upr = ci["upr"],
        n = length(deltas),
        noninferior = if (is.finite(ci["lwr"])) ci["lwr"] > (-opts$delta) else NA,
        stringsAsFactors = FALSE
      )
      bidx <- bidx + 1L
    }
  }
  ni_df <- rbind(ni_df, do.call(rbind, best_rows))

  ni_path <- "dev/ipca_vs_multiset_noninferiority.csv"
  utils::write.csv(ni_df, ni_path, row.names = FALSE)

  report_path <- "dev/ipca_vs_multiset_report.md"
  con <- file(report_path, open = "wt")
  on.exit(close(con), add = TRUE)
  writeLines("# iPCA vs Multiset Synthetic Battery", con)
  writeLines("", con)
  writeLines(sprintf("- Date: %s", as.character(Sys.time())), con)
  writeLines(sprintf("- Quick mode: %s", opts$quick), con)
  writeLines(sprintf("- Replicates per scenario: %d", reps), con)
  writeLines(sprintf("- Non-inferiority margin delta: %.4f", opts$delta), con)
  writeLines("", con)
  writeLines("## Scenario Summary", con)
  writeLines("", con)
  writeLines("```text", con)
  writeLines(capture.output(print(summary_df, row.names = FALSE)), con)
  writeLines("```", con)
  writeLines("", con)
  writeLines("## Non-inferiority Results", con)
  writeLines("", con)
  writeLines("```text", con)
  writeLines(capture.output(print(ni_df, row.names = FALSE)), con)
  writeLines("```", con)

  cat("\nCompleted synthetic battery.\n")
  cat(sprintf("Raw results: %s\n", raw_path))
  cat(sprintf("Summary: %s\n", summary_path))
  cat(sprintf("Non-inferiority: %s\n", ni_path))
  cat(sprintf("Report: %s\n", report_path))
}

main()

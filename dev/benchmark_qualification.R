# Shared, side-effect-free qualification helpers for development benchmarks.
bench_classify <- function(raw) {
  required <- if (is.null(raw$convergence_required)) rep(TRUE, nrow(raw)) else raw$convergence_required
  reason <- rep("eligible", nrow(raw))
  reason[is.na(required) | (required & is.na(raw$converged))] <- "unknown_convergence"
  reason[!is.na(required) & required & !is.na(raw$converged) & !raw$converged] <- "nonconverged"
  reason[is.na(raw$ok) | !raw$ok] <- "failed"
  reason
}

bench_summary <- function(raw, groups, metrics) {
  reason <- bench_classify(raw)
  keys <- unique(raw[groups])
  out <- lapply(seq_len(nrow(keys)), function(i) {
    take <- rep(TRUE, nrow(raw))
    for (key in groups) take <- take & raw[[key]] == keys[[key]][i]
    r <- reason[take]
    row <- keys[i, , drop = FALSE]
    row$attempted <- sum(take)
    for (status in c("failed", "nonconverged", "unknown_convergence", "eligible")) row[[status]] <- sum(r == status)
    for (metric in metrics) {
      x <- raw[[metric]][take]
      x <- x[r == "eligible" & is.finite(x)]
      row[[paste0(metric, "_n")]] <- length(x)
      stats <- if (length(x)) c(stats::median(x), stats::quantile(x, c(.25, .75), names = FALSE), min(x), max(x)) else rep(NA_real_, 5)
      for (j in seq_along(stats)) row[[paste0(metric, "_", c("median", "q1", "q3", "min", "max")[j])]] <- stats[j]
    }
    row
  })
  do.call(rbind, out)
}

bench_pairs <- function(raw, methods, metric, ids = c("scenario", "rep")) {
  stopifnot(length(methods) == 2L)
  raw$exclusion <- bench_classify(raw)
  select <- function(method) {
    out <- raw[raw$method == method, c(ids, metric, "exclusion"), drop = FALSE]
    if (anyDuplicated(out[ids])) stop("Duplicate benchmark replicate identities.")
    out
  }
  out <- merge(select(methods[1]), select(methods[2]), by = ids, all = TRUE, suffixes = c("_a", "_b"))
  for (side in c("a", "b")) {
    name <- paste0("exclusion_", side)
    out[[name]][is.na(out[[name]])] <- "missing_fit"
    out[[name]][out[[name]] == "eligible" & !is.finite(out[[paste0(metric, "_", side)]])] <- "nonfinite_metric"
  }
  out$eligible <- out$exclusion_a == "eligible" & out$exclusion_b == "eligible"
  out$delta <- ifelse(out$eligible, out[[paste0(metric, "_a")]] - out[[paste0(metric, "_b")]], NA_real_)
  out
}

bench_load_local <- function() {
  if (!file.exists("DESCRIPTION") || !requireNamespace("pkgload", quietly = TRUE)) stop("Run from the muscal checkout with pkgload installed.")
  pkgload::load_all(".", quiet = TRUE)
}

bench_provenance <- function(seed) {
  files <- sort(unique(c("DESCRIPTION", "NAMESPACE", list.files("R", full.names = TRUE),
    list.files("src", pattern = "\\.(cpp|h)$", full.names = TRUE),
    list.files("dev", pattern = "(benchmark.*R$|ipca.*R$)", full.names = TRUE))))
  list(seed = seed, utc = format(Sys.time(), tz = "UTC"),
    git_head = system2("git", c("rev-parse", "HEAD"), stdout = TRUE),
    git_status = system2("git", c("status", "--porcelain"), stdout = TRUE),
    sha256 = stats::setNames(vapply(files, digest::digest, character(1), algo = "sha256", file = TRUE), files),
    versions = vapply(c("muscal", "Matrix", "multivarious", "multidesign", "pkgload", "RSpectra"),
      function(p) if (requireNamespace(p, quietly = TRUE)) as.character(utils::packageVersion(p)) else NA_character_, character(1)),
    loaded_from = getNamespaceInfo(asNamespace("muscal"), "path"), session = utils::capture.output(utils::sessionInfo()))
}

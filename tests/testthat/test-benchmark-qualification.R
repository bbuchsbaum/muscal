.bench_helpers <- function() {
  e <- new.env(parent = baseenv())
  path <- testthat::test_path("..", "..", "dev", "benchmark_qualification.R")
  skip_if_not(file.exists(path), "Development benchmarks are excluded from source packages")
  sys.source(path, envir = e)
  e
}

test_that("benchmark summaries retain failure and convergence accounting", {
  h <- .bench_helpers()
  raw <- data.frame(method = "a", ok = c(TRUE, TRUE, FALSE, TRUE, TRUE),
    converged = c(TRUE, FALSE, NA, NA, TRUE), elapsed = c(2, 100, .1, 200, 4))
  out <- h$bench_summary(raw, "method", "elapsed")
  expect_equal(unlist(out[c("attempted", "failed", "nonconverged", "unknown_convergence", "eligible")], use.names = FALSE), c(5, 1, 1, 1, 2))
  expect_equal(out$elapsed_median, 3)
  expect_equal(out$elapsed_n, 2L)
  raw$convergence_required <- FALSE
  expect_equal(h$bench_summary(raw, "method", "elapsed")$eligible, 4L)
})

test_that("paired benchmarks use joint eligible identities with explicit exclusions", {
  h <- .bench_helpers()
  raw <- data.frame(scenario = "one", rep = rep(1:5, 2), method = rep(c("a", "b"), each = 5),
    ok = c(TRUE, TRUE, FALSE, TRUE, TRUE, rep(TRUE, 5)),
    converged = c(TRUE, FALSE, NA, NA, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE),
    metric = c(7, 8, 9, 10, 11, 2, 3, 4, 5, 6))
  out <- h$bench_pairs(raw, c("a", "b"), "metric")
  expect_equal(out$rep[out$eligible], 1L)
  expect_equal(out$delta[out$eligible], 5)
  expect_equal(out$exclusion_a, c("eligible", "nonconverged", "failed", "unknown_convergence", "eligible"))
  expect_equal(out$exclusion_b[5], "nonconverged")
  raw$metric[6] <- NA
  expect_false(any(h$bench_pairs(raw, c("a", "b"), "metric")$eligible))
  expect_error(h$bench_pairs(rbind(raw, raw[1, ]), c("a", "b"), "metric"), "Duplicate")
})


test_that("benchmark receipts bind the loaded checkout and source hashes", {
  h <- .bench_helpers()
  root <- normalizePath(testthat::test_path("..", ".."))
  old <- setwd(root)
  on.exit(setwd(old))
  p <- h$bench_provenance(117)
  expect_equal(p$seed, 117)
  expect_equal(normalizePath(p$loaded_from), root)
  expect_equal(p$sha256[["R/ipca.R"]], digest::digest(file = "R/ipca.R", algo = "sha256"))
  expect_true("multivarious" %in% names(p$versions))
})

test_that("Design-MFA profiling retains engine failures as raw rows", {
  .bench_helpers()
  expressions <- parse(testthat::test_path("..", "..", "dev", "benchmark_design_mfa.R"))
  env <- new.env(parent = globalenv())
  for (expr in expressions) {
    if (is.call(expr) && identical(expr[[1]], as.name("<-")) && is.symbol(expr[[2]]) &&
        as.character(expr[[2]]) %in% c("max_profiled_allocation", "profile_fit")) eval(expr, env)
  }
  out <- env$profile_fit("failed", stop("intentional fixture failure"), 4, 2)
  expect_null(out$fit)
  expect_false(out$row$ok)
  expect_match(out$row$error, "intentional fixture failure")
  expect_true(is.na(out$row$solver_converged))
})

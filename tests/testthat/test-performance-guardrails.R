library(testthat)
library(muscal)

skip_if_not_installed("genpca")
skip_if(
  tolower(Sys.getenv("MUSCAL_SKIP_PERF", unset = "false")) %in% c("1", "true", "yes"),
  "Performance guardrails disabled via MUSCAL_SKIP_PERF."
)

.perf_guardrail_budget <- function(default, strict = default * 0.6) {
  strict_flag <- tolower(Sys.getenv("MUSCAL_PERF_STRICT", unset = "false"))
  if (strict_flag %in% c("1", "true", "yes")) strict else default
}

.perf_guardrail_time <- function(fun, iterations, reps = 3L) {
  invisible(suppressMessages(fun()))
  times <- numeric(reps)

  for (r in seq_len(reps)) {
    gc()
    t0 <- proc.time()[["elapsed"]]
    for (i in seq_len(iterations)) {
      invisible(suppressMessages(fun()))
    }
    times[[r]] <- (proc.time()[["elapsed"]] - t0) / iterations
  }

  stats::median(times)
}

test_that("mfa representative fit stays within a soft runtime budget", {
  set.seed(601)
  X1 <- matrix(rnorm(300 * 30), 300, 30)
  X2 <- matrix(rnorm(300 * 25), 300, 25)

  elapsed <- .perf_guardrail_time(
    function() mfa(list(X1 = X1, X2 = X2), ncomp = 4),
    iterations = 5L
  )

  expect_lt(elapsed, .perf_guardrail_budget(0.20))
})

test_that("cv_muscal representative reconstruction workflow stays within a soft runtime budget", {
  set.seed(602)
  X1 <- matrix(rnorm(300 * 30), 300, 30)
  X2 <- matrix(rnorm(300 * 25), 300, 25)
  md <- multidesign::multidesign(
    cbind(X1, X2),
    data.frame(group = rep(1:3, each = 100))
  )
  folds <- multidesign::cv_rows(
    md,
    rows = list(1:40, 101:140, 201:240),
    preserve_row_ids = TRUE
  )

  elapsed <- .perf_guardrail_time(
    function() {
      cv_muscal(
        folds = folds,
        fit_fn = function(analysis) {
          Xa <- multidesign::xdata(analysis)
          mfa(
            list(
              X1 = Xa[, 1:ncol(X1), drop = FALSE],
              X2 = Xa[, (ncol(X1) + 1):ncol(Xa), drop = FALSE]
            ),
            ncomp = 4
          )
        },
        estimate_fn = function(model, assessment) {
          predict(model, multidesign::xdata(assessment), type = "reconstruction")
        },
        truth_fn = function(assessment) multidesign::xdata(assessment),
        metrics = "mse"
      )
    },
    iterations = 3L
  )

  expect_lt(elapsed, .perf_guardrail_budget(0.35))
})

test_that("infer_muscal representative bootstrap workflow stays within a soft runtime budget", {
  set.seed(603)
  X1 <- matrix(rnorm(300 * 30), 300, 30)
  X2 <- matrix(rnorm(300 * 25), 300, 25)
  Y <- matrix(rnorm(300 * 8), 300, 8)

  fit <- suppressMessages(
    anchored_mfa(
      Y = Y,
      X = list(X1 = X1, X2 = X2),
      row_index = list(X1 = 1:300, X2 = 1:300),
      ncomp = 3,
      max_iter = 20
    )
  )

  elapsed <- .perf_guardrail_time(
    function() infer_muscal(fit, method = "bootstrap", statistic = "sdev", nrep = 6, seed = 11),
    iterations = 3L
  )

  expect_lt(elapsed, .perf_guardrail_budget(1.20))
})

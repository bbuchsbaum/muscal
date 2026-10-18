evidence <- "dev/reviews/2026-10-18/evidence/repairs/F12-followup"
check_dir <- "/private/tmp/muscal-F12-followup-check-20261018"
stopifnot(!dir.exists(check_dir))
dir.create(check_dir)
result <- rcmdcheck::rcmdcheck(
  ".", args = "--no-manual", check_dir = check_dir, error_on = "never"
)
counts <- list(
  errors = length(result$errors), warnings = length(result$warnings),
  notes = length(result$notes), check_dir = check_dir,
  args = "--no-manual"
)
print(counts)
jsonlite::write_json(counts, file.path(evidence, "check-counts.json"),
                     auto_unbox = TRUE, pretty = TRUE)
saveRDS(result, file.path(check_dir, "rcmdcheck-result.rds"))
quit(status = as.integer(any(c(counts$errors, counts$warnings, counts$notes) > 0)))

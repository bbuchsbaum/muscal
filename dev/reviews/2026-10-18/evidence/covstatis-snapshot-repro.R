# Run from the repository root. This intentionally tests the reviewed Git
# snapshot, because another workspace actor repaired the live file mid-audit.
pkgload::load_all(".", compile = FALSE, quiet = TRUE)
snapshot <- "06eb85f1eee84eff16e1f0bcd23f54730ae18152"
src <- system2("git", c("show", paste0(snapshot, ":R/covstatis.R")), stdout = TRUE)
stopifnot(is.null(attr(src, "status")))
env <- new.env(parent = asNamespace("muscal"))
eval(parse(text = src), envir = env)
mk <- function(a) tcrossprod(cbind(a, a))
msg <- tryCatch({
  env$covstatis.list(list(mk(1:4), mk(2:5)), ncomp = 1)
  "<no error>"
}, error = conditionMessage)
stopifnot(grepl("non-conformable arguments", msg, fixed = TRUE))
cat("Snapshot:", snapshot, "\nConfirmed one-component error:", msg, "\n")

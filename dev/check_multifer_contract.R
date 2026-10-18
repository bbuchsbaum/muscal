#!/usr/bin/env Rscript
# Run after installing the pinned multifer revision in the workflow.
pkgload::load_all(".", quiet = TRUE)
if (!muscal:::.dmfa_multifer_effect_contract_available() ||
    !muscal:::.dmfa_multifer_stability_contract_available()) {
  stop("Unsupported multifer integration build: requires effect_significance, clustered list roles, bootstrap_fits, and statistic_stability_from_bootstrap with a statistic callback and probs/alpha. The supported source revision is pinned in .github/workflows/multifer-contract.yaml; public multifer 1.0.0 alone is not a sufficient version contract.", call. = FALSE)
}
source("tests/testthat/helper-design-mfa-sim.R")
for (p in c("test-design-mfa-multifer.R", "test-design-mfa-inference.R", "test-bada-multidesign.R")) {
  testthat::test_file(file.path("tests/testthat", p), reporter = "summary", stop_on_failure = TRUE)
}

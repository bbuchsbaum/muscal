# Review follow-up receipts

This follow-up addresses the zero-norm loading draw failure in F12 and the stale
COVSTATIS rationale. It preserves the original audit and repair-gate receipts.

`infer_muscal()` now validates nonzero loading directions before sampling the
observation and inside each replicate's statistic guard. Rejected draws retain
their original index and statistic-stage failure message. Valid draws remain in
their original loading units. All-invalid and single-valid bootstrap cases use
the existing valid-replicate guards. Numeric vector statistics may still be zero.

The COVSTATIS tag retains `research/experimental`; it acknowledges the repaired
scalar geometry and preserves the remaining subject-space scaling/rank caveats.

- `before.log`: focused inference suite with the new regressions, before the
  source fix; exit 1, reproducing the aggregation abort and sampling before an
  invalid observation is rejected.
- `focused.log` and metadata: inference, loading alignment and BaDA refit suites
  after the fix; 30 blocks, 205 expectations, no failures/errors/skips.
- `suite.log`, metadata and `suite-counts.json`: full local suite after the fix.
- `check.R`, `check.log`, metadata and `check-counts.json`: full package build,
  including vignette rebuilds, and `R CMD check --no-manual` after the fix.
- `tested-source-hashes.json`: SHA-256 hashes of the same 244 source paths as the
  original gate manifest. `source-delta.json` identifies the two changed source
  files and binds the original manifest. `verification.json` records outcomes
  and supplementary documentation hashes.

Commands run at the package root with `LC_ALL=en_US.UTF-8`,
`R_LIBS=/private/tmp/muscal-Q05-pinned-library`, and `R_MAKEVARS_USER` set to the
absolute path of `../GATE-apple-cxx17.Makevars`. Command arguments and exit status
are retained in the corresponding `.meta.json` files. The package check uses a
new task-owned directory at `/private/tmp/muscal-F12-followup-check-20261018`.
The full suite uses `testthat::test_local()`; focused suites use
`pkgload::load_all()` and `testthat::test_file()`.

These are local implementation checks. Q05's unpublished upstream multifer
revision still blocks its hosted CI prerequisite. No new scientific calibration,
publication or release qualification is claimed.

# Sequential repair status — 2026-10-18

All thirteen confirmed defects have locally verified repairs. Q01–Q04 and Q06 are
also locally complete. Q05's integration repair and exact pinned-engine smoke
pass; its CI prerequisite remains open because the supported multifer revision
has not been published upstream. Mote remains authoritative.

The original audit, plan, tags and adverse receipts remain a frozen snapshot.
[Current method tags](method-tags-current.csv) retain scientific qualifications;
passing these regressions does not promote experimental estimators to solid.
This catalog preserves the 79 original exported entries and adds three repaired
S3 methods explicitly. Concurrent DISTATIS additions remain outside the original
file-review scope.

| Finding | Verified local result | Main regression |
|---|---|---|
| F01 | Optional ICC requires explicit shared feature identities and valid nonzero directions. | test-component-retention-contract.R |
| F02 | MP support edge uses supplied noise variance, declared spectrum units, feature dimension and effective sample size. Unavailable evidence never passes. | test-rmt-edge.R |
| F03 | Feature counts derive from loading rows and cross-check stored indices. | test-component-retention-contract.R |
| F04 | Alpha affects enabled ICC only; a support edge is explicitly not an alpha-level test. | test-component-retention-contract.R |
| F05 | Whole-subject BaDA refits, fresh preprocessing, duplicate subject identities, permutation/sign alignment and failed-draw accounting replace the fixed-basis approximation. | test-bada-bootstrap-refit.R |
| F06 | Normalized Laplacians use adjacency degrees, with zero isolated rows. | test-clusterwise-graph.R |
| F07 | k-NN edges use undirected union, explicit self removal and unit-invariant tie handling. | test-clusterwise-graph.R |
| F08 | Stored fitted preprocessors support block projection and raw-unit reconstruction with schema checks. | test-clusterwise-projection.R |
| F09 | Common numerical rank is checked; unsupported rank errors or explicit global reduction replaces random padding. | test-clusterwise-rank.R |
| F10 | Empty public blocks and invalid k are rejected; zero/singleton/small graph cases are explicit. | test-clusterwise-graph.R |
| F11 | Concurrent COVSTATIS repair retained and verified against independent rank-one and rank-two geometry. | test-covstatis-scalar.R |
| F12 | Observed and replicate inference statistics have nonempty finite, compatible schemas and nonzero loading directions; invalid draws are counted failures before aggregation. | test-infer-muscal.R |
| F13 | Shared observation counts are checked before and after preprocessing. | test-penalized-mfa.R |

| Qualification | Result |
|---|---|
| Q01 | Benchmark summaries retain attempted/failed/nonconverged/unknown/eligible counts; paired comparisons use joint eligible identities. Raw failures, source hashes, dependency versions and seeds are retained. Training adjusted R-squared is labeled honestly. Synthetic aggregation tests pass; expensive benchmarks were not rerun. |
| Q02 | Exact one-to-one assignment and sign alignment; tied retained roots have unavailable coordinate intervals and subspace diagnostics. Independent assignment oracles pass. |
| Q03 | Frozen-reference conditional target, callback schema and held-out effect-recovery meaning are explicit. Existing subject and fold-isolation tests pass. No exchangeability/coverage/null-rejection calibration claim. |
| Q04 | Rebuilt native kernels validate dimensions, finite values, penalties and integer groups before indexing. Base-R numerical oracles and affected model suites pass. Generated Rcpp bindings were regenerated. |
| Q05 | Real scalar/multi-component permutation and stability routes pass, including failure accounting, on exact multifer `8fe24a8a211a4666ef93ccc9c079af477e235336` installed into an isolated library. Public main `b42afa827ebe4f38fa5236423e758026c6a6c551` fails the capability guard actionably. A manual CI lane pins the supported SHA, but it cannot fetch that unpublished revision. This issue stays open. |
| Q06 | Unique undiscovered specifications are consolidated into the existing discoverable suite; normalization, RNG and solver tolerance match the claimed equivalences. Duplicate helpers/cases removed; original preserved with a coverage map. |

Verification receipts are in [evidence/repairs](evidence/repairs/). Each completed
Mote child records focused checks and a source-hash manifest. The final integrated
source manifest for the original repair gate is `GATE-tested-source-hashes.json`.

The review follow-up repairs zero-norm loading validation inside the existing
observed/per-draw checks: invalid observations fail before sampling, degenerate
draws are excluded and counted, and valid draws retain their loading units.
The current COVSTATIS rationale now acknowledges the verified F11 scalar repair
while retaining its research/experimental classification and remaining
subject-space scaling/rank qualifications. The original gate receipts below
remain unchanged; supplementary receipts and the current 244-file source
manifest are in [F12-followup](evidence/repairs/F12-followup/).

- Follow-up full local suite: **59 files, 559 test blocks, 2,754 passing
  expectations, zero failures/errors, two pre-existing plotting skips**.
  Command exit 0; three dependency build-version warnings remain. The new
  regressions failed before the fix, and all 205 expectations across the focused
  inference, loading-alignment and BaDA-refit suites pass afterward.
- Follow-up `R CMD check --no-manual` with the same isolated Apple Clang/C++17
  configuration: **zero errors, zero warnings, zero notes**, exit 0. Examples,
  installed-package tests and all vignette rebuilds passed. The wrapper again
  printed the sandbox-limited Quarto version-query diagnostic described below;
  the package check itself reports `Status: OK`.

- Original gate full local suite: **59 files, 557 test blocks, 2,734 passing expectations,
  zero failures/errors, two pre-existing plotting skips**. Command exit 0.
  Three dependency build-version warnings remain. The first receipt exporter
  failed after tests completed because of a list column; the corrected exporter
  and suite rerun exit 0. Both logs are retained.
- Original gate first full R CMD check: **zero errors, one warning, one note**. All examples,
  tests and vignettes passed. The compiler-header warning was reproduced by
  installing untouched baseline `06eb85f1eee84eff16e1f0bcd23f54730ae18152`.
  The owned S3 documentation note was corrected by processing the actual package
  documentation blocks while writing only the affected page.
- Original gate final full R CMD check with isolated Apple Clang/C++17 configuration:
  **zero errors, zero warnings, zero notes**, exit 0. Examples, tests and all
  vignette rebuilds passed.
  The initial Apple compiler probe failed under its implicit C++98 default; the
  explicit C++17 probe passed. Both receipts are retained. User compiler settings
  were not changed, and warning checks were not suppressed. The rcmdcheck wrapper
  separately printed a sandbox-limited Quarto version query diagnostic; it did
  not affect package checks or the successful vignette rebuilds.

No hosted CI, push, merge, deployment or release is claimed. Positive-lambda
clusterwise stationarity/descent, general bootstrap coverage and scientific null
calibration remain research qualifications. The broader clusterwise recovery plan
still governs those separate milestones.

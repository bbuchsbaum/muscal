# File and method audit — 2026-10-18

Subsequent implementation: [sequential repair status](repair-status.md) and [current method tags](method-tags-current.csv). The audit and its original evidence below remain a frozen snapshot.

Two sub-agents individually reviewed all 51 computational source files: 41 R package files, two C++ files, seven development R scripts, and one build script. Relevant tests, manual pages, and vignettes were checked as supporting evidence. Generated documentation, binary build products, vendored skill material, and static assets were not individually audited. No implementation changes were made by this review.

The review started with a dirty tree at `2c5c1e83692814dc20502e6963573e57263c07c6`. Another workspace actor committed existing work while the audit ran; the source snapshot manifest records `06eb85f1eee84eff16e1f0bcd23f54730ae18152`. Findings refer to the file hashes in `source-manifest.json`, not an assumption that the working tree remains unchanged. The two reviewer appendices record their individual evidence and coverage.

**Concurrent edit qualification:** after that snapshot, another actor modified `R/covstatis.R`, `R/covstatis_multifer.R`, `R/plot_mfa.R`, and `NAMESPACE`, then added `R/distatis.R`, `R/plot_distatis.R`, and associated manual pages. Those additions and modified implementations are outside this audit's frozen scope. The current one-component COVSTATIS repro now succeeds after changes to the diagonal construction. F11 remains a confirmed finding about the recorded snapshot, with an apparent concurrent repair; this review does not certify that later implementation. New DISTATIS plotting dispatch was observed mid-edit and is not treated as a completed integration or a stable new defect. The full-suite receipt predates those edits. The consolidated tags below and in the CSV take precedence over earlier reviewer classifications.

## Confirmed findings

| ID | Severity | Location | Trigger and consequence |
|---|---|---|---|
| F01 | High | `R/utils.R:347-354` | `significant_components()` multiplies loadings from different feature spaces. Unequal lengths are recycled by R; an example with 3 and 6 features returns ICC=1 and p=0.0169 for an undefined cross-space cosine. It can retain components on invalid coherence evidence. Equal dimensions alone also do not establish feature correspondence. |
| F02 | High | `R/utils.R:297-312` | The purported Marchenko–Pastur edge divides a selected median eigenvalue by `(1+sqrt(gamma))^2` to estimate noise and immediately multiplies by the same factor. In every non-skipped branch the edge is just that median. This is a median-spectrum rule, not the advertised dimension-dependent MP edge. |
| F03 | Medium | `R/components_utils.R:50-58` | `estimate_components()` uses `ncol(V)` (retained components) for block feature counts instead of `nrow(V)`. With two 20-feature/3-component blocks and n=10, it evaluates gamma=0.6 instead of 4 and takes the wrong RMT applicability branch. F02 separately explains why gamma cancels inside the non-skipped branch. |
| F04 | Medium | `R/components_utils.R:37-59` | The documented `alpha` argument of `estimate_components()` is never read or forwarded, so changing the requested threshold has no effect. |
| F05 | High for inference | `R/bada.R:104-133` | `bootstrap.bada()` resamples barycenters but retains the original basis, scores, and variances. Its `boot_lds_*` summaries describe projections through a fixed fitted solution, rather than the documented sampling uncertainty of refitted component loadings. Relabel as conditional projected summaries or refit and align each replicate. This finding is established by control-flow inspection. |
| F06 | High | `R/penalized_mfa_clusterwise.R:748-751` | The normalized-Laplacian branch uses degrees from `rowSums(L+I)`. Since `L=D-A` has zero row sums, those values are all one: the normalization is a no-op. The default therefore applies unnormalized degree-dependent smoothing while claiming normalized smoothing. |
| F07 | High | `R/penalized_mfa_clusterwise.R:57-66` | `forceSymmetric(A, uplo="U")` discards lower-triangle-only directed neighbor edges. Feature reordering can change the graph and penalty after restoring feature identities. This is neither union nor mutual k-NN symmetrization. |
| F08 | High | `R/penalized_mfa_clusterwise.R:727-736,1079-1087` | The optimizer fits block centering/preprocessing but returns a projector with a new `pass()` preprocessor. Generic projections then use raw coordinates rather than the fitted coordinates. |
| F09 | Medium | `R/penalized_mfa_clusterwise.R:163-174,714-725,852-875` | With 2- and 5-feature blocks, `ncomp=3`, and positive lambda, local-rank reconstruction gradients are added to padded common-rank spatial gradients. The advertised variable-rank path crashes. |
| F10 | Medium | `R/penalized_mfa_clusterwise.R:57-61` | Fewer than seven pooled features with the default six neighbors creates NA sparse-matrix indices and errors. The neighbor count is not bounded by available features. |
| F11 | Medium | `R/covstatis.R:180-188` | A one-component fit passes a length-one numeric value to `diag()`, which R interprets as a matrix dimension. Valid one-dimensional fits can fail with `non-conformable arguments`; explicit diagonal dimensions or column scaling are needed. |
| F12 | Medium | `R/infer_muscal.R:51-53,97-131,166-171` | Missing observed `sdev` becomes `numeric(0)` and is accepted. A bootstrap can report every replicate successful while returning a zero-row summary. Validate a nonempty finite observed statistic and matching replicate shapes. |
| F13 | Low | `R/penalized_mfa.R:326-350` | Blocks with unequal observation counts are accepted until their scores are added, then error internally. This is a missing public-boundary check for shared observations, rather than evidence of silent incorrect results. |

Executable reproductions and raw receipts are in `evidence/`. The reviewer appendices provide further details and repair directions. The findings are flags, not fixes.

There are **13 snapshot findings** above: eight from the core-method review and five from the design/utilities review. Twelve remain unaddressed by the observed concurrent edits; F11 has the apparent concurrent repair noted above. Both reviewers' executable reproductions passed, and the parent independently reran their core/design checks. The explicit snapshot COVSTATIS reproducer also passed by detecting the historical error in the recorded Git version.

## Method tags

The complete catalog of the snapshot's 79 exported functions is `method-tags.csv`; the appendices also name S3 method groups and scientific helpers. Every snapshot export has a binary tag: **solid** or **research/experimental**, with scope and reason. A solid wrapper does not strengthen the scientific status of the model it wraps.

| Method family | Tag | Qualification |
|---|---|---|
| Complete-data MFA | **solid** | Descriptive factor fit and stored-preprocessing projection; iterative missing-data imputation remains research/experimental. |
| MCCA | **solid** | Descriptive regularized fit; this does not qualify component-selection or inference helpers. |
| BaDA point estimation | **solid** | The fixed-basis bootstrap uncertainty above is research/experimental. |
| COVSTATIS and its subject diagnostics | **research/experimental** | Classical descriptive construction, but the one-component implementation defect and subject-space rank/scaling qualifications prevent an unrestricted solid tag. It is not a qualified DISTATIS implementation. |
| iPCA and alpha tuning | **research/experimental** | Adaptive estimation and tuning need convergence/quality qualification for the intended data regime. |
| Penalized MFA, clusterwise PMFA | **research/experimental** | Optimization and graph/preprocessing defects above; avoid scientific claims from the affected paths. |
| Anchored/linked/aligned MFA, aligned/anchored MCCA | **research/experimental** | Extensions with study-specific alignment assumptions and iterative estimation. |
| Graph-anchored and coupled graph-anchored MFA | **research/experimental** | Graph, anchor, and optimization assumptions need application-specific validation. |
| Aligned RRR, aligned interbattery, response-aligned MFA | **research/experimental** | Supervised/common-space estimation and alignment assumptions remain conditional. |
| BaMFA, bilinear mixed models and tuning | **research/experimental** | Alternating latent models and heuristic/model-selection procedures. |
| Design-MFA, soft/spatial variants, tuning | **research/experimental** | Good fixed-design numerical oracle tests; the overall method's alignment, target-effect, tuning, and inferential claims remain bounded by the tested regimes. |
| Component retention, loading reliability, all infer_* and multifer adapters | **research/experimental** | Confirmed retention defects, incomplete component matching, exchangeability assumptions, and integration/calibration limits. |
| Coordinate maps, metric evaluation, orientation, plotting, dispatch and fit contracts | **solid** | Only their specified numerical/software behavior; all model-specific scientific limits still apply. COVSTATIS-derived diagnostics are qualified separately in the catalog. |
| Synthetic generators and research/benchmark scripts | **research/experimental** | Evidence for the specified synthetic regime, not a general scientific validation. The pkgdown build helper is solid software plumbing. |

## Verification and rejected flags

- The standard local suite completed: **487 test blocks across 48 files, 2,125 passing expectations, zero failures/errors, three skips**. Its log prints five warnings: three dependency build-version warnings and two insufficient-tail RMT warnings. The skipped multifer stability test means that particular engine contract was not exercised.
- Explicitly running `tests/testthat/aligned_models_synthetic_tests.R`, which is outside the standard test filename pattern, produced 18 passing expectations and two failures. This file labels some cases as specification tests. Their failures are reported separately: the MCCA comparison mixes default MFA normalization with explicit `None`, and the aligned-MFA comparison uses fresh random initializations with a threshold tighter than the default stopping tolerance. They do not, by themselves, establish implementation defects.
- The proposed `score_constraint="none"` defect was rejected: QR normalization is explicitly documented and its loading transformation preserves the fitted product.
- A fresh rank-one Design-MFA identity-permutation smoke test returned p=1 for both partial sums of squares and singular values. A historical error in another session's temporary log was not reproduced and is not a finding.
- The DISTATIS research script's independent geometry/rounded-target assertions passed, and it captured F11. Its distinction between raw interstructure products and normalized RV similarity is a methodological qualification, not a blanket assertion that every COVSTATIS normalization option is wrong.
- No full R CMD check, hosted CI, external reference-package parity study, or new calibration study was run. The test suite's green result does not negate the uncovered boundary cases.

Benchmark qualifications: repeated Design-MFA summaries omit convergence counts even though raw rows retain convergence; iPCA comparisons can include unconverged fits; the multiset script's adjusted R-squared is measured on training rows, and an installed package can be used without recording its source revision. Keep those limits attached to any performance/noninferiority claim. Details and locations are in the design reviewer appendix.

## Artifacts

- [Remediation and regression plan](remediation-plan.md): implementation order, contracts, regression fixtures, and completion criteria for F01–F13 and the additional audit qualifications.
- `core-review.md`: 21 package files, findings, method groups, supporting evidence, source hashes.
- `design-utilities-review.md`: 22 package/native files plus eight scripts, findings, method groups, evidence, source hashes.
- `method-tags.csv`: one entry for every NAMESPACE export, with tag, scope, reason, and defining file when local.
- `source-manifest.json`: all 51 computational source file hashes and the snapshot HEAD.
- `file-coverage.csv`: a row for each of those 51 files, its reviewer, source hash, and observed later drift.
- `evidence/`: raw test/reproducer logs and executable reproductions. Run from the repository root with Rscript and installed package dependencies.

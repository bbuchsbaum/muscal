# Remediation and regression plan

## Outcome and scope

Resolve findings F01–F13 in [the audit](README.md), with an independent regression for each finding and explicit limits on the resulting scientific claims. This document is an implementation plan; the repairs and new regression tests have not been implemented by this planning task.

The planning baseline is HEAD `06eb85f1eee84eff16e1f0bcd23f54730ae18152` plus the current dirty tree. DISTATIS/COVSTATIS work is in progress in the same workspace. Preserve it and coordinate ownership before editing overlapping files. The current COVSTATIS code uses `sweep()` for one-component scaling, and `tests/testthat/test-distatis.R` already contains a COVSTATIS geometry regression. F11 therefore needs verification and extension, not another independent implementation.

Planning verification: the current one-component, no-normalization `B=x x'` fixture described under F11 passed its score-shape and reconstruction assertions. Its receipt is `evidence/remediation-plan-f11-smoke.log`. The existing test was located by inspection; the full test suite and the proposed extended regression matrix were not run in this planning turn. The local `rcmdcheck` package and command signature below were also verified as available.

Keep the audit's source hashes and adverse receipts unchanged. Record each fix against the actual tested revision. A passing code regression closes its implementation defect; it does not automatically promote an entire statistical method to **solid**.

Use [the existing clusterwise recovery plan](../../penalized_mfa_clusterwise_recovery_plan.md) as the mathematical contract for F06–F10. This plan supplies concrete regression work for those findings; it does not authorize a second, incompatible graph or optimizer design. In particular, the existing recovery plan requires a common component rank and identifies further optimizer defects beyond these thirteen findings.

## Working decisions

1. A missing or invalid statistical test is **unavailable**, never an automatic pass. Do not silently fall back from requested RMT inference to retaining every positive component.
2. Cross-block loading coherence requires an explicit common feature domain. Equal feature counts do not prove correspondence. Make ICC optional and explicit rather than an unavoidable part of every RMT call.
3. A Marchenko–Pastur edge is not an alpha-level finite-sample hypothesis test. Do not make `alpha` change an arbitrary edge multiplier to manufacture an effect.
4. Clusterwise PMFA uses one rank `K` across all blocks, with `K <= min(p_s)` and a checked effective data rank. Default to a clear error for an unsupported rank; an explicit `rank_action="reduce"` may reduce every block together and must report the change.
5. The default graph symmetrization will be **union** of directed k-NN edges. Other policies must be explicit and separately tested. The normalized Laplacian uses adjacency degrees and zero rows for isolated vertices.
6. BaDA's loading bootstrap will refit the original estimation recipe. Keep any fixed-basis projected uncertainty as an explicitly named conditional target, not as a substitute for refitted loading uncertainty.
7. Preserve F11's concurrent repair. No implementation work in the new DISTATIS API is part of these thirteen repairs.

## Delivery sequence and ownership

Use one writer for each source scope. The component-retention changes share functions and should be delivered together rather than by competing writers. Each row below is a reviewable change set; the tests described later belong with the fix, not in a permanently red default branch.

| Change set | Findings | Source ownership | Prerequisite and completion evidence |
|---|---|---|---|
| 0. Freeze inputs and fixtures | All | This review folder; small synthetic test helpers | Refresh Git/dirty-state receipts; retain old reproductions; choose named test files and independent expected results. |
| 1. Validate inference statistics | F12 | `R/infer_muscal.R`, `test-infer-muscal.R` | Independent of other repairs; invalid observations fail before sampling; invalid replicates are accounted for. |
| 2. Validate shared observations | F13 | `R/penalized_mfa.R`, `test-penalized-mfa.R` | Independent; unequal row counts rejected before preprocessing/optimization. |
| 3. Close existing scalar-axis repair | F11 | Coordinate with current COVSTATIS owner; `test-covstatis.R` / existing test in `test-distatis.R` | Current fix retained; scalar and multi-axis geometry tests pass. |
| 4. Repair component retention | F01–F04 | `R/utils.R`, `R/components_utils.R`, relevant Rd files and retention tests | First lock the representation/noise contract; then implement ICC guards, dimensions, threshold, and alpha semantics together. |
| 5. Repair graph construction | F10 → F07 → F06 | Graph internals in `R/penalized_mfa_clusterwise.R`, graph tests | Small graphs first; then edge membership/permutation; then analytic Laplacians. |
| 6. Repair clusterwise rank and projection | F09, F08 | Same clusterwise writer; public wrapper and projection tests | Graph changes available for integration; common-rank validation before optimizer; fitted preprocessors retained. |
| 7. Repair BaDA bootstrap | F05 | `R/bada.R`, relevant helpers in `R/bada_multifer.R`, bootstrap tests | Reuse verified class/subject alignment and refit machinery where compatible; preserve the multifer target's declared semantics. |
| 8. Integrate, document, and qualify | All; Q01–Q06 below | Tests/docs/evidence; no competing source writers | Focused suites, full suite, package check, and explicit remaining scientific limits. |

Change sets 1–3 can proceed independently. Do not parallelize changes within sets 4–6 without disjoint files or isolated worktrees and an explicit integration owner. This is a proposed execution structure, not a record that agents have been launched.

## F01 — invalid cross-feature-space coherence

**Paths:** `R/utils.R::significant_components()`, `R/components_utils.R::estimate_components()`. Tests: new `tests/testthat/test-component-retention-contract.R`; update existing retention tests that currently assume implicit ICC.

**Implementation steps**

1. Separate RMT applicability and ICC applicability. Add an explicit ICC opt-in, with a declared common feature domain; do not compute ICC by default for arbitrary multiblock data.
2. Validate finite numeric loading matrices, equal component counts, and nonzero loading norms. For ICC, require common feature identities and an explicit correspondence assertion. Match named rows before computing any dot products. Reject duplicate/missing identities and incompatible spaces.
3. Permit RMT-only processing of heterogeneous blocks when its separate noise/representation contract is satisfied. Do not reject valid heterogeneous input merely because ICC cannot be defined.
4. Return explicit per-criterion status/reason. `keep` must describe only the requested, successfully evaluated criteria; disabled ICC cannot be presented as positive evidence. Keep the current approximate ICC null labeled experimental even after the dimensional bug is fixed.

**Regression tests**

- The audit's 3-by-1 / 6-by-1 loading fixture cannot return a finite ICC p-value or an ICC-based retained component. An explicit ICC request errors with a common-space message; an explicitly RMT-only request does not enter the ICC loop.
- Equal-size matrices with different declared feature IDs are rejected for ICC. This catches the more subtle version of the bug that a length check misses.
- Permuting one block's named feature rows, then aligning by identities, leaves the result unchanged.
- An explicitly matched two-block fixture has independently calculated squared cosine 1 for equal directions and 0 for orthogonal directions. Compare those numerical quantities directly, not two calls to the same implementation.
- Zero-norm, nonfinite, duplicated-ID, missing-ID, and inconsistent-component fixtures produce stable errors/statuses rather than recycled arithmetic.

**Done when:** incompatible spaces cannot produce coherence evidence, and RMT-only selection is demonstrably independent of ICC availability. This does not close null-calibration concerns about the approximate ICC p-value.

## F02 — cancelled MP-edge formula

**Paths:** `R/utils.R`, component-retention documentation. Tests: new `tests/testthat/test-rmt-edge.R` and the shared retention contract tests.

**Implementation steps**

1. First specify exactly what spectrum is being thresholded: covariance eigenvalues versus squared singular values, sample normalization, effective degrees of freedom after centering, total feature count, block weights, and the noise model. Do not infer an IID isotropic covariance spectrum from an arbitrary weighted MFA compromise.
2. Add a small pure numerical edge calculation for a **declared** isotropic covariance model with supplied positive finite noise variance. For covariance eigenvalues normalized by `n_eff`, use `gamma=p/n_eff` and edge `sigma2*(1+sqrt(gamma))^2`. Convert raw singular-value units explicitly if that input form is supported.
3. Remove the estimator that divides a selected median by the upper-edge factor. In the first repair, require independently supplied/justified noise variance. A later automatic noise estimator needs a specified full-bulk model and its own independent tests; do not estimate noise from only a handful of retained signal axes.
4. For unsupported fit/noise representations, return an unavailable diagnostic or a clear error for explicitly requested RMT. Never mark every component as passing. Preserve the old median rule only if deliberately exposed under a truthful heuristic name; do not silently retain it as RMT.
5. Keep finite-sample significance separate. An asymptotic support edge alone does not implement a test of size `alpha`.

**Regression tests**

- Analytic constants: supplied `sigma2=2` gives edge 4.5 at gamma 0.25 and 8 at gamma 1. These expectations are handwritten, independent of the production helper.
- With the same supplied noise variance, different positive gamma values produce the corresponding different edges. Use a fixed spectrum to expose the old cancellation.
- Put observed roots strictly below, equal to, and above the analytic edge; assert the documented comparison rule, including the equality boundary.
- Scaling data by `c` and noise variance by `c^2` scales covariance roots and the edge by `c^2`, preserving retention.
- Raw singular-value and covariance-eigenvalue representations, if both supported, agree after the stated `n_eff` conversion.
- A truncated weighted-MFA spectrum without a declared compatible noise model is explicitly unavailable, not an automatic success. Invalid variance, nonfinite roots, and invalid degrees of freedom are rejected.

**Done when:** the supplied-noise edge has independent numeric evidence, representation assumptions are enforced, and unsupported automatic RMT claims have been removed. Full automatic estimation remains experimental until independently specified and calibrated.

## F03 — component count substituted for feature count

**Paths:** `R/components_utils.R::estimate_components()`, relevant loading extraction and diagnostics. Tests: `test-component-retention-contract.R`.

**Implementation steps**

1. Resolve each `p_s` from the row count of its fitted loading matrix, consistent with post-preprocessing feature dimensions. Cross-check block-index metadata when available; disagreement is an error, not a reason to guess.
2. Derive sample/effective degrees of freedom from the representation contract in F02. A user-supplied valid `n` overrides only the documented sample-count field, not feature counts.
3. Return the resolved counts and criterion applicability in diagnostics so dimensional decisions can be inspected.

**Regression tests**

- Two 20-by-3 loading matrices with n=10 report feature counts `(20,20)`, total p=40, and the documented gamma of 4 for an uncentered covariance fixture. They must never report `(3,3)` or gamma 0.6.
- Heterogeneous 3-, 6-, and 9-feature blocks retain their distinct feature counts; this test explicitly disables ICC.
- Changing retained component columns while holding feature rows and sample metadata fixed does not change p or gamma.
- Mismatched block indices, zero-row loadings, and inconsistent metadata fail explicitly. A preprocessing fixture that drops a feature uses the post-fit dimension.

**Done when:** feature dimensionality is independently correct, regardless of whether RMT is applicable. Do not encode the current `gamma>1 -> pass everything` behavior as the new expected answer.

## F04 — ignored alpha

**Paths:** `R/components_utils.R`, `R/utils.R`, `man/estimate_components.Rd`, `man/significant_components.Rd`. Tests: `test-component-retention-contract.R`.

**Implementation steps**

1. Document `alpha` as a tail/significance probability, not a confidence level. Validate a finite scalar strictly between 0 and 1.
2. For explicitly requested, common-space ICC, forward `alpha` to the corrected decision procedure and report it with the decision. If ICC is combined with RMT, return per-criterion decisions as well as the actual combined criterion; the existing wrapper returns an RMT-only criterion alongside an ICC-filtered `keep`.
3. For an RMT-edge-only or variance-only call, an explicitly supplied `alpha` must warn/error as inapplicable rather than pretend to control significance. Do not invent a finite-sample RMT p-value. Use `missing(alpha)` to avoid noisy warnings when a caller did not supply the legacy default.

**Regression tests**

- A deterministic, matched five-feature ICC fixture with one component and identical directions has an approximate p-value between 0.001 and 0.01 under the currently documented approximation. Alpha 0.001 versus 0.01 must change the ICC decision and reported `keep`; record the approximation's experimental status.
- Wrapper and direct criterion calls agree for both alpha values, including the reported decision and resolved alpha.
- Invalid alpha values fail. Explicit alpha on an edge-only call produces the documented inapplicability message, rather than changing the edge.

**Done when:** alpha either governs the declared inferential decision or is explicitly rejected as irrelevant. Threading an unused argument through another helper is not sufficient.

## F05 — BaDA bootstrap uses a fixed fitted solution

**Paths:** `R/bada.R::bootstrap.bada()`, `summarize_boot()`, BaDA refit/alignment helpers; relevant Rd. Tests: new `tests/testthat/test-bada-bootstrap-refit.R` plus existing `test-bada-multidesign.R` and `test-infer-muscal.R` coverage.

**Implementation steps**

1. Preserve the subject as the resampling unit; sample whole subjects with replacement. Duplicate selections receive distinct bootstrap identities, so regrouping does not merge them and change their multiplicities.
2. Preserve the fitted class order and feature identities. Validate missing classes explicitly; do not substitute unlabeled row positions or silently fill missing class means.
3. Refit the original point-estimation recipe in each draw, including preprocessing and any residual-component options. Reuse compatible validated BaDA barycenter/refit helpers, but do not equate a conditional barycenter refit with a full raw-data refit without documenting what is frozen.
4. Extract refitted loadings/scores, then match component order and signs against the reference for separated roots. For tied/near-tied subspaces, report subspace stability or explicitly unavailable coordinatewise intervals; do not claim unique eigenvector intervals after arbitrary rotation.
5. Retain the public `boot_*` result shape where mathematically valid. If a fixed-basis projected mode remains useful, expose it with a distinct, explicit target and document conditional uncertainty. Preserve the separately declared `infer_bada()`/multifer projected-score target rather than silently changing that API's estimand.
6. Record refit failures, valid replicate count, and alignment/rank exclusions. Never summarize a failed or unidentifiable draw as successful loading inference.

**Regression tests**

- Build a small fixture with subject-specific loading directions and separated roots. For predetermined resample indices, independently average the selected subject class means, fit with base `svd()`, and align signs/order. Production replicate loadings and scores must agree with this independent oracle.
- The fixture must change the estimated basis across at least two draws; demonstrate that the old fixed-basis formula disagrees with the oracle. A check that bootstrap variance is merely nonzero would not catch the old bug.
- Verify duplicate subject multiplicity exactly, class-row permutation invariance after matching labels, and feature-name alignment.
- Explicit sign flips/component swaps of an otherwise identical fit leave aligned summaries unchanged. Tied-root fixtures exercise the documented subspace/unavailable policy.
- Cover ncomp=1, residual components, one failed draw, all failed draws, and seeded serial reproducibility. For residual options not yet refittable, require an explicit unsupported error rather than dropping those options.

**Done when:** each claimed loading replicate is a correctly weighted refit of the declared estimator and matches the independent fixture. Interval coverage is a separate seeded simulation study, not proven by these regressions.

## F06 — normalized Laplacian is a no-op

**Paths:** clusterwise graph construction/normalization. Tests: new `tests/testthat/test-clusterwise-graph.R` plus public clusterwise tests. Depends on the adjacency contract in F07/F10.

**Implementation steps**

1. Retain sparse adjacency `A` and compute `degree=rowSums(A)` before forming `L=D-A`.
2. Construct `D^(-1/2) L D^(-1/2)` using zero inverse degree for isolated vertices. Store the chosen graph type and degree vector; do not derive degrees from a Laplacian.
3. Route `normalized_laplacian=FALSE` to the combinatorial Laplacian and `TRUE` to the normalized one. Correct the default's documentation; expect numerical changes to spatial fits.

**Regression tests**

- Handwritten three-vertex path: combinatorial diagonal `(1,2,1)`; normalized diagonal `(1,1,1)` and both neighboring off-diagonals `-1/sqrt(2)`. Assert both complete matrices and that they differ.
- The path's normalized spectrum is `(0,1,2)` up to ordering; check it independently with base `eigen()` on the handwritten matrix.
- Add an isolated vertex: finite zero row/column; symmetry and nonnegative quadratic forms remain valid. Assert `L_norm %*% sqrt(degree)=0`, not the generally false zero-row-sum property for normalized L.
- A public fit returns the same requested Laplacian as the independent graph fixture. Scaling all edge weights by a positive constant preserves normalized L and scales combinatorial L.

**Done when:** graph normalization changes the intended metric, including at the public API, without NaNs for isolated vertices.

## F07 — triangle-dependent neighbor-edge loss

**Paths:** `spatial_constraints()`, `cluster_graph()` and graph internals. Tests: `test-clusterwise-graph.R`.

**Implementation steps**

1. Build directed edges with explicit vertex identities; remove self-edges by index, not by assuming self is the first distance tie.
2. Form union adjacency for the default, rather than copying one triangle. Keep zero diagonal and nonnegative symmetric weights.
3. Preserve feature/block identity through sparse construction. For ties at the kth-neighbor boundary, use stable feature IDs when available; without identity, include all exact boundary ties so row order cannot choose different scientific neighbors. Document that this can exceed k. Do not fabricate invariant identities from row positions.

**Regression tests**

- A directed lower-only edge must survive union symmetrization; a handwritten edge list supplies expected membership. If mutual/average options are exposed, test their distinct expected edge sets/weights separately.
- Use four non-equidistant points on a line with an independently enumerated directed one-neighbor list; verify the complete undirected edge set.
- Permute features within blocks, then blocks themselves, restore the original identities, and assert identical adjacency/Laplacian and penalty `sum(V*(L%*%V))` for a fixed loading matrix.
- Duplicate coordinates and kth-boundary distance ties preserve the documented policy and cannot create self-edges.

**Done when:** the graph and its penalty depend on declared coordinates/identities and policy, not storage order.

## F08 — fitted preprocessing is discarded

**Paths:** clusterwise wrapper/returned projector; shared preprocessing helpers. Tests: new `tests/testthat/test-clusterwise-projection.R` and existing preprocessing/fit-contract suites.

**Implementation steps**

1. Retain fitted block preprocessors from `prepare_block_preprocessors()` and construct the concatenated preprocessor using the established `multivarious::concat_pre_processors()` pattern in `penalized_mfa.list()`.
2. Preserve block names and feature indices. Define supported named-block and complete-list projection paths explicitly; do not rely on an inherited projector route that bypasses block preprocessing.
3. Transform raw new data exactly once. Preserve this method's existing `preproc=NULL` meaning (default centering), which differs from some other wrappers; changing that convention is not part of this fix.
4. Define score versus reconstruction output, including inverse preprocessing for raw-scale reconstruction if exposed. Serialize the fitted transformations rather than a specification that would refit on new data.

**Regression tests**

- Use strongly nonzero training means and a distinct held-out mean. For each block compare projection to the independently calculated `(X_new-training_mean) %*% V_s`.
- For standardization, use stored training center/scale in the oracle; ensure assessment means/scales are not fitted.
- Test default centering, explicit `pass()`, standardization, and a mixed per-block preprocessor list. Verify `NULL` follows the documented centering convention.
- Named block reordering preserves projected values after matching; unknown/missing blocks and feature-schema mismatch fail clearly.
- Saving/reloading the fit preserves projections. A spy/counter or equivalent independent transform check demonstrates exactly-once preprocessing.

**Done when:** raw-data projection matches explicit fitted-transform algebra for training and held-out data. This does not qualify the optimizer or spatial estimator by itself.

## F09 — incoherent variable-rank gradients

**Paths:** clusterwise argument/rank validation, initialization, loading combination, and spatial updates. Tests: `test-penalized-mfa-clusterwise.R` and new rank-contract cases.

**Implementation steps**

1. Enforce common `K <= min(p_s)` before graph construction or optimization. Remove padding/local-rank paths that imply a different objective from the stated blockwise Stiefel model.
2. Validate a positive integer K and each preprocessed block's effective rank using a scale-aware SVD threshold. Default to error; an explicit reduction policy reduces **all** blocks and reports requested/resolved rank.
3. Keep every `V_s`, gradient, optimizer state, and combined loading matrix at the same K. Do not fix the crash by recycling, trimming only one term, or appending random unsupported directions.
4. Update the existing equal-width test that expects silent ncomp=10 -> 5 truncation: assert a clear default error and test an explicit reduction option if implemented. This is an intentional contract correction, not a relaxation of the test.

**Regression tests**

- The audit's 2-/5-feature, K=3, lambda>0 case errors with a common-rank message before reaching the gradient. If reduction is requested, it uses K=2 in every returned block and records the change.
- K=2 on those same widths fits with coherent dimensions; cover both lambda=0 and positive lambda, and both exposed optimizer choices.
- Boundary K=min(p_s), K=1, rank-deficient centered data, constant columns, and fractional/nonfinite K obey the documented contract.
- For lambda=0, compare each returned loading subspace to an independent base SVD of the preprocessed block, using projector distance rather than signs.

**Done when:** the supported common-rank model either runs coherently or rejects unsupported ranks at the boundary. Full stationarity/descent qualification remains governed by the existing optimizer recovery milestone.

## F10 — k-NN indexing fails on small graphs

**Paths:** clusterwise neighbor construction and options validation. Tests: `test-clusterwise-graph.R`.

**Implementation steps**

1. Validate k as a finite nonnegative integer. Define k=0 as no k-NN edges; cap a requested positive k at the number of available non-self neighbors in the relevant graph-construction scope.
2. Use safe zero-length iteration/indexing rather than `1:n` or `2:(k+1)` for empty/singleton cases. Separate self-removal from neighbor truncation.
3. Reject empty public data blocks at the data boundary. An internal singleton graph has no edges and a zero Laplacian.

**Regression tests**

- P=1 internal graph and P=2 public two-block example; P=4 with the default k=6; k=0,1,P-1,P,and P+5. Assert valid complete expected edges or the documented empty graph, never NA indices.
- Invalid negative, fractional, NA, and infinite k fail clearly. Duplicate coordinates do not bypass self-removal.
- For a tie-free P=4 fixture, k>=P-1 produces a complete undirected graph and degree `(3,3,3,3)`.

**Done when:** every supported small input has explicit graph semantics, and no requested k creates invalid sparse indices.

## F11 — COVSTATIS scalar diagonal construction

**Status:** apparent repair already present, with an existing regression in `test-distatis.R`. Do not replace it with another implementation.

**Paths:** `R/covstatis.R` compromise scaling, `tests/testthat/test-covstatis.R`, existing `test-distatis.R` geometry test.

**Implementation/verification steps**

1. Retain the current columnwise `sweep()` scaling and explicit no-positive-component handling. Coordinate with the active COVSTATIS/DISTATIS writer.
2. Keep a core COVSTATIS one-component regression discoverable in its own test suite; reuse/move the existing case if appropriate instead of creating duplicate test definitions.
3. Extend coverage across normalization choices and actual rank-one inputs, not only truncation of a higher-rank Gram matrix.

**Regression tests**

- For centered `x=(-2,-1,1,2)`, use two copies of `B=x x'`, no double centering, and `norm_method="none"`. Scores/projection matrices are 4-by-1 and finite; score reconstruction equals B; the partial-score weighted barycenter equals global scores.
- Scale B to roots both below and above one. This catches length-one `diag(value)` dimension semantics that unit-root normalization can conceal.
- Exercise `none`, `frobenius`, and `mfa` normalization against an independent eigenvalue/normalization oracle; include K=2 as a non-regression case.
- A compromise with no positive retained root follows its documented explicit error policy.

**Done when:** these tests pass on the actual concurrent repair and the tested source hash is recorded. Do not require a red test on repaired code; the retained snapshot reproducer establishes that the original version failed.

## F12 — empty or incompatible inference statistics accepted

**Paths:** `R/infer_muscal.R` observation extraction, replicate extraction, accounting, and summary. Tests: `test-infer-muscal.R`.

**Implementation steps**

1. Validate the observed statistic before sampling: nonempty numeric finite vector for vector targets; nonempty finite numeric matrix for loading targets. Enforce the documented target type for custom extractors instead of flattening arbitrary lists/matrices silently.
2. Save the expected length/dimensions and, when present, component/feature identities. Every replicate must match; do not let `rbind()` recycle shorter vectors. Reorder named components only under an explicit validated identity contract; otherwise report mismatch.
3. Treat malformed replicate statistics as failures with replicate index/stage/reason. Only validated replicates contribute to `n_success` and summaries. If none survive, error. Include the effective replicate count in reported inference.
4. Validate nrep>=1 for permutation; define bootstrap uncertainty behavior for fewer than two valid draws explicitly (recommended: error). Sampling callback failures should use the same accounted failure path.
5. Retain and report failure rates. Omitting invalid draws does not establish calibrated inference when failures are selective; do not hide that limitation behind a finite summary.

**Regression tests**

- The audit's fit with absent `sdev` errors before its sampling/refit callback is invoked. A counter proves zero draws occurred.
- Observed `numeric(0)`, NULL, nonnumeric, NA/Inf, and a zero-dimension matrix fail with target-specific messages.
- Replicates of lengths 0,1,and3 against a length-2 observation are recorded as failures, not recycled. Test nonfinite vectors and matrices with mismatched dimensions/identities as well.
- A predetermined mixture of good/bad draws gives exact attempted/success/failure counts and exact means/quantiles based only on good draws. Compare to hand-calculated values.
- All-invalid replicates error; a valid scalar statistic retains a one-column matrix and finite summary; identity permutations produce p=1 with the correct valid-draw denominator.

**Done when:** successful inference always contains a nonempty validated target, and every reported success is accounted for consistently.

## F13 — unequal observation counts reach score addition

**Paths:** `R/penalized_mfa.R::penalized_mfa.list()` and delegating methods. Tests: `test-penalized-mfa.R`.

**Implementation steps**

1. Validate equal positive observation counts before preprocessing or optimization; include offending block names/counts in the message.
2. Assert preprocessing preserved the shared row layout before fitting. If future preprocessors are allowed to drop rows, require an explicit alignment contract rather than adding mismatched score matrices.
3. Apply the same validation through list/multiblock/multidesign entry points that reach this implementation. Keep heterogeneous **feature counts** valid where the selected penalty already supports them.

**Regression tests**

- Named blocks with 3 and 5 rows yield a clear shared-observation error, not `non-conformable arrays`. A preprocessor spy demonstrates rejection before preprocessing.
- Test both `compute_consensus=FALSE` and TRUE; the wrapper currently forms score summaries in either case.
- Equal-row blocks with different feature counts retain their existing supported behavior and warnings. Reuse existing fixtures so the row check cannot accidentally become an equal-width requirement.
- Empty-row input and a preprocessor that changes row counts are rejected explicitly; delegating wrappers preserve the message/contract.

**Done when:** invalid shared-observation inputs never enter optimization or consensus assembly, and valid heterogeneous widths continue to work.

## Additional audit qualifications: bounded follow-up work

These are separate from the thirteen confirmed defects. Addressing a documentation limit is not the same as proving an untested scientific claim.

| ID | Action and paths | Concrete tests / completion evidence |
|---|---|---|
| Q01. Benchmark qualification | `dev/run_design_mfa_benchmarks.R`, both iPCA scripts. Preserve raw failed/nonconverged rows; summarize attempted/failed/converged/unknown counts. Compute qualified timing/accuracy comparisons on explicitly eligible draws and report exclusions. Label `r2_y` as training adjusted R-squared; implement held-out fitting only if a predictive claim is retained. Record source hash/dirty state, package/dependency versions, and seed; avoid silently using an unrelated installed muscal. | Synthetic five-row benchmark table containing successes, failures, nonconvergence, and unknown convergence yields exact counts and eligible medians. Paired comparisons use the same jointly eligible replicate IDs and retain exclusion reasons. If held-out R2 is added, a fold-isolation fixture proves fitting/preprocessing never sees assessment rows. Extract only the small aggregation helper needed to test this; do not run expensive benchmarks as unit tests. |
| Q02. Loading reliability alignment | `loading_reliability()` and component orientation/alignment helpers. Match permutations/signs for separated components; diagnose near-ties and use subspace summaries or explicitly unavailable coordinate intervals. | Known sign flips/swaps leave aligned summaries unchanged; tied-root rotations preserve subspace measures without claiming uniquely identified coordinate intervals. A loading fixture with two near-equal roots exercises the gap policy. |
| Q03. Conditional target and exchangeability | Design-MFA inference/tuning docs and tests. State the frozen reference grid, conditional target, resampling unit, caller-owned null action, and the tuning score's effect-recovery interpretation. Keep the supplied-null callback requirement. | Invalid/mismatched callback outputs are rejected; subject-resampling and fold-local reference grids have direct identity tests. Calibrating a claimed rejection rate needs a separate protocol with a specified generator/null, replicate count, raw outcomes, binomial uncertainty, and locked code hashes. A callback contract test alone does not establish exchangeability. |
| Q04. Native argument contracts | `R/RcppExports.R`, native kernel inputs in `src/row_ops.cpp`. Put durable validation in the native implementations or stable R wrappers, not hand-edits to files regenerated by Rcpp. Validate list length, each system matrix dimension, RHS shape, graph dimension, finite inputs, and valid grouping indices. | Malformed shapes/types yield R errors before Eigen indexing. Valid small kernels agree with explicit base-R row-system operations, Laplacian multiplication, and `rowsum()`. No deliberate crash test in the main R process. |
| Q05. multifer integration | Existing adapter tests, `test-design-mfa-multifer.R`, supported engine-version checks. Declare supported capabilities; distinguish a missing engine from a broken compatible integration. | A pinned compatible multifer CI lane exercises scalar and multi-component inference plus the bootstrap stability route currently skipped locally. A deliberate unsupported-version fixture yields an actionable error. A live smoke remains separate from null calibration. |
| Q06. Undiscovered specification tests | Consolidate unique cases from `aligned_models_synthetic_tests.R` into discoverable `test-aligned-models-synthetic.R`; avoid duplicate helper/test definitions. Match normalization and control RNG/solver tolerance in true equivalence tests. | The standard suite discovers the unique cases. Equal-normalization MCCA and controlled-initialization named-index equivalence pass. A separate test can measure documented stochastic/convergence sensitivity without falsely asserting exact equality across unrelated initializations. |

Continue the broader clusterwise recovery plan's optimizer, finite-input, sparse-memory, and stationarity milestones before treating positive-lambda PMFA as solid. In particular, fixing the graph does not establish monotonic descent or stationarity of the existing fixed-step/Adam updates.

## Verification protocol

### Per repair

1. Capture the current source hash and run the smallest old-behavior reproducer, except where a concurrent repair already exists.
2. Add a regression expressing the **intended** contract. It must fail on the old version for the intended reason, not because of missing proposed API alone. Keep an independent numeric/output counterexample where relevant.
3. Make the bounded implementation/documentation change. Run the new regression plus relevant existing tests. Retain failures and actual exit status; do not infer success from a truncated log.
4. Record the new behavior, numerical tolerance rationale, changed public semantics, and limits. Update audit status/tag rationale only for the tested scope.

Canonical small floating-point fixtures should normally agree around `1e-10` to `1e-12`, adjusted for matrix norm, conditioning, and accumulated roundoff. Use exact equality for integer edge sets/counts and explicit error messages. Compare subspaces/reconstructions for sign/rotation-ambiguous decompositions. Do not loosen tolerances to hide a different objective or preprocessing path.

### Commands

The repository uses testthat/pkgload. From the repository root, a focused existing-file check is:

```sh
Rscript -e 'pkgload::load_all(".", compile=FALSE, quiet=TRUE); testthat::test_file("tests/testthat/test-infer-muscal.R", reporter="summary", stop_on_failure=TRUE)'
```

Substitute the affected existing/new test filename. `compile=FALSE` is appropriate only while using unchanged, current native code; native changes require recompilation and a new session before their tests.

After the integrated changes, run the standard suite once, saving raw output and the actual result:

```sh
Rscript -e 'testthat::test_local(".", reporter="summary", stop_on_failure=TRUE)'
```

Before a full build/check, inspect free space as required by the repository policy:

```sh
df -h /System/Volumes/Data
```

Use a task-owned temporary output directory and the repository's full package-check invocation (local example):

```sh
Rscript -e 'rcmdcheck::rcmdcheck(".", args="--no-manual", check_dir=tempfile("muscal-remediation-check-"), error_on="warning")'
```

Check `rcmdcheck` availability first. Missing dependencies/tooling are environmental blockers, not passes. Regenerate only affected Rd/NAMESPACE changes and inspect the diff; do not overwrite concurrent DISTATIS documentation or exports. R CMD check and local tests are distinct gates from hosted CI and publication. No merge, push, deployment, or release is part of this planning task.

## Closure record

For each F-ID, append to this plan or the existing issue record: implementation revision, regression names, before/after receipts, focused/full-check result, public behavior change, and remaining scientific qualification. Mark F11 closed only after the actual repair and its regression are verified together. A red snapshot reproducer and a passing current regression are different pieces of evidence and both should be retained.

Completion means every finding has either a verified repair or an explicitly unavailable/unsupported path that prevents the defective claim, with that restriction tested and documented. A paper-method name, a green generic suite, or a changed lifecycle tag is not a substitute for the specified regression.

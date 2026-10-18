# Muscal design / utility scientific audit

**Reviewer scope:** the 22 assigned source files only; read-only audit.

**Revision and source binding.** The review began at
`2c5c1e83692814dc20502e6963573e57263c07c6` and the source tree advanced
during review.  Final observations below apply to clean assigned sources at
`135e37742beec0f487777e0ba00c8d9e51f5124d` (the subsequent
`06eb85f` only changes Finder-metadata tracking).  The advance added the
Design-MFA multifer route and changed `R/design_mfa_inference.R` and
`R/multifer_feature_importance.R`; no conclusion silently refers to the
earlier dirty version.

**Post-review source drift.** After the review/reproducer completed,
`R/plot_mfa.R` became dirty in the shared worktree: reviewed hash
`42cfd1235b3a0aa48646b4a24edb0297dfa4a58e34abdee06a631cd1a452b067`,
current hash `52d97cfa1f1ed555f247408870ef422d059ee4d37c5f0176bbb8faedfbb637db`.
The one-line unreviewed addition at current line 690 dispatches a `distatis`
object to `plot_variance.distatis()`, but no definition of that function is
present under `R/`; this current dirty mutation will error for that class.
After `pkgload::load_all()`,
`plot_variance(structure(list(sdev=1), class="distatis"))` reproduced
`could not find function "plot_variance.distatis"`.
The manifest below intentionally binds the reviewed clean source, rather than
silently certifying this later change.

## Confirmed defects

### High — heterogeneous block dimensions produce invalid ICC and p-values

- **Location:** `R/utils.R:345-384`, especially 347-354.
- **Trigger:** `significant_components()` receives its documented
  `V_list` with blocks having different numbers of features, for example a
  3-by-1 loading matrix and a 6-by-1 loading matrix.
- **Mechanism:** the code evaluates `sum(vsc * vtc)`.  Those vectors refer
  to different feature spaces and have different lengths. R recycles the
  shorter vector (silently when the longer length is a multiple), then treats
  the result as a cosine similarity. There is no common-coordinate map,
  matching feature schema, or equal-length guard. The subsequent null mean
  `1/kbar` and z approximation therefore do not describe the computed
  quantity.
- **Consequence:** `keep`, `icc`, and `icc_pvalue` can claim
  significant cross-block component coherence for incomparable blocks. This
  is particularly serious because heterogeneous blocks are a normal
  multiblock use case.
- **Independent reproduction:** after
  `pkgload::load_all('.', compile = FALSE, quiet = TRUE)`, this fit
  `list(V_list=list(matrix(c(1,0,0),3,1), matrix(c(1,0,0,0,0,0),6,1)), sdev=1)`
  with `significant_components(..., n=20, k_vec=c(3,6), check_rmt=FALSE)`
  returned `icc=1`, `icc_pvalue=0.0169`, and `keep=1`. The two loading
  vectors cannot have a coordinatewise cosine because they have different
  feature domains.
- **Repair direction:** reject unequal feature spaces for ICC, or require an
  explicit validated feature correspondence/embedding and compute coherence
  in that common space. Do not retain an ICC p-value when neither exists.

### High — the RMT helper passes component count as feature count, and its MP edge cancels the aspect ratio

- **Location:** `R/components_utils.R:50-58`, especially line 51.
- **Trigger:** `estimate_components(fit, method = "rmt")` for any
  multiblock fit where block feature counts differ from the retained
  component count (the usual case).
- **Mechanism:** `k_vec <- rep(ncol(V_list[[1]]), length(V_list))` uses
  the number of component columns, not each loading matrix's feature rows.
  `significant_components()` uses `sum(k_vec)/n` to decide whether
  `gamma > 1` and should skip RMT at `R/utils.R:289-315`.
  Separately, lines 301/303 estimate
  `sigma2 = median(lambda)/(1+sqrt(gamma))^2`, then line 311 multiplies by
  exactly `(1+sqrt(gamma))^2`; the reported MP edge is therefore always
  the selected median eigenvalue for every non-skipped gamma. The aspect
  ratio has no numerical effect on that edge.
- **Consequence:** the wrong `k_vec` can choose the opposite branch:
  realistic high-dimensional blocks can be marked “RMT skipped” and pass all
  components under the real gamma, while `estimate_components()` uses a
  small component-count gamma and applies a threshold. The non-skipped
  threshold itself is not merely *mis-scaled* by wrong `k_vec`; it is
  algebraically independent of gamma, which is a separate defect in the
  intended MP calculation. It also feeds the dimension-incompatible ICC path
  above when blocks differ in row count.
- **Evidence:** `evidence/muscal-review-design-repro.R` asserts both
  branches: two 20-by-3 blocks with n=10 have true gamma=4 and skip RMT,
  while `estimate_components()` substitutes gamma=0.6 and applies the
  median-edge rule. The same script verifies equal reported edges for two
  different admissible gamma values.
- **Repair direction:** derive `k_vec <- vapply(V_list, nrow, integer(1))`
  after validating finite numeric matrices; then implement a stated
  Marchenko--Pastur noise estimator/edge without algebraic cancellation.
  Either remove/guard ICC as above or restrict it to a declared common
  feature space.

### Medium — `estimate_components(alpha=)` has no effect

- **Location:** `R/components_utils.R:37-59`.
- **Trigger:** any caller passes the documented `alpha` argument.
- **Mechanism:** `alpha` is accepted in the formal arguments but is never
  forwarded to `significant_components()` or otherwise read.
- **Consequence:** users can believe a different inferential threshold was
  applied when the fixed default alpha=0.05 was used internally.
- **Evidence:** static executable assertion in
  `evidence/muscal-review-design-repro.R`; the body of
  `estimate_components` contains no symbol `alpha`.
- **Repair direction:** either forward `alpha = alpha` to
  `significant_components()`, or remove it from the public signature and
  documentation until it governs a defined criterion.

### Medium — generic bootstrap accepts an absent `sdev` and reports successful empty inference

- **Location:** `R/infer_muscal.R:51-53, 97-131, 166-171, 192-204`.
- **Trigger:** a refit contract whose object and refits do not expose
  `sdev`, followed by the default `infer_muscal(..., statistic="sdev")`.
- **Mechanism:** `as.numeric(NULL)` is `numeric(0)`; this passes the
  numeric-vector check and each refit is counted as successful. The
  bootstrap summary becomes a 0-row tibble.
- **Consequence:** a result reports `n_success=nrep` despite containing no
  inferential statistic, which can be mistaken for valid uncertainty.
- **Reproduction:** a minimal valid refit specification returning
  `list(sdev=1)` while the observed fit has no `sdev` yielded
  `n_success=2` and a 0-row summary with `observed=numeric(0)`.
- **Repair direction:** require a nonempty finite observed vector and equal,
  nonempty finite replicate dimensions before resampling summaries.

## Important methodological concerns (not confirmed implementation defects)

1. **Experimental inference remains conditional.**
   `infer_design_mfa()` correctly refuses an implicit permutation
   (`R/design_mfa_inference.R:24-28,256-264`), but an arbitrary caller
   exchangeability callback is accepted and its preservation of the tested
   null is not checkable by the package. The cluster bootstrap also freezes
   the full-data reference grid and its weights
   (`R/design_mfa.R:471-515`), so intervals are conditional on that target;
   they do not include uncertainty from estimating a data-derived grid or
   observed design distribution. This is clearly closer to the documented
   conditional procedure than an error, but it is not population-level
   confirmatory inference.

2. **Loading intervals only correct signs, not rotations or component swaps.**
   `loading_reliability()` aligns every bootstrap basis independently by
   per-column signs (`R/components_utils.R:81-102`). Near-degenerate
   components can rotate or exchange order across bootstrap fits, in which
   case pointwise percentile intervals and sign consistency conflate
   nonidentifiability with instability. `infer_design_mfa()` does apply
   orthogonal Procrustes for its own loadings
   (`R/design_mfa_inference.R:36-57`); the generic routine should make the
   same distinction or explicitly limit claims to well-separated,
   ordered components.

3. **The tuning score is a proxy, not held-out response prediction.**
   `tune_design_mfa()` fits an assessment-only target effect and compares
   estimated effect fields (`R/design_mfa_tune.R:201-225`). Subject
   isolation and fold-local reference grid construction are good and no
   leakage was found in this path. Still, the target is itself noisy and
   rank-1/unpenalized by construction; its relative effect MSE is not a
   held-out trial/activation prediction loss. Report it as stability against
   an independently estimated effect, as the code does, rather than as
   predictive accuracy.

4. **MCCA/multifer adapter nulls are only conditionally valid.**
   The adapters in `R/multifer_feature_importance.R` and the new
   `R/design_mfa_multifer.R` declare conditional validity and caller-owned
   exchangeability. The live compatible effect-significance engine was smoke
   exercised here at rank one for `partial_ss` alone and for
   `c("partial_ss", "singular_values")`: both returned `engine="multifer"`
   with three successful draws and the latter produced its expected one-row
   singular-value table. This establishes only adapter execution, not
   contract-version integration or null calibration. The multifer stability
   capability remained unavailable in this environment, so its bootstrap path
   was not exercised. Component/feature p-values remain experimental until
   those integration and calibration checks exist.

5. **Native entry points rely on R-side contracts.**
   `src/row_ops.cpp:9-60` assumes `length(A_list)==nrow(S)`, every
   `A_i` is k-by-k, and a graph Laplacian is n-by-n. The internal callers
   construct these objects, but the generated R wrappers are callable and
   do not validate shape. This is a robustness boundary concern; it was not
   exercised with malformed native input to avoid risking a native crash.

6. **Seed restoration is a convention difference, not a defect.**
   `infer_muscal(seed=)` changes the caller RNG state
   (`R/infer_muscal.R:62-64`), whereas Design-MFA tuning/inference restores
   it. R's usual seeded-call convention permits this behavior and the generic
   function does not promise restoration, so this is only a reproducibility
   API consistency concern.

## Method classification

“Solid” below means the implementation has a clear deterministic contract and
the reviewed tests/oracles cover its numerical behavior. It is not a clinical,
causal, or population-inference endorsement.

| Source | Public/statistical methods reviewed | Tag and reason |
|---|---|---|
| `R/RcppExports.R`, `src/RcppExports.cpp`, `src/row_ops.cpp` | native row-system/operator and grouped sums | **Solid kernel, conditional interface.** Straight linear algebra and index checks for grouped sums; public wrappers need shape guards. |
| `R/components_utils.R` | `estimate_components` | **Experimental.** Wrong `k_vec`, cancelled MP edge, ignored `alpha`, and invalid heterogeneous-block ICC route above. |
|  | `loading_reliability` | **Experimental.** Bootstrap pointwise summaries lack rotational/permutation alignment. |
| `R/cv_muscal.R` | `cv_muscal` | **Solid wrapper, conditional evaluation.** It delegates splits and callbacks; it cannot prove caller preprocessing/selection stayed inside analysis folds. |
| `R/design_mfa.R`, `R/design_mfa_design.R`, `R/design_mfa_exact.R` | `design_mfa.default/list/multidesign`, `design_mfa_control`; `effects.design_mfa`, `components.design_mfa`, `effect_scores/loadings/map`, `estimability.design_mfa`, `alignment_coverage`, `subject_effects`, `predict.design_mfa`, `project.design_mfa_effect` | **Experimental research method.** Strong input/estimability/rank/metric guards and independent oracle/metamorphic tests exist, but partial-RDA/ASCA effects, soft alignment, ridge slopes, and spatial metric decomposition lack external calibration and should not be treated as confirmatory effect estimates. |
| `R/design_mfa_inference.R`, `R/design_mfa_multifer.R` | `infer_design_mfa`, `print.design_mfa_inference`, `as_multifer_adapter.design_mfa` | **Experimental/conditional.** Correctly requires explicit permutation exchangeability and clusters subjects. The live multifer permutation smoke path passed at rank one, but its null validity is caller-owned, target is frozen, and no calibration or multifer-stability evidence qualifies the route. |
| `R/design_mfa_maps.R` | `coordinate_feature_map(s)` | **Solid numerical mapping, conditional scientific use.** Nearest-neighbor/kernel construction and coverage guards are well specified; coordinate registration remains user supplied. |
| `R/design_mfa_tune.R` | `tune_design_mfa`, `nested_tune_design_mfa` | **Experimental selection/evaluation.** Subject-disjoint, fold-local construction is good; the score is an independently fitted effect proxy rather than prediction performance. |
| `R/fit_contract.R` | `muscal_fit_contract`, `predict.mfa`, `project/predict.anchored_mfa`, `project/predict.aligned_mfa` | **Solid software contract.** Stored-preprocessing projection checks are clear; statistical interpretation depends on the fitted method. |
| `R/infer_muscal.R` | `infer_muscal` | **Experimental.** Method-specific bootstrap/permutation callbacks own exchangeability; empty statistics are a confirmed defect. |
| `R/multifer_feature_importance.R` | `as_multifer_adapter.{mcca,aligned_mcca,anchored_mcca,aligned_rrr,bada,covstatis}` | **Experimental.** External-engine/API and null assumptions are conditional; no live-engine calibration evidence in this audit. |
| `R/muscal-package.R`, `R/tibble-imports.R` | package/import declarations | **Solid non-statistical plumbing.** |
| `R/orient_components.R` | `orient_components.{mfa,aligned_mfa,linked_mfa,aligned_interbattery,mcca,aligned_mcca,response_aligned_mfa,aligned_rrr,ipca,penalized_mfa,default}` | **Solid sign-orientation utility.** It preserves sign-equivalent reconstructions under its listed fields; it does not solve component matching/rotation. |
| `R/performance_metrics.R` | `metric_registry`, `default_metrics`, `performance_metrics` | **Solid metric computation.** Metrics have explicit task contracts; they are descriptive and do not create inferential validity. |
| `R/plot_mfa.R` | Exported helpers `plot_partial_scores`, `plot_block_weights`, `plot_block_similarity`, `plot_variance`, `plot_loadings`, `plot_convergence`, `plot_coverage`, `plot_feature_groups`, `plot_block_fit`, `plot_global_local`, `plot_barycenters`, and `plot_loading_alignment`; S3 methods `plot.mfa`, `plot.linked_mfa`, `plot.covstatis`, `plot.bamfa`, `plot.bada`, `plot.penalized_mfa`, `autoplot.mfa`, `autoplot.linked_mfa`, `autoplot.covstatis`, `autoplot.bamfa`, `autoplot.bada`, `autoplot.penalized_mfa`, `plot_partial_scores.{mfa,linked_mfa,covstatis,bamfa,bada,penalized_mfa}`, `plot_block_weights.{mfa,linked_mfa,covstatis}`, `plot_block_similarity.{mfa,covstatis}`, `plot_variance.{covstatis,bada,penalized_mfa}`, `plot_loadings.{mfa,linked_mfa,penalized_mfa}`, and `plot_convergence.{linked_mfa,bamfa,penalized_mfa}` | **Solid visualization utilities.** Plots inherit the scientific limits of their fitted inputs. |
| `R/synthdat.R` | `synthetic_multiblock` | **Experimental simulation fixture.** It is a useful low-rank Gaussian generator, not a realistic generative validation of alignment, dependence, or inference. |
| `R/utils.R` | `significant_components` | **Experimental and currently invalid for heterogeneous blocks.** See high-severity finding. |
|  | `prepare_block_preprocessors`, `within_class_scatter`, `between_class_scatter` and other non-inferential helpers | **Solid utilities** under their input contracts. |

## File coverage and evidence read

All assigned source files were individually read: `R/RcppExports.R`,
`R/components_utils.R`, `R/cv_muscal.R`, `R/design_mfa.R`,
`R/design_mfa_design.R`, `R/design_mfa_exact.R`,
`R/design_mfa_inference.R`, `R/design_mfa_maps.R`,
`R/design_mfa_multifer.R`, `R/design_mfa_tune.R`,
`R/fit_contract.R`, `R/infer_muscal.R`,
`R/multifer_feature_importance.R`, `R/muscal-package.R`,
`R/orient_components.R`, `R/performance_metrics.R`, `R/plot_mfa.R`,
`R/synthdat.R`, `R/tibble-imports.R`, `R/utils.R`,
`src/row_ops.cpp`, and `src/RcppExports.cpp`.

Related callers inspected included `R/linked_mfa.R:580-660`. Related tests
read were:
`test-components_utils.R`, `test-cv-muscal.R`,
`test-design-mfa-contract.R`, `test-design-mfa-coordinate-maps.R`,
`test-design-mfa-inference.R`, `test-design-mfa-metamorphic.R`,
`test-design-mfa-multifer.R`, `test-design-mfa-multilevel.R`,
`test-design-mfa-oracles.R`, `test-design-mfa-spatial.R`,
`test-design-mfa-spatial-recovery.R`, `test-design-mfa-tuning.R`,
`test-infer-muscal.R`, `test-multifer-feature-importance.R`,
`test-orient-components.R`, `test-performance-metrics.R`,
`test-summary-utils.R`, `test-synthdat.R`, `test-utils-extra.R`,
`test-utils-misc.R`, `test-performance-guardrails.R`, and
`helper-design-mfa-sim.R`.

Manual pages read: `man/design_mfa.Rd`, `design_mfa_control.Rd`,
`infer_design_mfa.Rd`, `tune_design_mfa.Rd`,
`nested_tune_design_mfa.Rd`, `coordinate_feature_map.Rd`,
`coordinate_feature_maps.Rd`, `infer_muscal.Rd`, `cv_muscal.Rd`,
`performance_metrics.Rd`, `estimate_components.Rd`,
`orient_components.Rd`, `as_multifer_adapter.Rd`,
`significant_components.Rd`, and `synthetic_multiblock.Rd`; relevant
sections of `vignettes/design_mfa.Rmd`, `model_evaluation.Rmd`,
`mfa.Rmd`, `mcca.Rmd`, `ipca.Rmd`, `aligned_mfa.Rmd`,
`linked_mfa.Rmd`, and `penalized_mfa.Rmd` were also checked.

Validation evidence: `pkgload::load_all('.', compile=FALSE, quiet=TRUE)`
passed. The parent’s contemporaneous full-suite receipt reports 487 blocks /
48 files, 2,125 passing expectations, zero failures/errors, and three skips;
raw log: `evidence/muscal-review-test.log`. That suite did not cover the
heterogeneous ICC case or absent-`sdev` generic inference boundary above.
`evidence/muscal-review-design-repro.R` also passed, including the live
rank-one multifer permutation call with both `partial_ss` and
`singular_values`; it dispatched at `R/design_mfa_inference.R:270-307` into
`R/design_mfa_multifer.R:291-385`, returned `engine="multifer"`, and had
three successful draws. The corresponding multifer stability capability was
false in this session, so no bootstrap-engine claim is made.

## Computational-script coverage

The following scripts were individually inspected after the source audit.
They are research/qualification artifacts and should not be used to upgrade a
method tag without their recorded raw outputs and independent review.

| Script | Tag | Findings |
|---|---|---|
| `dev/benchmark_design_mfa.R` | **Experimental benchmark.** | It retains per-fit `solver_converged`, residual, data size, allocation, environment, and streamed-vs-dense agreement fields. It has no warm-up/control distribution and is a single-process timing specimen, so it does not substantiate a general speed claim. |
| `dev/run_design_mfa_benchmarks.R` | **Experimental benchmark aggregation.** | **Qualification defect at lines 52-94:** `solver_converged` is omitted from `numeric_summary` and no nonnumeric convergence count is summarized. A summary can thus publish median timings from nonconverged fits with no visible qualification. Preserve/consolidate convergence counts and refuse a performance claim when any included fit is nonconverged. |
| `dev/calibrate_design_mfa_spatial.R` | **Experimental recovery calibration.** | It uses direct-coordinate, same-index synthetic fields and has transparent fixed thresholds. The null simulations at lines 107-110 only record effect energy; they do not calibrate a rejection rate, so this script supports recovery comparisons, not false-positive/inference qualification. |
| `dev/calibrate_design_mfa_spatial_confirmatory.R` | **Bounded confirmatory simulation, not general validation.** | Its fixed protocol/hash evidence, independent-coordinate and misspecification controls, and max-root familywise calculation (lines 364-381) are stronger than the developmental script. Its conclusions remain limited to this generator, the supplied within-probe shuffle, 60/100 simulations, and the frozen algorithm/hash set. |
| `dev/distatis_research.R` | **Research reproducer.** | Independent published-target and geometry assertions are useful. Current rank-one `covstatis` failures and assessor-scaling discrepancy are captured as observations, not silently passed. Parent's live receipt is `evidence/muscal-review-distatis.log`; it reports the rank-one non-conformable error and passed independent assertions. |
| `dev/ipca_benchmark_dense_vs_gram.R` | **Experimental benchmark.** | It records each fit's convergence and iteration counts, but still computes and prints speed/subspace comparisons whenever fits return, including nonconverged fits (lines 77-100). Filter or label such rows before interpreting equivalence or speed. |
| `dev/ipca_vs_multiset_synth.R` | **Experimental comparative simulation.** | **Qualification defect at lines 244-260 and 280-302:** `r2_y` is adjusted R-squared from a regression trained and evaluated on the same simulated rows; it is an association diagnostic, not prediction performance. The summary filters only `ok`, not `converged`, so iPCA rows with `converged=FALSE` can affect similarity, R2, time, and noninferiority results. It also loads an installed `muscal` when present without recording its commit/version, which can decouple results from the checkout. |
| `scripts/build-pkgdown.R` | **Solid build helper.** | In-process rendering with isolated temporary caches and fixed seeds is appropriate for this local workaround. It calls pkgdown internals, so it is version-sensitive; no statistical claims originate here. |

The nonstandard test discovery gap is a separate test-harness concern, not
evidence that a statistical method passed or failed. Parent's explicit receipt
for the omitted synthetic test file is
`evidence/muscal-review-extra-tests.log` (two failures), and should be
reported separately from algorithm diagnosis.

## SHA-256 source manifest (reviewed clean state)

```
067734453825964f8ca90182c1a17494331634eec93b8e2a30faf1b201b309cc  R/RcppExports.R
356cb29a55c93a6628ec675b3ce3da780e5401ed51bf195fd1e9ae60127ed0e5  R/components_utils.R
2aa3f5598c06fe30d34c257a869da6efa8776b19af83ca968ccc2e055b62fae2  R/cv_muscal.R
ab9f4a58e38abae812b743ea92574c9dac51a003df8945520814d0c8e3d0d68f  R/design_mfa.R
f5b152d3807c9600c0965cc33131e0ff3016c7a9b7936eb4b794b1c4376f88a0  R/design_mfa_design.R
f683d961a7ecb2bbee0e2f326a96e2b6068694ec34001a2c8933d172117a06a2  R/design_mfa_exact.R
1adcf61edee10ed57c17df9eb6390f11e95549693db850889b8fc86f44088ad9  R/design_mfa_inference.R
644eef818d6fc615f8ce617666d8b93bb2968761cf9f85f18a3ea6366aabc470  R/design_mfa_maps.R
eab1e33fc977514cce11b17f3294c2bbffdd5abbdd66860ec96d5fa4dbd01c41  R/design_mfa_multifer.R
3317c83cf6518841e5bba9c7fbd96268692ea5eee86b64286b0e78b2a3ebd4f6  R/design_mfa_tune.R
d9cd2f3d9ba96067afd5b8fd39bc2ab264db7432d5e7b83e9e3a033972d86b3c  R/fit_contract.R
83467073dbc0b3bc11ed1255e8e2f87443050b14f6aad4ad3992cb2e8c9a8e93  R/infer_muscal.R
78353d0a2ea49fc06d73eeaffe5c64025681707960a4a4f0c8088c71bf9278c6  R/multifer_feature_importance.R
dd2b40146e5fcc47edabc158b37243eeababf37c99ad5a99a5b48cc7f4a03d55  R/muscal-package.R
33f356ac1017cca76a556b5e3d8bcce65c815f9e4a01e64b1ffc9743df6d4a6d  R/orient_components.R
336ca96c9d0c108be60b928cbf4d9aaa12eec144c110557ea240b964a60caae9  R/performance_metrics.R
42cfd1235b3a0aa48646b4a24edb0297dfa4a58e34abdee06a631cd1a452b067  R/plot_mfa.R
571ff246733f0abde5379b9b3efae2b8046f2c5fd833b4e2fc6acb826615bf23  R/synthdat.R
16679d3665bb2a0b050acd5ae428df6e7762d7e313dc671903fe1017f9f456ce  R/tibble-imports.R
fe212daf1464301ffe8c496d0711fb0876e5a5fa162dd89c420289214ac55341  R/utils.R
2bf9bc7d6f0f31a8b54c7f8aad63b5d8ee6d6ff53685de1ec5e88eeeb9751f10  src/row_ops.cpp
7532c0e35850f12e002dd18045caf71d82fded857b3d1c2d0df536c8705e03bf  src/RcppExports.cpp
```

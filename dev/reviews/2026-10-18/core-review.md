# muscal core-method scientific/code audit

Reviewed at repository HEAD `2c5c1e83692814dc20502e6963573e57263c07c6` on 2026-10-18.  This was a read-only audit; no repository files were changed.  Source evidence is bound to the SHA-256 manifest at the end.  The original snapshot has **8 positive findings** (C1-C3, C6-C10); 7 remain reproducible in the later current checkout, while C6 was externally repaired during the audit. I separately ran `pkgload::load_all('.', compile = FALSE, quiet = TRUE)` successfully and the 21 relevant focused test files (all passed; only installed-under-R-4.5.2 dependency warnings).  Parent also reported the full current suite green.

## Confirmed findings

### C1 — default clusterwise PMFA k-NN graph fails for fewer than seven total clusters

* **Severity:** moderate — valid small spatial examples cannot be fit with defaults.
* **Path:** `R/penalized_mfa_clusterwise.R:41-68`, reached from `:744-745`.
* **Trigger:** total rows in `coords_list` is below `k_nn + 1` (default `k_nn = 6`).
* **Consequence:** `order(...)[2:(k_nn + 1)]` contains `NA` positions, which are passed as sparse-matrix column indices.  The default advertised method errors instead of using the available `min(k_nn, n_total - 1)` neighbours.
* **Evidence/reproducer:**
  ```r
  pkgload::load_all('.', compile = FALSE, quiet = TRUE)
  penalized_mfa_clusterwise(
    list(matrix(rnorm(4), 2, 2), matrix(rnorm(4), 2, 2)),
    list(matrix(c(0,0,0, 1,0,0), 2, 3, byrow=TRUE),
         matrix(c(0,1,0, 1,1,0), 2, 3, byrow=TRUE)),
    ncomp = 1, max_iter = 1)
  # Error: 'NA' indices are not (yet?) supported for sparse Matrices
  ```
  **Confidence: high.**

### C2 — `penalized_mfa()` accepts unequal row counts, then fails internally

* **Severity:** moderate — shared-observation prerequisite is not validated at the public boundary.
* **Path:** `R/penalized_mfa.R:326-350`.
* **Trigger:** at least two blocks have different `nrow` values.
* **Consequence:** each block is optimized, then the score matrices are added at line 349.  This gives an opaque `non-conformable arrays` error.  The method must either reject unequal row counts before fitting or define an alignment model; treating it as a conventional MFA consensus is scientifically invalid.
* **Evidence/reproducer:**
  ```r
  pkgload::load_all('.', compile = FALSE, quiet = TRUE)
  penalized_mfa(list(matrix(rnorm(12), 3, 4), matrix(rnorm(20), 5, 4)),
                ncomp = 1)
  # Error: non-conformable arrays
  ```
  **Confidence: high.**

### C3 — BaDA bootstrap intervals do not bootstrap component/loadings estimation

* **Severity:** high for inferential use — reported `boot_lds_*` intervals can be read as loading uncertainty, but the PCA basis is fixed.
* **Path:** `R/bada.R:104-138`, especially `:121-128`.
* **Trigger:** `bootstrap.bada(x, nboot, data = ...)` and interpretation of its documented component/loading uncertainty estimates.
* **Consequence:** each replicate resamples subject strata and creates a barycenter `Xc`, but then fixes `v_boot <- x$v`, `fs_boot <- x$fscores`, and the original variances.  It never refits PCA/discriminant directions (`pca_group` is commented out at 132-133) or aligns a refitted solution.  Thus `boot_j` is a projection through the original solution, not a bootstrap loading estimate; percentile summaries cannot quantify loading/component sampling variability claimed at `:53-74`.
* **Evidence:** direct control-flow inspection establishes fixed `x$v`, `x$fscores`, and `sdev(x)` for every replicate.  A valid loading bootstrap needs a per-replicate refit plus component matching/sign alignment before aggregation.  **Confidence: high.**

### C6 — historical snapshot: one-component COVSTATIS errored through scalar `diag()`

* **Severity:** moderate in the reviewed snapshot; apparently repaired in the current externally edited file.
* **Snapshot path:** `R/covstatis.R:180-188` at SHA-256 `caf44be8834e8f2876d093989616e2c14471785f7d6cbc19738eff235944168b`.
* **Trigger/evidence:** the frozen `dev/distatis_research.R` receipt at `evidence/muscal-review-distatis.log` recorded both rank-one and one-component MFA COVSTATIS failures as `non-conformable arguments`. The snapshot constructed score/projection matrices using scalar `diag(sqrt(fit$values[keep]))` forms.
* **Drift qualification:** during this audit an external writer changed `R/covstatis.R` to SHA-256 `314180145c74b4e628e995d1e389b213a04e4292f9bc1fd524708af295620671`, replacing those constructions with `sweep()` (current lines 206-207). The current simple rank-one check now succeeds. This report does not claim the repair is fully audited; it records the defect against the original review snapshot only.
* **Confidence: high for the historical defect and source drift; no claim about current broader COVSTATIS correctness.**

### C7 — `normalized_laplacian = TRUE` is algebraically identical to the unnormalized path

* **Severity:** high for spatial interpretation — the API/documentation claims degree-normalized smoothing but applies the same Laplacian.
* **Path:** `R/penalized_mfa_clusterwise.R:744-752`.
* **Trigger:** any call with default `normalized_laplacian = TRUE` and `lambda > 0`.
* **Consequence:** `cluster_graph()` has already returned `L = D - A`. The normalization computes `rowSums(L + I)`, which is exactly one because Laplacian rows sum to zero, then multiplies by the identity scaling. It does not compute `D_A^{-1/2}(D_A-A)D_A^{-1/2}` from adjacency degrees. Results labelled normalized retain degree dependence.
* **Evidence:** direct algebra from `cluster_graph()` at `:106-114` and normalization branch at `:748-751`; this also matches the independently recorded recovery plan. **Confidence: high.**

### C8 — graph symmetrization deletes directed neighbour edges

* **Severity:** high for spatial objective — the constructed graph depends on arbitrary feature ordering.
* **Path:** `R/penalized_mfa_clusterwise.R:57-66`.
* **Trigger:** a directed k-NN edge occurs only in the lower triangle of the adjacency matrix before symmetrization.
* **Consequence:** `Matrix::forceSymmetric(A, uplo = "U")` copies the upper triangle and discards lower-only edges; it is not a union or mutual-neighbour symmetrization. Reordering feature identities can move a one-way edge between triangles and change the graph/objective after mapping back.
* **Evidence/reproducer:**
  ```r
  A <- Matrix::Matrix(0, 3, 3, sparse = TRUE); A[3, 1] <- 1
  Matrix::forceSymmetric(A, uplo = "U")[3, 1]
  # 0; the lower-only edge was 1 before the call
  ```
  **Confidence: high.**

### C9 — clusterwise PMFA discards fitted preprocessing in its returned projector

* **Severity:** high for OOS projection/reconstruction — stored metadata applies raw, not training centering/scaling.
* **Path:** `R/penalized_mfa_clusterwise.R:727-736,1079-1087`.
* **Trigger:** default `preproc = NULL` (which fits `center()`) or any non-pass block preprocessor, then use the returned object for generic projection/transform.
* **Consequence:** fitting uses `Xp` from fitted block preprocessors, but the return object installs a newly fitted `pass()` preprocessor on a dummy matrix. An out-of-sample call therefore has a different coordinate system from the optimization and invalidates projection/reconstruction interpretation.
* **Evidence/reproducer:**
  ```r
  set.seed(2)
  xs <- list(matrix(rnorm(20, mean=10),5,4), matrix(rnorm(20, mean=-5),5,4))
  cs <- list(cbind(1:4,0,0), cbind(1:4,1,0))
  f <- penalized_mfa_clusterwise(xs, cs, ncomp=1, lambda=0, max_iter=1)
  z <- multivarious::transform(f$preproc, xs[[1]], colind=f$block_indices[[1]])
  max(abs(z - xs[[1]])) # 0, although fitting centered this block
  ```
  **Confidence: high.**

### C10 — variable-rank cluster blocks create incompatible spatial gradients

* **Severity:** high — documented variable-rank support crashes whenever spatial coupling is active.
* **Path:** `R/penalized_mfa_clusterwise.R:714-725,852-875,163-174`.
* **Trigger:** `ncomp` exceeds the feature count of one block but not another, with `lambda > 0`.
* **Consequence:** loadings are padded to the largest block rank in `Vbig_begin_iter`, while the block reconstruction gradient has its local rank. The spatial term has the padded column count and cannot be added to the local reconstruction gradient. This is not a coherent common-component objective.
* **Evidence/reproducer:**
  ```r
  penalized_mfa_clusterwise(
    list(matrix(rnorm(20),10,2), matrix(rnorm(50),10,5)),
    list(cbind(1:2,0,0), cbind(1:5,1,0)),
    ncomp=3, lambda=1, max_iter=1)
  # Error: number of rows are not compatible for +
  ```
  **Confidence: high.**

## File-by-file coverage, method status, and supporting evidence

Tags mean **solid** only for the stated computational scope (not inferential qualification); **experimental/research** means the method is an unqualified optimization, novel alignment/graph construction, heuristic tuner, or resampling/inference adapter.  Public S3 groups are named explicitly.  “No confirmed finding” does not assert proof beyond the code, focused tests, and checks stated above.

| Source file | Public statistical methods/S3 group reviewed and tag | Tests/docs read | Result |
|---|---|---|---|
| `R/aligned_interbattery.R` | `aligned_interbattery`; `project.aligned_interbattery`, `predict.aligned_interbattery` — **experimental** (joint subject coupling, whitening and graph/decorrelation penalties) | `test-aligned-interbattery.R`, `aligned_models_synthetic_tests.R`; `man/aligned_interbattery.Rd` | No confirmed defect. Boundary validation and OOS routing reviewed. |
| `R/aligned_mcca.R` | `aligned_mcca`, `anchored_mcca` — **experimental** (row-aligned/anchored MCCA); S3 orientation contract consumed through `orient_components` | `test-aligned-mcca.R`, `aligned_models_synthetic_tests.R`; `man/aligned_mcca.Rd` | No confirmed defect. The synthetic equivalence test omitted `normalization="None"` on its `aligned_mcca()` call; matched calls give subspace RMSE 0. |
| `R/aligned_mfa.R` | `aligned_mfa` — **experimental**; orientation/project/predict S3 methods are supplied in supporting code | `test-aligned-mfa.R`, `aligned_models_synthetic_tests.R`; `man/aligned_mfa.Rd` | No confirmed defect. The failed synthetic permutation test advances the RNG and retains default `tol=1e-6`; equal RNG with `tol=1e-10` gives RMSE 0. |
| `R/aligned_rrr.R` | `aligned_rrr`; `project.aligned_rrr`, `predict.aligned_rrr` — **experimental** supervised shared-response reduced-rank fit | `test-aligned-rrr.R`; `man/aligned_rrr.Rd` | No confirmed defect. |
| `R/aligned_rrr_multifer.R` | `infer_aligned_rrr` — **experimental/research inference adapter** | `test-infer-muscal.R`; `man/infer_aligned_rrr.Rd`, `man/infer_muscal.Rd` | No confirmed defect; inference remains contingent on multifer assumptions/resampling design. |
| `R/all_generic.R` | Generics `bada`, `covstatis`, `mfa`, `mcca`, `ipca`, `project_cov`, `project_subjects`, `project_covariate`, `project_table_covariate`, `project_feature_covariate` — **solid dispatch-only helpers** | corresponding method tests; no separate Rd | No confirmed dispatch defect. |
| `R/bada.R` | `bada.multidesign` — **solid point-estimation implementation**; `bootstrap.bada` — **experimental/unsafe for loading inference**; `reprocess.bada`, `project.bada` — **solid conditional projections** | `test-bada-multidesign.R`; `man/bada.Rd` | **C3** for bootstrap; no confirmed point-fit/projection defect. |
| `R/bada_multifer.R` | `infer_bada` — **experimental/research inference adapter** | `test-infer-muscal.R`; `man/infer_bada.Rd`, `man/infer_muscal.Rd` | No confirmed defect. |
| `R/bamfa.R` | `bamfa`, `bamfa.default`, `bamfa.list`, `bamfa.multiblock`, `bamfa.multidesign`, `predict.bamfa`, `print.bamfa` — **experimental** alternating global/local factor model | `test-bamfa.R`, `test-bamfa-recon.R`; `man/bamfa.Rd` | No confirmed defect. |
| `R/bilinear_mixed.R` | `bilinear_mixed`, `bilinear_mixed.default`, `bilinear_mixed.list`, `print.bilinear_mixed` — **experimental** repeated-measure latent model | `test-bilinear-mixed.R`; `man/bilinear_mixed.Rd` | No confirmed defect. |
| `R/bilinear_mixed_tune.R` | `bilinear_mixed_recommend`, `bilinear_mixed_easy`, `bilinear_mixed_tune`, `print.bilinear_mixed_tuning` — **experimental** heuristic/CV selection | `test-bilinear-mixed.R`; `man/bilinear_mixed_easy.Rd`, `man/bilinear_mixed_recommend.Rd`, `man/bilinear_mixed_tune.Rd` | No confirmed defect. Subject-level held-out folds are used; score inference uses held-out data only. |
| `R/covstatis.R` | `covstatis.list` — **solid descriptive COVSTATIS in snapshot apart from C6**; `project_cov.covstatis`, `reconstruct.covstatis`, `project_subjects.covstatis`, `project_covariate.covstatis`; `check_duality`, `roi_map`, `plot_subjects` — **solid helpers/S3 group** | `test-covstatis.R`, `test-duality.R`, `test-project-covariate.R`; `man/covstatis.Rd`; independent `dev/distatis_research.R` receipt | **C6 against snapshot only.** COVSTATIS is distinct from proposed RV-normalized DISTATIS assessor geometry. External edits during audit mean current repair is not fully reviewed. |
| `R/covstatis_multifer.R` | `infer_covstatis` — **experimental/research inference adapter** | `test-infer-muscal.R`; `man/infer_covstatis.Rd`, `man/infer_muscal.Rd` | No snapshot defect. Externally modified during audit; current SHA differs from manifest. |
| `R/graph_anchored_mfa.R` | `graph_anchored_mfa`, `coupled_graph_anchored_mfa`; `project.*`, `predict.*` for both classes — **experimental** graph/anchor-constrained factor models | `test-graph-anchored-mfa.R`, `test-graph-anchored-mfa-synthetic.R`, `test-coupled-graph-anchored-mfa.R`; `man/graph_anchored_mfa.Rd` | No confirmed defect. |
| `R/ipca.R` | `ipca.list`, `ipca.multiblock`, `project.ipca`, `predict.ipca` — **experimental implementation of iPCA**; `ipca_tune_alpha` — **experimental masked-entry tuner** | `test-ipca.R`; `man/ipca.Rd`, `man/ipca_tune_alpha.Rd` | No confirmed defect. The tuning preprocessor is fit to masked training values, then applied to evaluation values. |
| `R/linked_mfa.R` | `anchored_mfa`, `linked_mfa`; `project.anchored_mfa`, `predict.anchored_mfa` — **experimental** linked/anchor constrained MFA | `test-linked-mfa.R`; `man/linked_mfa.Rd` | No confirmed defect. |
| `R/mcca.R` | `mcca.list`, `mcca.multiblock`, `predict.mcca` — **solid descriptive MCCA implementation** (regularized and balanced modes are descriptive, not inferential) | `test-mcca.R`; `man/mcca.Rd` | No confirmed defect. |
| `R/mfa.R` | `mfa.list`, `mfa.multiblock` — **solid complete-data MFA**; missing-data `regularized`/`em` paths — **experimental** iterative imputation | tests found through `test-mcca.R`/general MFA callers and full-suite receipt; `man/mfa.Rd`, `man/mfa-plots.Rd` | No confirmed defect. `missing="nipals"` is honestly rejected as unimplemented. |
| `R/penalized_mfa.R` | `penalized_mfa`, `penalized_mfa.list`, `penalized_mfa.multiblock`, `penalized_mfa.multidesign`, `print.penalized_mfa`, `penalized_mfa_engine` — **experimental** manifold optimization with penalties | `test-penalized-mfa.R`, `test-penalized-mfa-opts.R`; `man/penalized_mfa.Rd`, `man/penalized_mfa_engine.Rd` | **C2**. |
| `R/penalized_mfa_clusterwise.R` | `penalized_mfa_clusterwise`, alias `pmfa_cluster`, `print.penalized_mfa_clusterwise`; helpers `spatial_constraints`, `cluster_graph` — **experimental** spatial graph-regularized optimization | `test-penalized-mfa-clusterwise.R`; `man/penalized_mfa_clusterwise.Rd`; `dev/penalized_mfa_clusterwise_recovery_plan.md` hypotheses independently checked | **C1, C7-C10**. |
| `R/response_aligned_mfa.R` | `response_aligned_mfa`; `project.response_aligned_mfa`, `predict.response_aligned_mfa` — **experimental** supervised/anchor response alignment | `test-response-aligned-mfa.R`, `aligned_models_synthetic_tests.R`; `man/response_aligned_mfa.Rd` | No confirmed defect. |

## Verification record and limits

* Focused files run: `test-aligned-interbattery.R`, `test-aligned-mcca.R`, `test-aligned-mfa.R`, `test-aligned-rrr.R`, `test-bada-multidesign.R`, `test-bamfa.R`, `test-bamfa-recon.R`, `test-bilinear-mixed.R`, `test-covstatis.R`, `test-duality.R`, `test-graph-anchored-mfa.R`, `test-graph-anchored-mfa-synthetic.R`, `test-coupled-graph-anchored-mfa.R`, `test-ipca.R`, `test-linked-mfa.R`, `test-mcca.R`, `test-penalized-mfa.R`, `test-penalized-mfa-opts.R`, `test-penalized-mfa-clusterwise.R`, `test-response-aligned-mfa.R`, and `test-infer-muscal.R`.
* The ordinary focused test files do not cover C1 small global graph size, C2 ragged rows, C6 snapshot rank-one COVSTATIS, or C7 normalization algebra. The separate synthetic specification file has two stale/sensitive failures: matched MCCA normalization and controlled-RNG/tighter-tolerance aligned-MFA runs both give exact subspace agreement. No test exercises a true BaDA bootstrap refit/alignment because the code does not perform one.
* **Source drift:** the manifest below is the original review snapshot. During audit, external edits changed `R/covstatis.R` from `caf44be...` to `314180145c74b4e628e995d1e389b213a04e4292f9bc1fd524708af295620671`, `R/covstatis_multifer.R` from `4e4e04ac...` to `a10e10e3f260aa4343976fd3ca7b769f50dedc523eb87d19b0ec034e871319dd`, and unassigned `R/plot_mfa.R` to `52d97cfa1f1ed555f247408870ef422d059ee4d37c5f0176bbb8faedfbb637db`. These current versions were not substituted for snapshot evidence.
* Method tags do not establish external scientific validity or calibrated uncertainty.  In particular, all `infer_*`, graph/anchor, linked/aligned, BaMFA, bilinear, penalized, and iPCA paths require study-specific simulation/validation before substantive claims.

## SHA-256 source manifest

```text
bbacd8c5d200db4a84b25a950343a7d3f886a0ddf843b90f055ecb0b4c9c99fc  R/aligned_interbattery.R
ed971857d2a9f718f57273f2e5fa6cad4ba3dfb1f48bbf975773b0011ac7e475  R/aligned_mcca.R
3b22efa5fdd45be7ca6b23f754670a497363e9bd5c2ea21241d5dade2b63d0a6  R/aligned_mfa.R
c778d9e5d587e07dada5b2d4d734db6346704e608b83f0e69363ec685d599c19  R/aligned_rrr.R
92d46b52b725c622717c7953fc8731a547db863e37ca87c41166636090c47a82  R/aligned_rrr_multifer.R
fea997c507b891aaef6b6ce7def0de5c58cf0507dc1d7747c918317a05b73941  R/all_generic.R
d3747727484f39930c633a919e013317ba9629204fb2edb12abf12d81713d8e8  R/bada.R
50c08282145eb41df94061d98059ebcfb3119977c11f02c3068ff7e67071adb0  R/bada_multifer.R
cf2d59d2c8eece5ceb0e424b6408ff833e7598781c4648fb17d43bb706a45b23  R/bamfa.R
0042a2cbc0baca877467280da1475335a119ad93e0c7a92fc8bbfa8a6982ffbf  R/bilinear_mixed.R
4c5fdfa2c3d74de33013765d5eb41e1ddeecfd58e07a9d62948ffa31f87e2a25  R/bilinear_mixed_tune.R
caf44be8834e8f2876d093989616e2c14471785f7d6cbc19738eff235944168b  R/covstatis.R
4e4e04acaf10181231f486d6a21fad436ed35142a0dbee997855ef131f5a9f03  R/covstatis_multifer.R
05439bb3937f0b3117ea9e94af02674546afd751833bbf3b4bc98a5436a605b7  R/graph_anchored_mfa.R
0412cd3b2e4af6a9d1722ba402a7c11b988f45deed7dc39c67f200ca99a0ee5b  R/ipca.R
1f07c0e35dc2380293e701c1eaf84e3d94d26a89ffebf18eb96f9e69a788738c  R/linked_mfa.R
833a28772f78fc6ba0e1ac8c9424bafbd023b41ec5f48d768f45f9e966f1c751  R/mcca.R
2dc46aa5190a00535eed37524962fab7431841b67ee7a1906848938588039296  R/mfa.R
a632277a48510d4661a6399923426f64a1e411b8b94df73120be08847e28e377  R/penalized_mfa.R
77339466b23913fa244313b06839b6b0ba42758f93d897e44e41162a898a45a1  R/penalized_mfa_clusterwise.R
d873deb1ed250106d4a71a2f691f27e38b90ea60c992397bcda5bc120d8a54ce  R/response_aligned_mfa.R
```

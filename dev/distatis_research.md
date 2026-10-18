# DISTATIS research and proposed integration

Research date: 2026-10-18. Checkout HEAD:
`2c5c1e83692814dc20502e6963573e57263c07c6`, with pre-existing changes.
This note records the research before implementation. The API is now implemented
in `R/distatis.R` and `R/plot_distatis.R`, with a worked example in
`vignettes/distatis.Rmd`. Proposed signatures and pre-refactor failure receipts
below describe the original research state; current validation is recorded in
`dev/distatis_validation.md`.

DISTATIS analyzes several distance matrices for the same objects. It produces
a consensus map of the objects, an assessor/table map, and partial object maps
showing how each assessor differs from the consensus. It fits muscal's existing
COVSTATIS interfaces, but a faithful implementation needs separate choices for
table normalization and between-table similarity. Wrapping the current MFA path
with distance conversion would give different results from the supplied paper.

## Sources and scope

- Abdi, Valentin, Chollet, and Chrea (2007), *Analyzing assessors and products in
  sorting tasks: DISTATIS, theory and applications*, Food Quality and Preference
  18, 627-640, [DOI](https://doi.org/10.1016/j.foodqual.2006.09.003).
  The supplied 16-page author manuscript is
  `/Users/bbuchsbaum/Downloads/Analyzing_assessors_and_products_in_sort.pdf`;
  its SHA256 is
  `5a14b7380fc974cdafa48bd999ddf13124a96e78bfb41625724b43dfe293089f`.
  Page references below refer to that PDF. Equations 1-14 define the method;
  Table 1 supplies eight beers by ten assessors, and Table 2 supplies scores.
  [Author-hosted copy](https://personal.utdallas.edu/~herve/Abdi-AVCC07_Pretty-tmp.pdf).
- The author's [DistatisR implementation](https://github.com/HerveAbdi/DistatisR/blob/master/R/distatis.R),
  and its [GetCmat](https://github.com/HerveAbdi/DistatisR/blob/master/R/GetCmat.R),
  [Dist2CP](https://github.com/HerveAbdi/DistatisR/blob/master/R/Dist2CP.R), and
  [DblCenterDist](https://github.com/HerveAbdi/DistatisR/blob/master/R/DblCenterDist.R)
  functions were inspected on the research date. The retrieved `distatis.R`
  SHA256 was `b252f75f331d3a883becc26aa01d268a182c93a9f34252aabfa75b268d2caaca`.
  These links track a mutable branch; this was source inspection, not an installed
  DistatisR parity run. DistatisR was not installed in the research environment.
- Follow-up reference for broader variants and later inference work:
  Abdi, Williams, Valentin, and Bennani-Dosse (2012), *STATIS and DISTATIS*,
  [DOI](https://doi.org/10.1002/wics.198),
  [author manuscript](https://personal.utdallas.edu/~herve/abdi_Wires_AWVB2012_Final.pdf).
  The present numerical reproduction targets the supplied 2007 paper.

## Algorithm to implement

Let there be K assessors and n objects, with aligned object order in every table.
For ordinary distances D[k], use uniform object masses and H = I - 11'/n.

1. Convert to a centered Gram matrix:
   `B[k] = -0.5 H (D[k]^2) H`, with elementwise squaring. If the caller supplies
   squared distances, omit that squaring. A centered covariance matrix is already
   in the corresponding Gram representation.
2. Default to the paper's MFA normalization:
   `S[k] = B[k] / largest_eigenvalue(B[k])`.
   Frobenius normalization and no normalization can be explicit alternatives.
3. Compute genuine RV similarities, independently of step 2:
   `C[k,l] = <S[k], S[l]>_F / (||S[k]||_F ||S[l]||_F)`.
   Compute this as a cross-product of unit-norm vectorized tables. The C diagonal
   is one because of this normalization, rather than an arbitrary diagonal edit.
4. Eigendecompose `C = P Theta P'`. Orient its leading eigenvector positively,
   and set `alpha = P[,1] / sum(P[,1])`. For PSD tables C is entrywise nonnegative.
5. Form the compromise `Splus = sum_k alpha[k] S[k]` and eigendecompose
   `Splus = V Lambda V'`.
6. Retain positive compromise roots and return object scores
   `F = V Lambda^(1/2)` and projection matrix `Q = V Lambda^(-1/2)`.
7. Return partial scores `F[k] = S[k] Q`. They satisfy
   `sum_k alpha[k] F[k] = F`. These are object-by-component matrices, not one
   coordinate per assessor.
8. Return assessor principal coordinates `G = P Theta^(1/2)` on positive roots in the separate
   interstructure space. Its rank and spectrum belong to C, independently of the
   object compromise rank. Squared column norms equal the eigenvalues of C.

Report compromise inertia as `lambda / sum(all positive compromise roots)`,
and interstructure inertia as `theta / sum(theta)`. The paper's compromise-quality
summary is `theta[1] / sum(theta)`; it describes agreement, not a p-value.

For sorting, each assessor supplies a partition. A helper can accept an
object-by-assessor data frame of group labels and compute
`D[k][i,j] = 1(group[i,k] != group[j,k])`. Group identifiers are categorical:
renaming groups must not change the result. These binary distances equal their
elementwise squares. With indicator matrix L, their Gram matrix is
`0.5 H L L' H`, hence PSD. General distances require separate PSD checks.

**Distance convention matters.** The inspected DistatisR `DblCenterDist()` computes
`-0.5 H Y H` without squaring Y. Its distance conversion therefore consumes
squared distances, despite the broad distance terminology. A muscal API should
default to ordinary R `dist` semantics and expose `squared = FALSE` explicitly.
For reference comparisons with ordinary distances, pass `D^2` to that DistatisR
path. Binary sorting examples alone cannot detect this distinction.

## What can be reused, and what needs changing

| Existing code | Integration consequence |
| --- | --- |
| `R/all_generic.R`, `covstatis()` | Use the same S3/projector conventions for `distatis()`. |
| `R/covstatis.R`, `double_center()` | Reuse centering after validating inputs; apply the distance square and -1/2 explicitly. |
| `covstatis.list()` normalization and compromise | Extract a shared Gram-fitting core with explicit table normalization and similarity choices. Preserve current COVSTATIS defaults. |
| Current `C = crossprod(table_matrix)` | This is RV only for unit Frobenius tables. Add RV normalization for DISTATIS without replacing the MFA-normalized tables used in the compromise. |
| `roi_scores`, `partial_scores`, `projmat`, `sdev`, `compromise_full` | Reuse their algebra and storage; expose object-oriented terminology for DISTATIS. |
| `rv_subject_scores()` | Currently scales interstructure eigenvectors by compromise `sdev`. DISTATIS must use `sqrt(eigC$values)` and its own axis count. |
| `check_duality()` | Equal object/assessor column norms are imposed by the current shared scaling; they are not a standard DISTATIS invariant. Test each spectrum and the partial-score barycenter instead. |
| `reconstruct.covstatis()` | Retain Gram reconstruction as the compatible default. A distance view needs a documented inverse transformation. |
| `project_cov.covstatis()`, `project_subjects.covstatis()` | Add distance preprocessing and correct RV-space supplementary projection; raw distances must never enter covariance preprocessing. |
| `project_covariate.covstatis()` | Object covariates can use object scores; assessor covariates need the new assessor coordinates. Existing cosine projections are uncentered, and should not be described as Pearson correlations. |
| `R/plot_mfa.R` COVSTATIS methods | Reuse plotting machinery but give object and assessor plots their own scores, labels, spectra, and rank checks. Partial maps and block weights transfer directly. |
| `R/covstatis_multifer.R` | Refit/null statistics currently recompute raw inner products. A future DISTATIS adapter must carry RV semantics through every refit. Inheritance alone is insufficient. |

Two additional implementation details deserve explicit tests. First,
`diag(sqrt(values[keep]))` interprets a length-one numeric value as a dimension,
which breaks valid one-component fits. Use `sweep()` or explicit diagonal sizes.
Second, current subject helpers index as many interstructure axes as object axes;
with fewer tables than retained object components this can exceed C's dimension.

Current supplementary summaries also need precise names: `distance_to_compromise`
actually computes a residual from the retained object subspace,
`||Snew - P Snew P||_F`, and its RV summary uses the truncated reconstruction.
For DISTATIS expose these as subspace diagnostics, or compute separately named
full-compromise RV and `||Snew - Splus||_F` summaries.

## Proposed public interface

These signatures are proposals, not available functions:

```r
fit <- distatis(distances, ncomp = 2, squared = FALSE,
                norm_method = "mfa")
fit <- distatis(sorting_distances(partitions), ncomp = 2)

stats::coef(fit)                  # object eigenvectors
multivarious::scores(fit)         # object scores, subject to generic integration
multivarious::reconstruct(fit)    # centered normalized compromise Gram matrix
project_distance(fit, new_dist)   # same objects, new assessor's partial map
project_subjects(fit, new_tables) # partial maps plus assessor diagnostics
plot(fit, type = "objects")
plot(fit, type = "subjects")
plot_partial_scores(fit)
plot_block_weights(fit)
```

Accept a list of `dist` objects or numeric square matrices, and an n-by-n-by-K
numeric array. An optional single-`dist` method can provide the K=1 MDS case.
Require an explicit `squared` setting for each fit; use it again for supplementary
inputs unless explicitly overridden. Do not infer it from the numeric values.

Use classes beginning with `distatis`, optionally followed by `covstatis` and the
existing projector classes. Keep compatible fields and add `object_scores`,
`subject_scores`, the complete object and interstructure spectra, preprocessing
settings, normalization constants, and diagnostics. If subclassing COVSTATIS,
override score helpers, projection, covariate dispatch and plots that assume the
old scaling. Check inherited generic dispatch rather than assuming it works.

In the shared core, use an explicit `interstructure = "inner_product"` default
for existing COVSTATIS behavior and request `"rv"` for DISTATIS. Storing this
setting makes downstream projection and inference reproducible. The first
DISTATIS API can fix RV behavior rather than exposing every STATIS variant.

For a supplementary assessor, normalize its centered Gram table as in training.
Let `r[k] = RV(Snew, S[k])`. Its principal assessor coordinates on positive
interstructure axes are `r' P Theta^(-1/2)`. A training table projected this way
must reproduce its row of G. Object partial coordinates remain `Snew Q`.
Both projections hold the fitted compromise and weights fixed. Supplementary
objects require a different out-of-sample MDS derivation and are outside v1.

If distance reconstruction is exposed, derive squared distances as
`diag(Splus) 1' + 1 diag(Splus)' - 2 Splus`. They describe the normalized
compromise geometry, not the original assessor-specific distance scales.

## Input and numerical decisions

- Require finite numeric square symmetric tables with a zero distance diagonal,
  nonnegative distances, common size, and consistent object identifiers/order.
  Reject conflicting names; do not silently reorder or fill missing distances.
  Define relative symmetry/diagonal tolerances and only repair roundoff within them.
- Use at least two objects. Allow K=1 for fitting if supported, but require more
  tables for assessor-comparison inference. Assessor and object axis limits are
  independent; a one-dimensional fit must still work.
- Default to an error for materially non-Euclidean Gram matrices. Report the
  table identifier, minimum root and negative inertia. Ignore only numerical
  negative roots within a scale-relative tolerance. Positive-part clipping or
  Cailliez/Lingoes corrections would be explicit future options, with their
  transformations recorded; they are not interchangeable with the exact method.
- Reject a zero Gram table: an all-in-one-group sorting has no contrast and makes
  RV undefined. All-singleton sorting is PSD and valid, though uninformative
  about preferred groupings; do not reject it as mathematically invalid.
- Use relative spectral thresholds and explicit component-count validation.
  Preserve full spectra for inertia denominators even when only two axes are saved.
- Detect a tied leading RV root rather than letting a solver pick arbitrary
  assessor weights. A first implementation can report an identifiability error;
  any deterministic tie rule must be specified and tested. Absolute values of a
  mixed-sign eigenvector are not a general repair.
- Avoid treating alpha as assessor accuracy: it measures agreement with this
  panel, whose size, composition and duplicated tables influence the compromise.

## Reproduction and adverse evidence

Run from the repository root:

```sh
LC_ALL=C Rscript --vanilla dev/distatis_research.R
```

The original base-R calculation uses explicit H matrices, independently of
muscal's centering helper. It asserts agreement with all 100 printed RV entries,
all ten printed weights, seven printed compromise roots, and 56 printed scores
after sign alignment. Tolerances reflect the source's two/three-decimal rounding.
It also checks the weighted partial-score identity and a non-binary Euclidean
distance example against centered coordinate cross-products. The current-checkout
comparison sources `R/covstatis.R`; it does not use an installed muscal version.

The run on R 4.5.1 passed those assertions. Full output is retained in
`dev/distatis_research_results.txt`.

| Quantity | Reproduced result |
| --- | --- |
| Leading RV roots | 5.982768, 1.404689, 0.7557803 |
| Compromise quality | 0.5982768 |
| Positive compromise roots | 0.663992, 0.490472, 0.404291, 0.338765, 0.228139, 0.130261, 0.050110 |
| Largest absolute weight difference from current COVSTATIS MFA | 0.03668511 |
| Current MFA C diagonal | Between 1.361111 and 2.640625, instead of all ones |
| Weighted partial-score barycenter error | 1.395e-15 |
| Current one-component MFA beer fit | `non-conformable arguments` |
| Current unnormalized rank-one sorting fit | `non-conformable arguments` |

The supplied paper's prose says that the first four compromise axes explain 84%
(PDF page 10). The reproduced roots give **82.28512%**, consistent with the
rounded Table 2 percentages 29 + 21 + 18 + 15 = 83 rather than 84. Use the
unrounded computed ratio, not that prose claim, as a future test target.
Other manuscript slips include a reference to ten beer samples in the procedure
despite the eight-row data table; the reproducer uses the actual table.

The run emitted an existing multivarious `prep()` deprecation warning. No package
build/check, installed DistatisR comparison, inference calibration, or performance
benchmark was performed. This is numerical research evidence, not release evidence.

## Implementation sequence and acceptance

1. Extract the shared Gram core and fix one-component scaling. Check that existing
   COVSTATIS defaults and tests keep their behavior, including its historical
   coordinate scaling unless deliberately addressed in a separate compatibility change.
2. Add DISTATIS distance conversion, validated list/array inputs, sorting helper,
   RV weights, separate spectra, projector methods and supplementary assessors.
3. Add object/assessor plots and the beer vignette, generated Rd/NAMESPACE entries,
   and pkgdown navigation. Existing unrelated checkout changes must be preserved.
4. Run focused COVSTATIS/DISTATIS tests and the repository's package checks.
   Validate installed DistatisR parity on a pinned version and non-binary distances
   with the appropriate squared-input convention.

Acceptance tests should cover the printed beer targets, direct Gram-versus-distance
equivalence, single-table classical MDS, alpha nonnegativity and sum-to-one,
weighted partial barycenters, supplementary training-table projections, full-rank
distance reconstruction, and independent object/assessor rank limits. Add object
and assessor permutation tests, within-table rescaling invariance under MFA or
Frobenius normalization, sorting-label renaming, one-component/rank-one cases,
degenerate leading roots, and invalid/non-Euclidean inputs. Compare subspaces or
Gram matrices when eigenvector signs or tied retained roots make raw coordinates
non-unique. Scale invariance is not expected with `norm_method = "none"`.

Inference should follow as a separately validated extension. The supplied paper's
partial-map hulls describe assessor dispersion; they are not confidence regions.
An independent assessor bootstrap with refitting and axis alignment would address
one uncertainty target, while repeated assessors require sampling at the assessor
cluster level. Permutation tests require an explicit null, such as independently
permuting object identities within each table on both rows and columns. A common
permutation of all tables leaves the fit invariant and is not that null. Reusing
the existing multifer adapter also requires deciding how residualization interacts
with RV renormalization; its current raw-inner-product refits do not qualify the
DISTATIS method. If DISTATIS inherits from COVSTATIS, guard `infer_covstatis()`
against accepting raw distance data for a DISTATIS object until that adapter exists.
Defer bootstrap regions and component-significance claims until
those contracts are specified and checked.

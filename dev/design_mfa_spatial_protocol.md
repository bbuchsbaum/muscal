# Design-MFA Spatial Recovery and Inference Protocol

## Status

This document defines the confirmatory protocol. A standard-mode result is
valid only when this file is tracked, unchanged, and committed before the run.
Quick-mode output is developmental smoke evidence and must never be relabelled
as confirmatory.

The earlier `design_mfa_spatial_calibration*.csv` files are retrospective,
exploratory artifacts. Their thresholds did not predate their data.

## Method under test

The method is the `design_mfa()` ASCA+/partial-RDA-style design decomposition.
It does not perform MFA block normalization. Spatial fits use explicit sparse
reference-to-subject maps. Responses are generated directly from continuous
coordinate fields and never through a fitted map.

## Fixed simulation sizes and seeds

Standard mode uses:

- 60 independent recovery replicates per non-null scenario;
- 100 independent null datasets and 100 planted-effect datasets;
- 199 permutations per inferential dataset;
- eight subjects, 36--50 trials per subject, and 30 reference features;
- recovery seeds `202609010 + replicate + 10000 * scenario_id`;
- inference seeds `202619010 + replicate + 10000 * null_indicator`.

Quick mode uses six recovery replicates, eight null datasets, eight
planted-effect datasets, and 19 permutations. Quick thresholds are diagnostic
only.

## Recovery estimands

For the saliency effect:

- weighted relative squared error of the complete centered effect matrix;
- absolute correlation between the first loading and the planted reference
  field;
- maximum absolute shared-versus-soft coefficient and effect differences in
  the identity-map reduction;
- coverage fractions and expected failure under deliberately incomplete
  reference support.

Loading signs are indeterminate, hence absolute correlations. Component
subspaces, rather than individual signed vectors, are primary when roots are
tied.

## Prespecified scenarios and comparators

1. `identity`: all subjects use the exact reference coordinates and identity
   maps. The shared and soft engines must agree.
2. `preserved_index`: each reference location is independently jittered for
   each subject while its index remains meaningful. Same-index shared fitting
   is the negative-control comparator and is expected to beat interpolation.
3. `independent_coordinates`: subject coordinates are independently sampled
   from the common spatial domain, so same-index identity is absent. Gaussian
   interpolation with `k = 3`, bandwidth `0.35`, is compared with a
   prespecified nearest-neighbor map.
4. `misspecified_map`: the independent-coordinate data are fitted with a
   deliberately left-right-reflected reference geometry and compared with the
   prespecified map.
5. `incomplete_coverage`: subject coordinates omit a reference region. The
   current fail-closed contract must reject unsupported reference features and
   report the absence of support; no recovery score is imputed.

No map parameter is selected from the assessment data. Any future data-driven
selection must occur entirely within training subjects and requires a new
protocol version.

## Recovery gates

- Identity-map shared and soft predictions, effect matrices, and coefficient
  matrices differ by at most `1e-7` in maximum absolute value.
- In `preserved_index`, mean same-index relative error is lower than mean soft
  relative error. This is a negative-control gate, not evidence against soft
  mapping where identity is absent.
- In `independent_coordinates`, mean Gaussian-map error is at most `1.25`
  times mean nearest-neighbor error and median absolute loading correlation is
  at least `0.90`.
- Mean error under the intentionally misspecified map is no lower than under
  the prespecified Gaussian map.
- Every incomplete-coverage replicate is rejected before fitting and the
  error identifies unsupported reference features.

## Exchangeability and inferential estimands

The only calibrated scheme in this protocol permutes randomized saliency
labels independently within subject and within probe stratum. It assumes trial
errors are exchangeable within those subject-by-probe strata under the
saliency null. Probe labels, responses, weights, subject membership, reference
grid, and reference weights remain fixed.

For each dataset the script records:

- the plus-one permutation p-value for term-level partial SS;
- component-root p-values reported by `infer_design_mfa()`;
- a max-root adjusted p-value for every retained component, using the maximum
  null singular value in each replicate;
- rejection at alpha `0.05` for the term and for the component family.

## Inferential gates

- Under the null, the 95% Wilson interval for the term-level rejection rate
  contains `0.05` and has upper endpoint at most `0.12`.
- Under the null, the 95% Wilson interval for familywise max-root rejection
  contains `0.05` and has upper endpoint at most `0.12`.
- Under the planted effect, term-level power is at least `0.80`.
- Under the planted effect, max-root familywise power is at least `0.70`.

Descriptive null effect energy is not an inferential calibration statistic.

## Provenance requirements

Every standard result row must include:

- the Git commit containing this protocol;
- SHA-256 hashes of this protocol, the calibration script, the simulation
  helper, and all `R/design_mfa*.R` implementation files;
- R, Matrix, BLAS, platform, and operating-system versions;
- the exact replicate seed, scenario, and permutation count.

The script must refuse standard mode if this protocol is untracked or differs
from its committed version. Results are written to new confirmatory filenames;
legacy exploratory CSVs are never overwritten.

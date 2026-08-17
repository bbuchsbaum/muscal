# `design_mfa()` Implementation Plan

## Status and objective

Implementation status (2026-08-15): the core method and its high-priority
extensions are implemented. The shared fixed-intercept engine accumulates
subject-streamed, feature-chunked sufficient statistics; the formula compiler,
estimability checks, equal-subject weighting, effect-specific thin SVD, query
API, and analytic/metamorphic tests cover Milestones 0--3. Validated extensions
include fixed-penalty centered subject slopes, whole-subject bootstrap, sparse
spatial metrics, coordinate-derived sparse soft-alignment maps, preconditioned
matrix-free fitting, and subject-disjoint held-out and nested tuning of rank and
all three penalties.

A post-implementation review hardening pass is also complete. Permutation
inference now has method-valid defaults; seeded inference and tuning preserve
the caller's RNG state; subject-keyed tuning inputs follow `Y` block identity
rather than an unrelated subject-label vector; resampling freezes the fitted
reference grid and weights; rank-truncated component fractions use total
effect energy; diagnostics distinguish term-contrast estimability from stricter
full-response predictability; and weak-pooling fits have a 5,000-iteration
default budget. Nonconverged fits still fail closed, while tuning retains and
reports the candidate/fold errors instead of replacing them with a generic
failure.

Current evidence, with the limits exposed, includes:

- exact streamed/dense coefficient and residual-SSE agreement. The existing
  20-subject, 100,000-feature shared-grid run took 22.33 seconds, but its peak R
  heap was 1.47 GB and its fitted object was 890 MB; the 20.8 MB figure is only
  the largest single profiled allocation. That run exercised neither the soft
  engine nor multi-subject coordinate maps at imaging scale
  (`dev/design_mfa_benchmark_results_imaging.csv`);
- a 20-subject, 4,000-feature benchmark covering shared, spatial, and soft
  engines in 1.48, 1.62, and 8.16 seconds on one run, with the soft solver
  converging in 28 iterations. These timings are descriptive, not a stable
  performance envelope (`dev/design_mfa_benchmark_results_standard.csv`);
- a replacement five-fresh-process 20-subject, 4,000-feature envelope with
  median times of 1.136, 1.236, and 4.690 seconds for shared, spatial, and soft
  engines. The soft engine used 27 operator applications in every replicate
  (`dev/design_mfa_benchmark_standard_repeated_summary.csv`);
- five fresh-process four-subject, 100,000-reference-feature runs covering all
  three engines and coordinate maps for every subject. Median soft time was
  28.762 seconds, median peak R heap was 983 MB, median fit size was 319 MB,
  and the largest allocation was 21.6 MB. These are descriptive measurements,
  not extrapolations (`dev/design_mfa_benchmark_imaging_soft_repeated_summary.csv`);
- a legacy 60-replicate spatial calibration that is invalid as comparative
  evidence: it permuted feature identity only for the same-index baseline and
  generated observations through the exact map given to the soft model. Its
  null-energy comparison measured smoothing, not type-I error. The files
  `dev/design_mfa_spatial_calibration*.csv` are retained only as superseded
  exploratory artifacts;
- a corrected 60-replicate, map-independent calibration with preserved feature
  identity. Mean relative error is 0.0120 for the same-index fit and 0.1314 for
  the soft fit (10.9 times larger), while the soft loading correlation remains
  0.9974. This supports soft-field recovery but reverses the earlier
  comparative claim (`dev/design_mfa_spatial_calibration_fair_summary.csv`);
- a confirmatory run against protocol commit `ab40f18`, comprising 60 recovery
  replicates per scenario, 100 null datasets, 100 planted-effect datasets, and
  199 permutations per inferential dataset. All 12 preregistered gates passed:
  term-level null rejection was 0.01 (95% Wilson interval 0.0018--0.0545),
  max-root familywise rejection was 0.04 (0.0157--0.0984), both power measures
  were 1.00, median independent-coordinate loading recovery was 0.9982, and
  the misspecified-map error ratio was 57.56. Raw rows, gates, provenance, and
  verified result hashes are in `dev/design_mfa_spatial_confirmatory_*.csv`
  and `dev/design_mfa_spatial_confirmatory_sha256.txt`;
- deterministic k-nearest coordinate-map tests against dense oracles and a
  100,000-point sparse construction benchmark; and
- a clean package build of the runnable saliency--probe--confidence vignette.

All nine Design-MFA test files pass, with 416 assertions. A fresh source-tree
full suite recorded 2,030 passes, five warnings, two skips, and one failure:
the pre-existing stochastic
`penalized_mfa_clusterwise()` objective-monotonicity test. That fix remains
scoped to `dev/penalized_mfa_clusterwise_recovery_plan.md`, not this method.
A fresh tarball check passed installation, code analysis, documentation,
examples including `--run-donttest`, the Design-MFA vignette, and vignette
rebuilding. Its installed suite recorded 2,026 passes, five warnings, three
skips, and only the same separately scoped clusterwise failure. Overall check
status was one error, two warnings, and zero notes: the test error, the existing
CRAN-incoming non-mainstream dependency warning, and the host
R-header/Homebrew-clang warning.

The remaining release gates are deliberately narrower: a variance-component
or REML estimator (the current subject deviations are explicitly fixed-penalty
ridge effects); canned permutation schemes and their scheme-specific type-I
calibration; broader comparison simulations against `aligned_rrr()`, BaDA, and
unsupervised MFA; and independent statistical and numerical review. The
current API does not imply those unfinished claims.

This plan introduces a new effect-oriented multilevel method for trialwise
multivariate fields. Its primary use case is repeated-measures neuroimaging in
which every trial has a beta volume and trial metadata describe ordered,
nominal, and response-related factors.

The method should answer:

> Which distributed spatial components capture the multivariate signal
> associated with each estimable design term at the group level?

The statistical core is multilevel partial reduced-rank decomposition: fit the
trial-level design first, then componentize each fitted design effect. This is
an ASCA+/partial-RDA-style decomposition with explicit subject structure. It is
not an unsupervised PCA whose components are interpreted after the fact.
Despite the exported `design_mfa()` API name, this method does not perform MFA
block normalization and should not be described as MFA in scientific claims.

The exact engine supports a shared feature grid. The spatial engine accepts
explicit sparse maps when subject fields merely correspond, and
`coordinate_feature_map()` / `coordinate_feature_maps()` construct those maps
deterministically from common-space coordinates without all-pairs distances.

## Relationship to existing muscal methods

- `aligned_rrr()` is the strongest supervised shared-response comparator and
  should be included in simulations, but it does not produce term-specific
  factorial decompositions.
- `response_aligned_mfa()` supplies useful weighting, naming, prediction, and
  fit-contract patterns, but its reconstruction objective is not the new
  method's statistical core.
- BaDA is a class-barycenter method and must not be used to fill or average
  sparse factorial cells.
- BAMFA is an unsupervised global/local decomposition and assumes a different
  row relationship.
- `penalized_mfa_clusterwise()` is not a backend. The exact-grid engine is
  independent of it; only validated sparse graph/map primitives may later be
  shared with the soft-alignment engine.

## Scientific contract

For subject `s`, let

```math
Y_s \in \mathbb{R}^{n_s \times p_s}
```

contain trialwise activation fields. Let `D_s` be a fixed-effect design matrix,
`Z_s` a subject/run design, and `W_s` trialwise precision weights. The model is

```math
Y_s = D_s B + Z_s U_s + E_s.
```

The first implementation uses subject intercepts as fixed nuisance effects and
subject-level resampling for population inference. The multilevel extension
adds penalized subject deviations or random slopes:

```math
B_{s,h} = B_{0,h} + U_{s,h}, \qquad
U_{s,h} \sim (0, \tau_h^2).
```

For term `h`, define a reduced model containing all terms except `h` that are
required by the chosen marginality rule. Residualize both the response and the
term design against that reduced model:

```math
\widetilde Y_h = R_{-h}Y, \qquad
\widetilde D_h = R_{-h}D_h.
```

Estimate the conditional term coefficient and evaluate it on an explicit
reference grid `G_h`:

```math
H_h = G_h \widehat B_h.
```

The effect components are obtained from

```math
H_h = U_h \Sigma_h V_h^\top.
```

The columns of `V_h` are group spatial loading fields; `U_h Sigma_h` describes
the effect trajectory across the reference design grid. Components maximize
fitted variance for term `h`, not total variance in `Y`.

### Interpretation rules

- Partial effects are conditional on the other specified terms.
- Partial sums of squares can overlap in an unbalanced design and are not
  presented as an additive variance pie.
- Sequential sums of squares are available only behind an explicit
  order-dependent option.
- Missing or rank-deficient cells produce non-estimable contrasts, not imputed
  cell means.
- Ordered factors use planned contrasts. For the motivating example, saliency
  has a primary linear trend and optional nonlinear polynomial trends.
- Confidence is an observed response category. Its components are labeled
  conditional associations, not causal experimental effects.
- Interactions obey model hierarchy unless the user explicitly overrides it.

### Term-hypothesis policy

The default partial test is hierarchical: it measures the extra fitted
variation contributed jointly by term `h` and all descendant interactions
containing it. This hypothesis is invariant to full-rank factor contrast
coding. `term_scope = "assigned"` is an explicit, coding-dependent option that
removes only the model-matrix columns assigned to `h` while higher-order
interactions remain in the conditioning set.

Reported effect trajectories are evaluated as contrasts on the reference grid,
not read directly from a treatment-coded coefficient. Reference-grid weighting
defaults to equal weighting of scientifically defined cells; observed-frequency
and custom weighting are explicit alternatives. The fit records both the term
scope and reference weights so that phrases such as "main effect" have an
auditable meaning.

## Proposed API

The generic should work with lists of subject matrices and aligned design
tables first. A `multidesign` method can adapt existing package structures.

```r
design_mfa(
  Y,
  design,
  formula,
  subject = NULL,
  nuisance = NULL,
  contrasts = NULL,
  terms = NULL,
  reference_grid = NULL,
  alignment = c("shared", "soft"),
  coords = NULL,
  reference_coords = NULL,
  feature_map = NULL,
  trial_weights = NULL,
  block_weight = "equal_subject",
  subject_model = c("fixed_intercept", "ridge_slopes"),
  rank = NULL,
  spatial_penalty = 0,
  alignment_penalty = 0,
  ss_type = c("partial", "sequential"),
  term_scope = NULL,
  reference_weight = c("equal_cells", "observed", "custom"),
  store = c("sufficient", "fitted", "residuals"),
  control = design_mfa_control(),
  ...
)
```

Example:

```r
fit <- design_mfa(
  Y = beta_volumes,
  design = trial_tables,
  formula = ~ saliency * probe * confidence,
  nuisance = ~ run,
  contrasts = list(
    saliency = "poly",
    probe = "sum",
    confidence = "poly"
  ),
  alignment = "shared",
  subject_model = "fixed_intercept"
)
```

The formula compiler must use ordinary R terms and contrast matrices rather
than special-case these factor names.

## Input contract

### Shared alignment

- `Y[[s]]` is `n_s x p` for every subject.
- Feature identity is established by identical, unique column names or an
  explicitly supplied feature schema. Equal column counts alone are
  insufficient.
- Each `design[[s]]` has exactly `n_s` rows.
- Subject identity is unique and stable.

### Soft alignment

- `Y[[s]]` may be `n_s x p_s`.
- `coords[[s]]` has `p_s` finite coordinate rows.
- `feature_map[[s]]` is a validated sparse `p_s x p_0` operator mapping a
  reference-grid field into subject `s`'s feature space.
- Maps record coverage, row sums, extrapolation, distances, and unmatched
  features.
- Automatically generated maps are deterministic and never use dense all-pairs
  distances.

### Design and weights

- Factors retain their declared level order; character columns are not silently
  alphabetized into scientific contrasts.
- The design compiler records model-matrix columns, term assignments, contrast
  matrices, rank, pivots, aliases, and estimability.
- Trial weights are finite and nonnegative. Full trial covariance is a later
  extension and must not be approximated silently by diagonal weights.
- Subject weighting defaults to equal population contribution, independently
  of the number of trials per subject.
- Under equal-subject weighting, precision weights are normalized within each
  subject and each subject receives total weight `1 / S`. Precision-weighted
  population estimation is a separate explicit option.
- Non-finite response values require an explicit missing-data policy. The first
  release rejects them rather than performing voxel-specific implicit row
  deletion.

## Fit object and query contract

`design_mfa` should be its own top-level S3 object. It should not be forced into
a single `multiblock_projector`, because every design term has a distinct
loading basis.

The object contains:

- parsed formula, contrast specification, term hierarchy, and reference grid;
- feature schema or spatial maps;
- subject, trial, and block weights;
- design rank and an estimability report;
- nuisance and subject-effect specifications;
- sufficient statistics needed for refit or diagnostics;
- one `design_mfa_effect` object per term;
- fit, convergence, and performance diagnostics;
- a `muscal_fit_contract` whose resampling unit is the subject.

Each effect object contains:

- term name and model-matrix columns;
- estimable degrees of freedom;
- group coefficient field and optional subject deviations;
- reference-grid effect matrix;
- singular values, design scores, spatial loadings, and explained fitted-effect
  variance;
- partial sums of squares and partial `R2` with their precise denominator;
- component orientation and stability metadata.

Required accessors:

```r
effects(fit)
components(fit, term)
effect_scores(fit, term, components = NULL)
effect_loadings(fit, term, components = NULL)
effect_map(fit, term, component = 1)
estimability(fit)
subject_effects(fit, term = NULL)
```

Required methods:

- `print.design_mfa()` gives design rank, subjects, features, terms, alignment,
  and warnings.
- `summary.design_mfa()` reports estimability, partial effect sizes, ranks, and
  stability.
- `plot.design_mfa()` defaults to effect trajectories, never a generic score
  scatterplot.
- `predict.design_mfa()` evaluates fitted design effects on new design rows; it
  is not a trial-class prediction method.
- `project.design_mfa_effect()` projects new fields onto one term's spatial
  basis after applying the stored feature/preprocessing contract.

## Computational architecture

### Exact shared-grid engine

Never form a `p x p` response covariance. With total design rank `q`, accumulate
feature-chunked sufficient statistics such as

```math
D_s^\top W_sD_s, \qquad D_s^\top W_sY_s.
```

This requires approximately `O(p q)` working memory and can stream subjects or
feature chunks. Effect rank is bounded by the term's estimable design rank, so
componentization uses a thin QR or a small Gram matrix.

The engine should have two equivalent paths:

- a transparent dense reference solver for small problems;
- an operator/chunked solver for imaging-scale data.

They must be differential-tested against one another.

### Multilevel engine

Implement in increasing complexity:

1. Subject fixed intercepts plus subject-clustered inference.
2. Subject-specific coefficient estimates combined by equal-subject or
   precision-weighted group estimation.
3. Ridge-partial-pooled subject deviations with an explicit penalty.
4. Optional empirical-Bayes or REML estimation only after its variance
   component behavior has its own oracle and coverage simulations.

Do not call a penalized fixed effect a random effect unless the variance model
and inferential target are explicit.

### Spatial generalized decomposition

For shared coordinates, allow a reference-grid penalty `L_0` and compute
effect components under `Q = I + lambda L_0`. Use sparse solves for the few
term directions rather than a voxel covariance eigenproblem.

### Soft-coordinate engine

Align coefficient fields before componentization. With a sparse map `P_s` from
reference to subject space, estimate

```math
B_{s,h} = B_{0,h}P_s^\top + U_{s,h}
```

by minimizing weighted trial reconstruction plus within-field roughness and
subject-deviation penalties. This is a convex quadratic in coefficient fields;
the effect SVD occurs afterward. It avoids trying to align non-identifiable
component columns during a nonconvex fit.

Use matrix-free conjugate gradients or another validated sparse linear solver.
Report map coverage and refuse a group map when reference support is
insufficient.

## Milestone 0: Specification and simulation harness

### Work

- Freeze terminology: design effect, conditional association, component,
  partial effect size, subject deviation, shared grid, and soft alignment.
- Implement a deterministic simulator that controls:
  - subjects and unequal trial counts;
  - ordered and nominal main effects and interactions;
  - absent cells and aliased contrasts;
  - subject intercepts and slopes;
  - smooth distributed spatial patterns;
  - noise covariance and trial precision;
  - exact and misregistered coordinate fields.
- Write a slow, small-problem reference implementation using explicit weighted
  least squares and base SVD.

### Acceptance gate

- The simulator returns the exact generating coefficient and effect fields.
- The reference implementation recovers noiseless estimable effects to
  numerical tolerance.
- Missing cells are classified correctly as estimable or non-estimable.

## Milestone 1: Design compiler and estimability

### Work

- Parse formulae, ordered contrasts, nuisance terms, subject identifiers, and
  term hierarchy.
- Preserve declared factor levels across subjects, including levels absent from
  one subject.
- Build full and reduced term-specific design operators.
- Use pivoted QR/SVD to determine numerical rank and aliased columns.
- Construct a balanced or user-provided reference grid independently of
  observed cell counts.
- Test contrast estimability before producing any effect matrix.
- Record whether a term is experimental, observational, or unspecified so
  inference cannot select an unsafe permutation automatically.

### Required tests

- Balanced and unbalanced factorial designs.
- Empty combinations, subject-specific missing levels, and collinear nuisance
  terms.
- Equivalent contrast bases yield the same fitted effect subspace.
- Row permutation and subject-list permutation equivariance.
- Ordered levels remain in declared scientific order.
- Full-model and partial-residualization calculations agree.

### Acceptance gate

- Every reported effect is demonstrably estimable.
- Reparameterizing a term without changing its column space leaves fitted
  values and effect components unchanged up to subspace tolerance.

## Milestone 2: Shared-grid fixed-effect engine

### Work

- Validate feature identity by schema, not dimensions alone.
- Implement weighted full and reduced model solves.
- Add subject intercepts and optional run/session nuisance terms.
- Implement sufficient-statistic and feature-chunked computation.
- Produce group coefficient and reference-grid effect matrices.
- Calculate partial sums of squares, residual sums of squares, and explicitly
  defined partial `R2`.
- Flag rank loss and ill conditioning with actionable diagnostics.

### Required tests

- Differential comparison with `lm.fit()` performed independently per feature
  on small fixtures.
- Weighted and unweighted analytic solutions.
- Dense, streamed-by-subject, and feature-chunked equivalence.
- Exact recovery of noiseless planted effects.
- Stability under extreme but finite weights and near-collinear designs.
- Equal-subject weighting is invariant to duplicating trials while adjusting
  their weights to preserve the subject contribution.

### Acceptance gate

- Coefficients, fitted values, residuals, and partial statistics agree with the
  reference solver under scale- and conditioning-aware tolerances.
- Imaging-scale paths never allocate a `p x p` matrix.

## Milestone 3: Effect componentization and user interface

### Work

- Implement thin SVD and small-Gram componentization of every effect matrix.
- Add deterministic orientation and explicit subspace handling for tied roots.
- Allow per-term rank selection bounded by estimable degrees of freedom.
- Implement fit objects, accessors, printing, summaries, effect trajectories,
  component maps, and subject overlays.
- Provide adapters from vectorized MNI fields back to image or surface metadata
  without making muscal depend on one neuroimaging container.
- Add a `muscal_fit_contract` and exact refit specification.

### Required tests

- Thin and dense SVD equivalence.
- Reconstruction identity for retained and full effect rank.
- Sign and order determinism for separated roots.
- Projector/subspace comparisons for tied roots.
- Query methods return the same values stored in the effect object.
- Reference-grid trajectories reproduce the fitted term matrix.
- Save/reload and raw-field projection contracts.

### Acceptance gate

- Every plotted quantity can be traced to a named fitted matrix.
- A complete motivating saliency × probe × confidence example can be queried
  without extracting undocumented slots.

## Milestone 4: Multilevel partial pooling

### Work

- Estimate subject-specific term coefficients and their sampling precision.
- Add group means and ridge-partial-pooled subject deviations.
- Support term-specific deviation penalties and optional random slopes.
- Separate population effect maps, subject deviations, and residual trial
  variation in the object model.
- Add a variance-component estimator only after the fixed-penalty model is
  validated.

### Required tests

- Dense closed-form ridge/BLUP oracle on tiny problems.
- Limits: infinite penalty approaches the common-effect model; zero penalty
  approaches independent subject coefficients when estimable.
- Subject permutation equivariance.
- Recovery across known between-subject variance regimes.
- Coverage and shrinkage simulations under unequal trials and noise.
- Honest failure when a requested random slope is not estimable.

### Acceptance gate

- Group coefficients are unbiased within Monte Carlo error under the canonical
  simulator.
- Subject deviations show the expected shrinkage limits and do not change when
  subject storage order changes.

## Milestone 5: Inference and tuning

### Work

- Bootstrap subjects as whole clusters, preserving all their trials.
- Align bootstrap components by signed permutation or Procrustes rotation and
  report subspace stability for tied components.
- Implement restricted/Freedman-Lane permutations only from an explicit
  exchangeability specification.
- Keep confidence-related inference conditional on subject and the other
  design terms; never present it as randomized-factor inference.
- Select rank plus spatial, partial-pooling, and soft-alignment penalties with
  nested held-out-subject validation.
- Define familywise term/component multiplicity policy before adding
  significance stars to plots.

### Required tests

- Seed reproducibility and subject-cluster preservation.
- Null type-I error simulations for every supported permutation scheme.
- Power simulations with planted effects.
- No trial leakage across training and validation subjects.
- Bootstrap orientation and subspace-stability oracles.
- Explicit refusal when no valid exchangeability scheme is supplied.

### Acceptance gate

- Null calibration and interval coverage meet predeclared Monte Carlo bounds.
- Tuning choices are reproducible and based entirely on held-out subjects.

## Milestone 6: Spatial penalty and soft alignment

### Work

- Reuse only the sparse graph/map primitives that have passed the clusterwise
  recovery plan's graph gates.
- Implement reference-grid smoothing in the shared engine.
- Implement sparse reference-to-subject interpolation maps and coverage
  diagnostics.
- Solve the soft coefficient-field objective with a matrix-free method.
- Compare exact alignment (`P_s = I`), softened identity alignment, and
  deliberately misregistered cases.

### Required tests

- `P_s = I` reduces to the shared-grid engine.
- Coordinate and feature permutations preserve restored results.
- Translation/rotation tests behave according to the declared distance model.
- Sparse solver agrees with a dense tiny-problem oracle.
- Objective and linear-system residuals meet tolerance.
- Increasing smoothing reduces graph roughness without silently changing the
  design estimand.
- Uncovered reference features are reported and excluded from interpretation.

### Retrospective corrected spatial validation

The corrected calibration used fields generated directly from independently
jittered coordinates. Feature order is preserved so the same-index comparator
retains the correspondence it is supposed to represent; the map handed to the
soft solver is never used to generate the responses. The following thresholds
were written down only after the exploratory results existed. They summarize
that run but are not prospective gates and cannot support confirmatory claims:

- exact and soft engines agree for identity maps;
- no graph or mapping path allocates all-pairs feature distances;
- mean same-index relative effect error is below 0.03 in the preserved-index
  negative control;
- mean same-index error is below mean soft-map error in that negative control,
  preventing a smoothing method from being advertised as superior when exact
  feature identity is already known;
- mean soft-map relative error remains below 0.20 and median absolute loading
  correlation exceeds 0.95, establishing recovery but not comparative benefit.

Null effect energy is descriptive only and is not a release gate. The
confirmatory protocol in `dev/design_mfa_spatial_protocol.md` must be committed
before its data are generated. A separate calibration must report rejection
rates from a scientifically justified, predeclared permutation scheme before
any type-I-error claim is made.

### Confirmatory spatial validation

The protocol was committed as `ab40f18` before the standard run. The runner
then completed 60 recovery replicates per scenario, 100 null datasets, 100
planted-effect datasets, and 199 permutations per inferential dataset. All 12
gates passed. Term-level null rejection was 0.01 with a 95% Wilson interval of
0.0018--0.0545; max-root familywise null rejection was 0.04 with interval
0.0157--0.0984; term and max-root power were both 1.00. Identity reduction,
preserved-index negative control, independent-coordinate recovery,
misspecification, and incomplete-support gates also passed. Each row records
the protocol commit, protocol/script/helper/implementation hashes, software
and BLAS versions, scenario, seed, and permutation count. A checked SHA-256
manifest covers the three result files.

## Milestone 7: Performance, documentation, and release

### Performance suite

Measure wall time, peak memory, feature throughput, and numerical agreement for:

- dense versus streamed exact fits;
- effect ranks from one-degree-of-freedom trends through moderate interactions;
- 20--40 subjects with unequal trials;
- ROI-scale and 100,000-feature shared-grid fields;
- soft maps with varying sparsity and coverage;
- bootstrap and held-out-subject refits.

Provisional performance gates, to be calibrated across repeated runs on a
named reference runner before enforcement:

- Shared-grid fitting has memory proportional to `p q`, never `p^2`.
- A 100,000-feature effect componentization uses the small term rank and does
  not call full base `svd()` on the trial-by-feature matrix.
- Streamed and dense engines agree numerically and the streamed engine remains
  within 1.5x the dense runtime when data already fit in memory.
- After baselining, PR benchmarks detect robust regressions; larger imaging and
  inference workloads run nightly or pre-release.

### Documentation

- Function reference with the precise estimand and non-additivity warning.
- A runnable saliency × probe × confidence vignette with unequal trial counts.
- A missing-cell/estimability example.
- A shared-grid versus soft-alignment example.
- An inference section distinguishing randomized effects from confidence
  associations.
- A comparison vignette against `aligned_rrr()`, BaDA, and unsupervised MFA.

### Release gate

- Full tests and `R CMD check` pass.
- All examples and vignettes run from a clean installation.
- Simulations establish recovery, null calibration, and coverage.
- Performance envelopes are recorded with hardware and dependency versions.
- An independent statistical review confirms the estimand, missing-cell
  behavior, and inference labels.
- An independent numerical review confirms the solver, invariants, and scaling.

## Suggested source and test layout

Source files:

- `R/design_mfa.R`: generic, methods, validation, and top-level fit.
- `R/design_mfa_design.R`: formula compilation, term hierarchy, contrasts, and
  estimability.
- `R/design_mfa_exact.R`: dense and streamed exact-grid solvers.
- `R/design_mfa_effects.R`: effect matrices, componentization, accessors, and
  orientation.
- `R/design_mfa_multilevel.R`: subject deviations and partial pooling.
- `R/design_mfa_spatial.R`: reference graphs, maps, and soft-coordinate solver.
- `R/design_mfa_inference.R`: subject bootstrap, restricted permutation, and
  tuning adapters.
- `R/plot_design_mfa.R`: effect trajectories, maps, and diagnostics.

Tests and development tools:

- `tests/testthat/helper-design-mfa-sim.R`
- `tests/testthat/test-design-mfa-design.R`
- `tests/testthat/test-design-mfa-exact.R`
- `tests/testthat/test-design-mfa-effects.R`
- `tests/testthat/test-design-mfa-multilevel.R`
- `tests/testthat/test-design-mfa-inference.R`
- `tests/testthat/test-design-mfa-spatial.R`
- `tests/testthat/test-design-mfa-api.R`
- `dev/benchmark_design_mfa.R`
- `vignettes/design_mfa.Rmd`

## Dependency and execution order

The critical path is:

```text
specification/simulator
  -> design compiler and estimability
  -> shared-grid fixed-effect engine
  -> effect componentization and query API
  -> multilevel partial pooling
  -> inference and tuning
```

The spatial branch is parallel after the shared engine:

```text
validated sparse graph/map primitives
  -> shared-grid spatial penalty
  -> soft coefficient-field alignment
  -> spatial simulations and imaging-scale benchmarks
```

The first useful release slice is Milestones 0--3: a correct, efficient,
shared-grid group design decomposition with subject intercepts and honest
estimability. Multilevel shrinkage, inference, and soft alignment should be
added only behind their respective validation gates.

## Completion definition

`design_mfa()` is complete when it:

1. handles unequal trial counts without cell averaging;
2. reports only estimable partial design effects;
3. yields effect-specific, group-interpretable spatial components;
4. treats subjects as the population resampling unit;
5. distinguishes randomized effects from observational confidence
   associations;
6. supports exact shared coordinates efficiently;
7. supports sparse soft correspondence without dense distances;
8. passes analytic, differential, metamorphic, simulation, inference, and
   performance gates; and
9. exposes all important scientific results through stable documented queries
   and visualizations.

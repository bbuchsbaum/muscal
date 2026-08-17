# Design-MFA Pre-Merge Remediation Plan

## Objective

Bring the trialwise design decomposition to a mergeable state by fixing every
confirmed correctness, diagnostics, repository-hygiene, test, documentation,
performance, and validation-evidence issue identified in the independent
reviews.

This plan treats the method scientifically as a multilevel
ASCA+/partial-RDA-style design decomposition. The `design_mfa()` API name may
remain as a package label, but no documentation or result may imply MFA block
normalization.

## Current disposition

The following core defects are already repaired in the working tree and must
remain protected by regression tests:

- direct residual SSE and subject-centered response energy;
- varying numeric reference grids;
- preservation of training `predvars` for `poly()`, `scale()`, and splines;
- hierarchical partial hypotheses by default;
- subject-keyed fold subsetting;
- iterative-solver-scaled model-comparison tolerances;
- method-valid permutation defaults, frozen bootstrap grids, RNG restoration,
  rank-truncation fractions, and expanded tuning diagnostics.

The work below closes the remaining gaps. No existing calibration artifact is
to be described as preregistered or prospective unless a committed protocol
demonstrably predates the run.

## Execution status (2026-08-16)

All Design-MFA work in P0--P3 is implemented and protected by tests:

- P0: warnings are independent of verbosity; subject-keyed tuning inputs fail
  early and canonicalize by `names(Y)`; prediction estimability is documented;
  package tarballs are ignored and removed from the index; effect functionals
  are exactly estimable and invariant to nuisance scaling; the unstable SSE
  subtraction fallback is gone.
- P1: bootstrap weights, Euclidean and spatial explained energy, all RNG exits,
  held-out and nested leakage, detailed solver failures, stable tuning failure
  categories, default plotting, every public example, and terminology are
  covered by discriminating regressions or independent oracles.
- P2: the soft engine applies maps matrix-free without storing `P_s^T P_s`,
  reuses map caches, honors feature chunks, and records one operator application
  per CG iteration. Five-fresh-process standard and four-subject,
  100,000-reference-feature envelopes are recorded.
- P3: legacy calibration is explicitly retrospective; protocol commit
  `ab40f18` predates the standard run; all 12 recovery, negative-control,
  misspecification, support, null-calibration, and power gates passed. The raw
  outputs have a verified SHA-256 manifest.

Verification currently records 416 passing Design-MFA assertions and 2,030
passes in the source suite. The sole source-suite failure is the independently
confirmed pre-existing `penalized_mfa_clusterwise()` objective-monotonicity
test, which remains scoped to its separate recovery plan. A temporary tarball
build passes with vignettes; installed-package documentation, examples,
`--run-donttest`, and vignette rebuilding pass. Installed tests record 2,026
passes, five warnings, three skips, and only the clusterwise failure. The
CRAN-style check status is one error, two warnings, and zero notes; it remains
non-clean only because of that separate clusterwise test, the repository's
existing CRAN-incoming dependency warning, and the host Homebrew-clang header
warning.

## P0: merge blockers

### 1. Preserve warnings independently of verbosity

Change `tune_design_mfa()` so candidate exclusion always emits one warning,
whether `verbose` is `TRUE` or `FALSE`.

- `verbose = TRUE` additionally emits per-fold progress messages.
- The warning remains the programmatic signal that at least one candidate was
  excluded.
- Successful results retain `fold_scores` and `failure_summary`.
- An all-candidate failure remains a detailed error rather than a warning.

Acceptance tests:

- partial candidate failure with `verbose = FALSE`: exactly one warning;
- partial candidate failure with `verbose = TRUE`: exactly one warning plus
  messages;
- no failures: no warning in either mode;
- all candidates fail: error includes candidate and fold information.

### 2. Validate all subject-keyed tuning inputs before folds

Add one upfront validator, called before fold construction, for
`trial_weights`, `feature_map`, and numeric `block_weight`.

- Named inputs must have unique, nonempty names equal to `names(Y)`.
- Error text lists missing and extra keys.
- Valid named inputs are canonicalized to `names(Y)` once.
- Unnamed inputs remain positional but must have exactly one entry per block.
- `.dmfa_subset_tuning_args()` retains defensive validation and subsets only
  canonicalized inputs by index before renaming to `names(out$Y)`.

Acceptance tests cover missing, extra, duplicated, empty, reversed, and valid
custom-subject-name cases for all three argument types. No malformed input may
reach a candidate fit or produce the generic all-candidates-failed error.

### 3. Document prediction estimability

Add a `predict.design_mfa()` details section explaining:

- full-response predictions require the complete reference functional to be
  estimable;
- a term's partial sum of squares can remain estimable when unrelated
  between-subject columns are aliased with subject intercepts;
- the current effect-matrix contract after the repair in P0.5;
- how to inspect `estimability(object)` and its effect/prediction diagnostics.

Regenerate Rd files and run the repository documentation-consistency gate.
Add a documentation test or source assertion that the generated Rd contains
the estimability contract.

### 4. Remove generated tarballs from version control

- Add a root `.gitignore` rule for package tarballs, narrowly scoped as
  `/muscal_*.tar.gz`.
- Run `git rm --cached -- muscal_0.1.0.tar.gz`; retain or rebuild local
  artifacts only outside the repository.
- Make build/check scripts write to a `mktemp -d` directory.
- Verify `git ls-files` contains no package tarball and a fresh build leaves
  `git status --short` unchanged.
- Correct any report or development note claiming the current working-tree
  tarball was not modified.

### 5. Make effect matrices exactly estimable and scale invariant

#### Contract

For a requested term or hierarchical term family:

- partial SS is the unique nested-model SSE/objective difference;
- the reported effect matrix is an estimable linear functional of the full
  model, not the difference of two non-unique reference predictions;
- effect matrices, scores, singular values, and loading subspaces are invariant
  to nonsingular rescaling or recoding of aliased nuisance/between-subject
  columns;
- full-response prediction remains strict and may be unavailable.

#### Numerically stable rank and estimability representation

Use weighted column equilibration for rank determination and estimability:

1. derive positive column scales from the weighted design Gram diagonal;
2. form the equilibrated Gram matrix and apply the rank threshold there;
3. solve unpenalized WLS in equilibrated coordinates and transform coefficients
   back to the original model-matrix coordinates;
4. store the equilibrated row-space basis plus the column scales;
5. transform reference functionals into the same coordinates before checking
   estimability.

Zero columns receive a safe unit scale and remain rank deficient. Scaling must
not change the scientific objective or the meaning of ridge/alignment
penalties. Penalized solvers either retain their original-coordinate objective
with a separate equilibrated rank/estimability representation or transform
their penalties algebraically; they must not silently change penalty meaning.

#### Effect construction

Replace

```r
L_full %*% beta_full - L_reduced %*% beta_reduced
```

with the estimable contribution of the columns removed by the hypothesis:

```r
L_removed %*% beta_full_removed
```

followed by the existing weighted reference-grid centering. For hierarchical
effects, `removed` includes the term and its descendant interactions; for
sequential effects it includes the newly added term columns.

Test estimability on a functional with zeros outside the removed columns. Do
not retain adjustment contributions from aliased "other" columns and do not
zero coefficients based on arbitrary generalized-inverse allocations.

Keep partial SS/objective comparisons based on the full and reduced fits.

#### Acceptance tests

Using `~ saliency + age`, where age is aliased with subject intercepts:

- partial SS is invariant for age multipliers `0.01`, `1`, and `1000`;
- effect matrices and component scores are invariant across those multipliers;
- singular values and loading projectors match within conditioning-aware
  tolerance;
- component rank and estimability status do not change;
- scores are identical across age values at fixed saliency;
- results match the `~ saliency` fit because age adds no column-space
  information beyond subject intercepts;
- strict full-response prediction remains nonestimable and is documented as
  such.

Add a dense oracle that constructs the removed-term functional directly from
an independently fitted weighted model. Repeat with centered and reordered
subjects and with an interaction hierarchy.

### 6. Remove the residual SSE cancellation footgun

Although current callers calculate SSE directly, `.dmfa_wls_from_stats()`
still exposes the unstable `response_ss - fitted_ss` fallback.

- Remove that fallback, or require an explicit opt-in used only by a documented
  centered-statistics caller.
- Prefer requiring direct residual SSE from every engine.
- Assert that no production call can obtain SSE by subtracting two large raw
  sums.
- Retain large-offset tests for shared, chunked, partial-pooling, and soft paths.

## P1: test and diagnostics hardening

### 7. Make bootstrap weight freezing discriminating

Replace the equal-cell weight test with an unbalanced categorical design fitted
using `reference_weight = "observed"`.

- Force a bootstrap sample whose subject multiplicities alter observed cell
  frequencies.
- Assert the refit uses exactly the original fitted grid and weights.
- Independently compute the weights that would have been re-derived and assert
  that they differ, proving the test would fail without freezing.

### 8. Complete explained-energy tests

For both Euclidean and spatial-metric decompositions:

- full retained rank: `sum(explained) == 1` within numerical tolerance;
- truncated rank: `sum(explained) < 1` and equals retained energy divided by
  total weighted effect energy;
- spatial total energy matches an independent
  `trace(H_w Q^{-1} H_w^T)` oracle;
- zero and numerically degenerate effects return an empty, finite result.

### 9. Complete RNG restoration tests

Test seeded inference, tuning, and nested tuning for all relevant exits:

- `.Random.seed` existed and the call succeeds;
- `.Random.seed` existed and the call errors;
- `.Random.seed` was absent and the call succeeds;
- `.Random.seed` was absent and the call errors.

Restore the test runner's original RNG state with `on.exit()`. Store the
validated, integer-rounded seed in result objects. Document RNG restoration for
both tuning functions and inference.

### 10. Strengthen leakage tests

Retain direct auditing of the subject identities in every analysis and target
fit. Add a metamorphic independence test:

- alter only assessment-subject responses with a large sentinel signal;
- analysis-fit coefficients, rank choice inputs, preprocessing state, and
  reference grid must remain unchanged;
- target estimates and held-out score must change;
- repeat inside one outer nested-tuning fold.

This test must fail under an intentionally injected all-subject training fit.

### 11. Improve solver failure diagnostics

- Soft full-model nonconvergence reports iterations, achieved relative
  residual, requested tolerance, and right-hand-side norm.
- Effect-level errors identify whether the containing or reduced fit failed and
  report the same solver diagnostics.
- Partial-pooling errors report final relative objective change and tolerance.
- Tests assert the diagnostic values, not only the phrase "did not converge".

### 12. Give tuning failures stable categories

Do not group `failure_summary` solely by numeric-bearing message strings.
Record stable fields such as:

- `stage`: analysis fit, target fit, or scoring;
- `reason`: nonconvergence, nonestimable effect, zero target energy, invalid
  candidate, or other fit error;
- `condition_class`;
- exact `message` retained at fold level.

Aggregate by stable category while preserving representative messages and the
candidate/fold rows.

### 13. Finish API and documentation coverage

- Let `plot(fit)` select the first estimable effect with retained components;
  preserve explicit `term` behavior.
- Cover default plotting, nonestimable terms, one-row grids, and absent
  components.
- Add small runnable examples for every new public Rd topic. Share a tiny
  fixture where appropriate; use `\donttest{}` only for genuinely expensive
  inference or tuning examples.
- Ensure examples demonstrate continuous predictors, estimability diagnostics,
  projection, and soft-map coverage without unsupported claims.

## P2: soft-alignment scaling and performance

### 14. Avoid materializing `P_s^T P_s`

Change the soft normal operator from

```r
G_s %*% B %*% crossprod(P_s)
```

to sparse operator composition equivalent to

```r
(G_s %*% B %*% t(P_s)) %*% P_s
```

and compute the Jacobi map diagonal from column sums of `P_s^2`. This removes
the stored reference-by-reference Gram matrix per subject.

Extend `feature_chunk_size` to cross-products involving `Y_s` and `P_s`, not
only final SSE evaluation. Cache map-only quantities once per top-level fit and
reuse them across full and reduced models.

Acceptance gates:

- bitwise-equivalent or tolerance-equivalent coefficients on deterministic
  fixtures;
- dense Kronecker oracle with genuinely non-diagonal `P_s^T P_s` still passes;
- an instrumented operator-count test proves one normal-operator application
  per CG iteration;
- peak memory scales with map nonzeros and `q * p`, not a stored
  `p_reference^2` object;
- `feature_chunk_size` changes peak workspace without changing results.

### 15. Establish a real performance envelope

- Run at least five fresh-process repetitions per configuration and report
  median, range/IQR, hardware, R, BLAS, Matrix, and compiler versions.
- Separate input bytes, peak heap, fit size, largest allocation, and temporary
  sparse-operator storage.
- Add a multi-subject 100,000-reference-feature soft-alignment run only after
  the Gram-removal work passes memory tests.
- Benchmark map construction for all subjects, not one subject.
- Keep shared, spatial, and soft results in separate rows and never extrapolate
  an unrun engine.

Performance tests should use soft budgets in ordinary CI and stricter repeated
benchmarks in a scheduled job.

## P3: spatial evidence and claim discipline

### 16. Correct calibration provenance

- Relabel existing legacy and fair CSVs as retrospective exploratory evidence.
- Remove the claim that the existing thresholds were prospectively declared.
- Create a versioned protocol file containing estimands, seeds, scenarios,
  comparators, thresholds, replicate count, and analysis code hash.
- Commit that protocol before generating confirmatory results.
- Store the protocol Git SHA and result-file hash in every confirmatory output.

### 17. Separate negative controls, recovery, and comparative benefit

The confirmatory spatial suite contains distinct scenarios:

1. known exact same-index identity: the shared model is the expected winner;
2. identity-map reduction: soft and shared engines agree;
3. independently sampled/jittered subject coordinates without exact index
   identity: compare soft interpolation with prespecified nearest-neighbor and
   other defensible baselines;
4. controlled map misspecification and incomplete coverage;
5. planted nulls for inferential calibration.

Responses must be generated from continuous coordinate fields, never through
the fitted interpolation operator. Tune map or penalty parameters inside
training subjects when they are selected from data.

### 18. Add actual inferential calibration

For every supported explicit permutation scheme:

- predeclare the exchangeability assumptions;
- report rejection rates and Monte Carlo confidence intervals under the null;
- report power under planted effects;
- distinguish term-level and component-level multiplicity;
- treat descriptive null energy only as a numerical diagnostic.

No type-I, coverage, or comparative-alignment claim ships before these gates
pass.

## Terminology decision

Retain `design_mfa()` for this merge only if all public prose describes it as
an ASCA+/partial-RDA-style design decomposition and explicitly states that it
does not perform MFA block normalization. Add a terminology regression check
over README, vignette, and Rd output. Because the API is not yet merged, a
final pre-merge decision may instead rename the primary API to
`design_asca()`; if chosen, update all exports, methods, docs, examples, and
tests atomically rather than adding two competing names.

## Execution order

1. Add failing regressions for P0.1--P0.6.
2. Implement quick contract/repository fixes P0.1--P0.4.
3. Implement equilibrated estimability and term-functional effects P0.5,
   followed by the SSE cleanup P0.6.
4. Complete P1 test and diagnostics work.
5. Refactor and benchmark the soft operator under P2.
6. Correct retrospective wording, commit the confirmatory protocol, and only
   then run P3 calibration.
7. Make the terminology decision and regenerate all documentation.
8. Run the final verification matrix from a clean temporary installation.

## Final verification matrix

- all Design-MFA contract, oracle, metamorphic, fuzz, and regression tests;
- deterministic large-offset and covariate-rescaling sweeps;
- shared/streamed/chunked equivalence;
- partial-pooling KKT oracle and penalty-limit tests;
- soft non-diagonal Kronecker oracle and operator-count instrumentation;
- bootstrap and permutation tests, including RNG restoration on errors;
- held-out and nested tuning leakage tests;
- documentation regeneration and runnable examples;
- vignette render from a clean installed package;
- `R CMD build` and `R CMD check` in a temporary directory;
- tracked-file cleanliness before and after build;
- repeated performance runs and memory audits;
- confirmatory calibration only against a protocol committed beforehand.

The known `penalized_mfa_clusterwise()` objective-monotonicity failure remains
separate from Design-MFA implementation work, but final package-check reporting
must continue to identify it precisely rather than calling the entire package
clean.

# DISTATIS implementation and validation

Validated locally on 2026-10-18 with R 4.5.1. The checkout HEAD was
`06eb85f1eee84eff16e1f0bcd23f54730ae18152`; these implementation changes are
uncommitted. This is local implementation evidence, not hosted CI or release evidence.

## Implemented behavior

`distatis()` accepts distance matrices, `dist` objects, lists, and three-dimensional
arrays. It defaults to ordinary distances, MFA normalization, and RV weighting.
`sorting_distances()` converts categorical partitions to binary distances.
The fit has separate object and assessor spectra, positive-component truncation,
partial object maps, normalized compromise reconstruction, and diagnostic metadata.

`project_distance()` and `project()` project a supplementary assessor's partial
object map. `project_subjects()` additionally projects its RV coordinates and
reports full-compromise similarity/distance separately from retained-subspace
residuals. Covariate projection, object/assessor plots, partial maps, weights,
RV heatmaps, and full-spectrum inertia plots are available. Plot axes use equal
distance units by default. Coincident assessor points can have overlapping labels:
for example, assessors A3 and A8 give identical partitions in the beer fixture.

The shared Gram-fitting core preserves COVSTATIS defaults. Comparing saved
pre-refactor fits with the refactored implementation at all three normalization
settings gave agreement within `1e-12` for C, weights, scores, partial scores,
projection matrices, and compromises. Single-component fitting was repaired.

Non-Euclidean inputs, zero Gram tables, missing values, conflicting object order,
and tied leading RV roots produce explicit errors. Inference is intentionally
unavailable: COVSTATIS inference and adapter entry points reject DISTATIS fits.

## Independent numerical evidence

The ordinary test suite checks all 100 printed RV entries, all ten weights,
seven compromise roots, and 56 sign-aligned product scores from the supplied
Abdi et al. (2007) beer example, using its printed rounding precision.
Coordinate-derived MDS cases test nonbinary distances, squared-input semantics,
rank-one/single-component behavior and distance reconstruction. Permutation,
rescaling, sorting-label renaming, partial barycenters, supplementary projection,
separate map ranks, full-compromise diagnostics, and invalid inputs are covered.

An additional installed-reference comparison used DistatisR 1.1.2, installed only
in a temporary library with prettyGraphs 2.2.0. On both the beer fixture and a
nonbinary Euclidean-distance fixture, direct comparisons with its `MFA` and
`NONE` options passed an absolute tolerance of `1e-10`. The largest observed
difference was `2.1316e-14`. Comparisons cover RV matrices, alpha, compromise,
full object spectra, sign-aligned object/partial coordinates, and the Gram matrix
of retained assessor coordinates. DistatisR received squared distance inputs.

Reproducer: `dev/validate_distatis_reference.R`; raw output:
`dev/distatis_reference_results.txt`. DistatisR is not a package dependency.
With DistatisR 1.1.2 available, run from the repository root:

```sh
LC_ALL=C Rscript --vanilla dev/validate_distatis_reference.R
```

An initial proposed comparison of Frobenius normalization with DistatisR's
`SUMPCA` failed. Inspecting the installed `SUMPCAnormCP()` established that
`SUMPCA` divides by trace, rather than Frobenius norm. That invalid option mapping
was removed; no direct installed-reference parity is claimed for Frobenius.
Its behavior is covered by the coordinate, scaling, and invariant tests.
The original failed receipt remains at
`/tmp/muscal-distatis-reference-first-attempt.log`.

Downloaded reference archive SHA256 values:

| Archive | SHA256 |
| --- | --- |
| DistatisR_1.1.2.tar.gz | `49972799ea518bbe4378ad330736dbee73347db7013a235275a57be3c8e3d6e9` |
| prettyGraphs_2.2.0.tar.gz | `26a80afa54c1dc62800b9ba7a3c3fa05f7b0b1aa7fb652e9ef0ed9f5304e466d` |

## Package checks

Roxygen Rd/NAMESPACE regeneration completed. The final package check used:

```r
rcmdcheck::rcmdcheck(
  path = ".", args = c("--as-cran", "--no-manual"),
  build_args = "--no-manual", quiet = TRUE, error_on = "never"
)
```

Final result: **0 errors, 2 warnings, 0 notes**. Examples, native installation,
S3 registrations, code/documentation checks, tests, and vignette build/rebuild
completed successfully. The test suite recorded **2,262 passed assertions,
0 failures, 5 warnings, and 4 skips**. Focused DISTATIS/COVSTATIS runs also passed;
they exposed dependency-version notices and the existing multivarious `prep()`
deprecation warning. Full-suite skips concern unavailable convergence information
(two), unavailable design-MFA multifer stability support (one), and unavailable
source Rd files in the installed test copy (one).

The two package-check warnings remain:

1. CRAN incoming feasibility flags the existing `Remotes` field and availability
   of the existing non-mainstream dependencies `multidesign` and `multifer`.
2. Installation reports an unknown compiler warning group
   `-Wfixed-enum-extension` from R's `R_ext/Boolean.h` with the local LLVM compiler.

Those metadata, dependencies, and compiler inputs were not changed by this work.
No PDF manual, hosted CI, publication, or inference qualification is claimed.

The final tested tarball is
`/tmp/muscal-distatis-check-final/muscal_0.1.0.tar.gz`, SHA256
`cc3d01f9e44c5f2877e360f40e0fff5fd7e82fc4fc4058f61e82297d85ac9b03`.
Its implementation files, NAMESPACE, DISTATIS tests, generated Rd, and vignette
were compared byte-for-byte with the checkout. R CMD build reformats DESCRIPTION
and adds its normal packaging metadata.

Raw evidence remains in `/tmp/muscal-distatis-check-final/muscal.Rcheck/`, with
the result object at `/tmp/muscal-distatis-check-final/result.rds`. The combined
check summary and process exit/timing metadata are at
`/tmp/muscal-distatis-check-final-run.log` and its `.meta.json` sidecar. The
earlier full check, before equal-unit plot adjustments, also recorded 0 errors,
2 warnings, and 0 notes in `/tmp/muscal-distatis-check/`.

The worked example is `vignettes/distatis.Rmd`. Both object and assessor maps
were rendered and visually inspected, and plot construction is tested through
`ggplot2::ggplot_build()`. `git diff --check` passed.

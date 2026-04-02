# Easy Front-End for Bilinear Mixed

Lightweight wrapper around \[bilinear_mixed()\] that uses
\[bilinear_mixed_recommend()\] defaults. Optionally performs compact
subject-block CV tuning via \[bilinear_mixed_tune()\].

## Usage

``` r
bilinear_mixed_easy(
  data,
  subject,
  z = NULL,
  y = NULL,
  row_design = NULL,
  mode = c("auto", "seed_axis", "seed_repeat", "both"),
  connectivity_type = c("auto", "cross", "symmetric"),
  profile = c("fast", "balanced", "adaptive"),
  tune = FALSE,
  tune_grid = NULL,
  n_folds = 3,
  metric = c("auto", "reconstruction", "trait_r2"),
  seed = 1,
  verbose = FALSE,
  ...
)
```

## Arguments

- data:

  A list of numeric connectivity matrices.

- subject:

  Subject identifier.

- z:

  Optional repeat-level design.

- y:

  Optional subject-level traits.

- row_design:

  Optional row-level covariates.

- mode:

  One of `"auto"`, `"seed_axis"`, `"seed_repeat"`, or `"both"`.

- connectivity_type:

  One of `"auto"`, `"cross"`, or `"symmetric"`.

- profile:

  One of `"fast"`, `"balanced"`, or `"adaptive"`.

- tune:

  Logical; if `TRUE`, run \[bilinear_mixed_tune()\].

- tune_grid:

  Optional grid for tuning; see \[bilinear_mixed_tune()\].

- n_folds:

  Number of subject-block CV folds for tuning.

- metric:

  One of `"auto"`, `"reconstruction"`, or `"trait_r2"`.

- seed:

  Integer random seed for fold assignment.

- verbose:

  Logical verbosity.

- ...:

  Optional overrides passed to fitting/tuning.

## Value

Fitted `"bilinear_mixed"` object if `tune = FALSE`, otherwise a
`"bilinear_mixed_tuning"` object.

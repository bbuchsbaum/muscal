# Tune Bilinear Mixed Hyperparameters via Subject-Block CV

Tunes a compact set of key parameters (\`r_seed\`, \`r_roi\`,
\`k_subject\`, optional \`lambda_y\`) using subject-blocked
cross-validation.

## Usage

``` r
bilinear_mixed_tune(
  data,
  subject,
  z = NULL,
  y = NULL,
  row_design = NULL,
  mode = c("auto", "seed_axis", "seed_repeat", "both"),
  connectivity_type = c("auto", "cross", "symmetric"),
  profile = c("fast", "balanced", "adaptive"),
  grid = NULL,
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

- grid:

  Optional candidate grid. Either:

  - a `data.frame` with columns among `r_seed`, `r_roi`, `k_subject`,
    `lambda_y`

  - or a named list of vectors for those fields.

- n_folds:

  Number of subject-block CV folds.

- metric:

  One of `"auto"`, `"reconstruction"`, or `"trait_r2"`.

- seed:

  Random seed for fold assignment.

- verbose:

  Logical verbosity.

- ...:

  Additional fixed arguments forwarded to \[bilinear_mixed()\].

## Value

A list of class `"bilinear_mixed_tuning"` with `best_params`, `results`,
and refit `fit`.

# Recommend Practical Defaults for Bilinear Mixed

Provides data-adaptive defaults for \[bilinear_mixed()\] so users can
avoid manual specification of many tuning parameters.

## Usage

``` r
bilinear_mixed_recommend(
  data,
  subject,
  z = NULL,
  y = NULL,
  row_design = NULL,
  mode = c("auto", "seed_axis", "seed_repeat", "both"),
  connectivity_type = c("auto", "cross", "symmetric"),
  profile = c("fast", "balanced", "adaptive"),
  max_r_seed = 6,
  max_r_roi = 20,
  max_k_subject = 5,
  sym_tol = 1e-08
)
```

## Arguments

- data:

  A list of numeric connectivity matrices.

- subject:

  Subject identifier with length `length(data)`.

- z:

  Optional repeat-level design.

- y:

  Optional subject-level traits.

- row_design:

  Optional row-level covariates (same number of rows as each
  connectivity matrix).

- mode:

  One of `"auto"`, `"seed_axis"`, `"seed_repeat"`, or `"both"`.

- connectivity_type:

  One of `"auto"`, `"cross"`, or `"symmetric"`.

- profile:

  One of `"fast"`, `"balanced"`, or `"adaptive"`.

- max_r_seed:

  Maximum suggested row rank.

- max_r_roi:

  Maximum suggested column rank.

- max_k_subject:

  Maximum suggested subject latent dimension.

- sym_tol:

  Numeric tolerance for symmetry detection when
  `connectivity_type = "auto"`.

## Value

Named list of recommended arguments for \[bilinear_mixed()\].

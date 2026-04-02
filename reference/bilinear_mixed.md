# Bilinear Repeated-Measures Mixed Model

Fits a matrix-aware latent model for repeated connectivity matrices with
optional subject-level supervision. The method supports two
complementary representations:

- `"seed_axis"`: treat row entities as the row axis of each matrix (\\s
  \times p\\).

- `"seed_repeat"`: treat row profiles as repeated observations nested
  within subject/repeat.

- `"both"`: fit both heads with a shared ROI basis.

## Usage

``` r
bilinear_mixed(
  data,
  subject,
  z = NULL,
  y = NULL,
  row_design = NULL,
  seed_design = NULL,
  mode = c("seed_axis", "seed_repeat", "both"),
  connectivity_type = c("auto", "cross", "symmetric"),
  symmetric = FALSE,
  sym_tol = 1e-08,
  complexity = c("manual", "fast", "balanced", "adaptive"),
  r_seed = 2,
  r_roi = 8,
  k_subject = 2,
  include_seed_interactions = TRUE,
  repeat_id = NULL,
  laplacian_repeat = NULL,
  laplacian_seed = NULL,
  laplacian_roi = NULL,
  laplacian_design = NULL,
  lambda_w = 0.01,
  lambda_t = 0.01,
  lambda_b = 0.001,
  lambda_a = 0.001,
  lambda_y = 0,
  lambda_repeat = 0,
  lambda_seed = 0,
  lambda_roi = 0,
  lambda_design_smooth = 0,
  max_iter = 50,
  tol = 1e-06,
  verbose = FALSE,
  ...
)

# Default S3 method
bilinear_mixed(data, ...)

# S3 method for class 'list'
bilinear_mixed(
  data,
  subject,
  z = NULL,
  y = NULL,
  row_design = NULL,
  seed_design = NULL,
  mode = c("seed_axis", "seed_repeat", "both"),
  connectivity_type = c("auto", "cross", "symmetric"),
  symmetric = FALSE,
  sym_tol = 1e-08,
  complexity = c("manual", "fast", "balanced", "adaptive"),
  r_seed = 2,
  r_roi = 8,
  k_subject = 2,
  include_seed_interactions = TRUE,
  repeat_id = NULL,
  laplacian_repeat = NULL,
  laplacian_seed = NULL,
  laplacian_roi = NULL,
  laplacian_design = NULL,
  lambda_w = 0.01,
  lambda_t = 0.01,
  lambda_b = 0.001,
  lambda_a = 0.001,
  lambda_y = 0,
  lambda_repeat = 0,
  lambda_seed = 0,
  lambda_roi = 0,
  lambda_design_smooth = 0,
  max_iter = 50,
  tol = 1e-06,
  verbose = FALSE,
  ...
)
```

## Arguments

- data:

  A list of numeric matrices. Each element is one observed connectivity
  matrix (typically row-entity \\\times\\ ROI, or symmetric ROI
  \\\times\\ ROI).

- subject:

  Subject identifier for each matrix in `data`.

- z:

  Optional repeat-level design matrix/data frame/vector (same number of
  rows as `length(data)`).

- y:

  Optional subject-level traits. Can be a named numeric vector, a
  matrix/data frame with one row per subject, or a vector ordered by
  subject level.

- row_design:

  Optional row-entity design (vector/matrix with one row per matrix
  row). Used by `mode = "seed_repeat"`.

- seed_design:

  Deprecated alias of `row_design`, retained for backward compatibility.

- mode:

  One of `"seed_axis"`, `"seed_repeat"`, or `"both"`.

- connectivity_type:

  One of `"auto"`, `"cross"`, or `"symmetric"`.

- symmetric:

  Deprecated logical alias for `connectivity_type = "symmetric"`.

- sym_tol:

  Numeric tolerance for automatic symmetry detection in
  `connectivity_type = "auto"`.

- complexity:

  One of `"manual"`, `"fast"`, `"balanced"`, or `"adaptive"`. If not
  `"manual"`, data-adaptive defaults are inferred for unspecified tuning
  parameters.

- r_seed:

  Rank for seed mode basis (row basis).

- r_roi:

  Rank for ROI mode basis (column basis).

- k_subject:

  Subject latent dimensionality.

- include_seed_interactions:

  Logical; if `TRUE`, includes \\z \otimes row\\design\\ interactions in
  the seed-repeat head.

- repeat_id:

  Optional repeat/condition index per matrix. Required when using
  `laplacian_repeat`.

- laplacian_repeat:

  Optional repeat Laplacian matrix.

- laplacian_seed:

  Optional seed-axis Laplacian.

- laplacian_roi:

  Optional ROI-axis Laplacian.

- laplacian_design:

  Optional design Laplacian on columns of `z`.

- lambda_w:

  Ridge penalty for `W`.

- lambda_t:

  Ridge penalty for subject scores `t`.

- lambda_b:

  Ridge penalty for design effects.

- lambda_a:

  Ridge penalty for supervision map `A`.

- lambda_y:

  Weight for supervised trait term.

- lambda_repeat:

  Smoothing strength for repeat smoothing.

- lambda_seed:

  Basis smoothing strength for seed basis.

- lambda_roi:

  Basis smoothing strength for ROI basis.

- lambda_design_smooth:

  Laplacian smoothing for design effects.

- max_iter:

  Maximum ALS iterations for the core mixed model.

- tol:

  Relative convergence tolerance for the core mixed model.

- verbose:

  Logical; if `TRUE`, prints iteration diagnostics.

- ...:

  Reserved for future extensions.

## Value

A list of class `"bilinear_mixed"` containing fitted head(s), bases,
core-model parameters, and backprojected effect maps.

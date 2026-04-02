# Estimate the number of components and loading reliability (shared helpers)

These helpers provide a common interface for component retention and
loading reliability across `bada`, `mfa`, `penalized_mfa`, and related
models. They are intentionally modular so alternative criteria (e.g.,
AIC, permutations) can be added later without changing call sites.

## Usage

``` r
estimate_components(
  fit,
  method = c("rmt", "variance"),
  sdev = NULL,
  V_list = NULL,
  n = NULL,
  tail_q = 0.2,
  alpha = 0.05
)

loading_reliability(
  fit,
  method = c("bootstrap"),
  boot_loadings = NULL,
  alpha = 0.05,
  V_ref = NULL
)
```

## Arguments

- fit:

  A fitted object (e.g., `bada`, `mfa`, `penalized_mfa`,
  `penalized_mfa_clusterwise`).

- method:

  Component or reliability method. For `estimate_components`, one of
  `"rmt"` (Marchenko–Pastur edge) or `"variance"` (keep all non-zero).
  For `loading_reliability`, currently `"bootstrap"`.

- sdev:

  Optional singular values; if `NULL`, taken from `fit$sdev` when
  available.

- V_list:

  Optional list of loading matrices; if `NULL`, the helper attempts to
  extract from `fit` (using `attr(fit, "V_list")` or splitting `fit$v`
  by `block_indices` when present).

- n:

  Optional number of observations (rows); defaults to `nrow` of scores
  if available.

- tail_q:

  Numeric; fraction of smallest eigenvalues used to estimate noise
  variance for the RMT test (default 0.2). Only used in
  `estimate_components` when `method = "rmt"`.

- alpha:

  Confidence level for hypothesis tests or bootstrap intervals.

- boot_loadings:

  For `loading_reliability`: a list of loading matrices from bootstrap
  refits (all same dimensions as the reference loadings).

- V_ref:

  For `loading_reliability`: optional reference loading matrix. If
  `NULL`, extracted from `fit` using the same logic as `V_list`.

## Value

- estimate_components:

  A list with `keep`, `criterion`, and `method_used`.

- loading_reliability:

  A list with matrices `mean`, `sd`, `lower`, `upper`, and
  `sign_consistency`.

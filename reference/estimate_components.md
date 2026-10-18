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
  alpha = 0.05,
  noise_variance = NULL,
  spectrum = c("unknown", "covariance", "singular"),
  noise_model = c("unknown", "isotropic"),
  n_eff = n,
  check_icc = FALSE,
  common_features = NULL
)

loading_reliability(
  fit,
  method = c("bootstrap"),
  boot_loadings = NULL,
  alpha = 0.05,
  V_ref = NULL,
  gap_tol = 1e-06
)
```

## Arguments

- fit:

  A fitted object (e.g., `bada`, `mfa`, `penalized_mfa`,
  `penalized_mfa_clusterwise`).

- method:

  Component or reliability method. For `estimate_components`, one of
  `"rmt"` (Marchenko\<e2\>\<80\>\<93\>Pastur edge) or `"variance"` (keep
  all non-zero). For `loading_reliability`, currently `"bootstrap"`.

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
  variance in the historical automatic RMT heuristic. Retained for
  compatibility; automatic noise estimation is now unavailable.

- alpha:

  Tail probability for bootstrap intervals in \`loading_reliability\`.
  In \`estimate_components\`, a significance threshold only for
  explicitly enabled experimental ICC. Explicit alpha on an edge-only or
  variance-only call warns that it is inapplicable; it does not change
  an MP edge.

- noise_variance, spectrum, noise_model, n_eff:

  Supplied-noise RMT model and spectral units; see
  \[significant_components()\]. No arbitrary weighted compromise
  spectrum is assumed to satisfy the isotropic covariance model.

- check_icc, common_features:

  Optional experimental common-feature ICC criterion and required
  feature identity contract; see \[significant_components()\].

- boot_loadings:

  For `loading_reliability`: a list of loading matrices from bootstrap
  refits (all same dimensions as the reference loadings).

- V_ref:

  For `loading_reliability`: optional reference loading matrix. If
  `NULL`, extracted from `fit` using the same logic as `V_list`.

- gap_tol:

  For loading reliability, a nonnegative tolerance for gaps between
  squared reference component scales, normalized by the largest scale
  squared. Reference \`fit\$sdev\` must follow the supplied reference
  column order.

## Value

- estimate_components:

  A list with `keep`, `criterion`, and `method_used`. RMT additionally
  returns `diagnostics` with feature dimensions, applicability and
  criterion details. Missing model or observation information yields
  `method_used="rmt_unavailable"`, NA decisions and no retained
  components, never a silent variance fallback.

- loading_reliability:

  A list with matrices `mean`, `sd`, `lower`, `upper`, and
  `sign_consistency`, plus `alignment` diagnostics. Exact one-to-one
  matching uses absolute normalized loading inner products, followed by
  sign alignment. Known tied/zero reference roots have unavailable (NA)
  coordinate summaries and geometric subspace distances. Without
  reference component scales, intervals are conditional on the supplied
  reference alignment and identifiability is reported as unknown. Only
  the retained reference gaps are checked; omitted-axis gaps and
  bootstrap coverage are not established by these summaries.

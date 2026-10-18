# Identify Significant Components via RMT and Inter-block Coherence Tests

Assess declared supplied-noise spectral edges and optional
common-feature loading coherence. Arbitrary multiblock loadings are not
assumed to share a feature space, and an asymptotic support edge is not
a calibrated p-value.

## Usage

``` r
significant_components(
  fit,
  n,
  k_vec = NULL,
  alpha = 0.05,
  check_rmt = TRUE,
  tail_frac = 0.3,
  noise_variance = NULL,
  spectrum = c("unknown", "covariance", "singular"),
  noise_model = c("unknown", "isotropic"),
  n_eff = n,
  check_icc = FALSE,
  common_features = NULL
)
```

## Arguments

- fit:

  A fitted multi-block model object containing at minimum a \`V_list\`
  attribute (list of loading matrices) and optionally \`sdev\` (singular
  values).

- n:

  Integer; the number of observations (rows) used to fit the model.

- k_vec:

  Optional integer vector of block dimensions (number of columns per
  original data block). If NULL, inferred from loading-matrix row
  counts. Supplied counts and block-index metadata must agree with those
  fitted feature counts, including any preprocessing that changes the
  feature dimension.

- alpha:

  Finite probability strictly between 0 and 1; the threshold for
  explicitly enabled experimental ICC (default 0.05). Explicit alpha
  without ICC warns that it is inapplicable; an MP support edge is not
  alpha-level.

- check_rmt:

  Logical; if TRUE (default), performs the Marchenko-Pastur edge test to
  check if eigenvalues exceed the noise threshold.

- tail_frac:

  Retained for compatibility. Automatic noise estimation from truncated
  fitted spectra is unavailable; this argument does not estimate noise.

- noise_variance:

  Positive finite noise variance, supplied independently in
  covariance-eigenvalue units. No automatic estimate is assumed.

- spectrum:

  Spectral units: \`"covariance"\` means \`sdev^2\` are covariance
  eigenvalues; \`"singular"\` means raw singular values, divided by
  \`n_eff\` after squaring. The default \`"unknown"\` leaves RMT
  unavailable.

- noise_model:

  \`"isotropic"\` explicitly declares the IID isotropic noise model. The
  default \`"unknown"\` leaves RMT unavailable, including for arbitrary
  weighted MFA compromise spectra.

- n_eff:

  Positive effective degrees of freedom of the covariance normalization,
  at most \`n\`; defaults to \`n\` for an uncentered model. For a
  centered unweighted sample covariance, supply \`n - 1\`.

- check_icc:

  Logical; explicitly request the experimental common-feature coherence
  criterion (default FALSE). It is independent of RMT applicability.

- common_features:

  Unique character feature identities declaring a common coordinate
  domain for ICC. Every loading matrix must have exactly these row
  identities; rows are matched by identity before comparison. Equal
  widths alone do not establish a common space.

## Value

A list with the following elements:

- keep:

  Component indices passing all enabled criteria. A requested
  unavailable criterion cannot retain components. With both disabled, no
  components are declared retained.

- rmt_pass:

  Logical edge decisions, or NA if RMT is unavailable/disabled.

- icc_pass:

  Experimental ICC decisions, or NA when disabled.

- icc:

  Numeric vector of inter-block coherence values per component.

- icc_pvalue:

  Numeric vector of p-values for the ICC test.

- mp_edge:

  The Marchenko-Pastur edge threshold (NA if RMT skipped).

- sigma2_est:

  Supplied noise variance (legacy field name; NA if unavailable).

- lambda:

  Eigenvalues (squared singular values).

- n, k_vec, alpha:

  Input parameters echoed back.

- rmt_status, rmt_reason:

  Whether the requested edge was evaluated, unavailable, or explicitly
  disabled, and why. Unavailable is never a pass.

- gamma, n_eff, spectrum, noise_model:

  Resolved aspect ratio and declared units/model.

- icc_status:

  \`"disabled"\` or \`"evaluated_experimental"\`.

## Details

The supplied-noise RMT criterion uses the asymptotic Marchenko-Pastur
upper support edge. This is not an alpha-level finite-sample
significance test. Its isotropic covariance representation must be
declared explicitly; no noise variance is inferred from a few retained
fitted components.

The ICC test measures the squared cosine similarity of loading vectors
across pairs of explicitly matched common-feature blocks. Its Gaussian
null approximation with Bonferroni correction remains experimental;
feature-domain validation does not establish statistical calibration.

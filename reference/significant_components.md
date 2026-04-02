# Identify Significant Components via RMT and Inter-block Coherence Tests

This function determines which components from a multi-block analysis
(e.g., penalized MFA) are statistically significant. It combines two
complementary tests: a Random Matrix Theory (RMT) test based on the
Marchenko-Pastur distribution, and an Inter-block Coherence (ICC) test
that assesses whether loadings are consistent across blocks.

## Usage

``` r
significant_components(
  fit,
  n,
  k_vec = NULL,
  alpha = 0.05,
  check_rmt = TRUE,
  tail_frac = 0.3
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
  block). If NULL, inferred from \`fit\` via \`block_indices\` or
  \`V_list\`.

- alpha:

  Numeric; significance level for hypothesis tests (default 0.05).

- check_rmt:

  Logical; if TRUE (default), performs the Marchenko-Pastur edge test to
  check if eigenvalues exceed the noise threshold.

- tail_frac:

  Numeric; fraction of smallest eigenvalues used to estimate noise
  variance for the RMT test (default 0.3).

## Value

A list with the following elements:

- keep:

  Integer vector of component indices that pass both tests.

- rmt_pass:

  Logical vector indicating which components pass the RMT test.

- icc_pass:

  Logical vector indicating which components pass the ICC test.

- icc:

  Numeric vector of inter-block coherence values per component.

- icc_pvalue:

  Numeric vector of p-values for the ICC test.

- mp_edge:

  The Marchenko-Pastur edge threshold (NA if RMT skipped).

- sigma2_est:

  Estimated noise variance (NA if RMT skipped).

- lambda:

  Eigenvalues (squared singular values).

- n, k_vec, alpha:

  Input parameters echoed back.

## Details

The RMT test uses the Marchenko-Pastur distribution to identify
eigenvalues that exceed the expected bulk edge under the null hypothesis
of pure noise. Noise variance is estimated robustly from the tail of the
eigenvalue distribution.

The ICC test measures the squared cosine similarity of loading vectors
across all pairs of blocks. Under the null hypothesis of random
loadings, the expected value is 1/k (where k is the harmonic mean of
block dimensions). A z-test with Bonferroni correction is applied.

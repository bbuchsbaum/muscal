# Aligned Multiple Factor Analysis (Aligned MFA)

Aligned MFA estimates a shared score matrix for a set of \*latent
reference rows\* when multiple blocks have different row sets. Each
block \`X\[\[k\]\]\` is linked to the shared rows via an integer index
vector \`row_index\[\[k\]\]\`, but \*\*no single block is privileged\*\*
as an anchor.

This is the symmetric counterpart to \[anchored_mfa()\], which uses an
explicit reference block \`Y\` and estimates scores in the row space of
\`Y\`.

## Usage

``` r
aligned_mfa(
  X,
  row_index,
  N = NULL,
  preproc = multivarious::center(),
  ncomp = 2,
  normalization = c("MFA", "None", "custom"),
  alpha = NULL,
  feature_groups = NULL,
  feature_lambda = 0,
  max_iter = 50,
  tol = 1e-06,
  ridge = 1e-08,
  verbose = FALSE,
  ...
)
```

## Arguments

- X:

  A list of numeric matrices/data.frames. Each element \`X\[\[k\]\]\` is
  \`n_k × p_k\`.

- row_index:

  A list of integer vectors. \`row_index\[\[k\]\]\` has length \`n_k\`
  and maps rows of \`X\[\[k\]\]\` to shared rows in \`1..N\`.

- N:

  Optional integer specifying the number of shared rows. If \`NULL\`
  (default), \`N\` is inferred as \`max(unlist(row_index))\`.

- preproc:

  A \`multivarious\` preprocessing pipeline (a
  \`pre_processor\`/\`prepper\`) or a list of them. If a list, it must
  have length \`length(X)\` and will be applied to \`X\` in that order.

- ncomp:

  Integer number of components to extract.

- normalization:

  Block weighting scheme. \`"MFA"\` uses inverse squared first singular
  value per block; \`"None"\` uses uniform weights; \`"custom"\` uses
  \`alpha\`.

- alpha:

  Optional numeric vector of per-block weights (length \`length(X)\`),
  used when \`normalization = "custom"\`.

- feature_groups:

  Feature prior specification. One of: \* \`NULL\` (no feature prior),
  \* \`"colnames"\` to group features with identical column names across
  blocks, \* a \`data.frame\` with columns \`block\`, \`feature\`,
  \`group\` and optional \`weight\`.

- feature_lambda:

  Non-negative scalar controlling strength of the feature prior.

- max_iter:

  Maximum number of alternating least-squares iterations.

- tol:

  Relative tolerance on the objective for convergence.

- ridge:

  Non-negative ridge stabilization added to normal equations.

- verbose:

  Logical; if \`TRUE\`, prints iteration progress.

- ...:

  Unused (reserved for future extensions).

## Value

An object inheriting from \`multivarious::multiblock_biprojector\` with
additional class \`"aligned_mfa"\`. The object contains global scores in
\`s\`, concatenated loadings in \`v\`, and block mappings in
\`block_indices\`. Additional fields include \`V_list\`, \`row_index\`,
\`alpha_blocks\`, and \`objective_trace\`.

## Details

\## Model The fitted model has the form: \$\$X_k \approx
S\[\mathrm{idx}\_k,\] V_k^\top\$\$ where \`S\` is \`N × ncomp\`, each
\`V_k\` is \`p_k × ncomp\`, and \`idx_k\` maps rows of \`X_k\` to
\`1..N\`. Repeated indices are allowed (e.g., repeated measures), and
contributions are aggregated in the score updates.

\## Feature similarity prior When \`feature_lambda \> 0\` and
\`feature_groups\` is supplied, Aligned MFA applies the same
group-shrinkage penalty as \[anchored_mfa()\], pulling loading vectors
of grouped features toward a shared group center.

## Examples

``` r
# \donttest{
set.seed(1)
N <- 30
X1 <- matrix(rnorm(20 * 10), 20, 10)
X2 <- matrix(rnorm(15 * 8), 15, 8)
idx1 <- sample.int(N, nrow(X1), replace = FALSE)
idx2 <- sample.int(N, nrow(X2), replace = FALSE)

fit <- aligned_mfa(list(X1 = X1, X2 = X2), list(X1 = idx1, X2 = idx2), ncomp = 2)
#> Applying the same preprocessor definition independently to each block.
stopifnot(nrow(multivarious::scores(fit)) == N)
# }
```

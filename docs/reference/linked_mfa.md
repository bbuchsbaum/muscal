# Anchored Multiple Factor Analysis (Anchored MFA)

Anchored MFA generalizes Multiple Factor Analysis (MFA) to the case
where a single reference block \`Y\` (with \`N\` rows) is linked to
multiple blocks \`X_k\` that may have different numbers of rows. Each
row of \`X_k\` is mapped to a row of \`Y\` via an index vector
\`row_index\[\[k\]\]\`. The model estimates a shared score matrix \`S\`
for the rows of \`Y\` and block-specific loading matrices for \`Y\` and
each \`X_k\`.

Optionally, a feature prior can be provided to encourage corresponding
or similar features across different \`X_k\` blocks to have similar
loading vectors.

## Usage

``` r
anchored_mfa(
  Y,
  X,
  row_index,
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

linked_mfa(
  Y,
  X,
  row_index,
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

- Y:

  Numeric matrix/data.frame (\`N × q\`) serving as the reference block.

- X:

  A list of numeric matrices/data.frames. Each element \`X\[\[k\]\]\` is
  \`n_k × p_k\`.

- row_index:

  A list of integer vectors. \`row_index\[\[k\]\]\` has length \`n_k\`
  and maps rows of \`X\[\[k\]\]\` to rows of \`Y\` (values in \`1..N\`).

- preproc:

  A \`multivarious\` preprocessing pipeline (a
  \`pre_processor\`/\`prepper\`) or a list of them. If a list, it must
  have length \`1 + length(X)\` and will be applied to \`c(list(Y), X)\`
  in that order.

- ncomp:

  Integer number of components to extract.

- normalization:

  Block weighting scheme. \`"MFA"\` uses inverse squared first singular
  value per block; \`"None"\` uses uniform weights; \`"custom"\` uses
  \`alpha\`.

- alpha:

  Optional numeric vector of per-block weights (length \`1 +
  length(X)\`), used when \`normalization = "custom"\`. The first weight
  corresponds to \`Y\`.

- feature_groups:

  Feature prior specification. One of: \* \`NULL\` (no feature prior),
  \* \`"colnames"\` to group X-features with identical column names
  across blocks, \* a \`data.frame\` with columns \`block\`,
  \`feature\`, \`group\` and optional \`weight\`. \`block\` refers to a
  name or index in \`X\` (not including \`Y\`), and \`feature\` is a
  column name or index within that block.

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
additional classes \`"anchored_mfa"\` and \`"linked_mfa"\`. The object
contains global scores in \`s\`, concatenated loadings in \`v\`, and
block mappings in \`block_indices\`. Additional fields include
\`V_list\`, \`B\`, \`row_index\`, \`alpha_blocks\`, and
\`objective_trace\`.

## Details

\## Model The fitted model has the form: \$\$Y \approx S B^\top\$\$
\$\$X_k \approx S\[\mathrm{idx}\_k,\] V_k^\top\$\$ where \`S\` is \`N ×
ncomp\`, \`B\` is \`q × ncomp\`, and each \`V_k\` is \`p_k × ncomp\`.

\## Feature similarity prior (v1) When \`feature_lambda \> 0\` and
\`feature_groups\` is supplied, Anchored MFA applies a group-shrinkage
penalty that pulls the loading vectors of features in the same group
toward a shared group center.

\`linked_mfa()\` is a legacy alias for \[anchored_mfa()\] retained for
backward compatibility.

## Examples

``` r
# \donttest{
set.seed(1)
N <- 30
Y <- matrix(rnorm(N * 5), N, 5)
X1 <- matrix(rnorm(20 * 10), 20, 10)
X2 <- matrix(rnorm(15 * 8), 15, 8)
idx1 <- sample.int(N, nrow(X1), replace = FALSE)
idx2 <- sample.int(N, nrow(X2), replace = FALSE)

fit <- anchored_mfa(Y, list(X1 = X1, X2 = X2), list(X1 = idx1, X2 = idx2), ncomp = 2)
#> Applying the same preprocessor definition independently to each block.
stopifnot(nrow(multivarious::scores(fit)) == N)
# }
```

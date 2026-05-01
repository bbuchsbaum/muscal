# Response-Aligned Multiple Factor Analysis

\`response_aligned_mfa()\` learns a shared response-informed latent
space for a set of multiblock predictors \`X\` paired with blockwise
multivariate responses \`Y\`, without assuming that rows of \`Y\` are
shared across blocks.

## Usage

``` r
response_aligned_mfa(
  Y,
  X,
  preproc = multivarious::center(),
  response_preproc = multivarious::center(),
  ncomp = 2,
  normalization = c("MFA", "None", "custom"),
  alpha = NULL,
  response_alpha = 1,
  response_weights = NULL,
  anchor_response = NULL,
  anchor_response_alpha = 1,
  anchor_response_weights = NULL,
  anchor_map = NULL,
  anchor_weight = NULL,
  coupling_lambda = 0,
  feature_groups = NULL,
  feature_graph = NULL,
  feature_lambda = 0,
  graph_lambda = 0,
  graph_form = c("laplacian", "adjacency", "normalized_laplacian"),
  max_iter = 50,
  tol = 1e-06,
  ridge = 1e-08,
  verbose = FALSE,
  use_future = FALSE,
  ...
)
```

## Arguments

- Y:

  A list of numeric matrices/data.frames. Each element \`Y\[\[k\]\]\` is
  \`n_k x q\` and all blocks must share the same response column
  dimension.

- X:

  A list of numeric matrices/data.frames. Each element \`X\[\[k\]\]\` is
  \`n_k x p_k\`.

- preproc:

  A \`multivarious\` preprocessing pipeline (a \`pre_processor\` or
  \`prepper\`) or a list of them for the \`X\` blocks.

- response_preproc:

  A single shared \`multivarious\` preprocessing pipeline used to define
  the common response space.

- ncomp:

  Integer number of latent components.

- normalization:

  Block weighting scheme for the \`X\` blocks. \`"MFA"\` uses inverse
  squared first singular values; \`"None"\` uses uniform weights;
  \`"custom"\` uses \`alpha\`.

- alpha:

  Optional numeric vector of per-block \`X\` weights used when
  \`normalization = "custom"\`.

- response_alpha:

  Optional numeric scalar or vector of per-block response weights.

- response_weights:

  Optional rowwise response weights. May be \`NULL\`, a single
  non-negative scalar, a length-\`length(X)\` vector of blockwise
  constants, or a list mirroring \`X\` with one non-negative weight
  vector per block row.

- anchor_response:

  Optional anchor-level multivariate response table with one row per
  anchor state. Requires active \`anchor_map\` information and must have
  the same number of columns as the blockwise responses.

- anchor_response_alpha:

  Optional non-negative scalar controlling the contribution of
  \`anchor_response\` to the shared response loading update.

- anchor_response_weights:

  Optional rowwise weights for \`anchor_response\`. May be \`NULL\`, a
  scalar, or a non-negative vector of length \`nrow(anchor_response)\`.

- anchor_map:

  Optional anchor structure. May be \`NULL\`, a list of integer vectors
  whose values identify anchor states and whose \`NA\` values indicate
  unanchored rows, or a list of non-negative row-stochastic matrices
  encoding soft anchor assignments.

- anchor_weight:

  Optional rowwise anchor weights. May be \`NULL\`, a scalar, a
  length-\`length(X)\` vector of blockwise constants, or a list
  mirroring \`X\` with one non-negative weight vector per block row.

- coupling_lambda:

  Optional non-negative scalar or vector of per-block anchor-coupling
  strengths.

- feature_groups:

  Feature prior specification for the \`X\` loadings. One of \`NULL\`,
  \`"colnames"\`, or a \`data.frame\` with columns \`block\`,
  \`feature\`, \`group\`, and optional \`weight\`.

- feature_graph:

  Optional graph specification for the \`X\` loadings. One of \`NULL\`,
  \`"colnames"\`, a \`data.frame\` with columns \`block1\`,
  \`feature1\`, \`block2\`, \`feature2\` and optional \`weight\`, or a
  square matrix interpreted according to \`graph_form\`. This is an
  alternative to \`feature_groups\`, not an additional penalty.

- feature_lambda:

  Non-negative scalar controlling the strength of the feature-group
  shrinkage.

- graph_lambda:

  Non-negative scalar controlling the strength of the feature-graph
  penalty.

- graph_form:

  Interpretation of \`feature_graph\` when it is matrix-like, or the
  Laplacian construction used for edge-based inputs.

- max_iter:

  Maximum ALS iterations.

- tol:

  Relative convergence tolerance on the objective.

- ridge:

  Non-negative ridge stabilization applied to scores and loadings.

- verbose:

  Logical; if \`TRUE\`, prints iteration diagnostics.

- use_future:

  Logical; if \`TRUE\`, block-wise computations are performed via
  \`furrr::future_map()\` when available.

- ...:

  Unused (reserved for future extensions).

## Value

An object inheriting from \`multivarious::multiblock_biprojector\` with
additional class \`"response_aligned_mfa"\`. The object stores
block-specific scores in \`Z_list\`, shared response loadings in \`B\`,
block loadings in \`V_list\`, rowwise response weights in
\`response_weights\`, optional anchor scores in \`S\`, optional
anchor-level responses in \`anchor_response\`, parsed anchor maps in
\`anchor_map\`, rowwise anchor weights in \`anchor_weight\`, optional
\`anchor_response_fit\` diagnostics, pooled score normalization weights
in \`score_weights\`, and optional feature-graph metadata in
\`graph_laplacian\`, \`graph_adjacency\`, \`graph_lambda\`, and
\`graph_form\`. For compatibility with generic score accessors, the
\`s\` slot stores a concatenated stack of \`Z_list\`, \`score_index\`
maps those rows back to blocks, and \`score_representation =
"stacked_block_scores"\` records that this is not a shared
observation-level score table.

## Details

The fitted model minimizes \$\$ \sum_k \alpha_k \\X_k - Z_k
V_k^\top\\\_F^2 + \sum_k \eta_k \\R_k^{1/2}(Y_k - Z_k B^\top)\\\_F^2 +
\eta_0 \\R_0^{1/2}(Y^{(a)} - S B^\top)\\\_F^2 + \sum_k \mu_k
\\M_k^{1/2}(Z_k - A_k S)\\\_F^2 \$\$ plus optional feature-side coupling
on the rows of the block loading matrices \`V_k\` and ridge penalties on
\`Z_k\`, \`V_k\`, \`B\`, and the optional anchor score matrix \`S\`. In
v1, feature coupling may be supplied either through
\`feature_groups\`/\`feature_lambda\` or through
\`feature_graph\`/\`graph_lambda\`, but not both simultaneously.

Predictor blocks are preprocessed blockwise, while the response blocks
are preprocessed in a shared pooled response space that also includes
\`anchor_response\` when present. Identifiability is enforced on the
score side by orthonormalizing a weighted pooled stack of block-specific
scores after each iteration. A deterministic sign convention is then
applied componentwise so that refits do not differ only by arbitrary
sign flips.

By default, out-of-sample prediction is pure with respect to the
response target: \`x_new -\> z_hat -\> y_hat\`, optionally refined by
\`new_anchor_map\` when test rows have anchor-state information.
Supplying \`new_response\` requires \`conditional = TRUE\` in
\[project()\] or \[predict()\] and switches to explicit conditional
completion rather than pure prediction.

## Examples

``` r
# \donttest{
set.seed(1)
X1 <- matrix(rnorm(30 * 6), 30, 6)
X2 <- matrix(rnorm(28 * 5), 28, 5)
Y1 <- matrix(rnorm(30 * 3), 30, 3)
Y2 <- matrix(rnorm(28 * 3), 28, 3)

fit <- response_aligned_mfa(
  Y = list(X1 = Y1, X2 = Y2),
  X = list(X1 = X1, X2 = X2),
  ncomp = 2
)
#> Applying the same preprocessor definition independently to each block.

pred <- predict(fit, X1[1:5, , drop = FALSE], block = "X1", type = "response")
stopifnot(nrow(pred) == 5)
# }
```

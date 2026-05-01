# Coupled Graph-Anchored Multiple Factor Analysis

\`coupled_graph_anchored_mfa()\` extends \[graph_anchored_mfa()\] by
introducing block-specific row scores \`Z_k\` that are softly coupled to
a shared stimulus-level score matrix \`S\`. This is useful when each
auxiliary block is expected to express the common stimulus structure
with subject- or domain-specific deviations.

## Usage

``` r
coupled_graph_anchored_mfa(
  Y,
  X,
  row_index,
  block_info = NULL,
  preproc = multivarious::center(),
  ncomp = 2,
  normalization = c("MFA", "None", "custom"),
  alpha = NULL,
  score_constraint = c("none", "orthonormal"),
  feature_graph = NULL,
  graph_lambda = 0,
  graph_form = c("laplacian", "adjacency", "normalized_laplacian"),
  score_graph = NULL,
  score_graph_lambda = 0,
  score_graph_form = c("laplacian", "adjacency", "normalized_laplacian"),
  score_graph_k = 10,
  score_graph_weight_mode = c("heat", "binary"),
  score_graph_sigma = NULL,
  coupling_lambda = 1,
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

  Numeric matrix/data.frame (\`N × q\`) serving as the anchored target
  space.

- X:

  Auxiliary blocks. Either a flat named list of matrices/data frames, or
  a nested list \`X\[\[subject\]\]\[\[domain\]\]\`.

- row_index:

  A structure mirroring \`X\`. Each vector maps rows of the
  corresponding auxiliary block to rows of \`Y\`.

- block_info:

  Optional data frame describing flattened blocks. If supplied, it must
  have one row per flattened block. Recommended columns are \`block\`,
  \`subject\`, and \`domain\`.

- preproc:

  A \`multivarious\` preprocessing pipeline (a \`pre_processor\` or
  \`prepper\`) or a list of them. If a list, it must have length \`1 +
  length(flattened_X)\` and will be applied to \`c(list(Y),
  flattened_X)\`.

- ncomp:

  Integer number of latent components.

- normalization:

  Block weighting scheme. \`"MFA"\` uses inverse squared first singular
  values; \`"None"\` uses uniform weights; \`"custom"\` uses \`alpha\`.

- alpha:

  Optional numeric vector of per-block weights. When \`normalization =
  "custom"\`, it must have length \`1 + length(flattened_X)\`, with the
  first weight corresponding to \`Y\`.

- score_constraint:

  Identification strategy for the anchored score matrix. \`"none"\` uses
  the historical unconstrained update followed by QR normalization
  inside each ALS iteration. \`"orthonormal"\` enforces \`S transpose S
  = I\` with a constrained majorization/polar update.

- feature_graph:

  Feature-graph specification; see Details.

- graph_lambda:

  Non-negative scalar controlling graph-penalty strength.

- graph_form:

  Interpretation of \`feature_graph\` when it is matrix-like, or the
  Laplacian construction used for edge-based inputs.

- score_graph:

  Optional score-graph specification; see Details.

- score_graph_lambda:

  Non-negative scalar controlling row/score-graph smoothing strength.

- score_graph_form:

  Interpretation of \`score_graph\` when it is matrix-like, or the
  Laplacian construction used for edge-based inputs.

- score_graph_k:

  Integer number of neighbors used when \`score_graph = "knn"\`.

- score_graph_weight_mode:

  Weighting applied when \`score_graph = "knn"\`. \`"heat"\` uses a
  Gaussian similarity kernel on neighbor distances; \`"binary"\` assigns
  weight 1 to every retained neighbor edge.

- score_graph_sigma:

  Optional positive bandwidth used when \`score_graph = "knn"\` and
  \`score_graph_weight_mode = "heat"\`. If \`NULL\`, a robust value is
  estimated from the k-th neighbor distances.

- coupling_lambda:

  Non-negative scalar controlling the strength of the consensus penalty
  tying each block-specific score matrix \`Z_k\` to the shared anchor
  scores \`S\[row_index\[\[k\]\], \]\`.

- max_iter:

  Maximum ALS iterations.

- tol:

  Relative convergence tolerance on the objective.

- ridge:

  Non-negative ridge stabilization applied to loading and score updates.

- verbose:

  Logical; if \`TRUE\`, prints iteration diagnostics.

- use_future:

  Logical; if \`TRUE\`, block-wise computations that do not depend on
  one another are performed via \`furrr::future_map()\` when available.

- ...:

  Unused (reserved for future extensions).

## Value

An object inheriting from \`multivarious::multiblock_biprojector\` with
additional classes \`coupled_graph_anchored_mfa\` and \`linked_mfa\`.
The object contains shared anchor scores in \`s\`, block-specific row
scores in \`Z_list\`, auxiliary loadings in \`V_list\`, and anchor
loadings in \`B\`.

## Details

The fitted model has the form: \$\$Y \approx S B^\top\$\$ \$\$X_k
\approx Z_k V_k^\top\$\$ with a consensus penalty \$\$\mu \sum_k \\Z_k -
S\[\mathrm{idx}\_k,\]\\\_F^2.\$\$

Feature-side and score-side graph penalties are inherited from
\[graph_anchored_mfa()\]: \$\$\lambda_G \mathrm{tr}(V^\top L_G V) +
\lambda_S \mathrm{tr}(S^\top L_S S).\$\$ As in \[graph_anchored_mfa()\],
\`score_constraint = "none"\` uses the historical unconstrained/QR score
update, while \`score_constraint = "orthonormal"\` enforces \`S
transpose S = I\` directly during fitting.

When \`coupling_lambda\` is large, block-specific scores are strongly
tied to the shared anchor score space; when it is small, blocks can
deviate more freely while still being linked through the common response
target \`Y\`.

# Aligned Multiblock Canonical Correlation Analysis (Aligned MCCA)

Aligned MCCA extends MAXVAR generalized CCA (GCCA) to the setting where
each block has its own row set, linked to a common set of latent
reference rows via an integer index vector. No block is privileged: the
shared scores live in the reference row space \`1..N\`, and each block
contributes through its row mapping.

This is the correlation-based counterpart to \[aligned_mfa()\],
analogous to how \[mcca()\] relates to \[mfa()\] when all blocks share
rows.

## Usage

``` r
aligned_mcca(
  X,
  row_index,
  N = NULL,
  preproc = multivarious::center(),
  ncomp = 2,
  normalization = c("MFA", "balanced", "None", "custom"),
  alpha = NULL,
  ridge = 1e-06,
  max_iter = 50,
  tol = 1e-06,
  verbose = FALSE,
  use_future = FALSE,
  block_weights = NULL,
  ...
)
```

## Arguments

- X:

  A list of numeric matrices/data.frames. Each element \`X\[\[k\]\]\` is
  \`n_k × p_k\`.

- row_index:

  A list of integer vectors. \`row_index\[\[k\]\]\` has length \`n_k\`
  and maps rows of \`X\[\[k\]\]\` to shared rows in \`1..N\`. Repeated
  indices are allowed.

- N:

  Optional integer specifying the number of shared rows. If \`NULL\`
  (default), \`N\` is inferred as \`max(unlist(row_index))\`.

- preproc:

  A \`multivarious\` preprocessing pipeline (a
  \`pre_processor\`/\`prepper\`) or a list of them. If a list, it must
  have length \`length(X)\` and will be applied to \`X\` in that order.

- ncomp:

  Integer; number of canonical components to compute.

- normalization:

  Block weighting scheme. One of \`"MFA"\` (default), \`"balanced"\`,
  \`"None"\`, or \`"custom"\`. See "Block-weighting schemes".

- alpha:

  Optional numeric vector of non-negative weights of length
  \`length(X)\`. Required for \`normalization = "custom"\`; used as the
  IRLS target when \`normalization = "balanced"\`.

- ridge:

  Non-negative numeric scalar (or vector of length \`length(X)\`)
  controlling ridge stabilisation. The effective ridge used per block is
  scaled by the block's overall energy.

- max_iter:

  Maximum number of iterations for \`normalization = "balanced"\`.
  Ignored otherwise.

- tol:

  Relative tolerance on the gradient norm / objective change used to
  declare IRLS convergence for \`normalization = "balanced"\`.

- verbose:

  Logical; if \`TRUE\`, prints IRLS iteration progress.

- use_future:

  Logical; if \`TRUE\`, block-wise projector fitting is performed via
  \`furrr::future_map()\` when available.

- block_weights:

  Deprecated alias for \`alpha\`. When supplied, implicitly sets
  \`normalization = "custom"\` unless \`normalization\` is already
  explicit.

- ...:

  Unused (reserved for future extensions).

## Value

An object inheriting from \`multivarious::multiblock_biprojector\` with
additional class \`"aligned_mcca"\`. Relevant fields include \`s\` (\`N
× ncomp\` reference-space scores), \`v\` (concatenated canonical
directions), \`sdev\`, \`block_indices\`, \`alpha_blocks\` (block
weights used), \`block_weights\` (alias), \`normalization\`,
\`canonical_weights\`, \`partial_scores\`, \`cor_loadings\` /
\`scaled_loadings\` (feature-score correlations),
\`scaled_loadings_by_block\`, \`block_contribs\` (\`ncomp × B\` matrix
of \`s_k^T M_b s_k\`), and (for \`"balanced"\`) \`alpha_per_component\`,
\`balance_trace\`, \`balance_converged\`, \`balance_iters\`.

## Details

\## Model (MAXVAR GCCA with row alignment) Let \`S\` be the shared score
matrix for the \`N\` reference rows. For each block \`X_k\` with \`n_k\`
rows, \`row_index\[\[k\]\]\` maps its rows to reference rows: \$\$X_k \\
\mathrm{linked\\ to}\\ S\[\mathrm{idx}\_k, \]\$\$ Aligned MCCA computes
\`S\` as the leading eigenvectors of a weighted sum of block-wise ridge
projection operators lifted into the reference space: \$\$H = \sum_b
alpha_b M_b\$\$ with \\M_b = lift(P_b)\\ and \\P_b = X_b (X_b^T X_b +
kappa_b I)^{-1} X_b^T\\.

\## Block-weighting schemes The block-weight vector \`alpha_blocks\` is
formed according to \`normalization\`: \* \`"MFA"\` (default) —
\`alpha_b = 1 / (sigma_1(X_b)^2)\` (the CCA analogue of the MFA
normalisation; prevents a block from dominating purely by scale). \*
\`"balanced"\` — per-component projected gradient ascent on the
geometric mean criterion \`Σ_b β_b log(g^T M_b g)\` so every block
contributes meaningfully to every component (uses \`target = alpha\` if
supplied, else uniform). \* \`"None"\` — uniform weights \`alpha_b =
1\`. \* \`"custom"\` — use \`alpha\` as supplied (length \`length(X)\`).

\## High-dimensional stability (p \> n) The implementation works in
observation space and uses ridge-stabilised solves on \`n_b × n_b\`
systems, so \`p_b \>\> n_b\` blocks are handled without special casing.
Ridge stabilisation (\`kappa_b = ridge_b \* mean(diag(K_b))\`, with an
automatic ridge floor if \`ridge = 0\` yields a singular system) keeps
the Cholesky well-posed.

## Examples

``` r
# \donttest{
set.seed(1)
N <- 30
X1 <- matrix(rnorm(20 * 10), 20, 10)
X2 <- matrix(rnorm(15 * 8), 15, 8)
idx1 <- sample.int(N, nrow(X1), replace = TRUE)
idx2 <- sample.int(N, nrow(X2), replace = TRUE)

fit <- aligned_mcca(list(X1 = X1, X2 = X2), list(X1 = idx1, X2 = idx2), N = N, ncomp = 2)
#> Applying the same preprocessor definition independently to each block.
stopifnot(nrow(multivarious::scores(fit)) == N)
# }
```

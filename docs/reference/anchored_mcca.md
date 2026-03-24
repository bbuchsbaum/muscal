# Anchored Multiblock Canonical Correlation Analysis (Anchored MCCA)

Convenience wrapper around \[aligned_mcca()\] for the common case with a
fully observed reference block \`Y\` (e.g., stimulus features) and one
or more linked blocks \`X_k\` (e.g., brain data). Rows of each \`X_k\`
are linked to rows of \`Y\` via \`row_index\[\[k\]\]\`.

Compared to the legacy name "linked", "anchored" makes explicit that the
shared scores live in the row space of \`Y\` (with \`N = nrow(Y)\`).

## Usage

``` r
anchored_mcca(
  Y,
  X,
  row_index,
  preproc = multivarious::center(),
  ncomp = 2,
  ridge = 1e-06,
  block_weights = NULL,
  use_future = FALSE,
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

  Integer; number of canonical components to compute.

- ridge:

  Non-negative numeric scalar (or vector of length equal to the number
  of blocks) controlling ridge stabilization. The effective ridge used
  per block is scaled by the block's overall energy so the default works
  across variable scalings.

- block_weights:

  Optional numeric vector of non-negative weights. If unnamed and of
  length \`length(X)\`, it is interpreted as weights for \`X\` blocks
  and \`Y\` is given weight 1. If of length \`1 + length(X)\`, the first
  weight corresponds to \`Y\`.

- use_future:

  Logical; if \`TRUE\`, block-wise computations are performed via
  \`furrr::future_map()\` when available.

- ...:

  Additional arguments forwarded to \[aligned_mcca()\] (currently
  unused).

## Value

An object of class \`"anchored_mcca"\` (and \`"aligned_mcca"\`).

## Examples

``` r
# \donttest{
set.seed(1)
N <- 30
Y <- matrix(rnorm(N * 5), N, 5)
X1 <- matrix(rnorm(20 * 10), 20, 10)
X2 <- matrix(rnorm(15 * 8), 15, 8)
idx1 <- sample.int(N, nrow(X1), replace = TRUE)
idx2 <- sample.int(N, nrow(X2), replace = TRUE)

fit <- anchored_mcca(Y, list(X1 = X1, X2 = X2), list(X1 = idx1, X2 = idx2), ncomp = 2)
#> Applying the same preprocessor definition independently to each block.
stopifnot(nrow(multivarious::scores(fit)) == N)
# }
```

# Barycentric Multiple Factor Analysis (BaMFA)

Performs Barycentric Multiple Factor Analysis (BaMFA) using an
alternating optimization approach based on a two-level factor model. It
decomposes multi-block data (e.g., multiple subjects) into a shared
global subspace (`G`) and block-specific subspaces (`B_i`).

## Usage

``` r
bamfa(
  data,
  k_g = 2,
  k_l = 2,
  niter = 10,
  preproc = multivarious::center(),
  lambda_l = 0,
  tol = 1e-05,
  subject = NULL,
  ...
)

# Default S3 method
bamfa(
  data,
  k_g = 2,
  k_l = 2,
  niter = 10,
  preproc = multivarious::center(),
  lambda_l = 0,
  tol = 1e-05,
  ...
)

# S3 method for class 'list'
bamfa(
  data,
  k_g = 2,
  k_l = 2,
  niter = 10,
  preproc = multivarious::center(),
  lambda_l = 0,
  tol = 1e-05,
  ...
)

# S3 method for class 'multiblock'
bamfa(
  data,
  k_g = 2,
  k_l = 2,
  niter = 10,
  preproc = multivarious::center(),
  lambda_l = 0,
  tol = 1e-05,
  ...
)

# S3 method for class 'multidesign'
bamfa(
  data,
  k_g = 2,
  k_l = 2,
  niter = 10,
  preproc = multivarious::center(),
  lambda_l = 0,
  tol = 1e-05,
  subject,
  ...
)
```

## Arguments

- data:

  A `list` of matrices (each a block, e.g., subject), a `multiblock`
  object, or a `multidesign` object. If a list or multiblock, all blocks
  must have the same number of rows (observations, e.g., conditions). If
  a multidesign, the `subject` parameter must be specified to indicate
  how to split the data.

- k_g:

  Integer: Number of global components to extract (shared across
  blocks).

- k_l:

  Integer: Number of local components to extract (specific to each
  block).

- niter:

  Integer: Maximum number of iterations (default: 10).

- preproc:

  A preprocessing pipeline object from `multivarious` (default:
  [`multivarious::center()`](https://bbuchsbaum.github.io/multivarious/reference/center.html)).

- lambda_l:

  Numeric: Non-negative ridge penalty applied when estimating the local
  scores `U_i` using the internal `ls_ridge` function. Default: 0 (no
  penalty).

- tol:

  Numeric: Convergence tolerance. Stops if the relative change in the
  objective function (Mean Squared Error) and/or the global loadings `G`
  is less than `tol`. Default: 1e-5.

- subject:

  Optional: A variable name identifying the blocking/subject variable
  when using a multidesign object. Required only for the multidesign
  method.

- ...:

  Additional arguments (currently unused).

## Value

A `multiblock_projector` object with class `"bamfa"`. The base projector
stores the global loading matrix in `v`, the concatenated preprocessor
in `preproc`, and block mappings in `block_indices`. Additional named
elements passed via `...` include:

- `B_list` – block-specific local loading matrices.

- `S_list` – block-specific global score matrices.

- `U_list` – block-specific local score matrices.

- `k_g` – number of global components in the final model.

- `k_l` – requested number of local components.

- `lambda_l` – regularization parameter used for local scores `U_i`.

- `niter_actual` – actual number of iterations performed.

- `objective_trace` – objective function value at each iteration.

- `g_change_trace` – relative change in global loadings at each
  iteration.

- `data_names` – names of the input blocks/subjects.

- `proclist_fitted` – list of fitted preprocessors for each block.

## Details

The algorithm models each data block \\X_i\\ as: \$\$X_i = S_i G^T + U_i
B_i^T + E_i\$\$ where \\G\\ represents shared global loadings, \\B_i\\
represents block-specific local loadings, \\S_i\\ and \\U_i\\ are the
corresponding scores, and \\E_i\\ is noise. Loadings are constrained to
be orthonormal (\\G^T G = I, B_i^T B_i = I\\).

The algorithm aims to minimize the total reconstruction error:
\$\$\sum\_{i=1}^{m} \|\|X_i - S_i G^T - U_i B_i^T\|\|\_F^2\$\$ using an
iterative alternating optimization strategy (similar to
Expectation-Maximization):

1.  **Preprocessing:** Each block is preprocessed using the provided
    `preproc` pipeline.

2.  **Initialization:** Initialize global loadings `G` (via SVD on the
    mean block).

3.  **Iterate (E-step like):** For each block `i`, holding `G` fixed,
    update scores `S_i` (global projection), calculate residuals, find
    the local basis `B_i` from residuals (via SVD), and estimate local
    scores `U_i` (via projection, potentially regularized by
    `lambda_l`).

4.  **Iterate (M-step like):** Holding scores `S_i` fixed, update the
    global loadings `G` by performing SVD on an aggregated cross-product
    matrix `sum(X_i^T S_i)`.

5.  **Repeat steps 3-4** for `niter` iterations or until convergence
    based on `tol`.

## Caveats and Limitations

- **Model Choice:** Assumes a linear factor model with orthogonal global
  and local components.

- **Initialization:** Results can be sensitive to the initialization of
  `G`.

- **Local Minima:** The alternating optimization algorithm may converge
  to a local minimum.

- **Interpretation:** The separation into global and local components
  depends on the model fit and ranks (`k_g`, `k_l`).

- **Regularization (`lambda_l`):** Penalizes the squared Frobenius norm
  of the local scores `U_i` via ridge regression.

- **Inference:** The method does not provide p-values or confidence
  intervals.

## See also

[`pca`](https://bbuchsbaum.github.io/multivarious/reference/pca.html),
[`multiblock`](https://rdrr.io/pkg/multidesign/man/multiblock.html)

## Examples

``` r
# Generate example multi-block data (e.g., 3 subjects, 10 conditions, 50 features)
set.seed(123)
n_obs <- 10
n_features <- 50
n_subjects <- 3
data_list <- lapply(1:n_subjects, function(i) {
  matrix(rnorm(n_obs * n_features), n_obs, n_features) +
  matrix(rnorm(n_obs * 1, mean=i), n_obs, n_features) # Add subject offset
})
names(data_list) <- paste0("Subject_", 1:n_subjects)

# Run BaMFA with k_g=3 global, k_l=2 local components
result <- bamfa(data_list, k_g = 3, k_l = 2, niter = 10)
print(result)
#> Barycentric Multiple Factor Analysis (BaMFA) 
#> 
#> Model Structure: 
#>   Global components (k_g): 3 
#>   Local components (k_l requested): 2 
#>   Convergence after 10 iterations
#>   Local score regularization (lambda_l): 0 
#> 
#> Block Information: 
#>   Block 1 ( Subject_1 ): 
#>     Features: 50 
#>     Local components retained: 2 
#>   Block 2 ( Subject_2 ): 
#>     Features: 50 
#>     Local components retained: 2 
#>   Block 3 ( Subject_3 ): 
#>     Features: 50 
#>     Local components retained: 2 
#> 
#> Final reconstruction error (per feature): 0.458139 
```

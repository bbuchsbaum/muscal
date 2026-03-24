# Generate synthetic multiblock data

This helper function creates several blocks of multivariate data that
share a common set of latent factor scores. Optionally the variables can
be placed on the unit sphere to yield spatial coordinates and a sparse
k-nearest-neighbour graph. The returned object also contains the
ground-truth loadings and scores used for simulation.

## Usage

``` r
synthetic_multiblock(
  S = 5,
  n = 100,
  p = 200,
  r = 3,
  sigma = 0.1,
  sphere = FALSE,
  k_nn = 6,
  seed = 1
)
```

## Arguments

- S:

  Number of subjects/blocks to generate (default 5).

- n:

  Number of observations (rows) per block (default 100).

- p:

  Number of variables (columns) per block, or a vector of length S
  specifying different dimensions per block (default 200).

- r:

  The rank of the shared component structure (default 3).

- sigma:

  The standard deviation of the noise added to the data (default 0.1).

- sphere:

  Logical; if TRUE, variables are placed on the unit sphere and a
  k-nearest-neighbour graph is computed (default FALSE).

- k_nn:

  Number of nearest neighbors for spatial correlations when \`sphere =
  TRUE\` (default 6).

- seed:

  Random seed for reproducibility (default 1).

## Value

A list containing:

- data_list:

  A list of data matrices, one per block.

- coords_list:

  Coordinates on the unit sphere (if \`sphere = TRUE\`), otherwise NULL.

- V_true:

  List of ground-truth loading matrices.

- F_true:

  Ground-truth factor score matrix.

- Sadj:

  Spatial adjacency Laplacian matrix (if \`sphere = TRUE\`), otherwise
  NULL.

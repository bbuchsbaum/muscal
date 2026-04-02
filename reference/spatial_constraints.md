# Simple spatial constraints function (creates k-NN adjacency matrix)

Simple spatial constraints function (creates k-NN adjacency matrix)

## Usage

``` r
spatial_constraints(coords_list, k_nn = 6, nblocks = NULL, ...)
```

## Arguments

- coords_list:

  List of coordinate matrices (each k_s x 3)

- k_nn:

  Number of nearest neighbors (default 6)

- nblocks:

  Number of blocks (ignored, for compatibility)

## Value

Sparse adjacency matrix

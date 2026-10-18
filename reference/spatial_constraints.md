# Pooled-feature union k-NN adjacency

Pooled-feature union k-NN adjacency

## Usage

``` r
spatial_constraints(coords_list, k_nn = 6, nblocks = NULL, ...)
```

## Arguments

- coords_list:

  List of coordinate matrices (each k_s x 3)

- k_nn:

  Finite nonnegative integer neighbor count (default 6), capped at
  available non-self vertices. Zero creates no edges. All boundary ties
  within floating-point tolerance are included, so actual degree may
  exceed k_nn.

- nblocks:

  Number of blocks (ignored, for compatibility)

## Value

Sparse adjacency matrix

## Details

Features are pooled across the supplied coordinate blocks. Directed
neighbors are symmetrized by union with unit edge weights. Self removal
uses vertex identity, including duplicate coordinates. Boundary ties are
included regardless of input order. The graph does not distinguish
within-block smoothing from between-block correspondence; that
scientific distinction remains part of the broader experimental recovery
plan.

# Build a cluster graph from spatial coordinates

Build a cluster graph from spatial coordinates

## Usage

``` r
cluster_graph(coords_list, adjacency_opts = list(), return_graph = FALSE)
```

## Arguments

- coords_list:

  List of coordinate matrices (each k_s x 3)

- adjacency_opts:

  Options for adjacency construction

- return_graph:

  Return adjacency, adjacency degrees, and combinatorial Laplacian
  together instead of only the Laplacian (default FALSE).

## Value

A sparse combinatorial Laplacian, or graph data if requested.

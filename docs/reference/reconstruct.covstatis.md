# Low-rank reconstruction of the compromise matrix

Recreates the compromise matrix using only the selected components.

## Usage

``` r
# S3 method for class 'covstatis'
reconstruct(x, comp = 1:multivarious::ncomp(x), ...)
```

## Arguments

- x:

  A fitted \`covstatis\` model

- comp:

  Integer vector of component indices to retain

- ...:

  Additional arguments (currently unused).

## Value

A matrix of the same size as each input matrix, representing a low-rank
approximation of the compromise matrix in the preprocessed space (after
double-centering and normalization, if they were applied)

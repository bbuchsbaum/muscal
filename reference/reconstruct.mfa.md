# Reconstruct Data from an MFA Fit

MFA loadings are expressed in the block-weighted column metric
(\`components(x) = A \`scores weights. These methods back-project with
the metric-orthonormal axes \`ov\` instead, giving the
rank-\`length(comp)\` generalized SVD approximation of the preprocessed
data, and then undo the preprocessing.

## Usage

``` r
# S3 method for class 'mfa'
reconstruct(
  x,
  comp = 1:multivarious::ncomp(x),
  rowind = 1:nrow(multivarious::scores(x)),
  colind = 1:nrow(multivarious::components(x)),
  ...
)

# S3 method for class 'mfa'
reconstruct_new(
  x,
  new_data,
  comp = 1:multivarious::ncomp(x),
  colind = 1:nrow(multivarious::components(x)),
  ...
)
```

## Arguments

- x:

  A fitted \`mfa\` object.

- comp:

  Integer vector of components to use.

- rowind:

  Training rows to reconstruct.

- colind:

  Concatenated variables (columns) to reconstruct.

- ...:

  Unused.

- new_data:

  Matrix of new observations: all concatenated variables, or only the
  variables in \`colind\`.

## Value

A numeric matrix on the original data scale.

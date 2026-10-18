# Report Reference-Feature Alignment Coverage

Report Reference-Feature Alignment Coverage

## Usage

``` r
alignment_coverage(x)
```

## Arguments

- x:

  A fitted \`design_mfa\` object.

## Value

A data frame with the number and fraction of subjects supporting each
reference feature and the aggregate sparse-map weight. Shared alignment
reports complete support.

## Examples

``` r
example_data <- design_mfa_example_data()
identity_maps <- lapply(example_data$Y, function(y) {
  P <- Matrix::Diagonal(ncol(y))
  dimnames(P) <- list(colnames(y), colnames(y))
  P
})
names(identity_maps) <- names(example_data$Y)
soft_fit <- design_mfa(
  example_data$Y, example_data$design, ~ saliency,
  alignment = "soft", feature_map = identity_maps
)
head(alignment_coverage(soft_fit))
#>   feature subject_support coverage_fraction map_weight
#> 1 voxel_1               6                 1          6
#> 2 voxel_2               6                 1          6
#> 3 voxel_3               6                 1          6
#> 4 voxel_4               6                 1          6
#> 5 voxel_5               6                 1          6
#> 6 voxel_6               6                 1          6
```

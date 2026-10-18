# Build Coordinate Maps for Several Subjects

Build Coordinate Maps for Several Subjects

## Usage

``` r
coordinate_feature_maps(
  reference_coords,
  subject_coords,
  require_full_reference = TRUE,
  ...
)
```

## Arguments

- reference_coords:

  Numeric reference-feature-by-coordinate matrix.

- subject_coords:

  A list of subject coordinate matrices.

- require_full_reference:

  Require every reference feature to receive support from at least one
  subject map.

- ...:

  Arguments passed to \[coordinate_feature_map()\].

## Value

A named list of sparse feature maps. Its \`coverage\` attribute reports
subject support and aggregate map weight for every reference feature.

## Examples

``` r
reference <- cbind(x = 0:3, y = c(0, 0, 1, 1), z = 0)
rownames(reference) <- paste0("r", 1:4)
subjects <- list(
  S1 = structure(reference + 0.03, dimnames = list(paste0("S1_", 1:4), colnames(reference))),
  S2 = structure(reference - 0.03, dimnames = list(paste0("S2_", 1:4), colnames(reference)))
)
maps <- coordinate_feature_maps(subject_coords = subjects,
                                reference_coords = reference,
                                k = 2, bandwidth = 0.2)
attr(maps, "coverage")
#>   feature subject_support coverage_fraction map_weight
#> 1      r1               2                 1          2
#> 2      r2               2                 1          2
#> 3      r3               2                 1          2
#> 4      r4               2                 1          2
```

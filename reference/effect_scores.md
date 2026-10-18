# Extract Effect Scores or Spatial Loadings

Extract Effect Scores or Spatial Loadings

## Usage

``` r
effect_scores(x, term, components = NULL)

effect_loadings(x, term, components = NULL)
```

## Arguments

- x:

  A fitted \`design_mfa\` object.

- term:

  A fitted design term.

- components:

  Optional component indices.

## Value

A numeric matrix.

## Examples

``` r
example_data <- design_mfa_example_data()
example_fit <- design_mfa(example_data$Y, example_data$design, ~ saliency)
head(effect_scores(example_fit, "saliency"))
#>        Component_1   Component_2
#> [1,] -3.703423e+00 -0.0003778003
#> [2,] -1.851782e+00  0.0001617018
#> [3,] -7.656625e-06  0.0004280810
#> [4,]  1.851776e+00  0.0001699288
#> [5,]  3.703437e+00 -0.0003819114
head(effect_loadings(example_fit, "saliency"))
#>         Component_1  Component_2
#> voxel_1 -0.54000494  0.418006209
#> voxel_2 -0.38578208  0.383434818
#> voxel_3 -0.23152065 -0.009654883
#> voxel_4 -0.07718883 -0.397461994
#> voxel_5  0.07718051 -0.421042401
#> voxel_6  0.23151998 -0.056320346
```

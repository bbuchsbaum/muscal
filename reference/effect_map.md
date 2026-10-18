# Extract One Effect Loading Map

Extract One Effect Loading Map

## Usage

``` r
effect_map(x, term, component = 1L)
```

## Arguments

- x:

  A fitted \`design_mfa\` object.

- term:

  A fitted design term.

- component:

  A single component index.

## Value

A named numeric loading vector.

## Examples

``` r
example_data <- design_mfa_example_data()
example_fit <- design_mfa(example_data$Y, example_data$design, ~ saliency)
effect_map(example_fit, "saliency", component = 1L)
#>     voxel_1     voxel_2     voxel_3     voxel_4     voxel_5     voxel_6 
#> -0.54000494 -0.38578208 -0.23152065 -0.07718883  0.07718051  0.23151998 
#>     voxel_7     voxel_8 
#>  0.38578967  0.54001382 
```

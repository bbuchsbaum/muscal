# Predict Group-Average Activation Fields

Predict Group-Average Activation Fields

## Usage

``` r
# S3 method for class 'design_mfa'
predict(object, new_design = NULL, type = c("response"), ...)
```

## Arguments

- object:

  A fitted \`design_mfa\` object.

- new_design:

  Optional design rows. Defaults to the fitted reference grid.

- type:

  Currently \`"response"\`.

- ...:

  Unused.

## Value

A matrix of group-average fitted activation fields.

## Details

Full-response prediction is stricter than extracting a fitted design
effect. Every coefficient functional needed for each requested
prediction row must be estimable under the complete fitted model. A
within-subject term effect and its partial sum of squares can therefore
remain estimable when an unrelated between-subject covariate is aliased
with the subject intercepts, even though \`predict()\` refuses the
corresponding full response. Reported term effects use only the
estimable contribution of the columns removed by that term's hypothesis.
Use \[estimability()\] to compare the \`effect_nonestimable_grid_rows\`
and \`prediction_nonestimable_grid_rows\` diagnostics before predicting.

## Examples

``` r
example_data <- design_mfa_example_data()
example_fit <- design_mfa(example_data$Y, example_data$design, ~ saliency + probe)
predict(example_fit, new_design = data.frame(
  saliency = ordered(c(20, 100), levels = c(20, 40, 60, 80, 100)),
  probe = factor(c("old", "lure"), levels = c("old", "new", "lure"))
))
#>     voxel_1    voxel_2    voxel_3    voxel_4    voxel_5    voxel_6   voxel_7
#> 1  1.324822  0.9917741  0.8935248  0.7113824  0.1898471 -0.6709828 -1.615502
#> 2 -1.675080 -1.2420106 -1.0435223 -0.7611432 -0.1395912  0.8210203  1.865286
#>    voxel_8
#> 1 -2.32510
#> 2  2.67483
```

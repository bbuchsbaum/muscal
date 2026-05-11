# Bidirectional Response Alignment with aligned_interbattery()

## Why use `aligned_interbattery()`?

Use
[`aligned_interbattery()`](https://bbuchsbaum.github.io/muscal/reference/aligned_interbattery.md)
when you have two paired block families for the same subjects and
neither side is the privileged response surface. You want a shared
common space, but you also want to predict from `X` to `Y`, from `Y` to
`X`, or complete a partially observed bundle using both sides together.

The key contract is symmetric: each subject gets common scores plus one
score system for the `X` side and one for the `Y` side. That makes this
a better fit than forcing one side to play the role of “predictors” and
the other to play the role of “responses.”

## What does a quick fit look like?

``` r

library(muscal)
library(knitr)
```

``` r

fit <- aligned_interbattery(
  X = train$X,
  Y = train$Y,
  ncomp = 2,
  max_iter = 10
)
```

``` r

kable(quick_summary, align = c("l", "r"))
```

| path                 | mean_test_mse |
|:---------------------|--------------:|
| x:D1 + x:D2 -\> y:E1 |        0.0038 |
| y:E1 -\> x:D1        |        0.0159 |
| x:D1 + y:E1 -\> x:D2 |        0.0039 |

The important point is not which direction wins on this toy simulation.
It is that one fitted object supports all three contracts: `X -> Y`,
`Y -> X`, and joint completion from partially observed bundles.

## What do the inputs look like?

You can supply nested `subject -> domain` lists on both sides. Blocks
with the same domain name and column schema share loadings across
subjects, which is the main way the model learns a common geometry
rather than fitting each subject independently.

``` r

c(
  n_subjects = length(train$X),
  x_domains_per_subject = length(train$X$S1),
  y_domains_per_subject = length(train$Y$S1)
)
#>            n_subjects x_domains_per_subject y_domains_per_subject 
#>                     3                     2                     1
```

Here each subject has two `X`-side domains (`D1`, `D2`) and one `Y`-side
domain (`E1`). In a real study those domains might be assays, feature
families, or paired modalities collected on the same rows.

## How do you project from one side or both?

If you observe both sides for the same new rows,
[`project()`](https://bbuchsbaum.github.io/muscal/reference/project.md)
can return the common scores together with the side-specific score
systems.

``` r

common_scores <- project(
  fit,
  partial_both,
  from = "both",
  side = "common"
)

dim(common_scores)
#> [1] 9 2
```

If you only want to reconstruct the blocks you actually observed, use
`type = "reconstruction"`.

``` r

reconstructed <- predict(
  fit,
  partial_both,
  from = "both",
  type = "reconstruction"
)

names(reconstructed)
#> [1] "x" "y"
```

That separation matters. Reconstruction is for the observed side bundle.
Prediction is for completing the opposite side or an unobserved block
given the common space you inferred from what you have.

## Where next?

- [`vignette("response_alignment")`](https://bbuchsbaum.github.io/muscal/articles/response_alignment.md)
  when one side really is the privileged response surface
- [`?aligned_interbattery`](https://bbuchsbaum.github.io/muscal/reference/aligned_interbattery.md)
  for the fitting contract and supported input layouts
- [`?project.aligned_interbattery`](https://bbuchsbaum.github.io/muscal/reference/project.aligned_interbattery.md)
  and
  [`?predict.aligned_interbattery`](https://bbuchsbaum.github.io/muscal/reference/predict.aligned_interbattery.md)
  for the out-of-sample projection and completion rules

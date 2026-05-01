# Choosing a Supervised Common-Space Model

## Why this comparison matters

Sometimes the response table really defines a small set of shared states
across blocks. Sometimes every block has its own paired multivariate
response with no exact row sharing. Sometimes you are in between: some
rows map back to known anchor states, while others are novel.

This vignette is about the cases where one side is the privileged
response surface. Within that scope, `muscal` has three honest answers:

- [`anchored_mfa()`](https://bbuchsbaum.github.io/muscal/reference/linked_mfa.md)
  when the response lives at the anchor-state level
- [`response_aligned_mfa()`](https://bbuchsbaum.github.io/muscal/reference/response_aligned_mfa.md)
  when each block has paired multivariate `Y`, with or without partial
  anchor structure
- [`aligned_rrr()`](https://bbuchsbaum.github.io/muscal/reference/aligned_rrr.md)
  when you want the clean reduced-rank regression baseline

This vignette shows the same prediction problem under each regime and
focuses on the contract you actually need, not on forcing every case
through repeated rows of `Y`.

If neither side is privileged and you want symmetric prediction or
completion between paired `X` and `Y` block families, use
[`aligned_interbattery()`](https://bbuchsbaum.github.io/muscal/reference/aligned_interbattery.md)
instead; see
[`vignette("aligned_interbattery")`](https://bbuchsbaum.github.io/muscal/articles/aligned_interbattery.md).

## Load the package

``` r

library(muscal)
library(multivarious)
library(knitr)
```

## Quick start: paired responses, no exact anchor assumption

If every block has its own paired multivariate response, start with
[`response_aligned_mfa()`](https://bbuchsbaum.github.io/muscal/reference/response_aligned_mfa.md).
If you want the cleaner supervised baseline with no predictor
reconstruction term, fit
[`aligned_rrr()`](https://bbuchsbaum.github.io/muscal/reference/aligned_rrr.md)
alongside it.

``` r

fit_response <- response_aligned_mfa(
  Y = blockwise$train$Y,
  X = blockwise$train$X,
  ncomp = 2,
  preproc = multivarious::pass(),
  response_preproc = multivarious::pass(),
  normalization = "None",
  ridge = 1e-8,
  max_iter = 80,
  tol = 1e-9
)

fit_rrr <- aligned_rrr(
  Y = blockwise$train$Y,
  X = blockwise$train$X,
  ncomp = 2,
  preproc = multivarious::pass(),
  response_preproc = multivarious::pass(),
  ridge = 1e-8,
  max_iter = 80,
  tol = 1e-9
)
```

``` r

kable(quick_summary, align = c("l", "r"))
```

| model                  | mean_test_mse |
|:-----------------------|--------------:|
| response_aligned_mfa() |        0.0093 |
| aligned_rrr()          |        0.0092 |

Both models predict from block-specific `X` into a shared multivariate
response space. The practical difference is structural:
[`response_aligned_mfa()`](https://bbuchsbaum.github.io/muscal/reference/response_aligned_mfa.md)
also models the predictor blocks through a shared latent space, while
[`aligned_rrr()`](https://bbuchsbaum.github.io/muscal/reference/aligned_rrr.md)
is the lean reduced-rank baseline.

## Which data contract do you actually have?

``` r

chooser <- data.frame(
  Data_regime = c(
    "Exact shared anchor states",
    "Some anchored rows, some novel rows",
    "No anchor states, paired blockwise responses"
  ),
  What_is_observed = c(
    "One anchor-level Y table plus row maps from each X block",
    "Blockwise Y_k plus optional anchor_map for the rows you can link",
    "Only blockwise paired (X_k, Y_k)"
  ),
  Recommended_fit = c(
    "anchored_mfa()",
    "response_aligned_mfa()",
    "response_aligned_mfa() plus aligned_rrr() baseline"
  )
)

kable(chooser, align = "l")
```

| Data_regime | What_is_observed | Recommended_fit |
|:---|:---|:---|
| Exact shared anchor states | One anchor-level Y table plus row maps from each X block | anchored_mfa() |
| Some anchored rows, some novel rows | Blockwise Y_k plus optional anchor_map for the rows you can link | response_aligned_mfa() |
| No anchor states, paired blockwise responses | Only blockwise paired (X_k, Y_k) | response_aligned_mfa() plus aligned_rrr() baseline |

The model choice should follow that contract directly. The point is not
to force every problem into repeated rows of `Y`, but to keep the
common-space assumption honest.

## When the response really is an anchor table

Use
[`anchored_mfa()`](https://bbuchsbaum.github.io/muscal/reference/linked_mfa.md)
when the scientifically meaningful response lives at the anchor-state
level and each predictor block only tells you which anchor rows it
touches.

``` r

fit_anchor <- anchored_mfa(
  Y = exact_anchor$Y_anchor,
  X = exact_anchor$train$X,
  row_index = exact_anchor$train$row_index,
  ncomp = 2,
  preproc = multivarious::pass(),
  normalization = "None",
  ridge = 1e-8,
  max_iter = 80,
  tol = 1e-9
)
```

``` r

kable(anchor_summary, align = c("l", "r", "r", "r"))
```

| model          | mean_test_mse | n_anchor_states | learned_score_rows |
|:---------------|--------------:|----------------:|-------------------:|
| anchored_mfa() |         0.116 |              24 |                 24 |

That is the clean anchored contract: the shared score table lives over
anchor states, not over block rows, and new `X` rows are projected back
to that anchor-level response surface.

## When some rows are anchored but others are novel

This is the regime
[`response_aligned_mfa()`](https://bbuchsbaum.github.io/muscal/reference/response_aligned_mfa.md)
was generalized to handle. You keep the blockwise responses `Y_k`, pass
anchor information only where it is real, and let the novel rows stay
free.

``` r

fit_hybrid <- response_aligned_mfa(
  Y = blockwise$train$Y,
  X = blockwise$train$X,
  ncomp = 2,
  preproc = multivarious::pass(),
  response_preproc = multivarious::pass(),
  normalization = "None",
  anchor_response = blockwise$anchor_response,
  anchor_map = blockwise$train$anchor_map,
  coupling_lambda = 2,
  ridge = 1e-8,
  max_iter = 80,
  tol = 1e-9
)
```

``` r

kable(hybrid_summary, align = c("l", "r"))
```

| prediction_path          | mean_test_mse |
|:-------------------------|--------------:|
| X only                   |        0.0093 |
| X + test-time anchor_map |        0.0093 |

The important point is not just that the fit stores both `Z_list` and
`S`. It is that prediction stays honest: default prediction uses only
`X`, while known test-time anchor information can refine the score solve
when you truly have it.

## How should you choose?

- Use
  [`anchored_mfa()`](https://bbuchsbaum.github.io/muscal/reference/linked_mfa.md)
  when the shared response really is an anchor-state object and your
  blocks connect to it through row maps.
- Use
  [`response_aligned_mfa()`](https://bbuchsbaum.github.io/muscal/reference/response_aligned_mfa.md)
  when each block has its own paired multivariate response and you want
  a shared latent geometry for the predictor blocks.
- Add `anchor_map` and `anchor_response` to
  [`response_aligned_mfa()`](https://bbuchsbaum.github.io/muscal/reference/response_aligned_mfa.md)
  only when those anchors are genuine; missing anchor rows should stay
  missing.
- Fit
  [`aligned_rrr()`](https://bbuchsbaum.github.io/muscal/reference/aligned_rrr.md)
  alongside
  [`response_aligned_mfa()`](https://bbuchsbaum.github.io/muscal/reference/response_aligned_mfa.md)
  when prediction is the main target and you want to know whether the
  extra predictor-reconstruction structure is actually buying anything.

## Where next?

- [`vignette("linked_mfa")`](https://bbuchsbaum.github.io/muscal/articles/linked_mfa.md)
  for the anchored row-mapping workflow in more depth
- [`vignette("aligned_mfa")`](https://bbuchsbaum.github.io/muscal/articles/aligned_mfa.md)
  for the unsupervised same-latent-row analogue
- [`vignette("aligned_interbattery")`](https://bbuchsbaum.github.io/muscal/articles/aligned_interbattery.md)
  when you want symmetric `X`/`Y` prediction rather than a privileged
  response surface
- [`vignette("model_evaluation")`](https://bbuchsbaum.github.io/muscal/articles/model_evaluation.md)
  for bootstrap, permutation, and held-out evaluation workflows across
  supported methods
- [`?response_aligned_mfa`](https://bbuchsbaum.github.io/muscal/reference/response_aligned_mfa.md)
  and
  [`?aligned_rrr`](https://bbuchsbaum.github.io/muscal/reference/aligned_rrr.md)
  for the full argument contracts

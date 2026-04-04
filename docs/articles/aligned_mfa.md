# Aligned Multiple Factor Analysis

## Why aligned MFA?

Standard MFA assumes every block has exactly the same rows.
[`aligned_mfa()`](https://bbuchsbaum.github.io/muscal/reference/aligned_mfa.md)
relaxes that assumption. You use it when blocks point back to the same
latent population, but each block observes only a subset of rows, or
observes them in different patterns.

The key difference from
[`anchored_mfa()`](https://bbuchsbaum.github.io/muscal/reference/linked_mfa.md)
is that no single block is treated as the privileged reference. The
shared scores live in a latent reference row space, and each block
contributes through its own `row_index` mapping.

## What does a quick fit look like?

``` r
library(muscal)
library(multivarious)
library(ggplot2)
```

``` r
c(
  n_expression_rows = nrow(blocks$expression),
  n_imaging_rows = nrow(blocks$imaging),
  n_behavior_rows = nrow(blocks$behavior)
)
#> n_expression_rows    n_imaging_rows   n_behavior_rows 
#>                60                55                55
```

``` r
fit <- aligned_mfa(blocks, row_index, N = N, ncomp = 2, normalization = "None")
S <- multivarious::scores(fit)

stopifnot(all(dim(S) == c(N, 2)))
stopifnot(all(is.finite(S)))
stopifnot(length(fit$objective_trace) >= 1)
```

``` r
cc <- cancor(
  scale(S_true, center = TRUE, scale = FALSE),
  scale(S, center = TRUE, scale = FALSE)
)$cor[1:2]

stopifnot(all(is.finite(cc)))
stopifnot(mean(cc) > 0.77)

round(cc, 3)
#> [1] 0.913 0.779
```

![Aligned MFA estimates one score vector for each latent reference row.
Here the points are colored by how many blocks cover that
row.](aligned_mfa_files/figure-html/quick-score-plot-1.png)

Aligned MFA estimates one score vector for each latent reference row.
Here the points are colored by how many blocks cover that row.

Rows observed in two blocks are better constrained than rows observed
once, but the model still estimates a single coherent score space over
all `N` latent rows.

## How do row mappings work?

``` r
lapply(row_index, head, 8)
#> $expression
#> [1]  1  4  5  6  7  9 10 11
#> 
#> $imaging
#> [1]  1  2  3  4  5  7  8 10
#> 
#> $behavior
#> [1]  1  2  3  4  5  6  9 10
```

Each integer vector says which latent reference row each observed row
belongs to. This is the entire alignment mechanism: if row 7 in
`expression` and row 3 in `imaging` both map to latent row 21, they
contribute to the same shared score vector.

## What happens when all rows line up?

When every block observes the same rows in the same order,
[`aligned_mfa()`](https://bbuchsbaum.github.io/muscal/reference/aligned_mfa.md)
should agree with ordinary
[`mfa()`](https://bbuchsbaum.github.io/muscal/reference/mfa.md).

``` r
fit_aligned_full <- aligned_mfa(
  full_blocks,
  list(expression = seq_len(N), imaging = seq_len(N), behavior = seq_len(N)),
  N = N,
  ncomp = 2,
  normalization = "None",
  ridge = 1e-10,
  max_iter = 200,
  tol = 1e-10
)
fit_mfa <- mfa(full_blocks, ncomp = 2, normalization = "None")

P1 <- multivarious::scores(fit_aligned_full) %*%
  solve(crossprod(multivarious::scores(fit_aligned_full)), t(multivarious::scores(fit_aligned_full)))
P2 <- multivarious::scores(fit_mfa) %*%
  solve(crossprod(multivarious::scores(fit_mfa)), t(multivarious::scores(fit_mfa)))
rel <- norm(P1 - P2, type = "F") / (norm(P2, type = "F") + 1e-12)

stopifnot(is.finite(rel))
stopifnot(rel < 0.05)

rel
#> [1] 1.523229e-05
```

That reduction check is important: it confirms that you are getting a
genuine extension of MFA rather than a separate, incompatible factor
model.

## What should you inspect while fitting?

``` r
sapply(fit$V_list, dim)
#>      expression imaging behavior
#> [1,]         20      18       16
#> [2,]          2       2        2
```

![The aligned-MFA objective should fall quickly and then
flatten.](aligned_mfa_files/figure-html/objective-trace-1.png)

The aligned-MFA objective should fall quickly and then flatten.

The loading dimensions tell you how much information each block carries.
The objective trace tells you whether the alternating updates have
settled cleanly.

## Where next?

- [`vignette("linked_mfa")`](https://bbuchsbaum.github.io/muscal/articles/linked_mfa.md)
  for the anchored version with an explicit reference block
- [`vignette("mfa")`](https://bbuchsbaum.github.io/muscal/articles/mfa.md)
  for the classical same-row setting
- [`?aligned_mfa`](https://bbuchsbaum.github.io/muscal/reference/aligned_mfa.md)
  for `feature_groups`, `feature_lambda`, and weighting options

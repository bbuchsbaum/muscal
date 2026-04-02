# Barycentric MFA (BaMFA)

## Why BaMFA?

Multi-subject data often contains two kinds of structure: patterns
shared across all subjects (*global* components) and patterns unique to
individuals (*local* components). Standard MFA captures the global
structure but discards subject-specific variation. If you care about
individual differences — which subjects deviate, and how — you need a
method that explicitly models both.

BaMFA (Barycentric Multiple Factor Analysis) decomposes multi-block data
into **global** and **local** components simultaneously. Global
components capture the consensus structure; local components capture
what makes each subject unique.

## Quick start

``` r
library(muscal)
library(multivarious)
library(ggplot2)
```

We simulate 6 subjects with shared structure (3 global factors) plus
subject-specific structure (2 local factors per subject).

``` r
length(data_list)
#> [1] 6
sapply(data_list, dim)
#>      Subject_1 Subject_2 Subject_3 Subject_4 Subject_5 Subject_6
#> [1,]        50        50        50        50        50        50
#> [2,]        20        20        20        20        20        20
```

Fit BaMFA with 3 global and 2 local components:

``` r
fit <- bamfa(data_list, k_g = 3, k_l = 2, niter = 20)
fit
#> Barycentric Multiple Factor Analysis (BaMFA) 
#> 
#> Model Structure: 
#>   Global components (k_g): 3 
#>   Local components (k_l requested): 2 
#>   Convergence after 20 iterations
#>   Local score regularization (lambda_l): 0 
#> 
#> Block Information: 
#>   Block 1 ( Subject_1 ): 
#>     Features: 20 
#>     Local components retained: 2 
#>   Block 2 ( Subject_2 ): 
#>     Features: 20 
#>     Local components retained: 2 
#>   Block 3 ( Subject_3 ): 
#>     Features: 20 
#>     Local components retained: 2 
#>   Block 4 ( Subject_4 ): 
#>     Features: 20 
#>     Local components retained: 2 
#>   Block 5 ( Subject_5 ): 
#>     Features: 20 
#>     Local components retained: 2 
#>   Block 6 ( Subject_6 ): 
#>     Features: 20 
#>     Local components retained: 2 
#> 
#> Final reconstruction error (per feature): 0.0733104
```

![Mean global scores across subjects. These capture the consensus
structure shared by all
blocks.](bamfa_files/figure-html/quick-score-plot-1.png)

Mean global scores across subjects. These capture the consensus
structure shared by all blocks.

## Global vs. local components

The key idea: BaMFA decomposes each block’s data as:

``` math
\mathbf{X}_i \approx \mathbf{S}_i \mathbf{B}_i^\top + \mathbf{U}_i \mathbf{V}_i^\top
```

where:

- $`\mathbf{S}_i`$ are the **global scores** — similar across subjects
- $`\mathbf{B}_i`$ are the **global loadings**
- $`\mathbf{U}_i`$ are the **local scores** — unique to subject $`i`$
- $`\mathbf{V}_i`$ are the **local loadings**

The global loadings $`\mathbf{B}_i`$ are encouraged to be similar across
subjects, while the local loadings $`\mathbf{V}_i`$ are free to differ.

## Examining the components

``` r
# Global loadings: one per block, should be similar
sapply(fit$B_list, dim)
#>      [,1] [,2] [,3] [,4] [,5] [,6]
#> [1,]   20   20   20   20   20   20
#> [2,]    2    2    2    2    2    2

# Local loadings: one per block, capture individual structure
sapply(fit$V_list, dim)
#> list()
```

## Convergence

``` r
plot_convergence(fit)
```

![BaMFA convergence trace. The alternating least squares algorithm
typically converges within 10-20
iterations.](bamfa_files/figure-html/convergence-1.png)

BaMFA convergence trace. The alternating least squares algorithm
typically converges within 10-20 iterations.

## Per-subject score plots

View the global scores for a specific subject:

``` r
autoplot(fit, block = 1)
```

![Global scores for Subject
1.](bamfa_files/figure-html/subject-scores-1.png)

Global scores for Subject 1.

## Partial factor scores

Compare how the global structure looks across subjects:

``` r
plot_partial_scores(fit, connect = TRUE, show_consensus = TRUE)
```

![Partial factor scores: each subject's global scores overlaid.
Agreement indicates consistent global
structure.](bamfa_files/figure-html/partial-scores-1.png)

Partial factor scores: each subject’s global scores overlaid. Agreement
indicates consistent global structure.

Short connecting lines mean subjects agree about the global structure of
an observation. Large deviations indicate observations where subjects
diverge.

## The local penalty (lambda_l)

The `lambda_l` parameter controls regularization on local components. By
default it is 0 (no regularization). Increasing it shrinks local
components, pushing more structure into the global subspace:

``` r
# Strong local regularization — most structure forced into global
fit_sparse <- bamfa(data_list, k_g = 3, k_l = 2, niter = 20, lambda_l = 1)
```

## Choosing k_g and k_l

- **`k_g`** (global components): Set this to the number of shared
  factors you expect. Start with the number of dominant eigenvalues from
  pooled PCA.
- **`k_l`** (local components): Set this to capture meaningful
  individual variation. Too many local components can absorb noise.

A useful diagnostic: if partial factor scores are tight (short lines in
the plot), `k_g` is capturing the shared structure well. If they’re
spread out, you may need more global components.

## When to use BaMFA

BaMFA is appropriate when:

- You have multi-subject or multi-condition data blocks
- You want to separate shared (global) from individual (local) structure
- Individual differences are scientifically interesting, not just noise
- You need interpretable per-subject components

BaMFA is *not* appropriate when:

- You only care about the consensus (use
  [`mfa()`](https://bbuchsbaum.github.io/muscal/reference/mfa.md) —
  simpler and faster)
- Blocks have different observations (use
  [`anchored_mfa()`](https://bbuchsbaum.github.io/muscal/reference/linked_mfa.md))
- Your data are covariance matrices (use
  [`covstatis()`](https://bbuchsbaum.github.io/muscal/reference/covstatis.md))

## Next steps

- [`vignette("mfa")`](https://bbuchsbaum.github.io/muscal/articles/mfa.md)
  — Standard MFA for consensus-only analysis
- [`vignette("penalized_mfa")`](https://bbuchsbaum.github.io/muscal/articles/penalized_mfa.md)
  — Penalized MFA for loading similarity
- [`?bamfa`](https://bbuchsbaum.github.io/muscal/reference/bamfa.md) —
  Full parameter documentation

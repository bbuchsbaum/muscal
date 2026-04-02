# muscal

`muscal` is an R package for multiblock and matrix-based multivariate analysis.
It focuses on methods for integrating several tables measured on the same
observations, linked designs with partial overlap, and repeated
covariance/connectivity matrices.

The package includes classical "French school" methods such as Multiple Factor
Analysis (MFA) and STATIS-style workflows, along with newer extensions for
adaptive weighting, regularization, and repeated-matrix modeling. Most fitted
objects work naturally with the `multivarious` ecosystem for scoring,
projection, and preprocessing.

## What to use when

- `mfa()` for several data blocks measured on the same observations.
- `ipca()` when you want adaptive block weighting instead of fixed MFA
  normalization.
- `penalized_mfa()` when you want MFA with loading-similarity penalties.
- `anchored_mfa()` for linked blocks that map back to a reference table but do
  not share identical rows.
- `aligned_mfa()` and `aligned_mcca()` for aligned multiblock integration and
  correlation analysis.
- `covstatis()` for lists of covariance or connectivity matrices.
- `bada()` and `bamfa()` for barycentric and discriminant variants.
- `bilinear_mixed()` for repeated connectivity matrices with subject- and
  design-level structure.

## Installation

Install the package from r-universe:

```r
install.packages(
  "muscal",
  repos = c("https://bbuchsbaum.r-universe.dev", "https://cloud.r-project.org")
)
```

Or install the current GitHub version with `pak`:

```r
install.packages("pak")
pak::pak("bbuchsbaum/muscal")
```

## Quick start

```r
library(muscal)

sim <- synthetic_multiblock(
  S = 4,
  n = 60,
  p = c(20, 30, 15, 25),
  r = 3,
  sigma = 0.3,
  seed = 42
)

fit <- mfa(sim$data_list, ncomp = 3)

dim(multivarious::scores(fit))
#> [1] 60  3

plot_variance(fit, type = "bar")
ggplot2::autoplot(fit)
```

That example covers the most common case: several feature blocks measured on
the same subjects. If your data are linked across uneven row sets, covariance
matrices, or repeated connectivity observations, jump to the method-specific
articles below.

## Where to start

- Package site: <https://bbuchsbaum.github.io/muscal/>
- Package overview: [docs/index.md](docs/index.md)
- Articles index: [docs/articles/index.md](docs/articles/index.md)
- Reference index: [docs/reference/index.md](docs/reference/index.md)
- MFA article: [docs/articles/mfa.md](docs/articles/mfa.md)
- MCCA intro: [vignettes/mcca.Rmd](vignettes/mcca.Rmd)
- Integrative PCA article: [docs/articles/ipca.md](docs/articles/ipca.md)
- Anchored MFA article: [docs/articles/linked_mfa.md](docs/articles/linked_mfa.md)
- Aligned MFA intro: [vignettes/aligned_mfa.Rmd](vignettes/aligned_mfa.Rmd)
- Aligned MCCA intro: [vignettes/aligned_mcca.Rmd](vignettes/aligned_mcca.Rmd)
- COVSTATIS article: [docs/articles/covstatis.md](docs/articles/covstatis.md)
- BaDA intro: [vignettes/bada.Rmd](vignettes/bada.Rmd)
- Bilinear mixed-model article: [docs/articles/bilinear_mixed.md](docs/articles/bilinear_mixed.md)

## Package scope

Beyond the core factor methods, `muscal` also includes helpers for simulation
(`synthetic_multiblock()`), component estimation
(`estimate_components()`, `significant_components()`), projection of new data,
and a fairly broad set of plotting helpers for scores, loadings, block weights,
partial scores, convergence traces, and coverage diagnostics.

## Development

For local site builds, use `Rscript scripts/build-pkgdown.R`. On this setup,
the default `pkgdown::build_site()` call can fail in its `callr` subprocess
even when the package itself installs and the GitHub Pages build path works.

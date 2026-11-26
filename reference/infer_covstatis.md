# Test leading-root concentration in COVSTATIS

Test whether the leading eigenvalue accounts for an unusually large
share of the spectrum under independent ROI relabeling across tables.

## Usage

``` r
infer_covstatis(
  x,
  data,
  axes = c("compromise", "interstructure", "both"),
  B = 1000L,
  alpha = 0.05,
  seed = NULL,
  strict = TRUE,
  ...,
  max_steps = 1L,
  auto_subspace = TRUE,
  tie_threshold = 0.01
)
```

## Arguments

- x:

  A fitted `covstatis` object.

- data:

  Original list of covariance/correlation matrices used to fit `x`.

- axes:

  Matrix to test. `"interstructure"` uses the table-by-table matrix of
  inner products, `"compromise"` uses the weighted compromise, and
  `"both"` returns one test for each.

- B:

  Maximum number of Monte Carlo draws per test.

- alpha:

  Significance threshold passed to
  [`multifer::infer()`](https://rdrr.io/pkg/multifer/man/infer.html).

- seed:

  Optional random seed.

- strict:

  Logical; passed to
  [`multifer::infer()`](https://rdrr.io/pkg/multifer/man/infer.html).

- ...:

  Additional arguments passed to
  [`multifer::infer()`](https://rdrr.io/pkg/multifer/man/infer.html).

- max_steps:

  Only the leading root can be tested (`1`, the default).

- auto_subspace:

  Group nearly tied leading roots in the reported unit. Roots tied to
  numerical precision always share a unit.

- tie_threshold:

  Relative eigenvalue gap used to group near ties.

## Value

A
[`multifer::infer_result`](https://rdrr.io/pkg/multifer/man/infer_result.html)
object, or a named list of two such objects when `axes = "both"`.

## Details

The statistic is the largest positive eigenvalue divided by the sum of
**all** positive eigenvalues of the matrix being tested, or zero when
there is no positive inertia. Both the observed statistic and every null
draw use the full spectrum, regardless of the fitted `ncomp`.

Conditional on each preprocessed table up to ROI labeling, the null
assumes the ROI labels are independently exchangeable within each table.
Each draw applies the same permutation to rows and columns of a table
and recomputes the interstructure and, for the compromise, its table
weights. With unequal variable weights this assumption concerns the
weighted tables.

The p-value tests concentration under this global correspondence null.
The unit's `selected` flag records rejection of this global null. When
leading roots are tied, the reported unit spans the tied eigenspace and
the p-value still concerns the largest root alone. It supplies no
estimate of how many dimensions to retain. Subsequent axes are
descriptive. Independent ROI permutations of deflated tables can move
their shared removed subspace. Conditional dimension selection would
require a separate residual null.

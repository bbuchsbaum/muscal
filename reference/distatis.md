# DISTATIS for Multiple Distance Matrices

Fits a consensus object map, an assessor map based on RV similarities,
and partial object maps for distances measured on the same objects.

## Usage

``` r
distatis(
  data,
  ncomp = 2,
  squared = FALSE,
  norm_method = c("mfa", "frobenius", "none"),
  labels = NULL,
  tol = 1e-08,
  ...
)

# Default S3 method
distatis(
  data,
  ncomp = 2,
  squared = FALSE,
  norm_method = c("mfa", "frobenius", "none"),
  labels = NULL,
  tol = 1e-08,
  ...
)
```

## Arguments

- data:

  A list of numeric square distance matrices or
  [`stats::dist()`](https://rdrr.io/r/stats/dist.html) objects, an
  object-by-object-by-assessor numeric array, or a single distance
  matrix or `dist` object. Object names and order must agree across
  tables.

- ncomp:

  Positive integer; maximum number of axes retained in each map. Object
  and assessor ranks are determined independently.

- squared:

  Logical; are the supplied distances already squared? Defaults to
  `FALSE`, including for `dist` objects.

- norm_method:

  Table normalization: `"mfa"` (default) divides each centered Gram
  matrix by its largest eigenvalue; `"frobenius"` gives unit Frobenius
  norm; `"none"` preserves scale. RV similarities are used in all three
  cases.

- labels:

  Optional unique object labels. These override display labels, but do
  not bypass checks of input object names and order.

- tol:

  Relative numerical tolerance, strictly between zero and one, for input
  symmetry, Euclidean geometry, spectral rank, and leading-root ties.

- ...:

  Reserved for methods. Unused arguments are rejected.

## Value

A `distatis` object inheriting from `covstatis` and `projector`.
`object_scores` (also `roi_scores`) and `subject_scores` contain the two
maps; `partial_scores` contains one object map per assessor; `alpha`
contains weights; `C` is the RV matrix; `compromise_full` is the
normalized compromise Gram matrix. `eigenvalues` and
`interstructure_eigenvalues` retain the full spectra, while `inertia`
and `subject_inertia` use their complete positive inertia as
denominators. `normalization_constants`, `input_labels`, `squared`,
`norm_method`, `tol`, and `diagnostics` record preprocessing.

## Details

For ordinary distance matrices \\D_k\\, DISTATIS computes
\\B_k=-HD_k^2H/2\\, with elementwise squaring and \\H=I-11^T/n\\. The
normalized Gram matrices \\S_k\\ are compared using RV coefficients. The
leading eigenvector of the RV matrix, oriented nonnegatively and
rescaled to sum to one, supplies the compromise weights. The weighted
compromise has object coordinates \\V\Lambda^{1/2}\\; each partial map
is \\S_kV\Lambda^{-1/2}\\. Their weighted mean equals the object map.
Assessor coordinates use the RV matrix's own spectrum.

Non-Euclidean distances, zero Gram matrices, missing values, and a tied
leading RV eigenvalue are rejected. Small negative eigenvalues within
`tol` are treated as numerical zero for rank and inertia calculations.
No distance correction or imputation is performed. Weights describe
agreement with the fitted panel, rather than assessor accuracy.
Inference is not currently implemented.

## References

Abdi, H., Valentin, D., Chollet, S., & Chrea, C. (2007). Analyzing
assessors and products in sorting tasks: DISTATIS, theory and
applications. Food Quality and Preference, 18, 627-640.
[doi:10.1016/j.foodqual.2006.09.003](https://doi.org/10.1016/j.foodqual.2006.09.003)
.

## See also

[`covstatis()`](https://bbuchsbaum.github.io/muscal/reference/covstatis.md),
[`sorting_distances()`](https://bbuchsbaum.github.io/muscal/reference/sorting_distances.md),
[`project_distance()`](https://bbuchsbaum.github.io/muscal/reference/project_distance.md),
[`project_subjects()`](https://bbuchsbaum.github.io/muscal/reference/project_subjects.md),
[`reconstruct.distatis()`](https://bbuchsbaum.github.io/muscal/reference/reconstruct.distatis.md),
[`plot.distatis()`](https://bbuchsbaum.github.io/muscal/reference/plot.distatis.md)

## Examples

``` r
partitions <- data.frame(A = c(1, 1, 2, 3), B = c(1, 1, 2, 2),
                         C = c(1, 2, 3, 3))
distances <- sorting_distances(partitions)
fit <- distatis(distances)
fit$alpha
#> [1] 0.3269993 0.3460015 0.3269993
multivarious::scores(fit)
#>   Dim1          Dim2
#> 1 -0.5  0.000000e+00
#> 2 -0.5  3.723435e-17
#> 3  0.5 -3.301511e-01
#> 4  0.5  3.301511e-01
project_distance(fit, distances[[1]])
#>   Dim1          Dim2
#> 1 -0.5 -1.347650e-17
#> 2 -0.5 -1.347650e-17
#> 3  0.5 -1.009639e+00
#> 4  0.5  1.009639e+00
```

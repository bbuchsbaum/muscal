# Supplementary Assessor Maps and Diagnostics for DISTATIS

Supplementary Assessor Maps and Diagnostics for DISTATIS

## Usage

``` r
# S3 method for class 'distatis'
project_subjects(x, new_data, subject_ids = NULL, squared = x$squared, ...)
```

## Arguments

- x:

  A fitted `distatis` model.

- new_data:

  A distance matrix or `dist` object on the fitted objects in their
  original order. This projects a new assessor, not new objects.

- subject_ids:

  Optional unique names for supplementary assessors.

- squared:

  Whether new distances are already squared. Defaults to the setting
  used when fitting `x`.

- ...:

  Unused arguments, rejected.

## Value

A list with `subject_scores` in the fitted RV space, `subject_cosines`
(cosines of those retained coordinates), `object_scores` (also
`roi_scores`) containing partial object maps, and `scalar_summaries`.
The latter reports RV against the full compromise, Frobenius
`distance_to_compromise`, and `distance_to_subspace` for the retained
object subspace. `new_data` accepts the same representations as
[`distatis()`](https://bbuchsbaum.github.io/muscal/reference/distatis.md).
A training table reproduces its fitted assessor and partial object
coordinates.

# Project and Summarize New Subject Data (generic)

Generic for projecting one or more new subject-level data instances
(e.g., matrices) into the compromise space of a fitted model. This
function provides a comprehensive analysis by computing several key
metrics for each new instance.

This function provides a comprehensive analysis for one or more new
subject matrices by projecting them into the compromise space of a
fitted `covstatis` model. It computes several key metrics for each new
matrix, following the projection logic of `DISTATIS`.

## Usage

``` r
project_subjects(x, new_data, ...)

# S3 method for class 'covstatis'
project_subjects(x, new_data, subject_ids = NULL, ...)
```

## Arguments

- x:

  A fitted model object (e.g., \`covstatis\`).

- new_data:

  A single data instance or a list of instances to project. Each
  instance must have dimensions compatible with the training data of the
  model.

- ...:

  Other arguments passed to methods.

- subject_ids:

  Optional character vector of identifiers for the new subjects. If not
  provided, names will be taken from `new_data` or generated
  automatically.

## Value

A list containing scores, projections, and summary statistics. The exact
contents depend on the specific method.

A list containing the following elements:

- subject_scores:

  A matrix (`n_subjects` × `ncomp`) of subject-level coordinates in the
  compromise space.

- subject_cosines:

  A matrix (`n_subjects` × `ncomp`) of cosines, indicating alignment
  with each dimension.

- scalar_summaries:

  A `data.frame` with one row per subject, containing the
  `rv_coefficient` and `distance_to_compromise` (Frobenius distance to
  the full compromise).

- roi_scores:

  A list of matrices, where each element contains the ROI-level scores
  for a subject.

## Details

The function performs the following steps for each new matrix:

1.  Applies the same pre-processing (double-centering, normalization) as
    the original model.

2.  Calculates ROI-level factor scores (`roi_scores`), representing the
    coordinates of each ROI in the compromise space.

3.  Calculates subject-level scores (`subject_scores`) in the RV/table
    space by projecting the new table against the fitted interstructure
    axes.

4.  Calculates `subject_cosines` indicating the alignment of the
    subject's scores with each compromise dimension.

5.  Calculates a global `rv_coefficient` measuring the overall
    similarity between the new matrix and the group compromise matrix.

6.  Calculates the `distance_to_compromise`, the Frobenius distance
    between the new matrix and the full compromise matrix.

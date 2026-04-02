# Project a Subject-Level Covariate (generic)

Generic for projecting a subject-level covariate into the space of a
fitted model, treating it as a supplementary variable without re-fitting
the model.

This function projects a subject-level covariate into the space of a
fitted \`covstatis\` model, treating it as a supplementary variable
without re-fitting the model.

## Usage

``` r
project_covariate(x, y, ...)

# S3 method for class 'covstatis'
project_covariate(
  x,
  y,
  what = c("dimension", "observation"),
  scale = c("cosine", "beta"),
  ...
)
```

## Arguments

- x:

  A fitted model object (e.g., \`covstatis\`).

- y:

  Numeric vector representing the covariate, with length equal to the
  number of subjects in the model.

- ...:

  other arguments passed to methods.

- what:

  \`"dimension"\` (default) returns cosine / β per compromise dimension;
  \`"observation"\` returns an ROI-wise pattern.

- scale:

  \`"cosine"\` (−1…1) or \`"beta"\` (regression coeff.).

## Value

The return value depends on the method, but is typically a projection of
the covariate onto the model's components or a spatial pattern.

If \`what = "dimension"\`, a named numeric vector of length
\`ncomp(x)\`. If \`what = "observation"\`, an \`R × ncomp\` matrix.

## Details

The interpretation of the output depends on the \`what\` parameter: -
\*\*Dimension-wise cosine\*\*: This is the correlation between the
covariate \`y\` and the subject G-scores. It indicates which compromise
dimensions are most strongly associated with the covariate. -
\*\*Observation matrix\*\*: This is a weighted sum of each subject's
partial factor scores (\`Partial-F\`). It represents the spatial
signature in the compromise space associated with a one-unit change in
the covariate.

## Examples

``` r
# Create a list of correlation matrices
Xs <- lapply(1:5, function(i) matrix(rnorm(10*10), 10, 10))
Xs <- lapply(Xs, cor)

# Apply COVSTATIS
res <- covstatis(Xs, ncomp=3)

# Create a random covariate vector (e.g., episodic memory scores)
y <- rnorm(length(Xs))

# Project the covariate to get dimension-wise coordinates
dim_cos <- project_covariate(res, y, what = "dimension", scale = "cosine")

# Project the covariate to get an ROI-wise pattern
roi_beta <- project_covariate(res, y, what = "observation", scale = "beta")
```

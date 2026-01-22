#' Barycentric Discriminant Analysis (generic)
#'
#' A generic function for fitting a **Barycentric Discriminant Analysis** (BaDA) model
#' to multi–subject multivariate data.  Concrete methods (e.g. for objects of class
#' `multidesign`) implement the actual estimation procedure and return an object that
#' inherits from class `\code{bada}`.
#'
#' @section Documentation strategy:
#'   * This block documents **all** BaDA S3 methods via the shared
#'     `@rdname bada` tag.  Method‐specific files should therefore **omit** the
#'     `@return` field and simply add `@inheritParams bada` (or `@inherit bada`)
#'     plus any method‑specific `@details` or `@note` sections.
#'
#' @param data    A data object for which a BaDA method is defined (see individual
#'                methods for details).
#' @param y       A factor or bare column name giving the grouping variable.
#' @param subject A factor or bare column name identifying repeated‑measures, i.e.
#'                the subject/block variable.
#' @param preproc A preprocessing pipeline from the multivarious package.  Defaults to
#'                `multivarious::center()` in the concrete methods.
#' @param ncomp   Integer; the number of discriminant components to compute.
#' @param ...     Additional arguments passed to the underlying method.
#'
#' @return An object that inherits from class `bada`.  See the help for the
#'   corresponding method (e.g. [\code{bada.multidesign}]) for the exact
#'   structure.
#'
#' @references
#' Abdi, H., Williams, L. J., & Bera, M. (2017). *Barycentric discriminant
#' analysis*. In **Encyclopedia of Social Network Analysis and Mining** (pp.
#' 1–20).
#'
#' @export
#' @rdname bada
bada <- function(data, y, subject, preproc, ncomp, ...) {
  UseMethod("bada")
}


#' STATIS for Covariance Matrices (generic)
#'
#' A generic function implementing the STATIS approach for a collection of
#' covariance matrices.  Method implementations should return an object of class
#' `\code{covstatis}`.
#'
#' @inheritParams bada
#' @param ncomp     Integer; number of components to compute.
#' @param normalize Logical; if `TRUE` (default) each covariance matrix is scaled
#'                  to have unit Frobenius norm before analysis.
#' @param dcenter   Logical; if `TRUE` (default) each covariance matrix is double
#'                  centred (Gower transformation) prior to analysis.
#' @param ...       Additional arguments passed to the method.
#'
#' @return An object inheriting from class `covstatis`.
#'
#' @seealso [covstatis.list] for the reference implementation.
#'
#' @export
#' @rdname covstatis
covstatis <- function(data, ncomp = 2, normalize = TRUE, dcenter = TRUE, ...) {
  UseMethod("covstatis")
}


#' Multiple Factor Analysis (generic)
#'
#' A generic front–end for **Multiple Factor Analysis** (MFA).  Concrete methods
#' should perform the estimation and return an object inheriting from class
#' `\code{mfa}`.
#'
#' @inheritParams bada
#' @param normalization Character string specifying the block‑weighting scheme
#'   (see [mfa.multiblock]).
#' @param A,M Optional user‑supplied column/row weight matrices used when
#'   `normalization = "custom"`.
#' @param ... Additional arguments passed to the method.
#'
#' @return An object of class `mfa`.
#'
#' @references
#' Abdi, H., Williams, L. J., & Valentin, D. (2013). *Multiple factor analysis:
#' principal component analysis for multi‑table and multi‑block data sets*.
#' **Wiley Interdisciplinary Reviews: Computational Statistics, 5**(2), 149–179.
#'
#' @export
#' @rdname mfa
mfa <- function(data, preproc, ncomp = 2, normalization = "MFA", A = NULL, M = NULL, ...) {
  UseMethod("mfa")
}



#' Project a Covariance Matrix (generic)
#'
#' Generic for projecting a new covariance matrix onto a previously fitted
#' `covstatis` model.
#'
#' @param x        A model object for which a `project_cov` method exists (e.g.
#'                 an object of class `covstatis`).
#' @param new_data A symmetric numeric matrix to be projected.
#' @param ...      Additional arguments passed to methods.
#'
#' @return A numeric matrix of projected scores.
#'
#' @export
#' @rdname project_cov
project_cov <- function(x, new_data, ...) {
  UseMethod("project_cov")
}


#' Project and Summarize New Subject Data (generic)
#'
#' @description
#' Generic for projecting one or more new subject-level data instances (e.g., matrices)
#' into the compromise space of a fitted model. This function provides a comprehensive
#' analysis by computing several key metrics for each new instance.
#'
#' @param x A fitted model object (e.g., `covstatis`).
#' @param new_data A single data instance or a list of instances to project. Each instance
#'   must have dimensions compatible with the training data of the model.
#' @param ... Other arguments passed to methods.
#'
#' @return A list containing scores, projections, and summary statistics. The exact
#'   contents depend on the specific method.
#'
#' @export
#' @rdname project_subjects
project_subjects <- function(x, new_data, ...) {
  UseMethod("project_subjects")
}


#' Project a Subject-Level Covariate (generic)
#'
#' @description
#' Generic for projecting a subject-level covariate into the space of a fitted model,
#' treating it as a supplementary variable without re-fitting the model.
#'
#' @param x A fitted model object (e.g., `covstatis`).
#' @param y Numeric vector representing the covariate, with length equal to the number of subjects in the model.
#' @param ... other arguments passed to methods.
#'
#' @return The return value depends on the method, but is typically a projection
#'   of the covariate onto the model's components or a spatial pattern.
#'
#' @export
#' @rdname project_covariate
project_covariate <- function(x, y, ...) {
  UseMethod("project_covariate")
}


#' Re‑export selected generics from **multivarious**
#'
#' The package extends the `project()` and `reprocess()` generics defined in
#' \pkg{multivarious}.  Re‑exporting them avoids forcing users to load
#' \pkg{multivarious} explicitly.
#'
#' @importFrom multivarious project
#' @export
#' @rdname project
multivarious::project

#' @importFrom multivarious reprocess
#' @export
#' @rdname project
multivarious::reprocess
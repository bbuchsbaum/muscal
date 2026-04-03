#' Double center a matrix (Gower transformation)
#'
#' This function applies a double centering operation to a matrix, which is a necessary
#' preprocessing step for many matrix correlation methods. For covariance matrices,
#' the standard STATIS approach uses a Gower-like transformation.
#'
#' @param x A numeric matrix to be double-centered
#' @return A double-centered matrix with the same dimensions as x
#' @references Abdi, H., Williams, L. J., Valentin, D., & Bennani-Dosse, M. (2012).
#'            STATIS and DISTATIS: optimum multitable principal component analysis and
#'            three way metric multidimensional scaling. WIREs Computational Statistics, 4(2), 124–167.
#' @noRd
#' @keywords internal
double_center <- function(x) {
  # numerically cheaper than forming H %*% x %*% H
  r <- rowMeans(x)
  c <- colMeans(x)
  m <- mean(x)
  x - r - matrix(c, nrow(x), ncol(x), byrow = TRUE) + m
}

#' Normalize a matrix to unit Frobenius norm
#'
#' Scales a matrix so that its Frobenius norm (sqrt of sum of squared elements) is 1.
#'
#' @param S A numeric matrix to normalize
#' @return The normalized matrix with the same dimensions as S
#' @noRd
#' @keywords internal
norm_crossprod <- function(S) {
  if ((ss <- c(crossprod(c(S)))) < 1e-12) {
    stop("Matrix has a zero Frobenius norm, cannot normalize.")
  }
  S / sqrt(ss)     # ~15% quicker than sum(S^2)
}

#' Compute a matrix of pairwise inner products
#'
#' Creates a symmetric matrix where each element [i,j] is the inner product
#' between matrices S[[i]] and S[[j]]. The diagonal is set to 1 if normalized.
#'
#' @param S A list of numeric matrices of the same dimensions
#' @param fast Logical; if TRUE (default), uses a faster vectorized implementation
#' @param normalize Logical; if TRUE, sets diagonal to 1.
#' @return A symmetric matrix of inner products with dimensions length(S) × length(S)
#' @noRd
#' @keywords internal
compute_prodmat <- function(S, fast = TRUE, normalize = TRUE) {
  if (fast) {
    # Vectorize every matrix → columns of M
    M <- vapply(S, c, numeric(length(S[[1]])), USE.NAMES = FALSE)
    C <- crossprod(M)               # all inner products in one go
    if (normalize) {
      diag(C) <- 1                     # enforce exact ones
    }
    return(C)
  }
  
  # Fallback – simple but O(p²·n²)
  C <- matrix(0, length(S), length(S))
  for (i in seq_along(S)) {
    for (j in i:length(S)) {
      C[i, j] <- C[j, i] <- sum(S[[i]] * S[[j]])
    }
  }
  if (normalize) {
    diag(C) <- 1
  }
  C
}


#' @md
#' @rdname covstatis
#' @details
#' The `covstatis.list` method implements STATIS analysis for a list of covariance matrices.
#' 
#' The method operates as follows:
#' 1. Optionally double-centers each matrix (Gower transformation)
#' 2. Optionally normalizes each matrix to unit Frobenius norm
#' 3. Computes an RV matrix of inner products between matrices
#' 4. Determines optimal weights (`alpha`) via first eigenvector of the RV matrix. These weights are accessible in the output and can be useful for outlier diagnostics.
#' 5. Creates a compromise matrix as weighted sum of normalized matrices
#' 6. Performs eigendecomposition on the compromise matrix
#' 
#' The `partial_scores` in the returned object are a list of matrices (one for each subject), with each matrix containing the ROI-level factor scores (ROI x D).
#' 
#' Additional arguments can be passed via `...`:
#' \describe{
#'   \item{labels}{Optional character vector of labels for the rows/columns of the
#'                 covariance matrices. If NULL, tries to use row names from the
#'                 first matrix, or generates sequential labels.}
#'   \item{norm_method}{The normalization method to apply to each matrix. One of
#'                 `"frobenius"` (default when `normalize=TRUE`), `"mfa"`, or `"none"`.
#'                 Frobenius norm scales each matrix to have a total sum-of-squares of 1.
#'                 MFA (multiple factor analysis) scales each matrix so its first
#'                 eigenvalue is 1. `"none"` applies no normalization.}
#' }
#'
#' @examples
#' # Create a list of correlation matrices
#' Xs <- lapply(1:5, function(i) matrix(rnorm(10*10), 10, 10))
#' Xs <- lapply(Xs, cor)
#' 
#' # Apply COVSTATIS
#' res <- covstatis(Xs, ncomp=3)
#' 
#' # Project a new correlation matrix
#' new_mat <- cor(matrix(rnorm(10*10), 10, 10))
#' proj <- project_cov(res, new_mat)
#' 
#' @importFrom stats coefficients
#' @export
covstatis.list <- function(data, ncomp=2, normalize=TRUE, dcenter=TRUE, ...) {
  # Extract additional arguments from ...
  dots <- list(...)
  norm_method <- dots$norm_method %||% (if (isTRUE(normalize)) "frobenius" else "none")
  labels <- dots$labels

  norm_method <- match.arg(norm_method, c("frobenius", "mfa", "none"))
  nr <- sapply(data, nrow)
  nc <- sapply(data, ncol)
  
  assertthat::assert_that(all(nr == nr[1]))
  assertthat::assert_that(all(nc == nc[1]))
  assertthat::assert_that(all(nc[1] == nr[1]))
  
  block_labels <- names(data) %||% paste0("B_", seq_along(data))
  
  if (is.null(labels)) {
    labels <- row.names(data[[1]])
    if (is.null(labels)) {
      labels <- paste0("Obs_", 1:nr[1])
    }
  } else {
    assertthat::assert_that(length(labels) == nr[1])
  }
  
  # Preprocessing
  S <- data
  
  if (dcenter) {
    S <- lapply(S, double_center)
  }
  
  if (norm_method == "frobenius") {
    S <- lapply(S, norm_crossprod)
  } else if (norm_method == "mfa") {
    mfa_norm <- function(mat) {
      eig1 <- eigen(mat, symmetric = TRUE, only.values = TRUE)$values[1]
      if (eig1 <= 1e-10) {
        stop("First eigenvalue too small or non-positive for MFA normalization.")
      }
      mat / eig1
    }
    S <- lapply(S, mfa_norm)
  }
  
  # Cache vectorized preprocessed tables for RV-space diagnostics and
  # supplementary subject projection.
  table_matrix <- vapply(S, c, numeric(length(S[[1]])), USE.NAMES = FALSE)

  # RV-matrix and weights
  C <- crossprod(table_matrix)
  if (norm_method == "frobenius") {
    diag(C) <- 1
  }
  eigC <- eigen(C, symmetric = TRUE)
  alpha_vec <- eigC$vectors[,1]
  alpha_vec <- alpha_vec * sign(sum(alpha_vec))
  if (abs(sum(alpha_vec)) < 1e-12) {
    alpha_vec <- abs(alpha_vec)            # fallback in rare degenerate case
  }
  alpha <- alpha_vec / sum(alpha_vec)
  
  # Compromise and eigen-decomposition
  Sall <- Reduce("+", Map(`*`, S, alpha))
  fit <- eigen(Sall, symmetric = TRUE)
  
  # Component retention - take the first ncomp whose eigenvalues exceed tolerance
  tol <- max(fit$values) * 1e-8
  valid_comps <- which(fit$values > tol)
  keep <- valid_comps[seq_len(min(ncomp, length(valid_comps)))]
  
  scores <- fit$vectors[,keep,drop=FALSE] %*% diag(sqrt(fit$values[keep]))
  projmat <- fit$vectors[,keep,drop=FALSE] %*% diag(1/sqrt(fit$values[keep]))
  
  partial_scores <- lapply(S, function(mat) mat %*% projmat)
  names(partial_scores) <- block_labels
  
  # Create result object using projector - using 'basis' instead of 'v' for clarity
  ret <- multivarious::projector(fit$vectors[,keep,drop=FALSE], classes="covstatis", 
                                 roi_scores=scores,
                                 projmat=projmat,
                                 sdev=sqrt(fit$values[keep]),
                                 norm_method=norm_method, 
                                 dcenter=dcenter, 
                                 alpha=alpha, 
                                 C=C,
                                  eigC=eigC,
                                 compromise_full=Sall,
                                 partial_scores=partial_scores,
                                 table_matrix=table_matrix,
                                 block_labels=block_labels, 
                                 labels=labels)
  ret
}

# Internal helper to apply the same preprocessing steps to a new matrix
.pre_process_new_cov <- function(x, new_data) {
  assertthat::assert_that(is.matrix(new_data), msg = "`new_data` must be a matrix.")
  assertthat::assert_that(nrow(new_data) == nrow(coefficients(x)), 
                          msg="`new_data` must have the same dimensions as training data")
  assertthat::assert_that(isSymmetric(new_data), msg="`new_data` must be symmetric")
  
  if (x$dcenter) {
    new_data <- double_center(new_data)
  }
  
  if (x$norm_method == "frobenius") {
    new_data <- norm_crossprod(new_data)
  } else if (x$norm_method == "mfa") {
    eig1 <- eigen(new_data, symmetric = TRUE, only.values = TRUE)$values[1]
    if (eig1 <= 1e-10) {
      stop("First eigenvalue of a new matrix is too small or non-positive for MFA normalization.")
    }
    new_data <- new_data / eig1
  }
  
  new_data
}

#' @md
#' @rdname project_cov
#' @details
#' For `covstatis` objects, a new covariance/correlation matrix is transformed using 
#' the same preprocessing steps (centering and normalization) as were applied during 
#' model fitting, then projected onto the model space. This returns the "partial factor scores"
#' or ROI-level coordinates for the new subject.
#' 
#' @return A numeric matrix of projected scores with dimensions `nrow(new_data)` × `ncomp`.
#' @export
project_cov.covstatis <- function(x, new_data, ...) {
  processed_data <- .pre_process_new_cov(x, new_data)
  processed_data %*% x$projmat
}

#' Low-rank reconstruction of the compromise matrix
#'
#' Recreates the compromise matrix using only the selected components.
#'
#' @param x A fitted `covstatis` model
#' @param comp Integer vector of component indices to retain
#' @param ... Additional arguments (currently unused).
#' @return A matrix of the same size as each input matrix, representing
#'         a low-rank approximation of the compromise matrix in the preprocessed space
#'         (after double-centering and normalization, if they were applied)
#' @importFrom multivarious ncomp reconstruct
#' @export
reconstruct.covstatis <- function(x, comp = 1:multivarious::ncomp(x), ...) {
  chk::chk_numeric(comp)
  chk::chk_true(max(comp) <= multivarious::ncomp(x))

  V <- coefficients(x)[, comp, drop = FALSE]  # eigen-vectors
  lambda <- as.numeric(x$sdev[comp]^2)        # eigen-values
  V %*% (lambda * t(V))                       # fast diag trick
}

#' @rdname project_subjects
#' @description 
#' This function provides a comprehensive analysis for one or more new subject matrices 
#' by projecting them into the compromise space of a fitted `covstatis` model. It
#' computes several key metrics for each new matrix, following the projection logic
#' of `DISTATIS`.
#' @param subject_ids Optional character vector of identifiers for the new subjects.
#'   If not provided, names will be taken from `new_data` or generated automatically.
#' @details
#' The function performs the following steps for each new matrix:
#' 1.  Applies the same pre-processing (double-centering, normalization) as the original model.
#' 2.  Calculates ROI-level factor scores (`roi_scores`), representing the coordinates of each ROI in the compromise space.
#' 3.  Calculates subject-level scores (`subject_scores`) in the RV/table space by
#'     projecting the new table against the fitted interstructure axes.
#' 4.  Calculates `subject_cosines` indicating the alignment of the subject's scores with each compromise dimension.
#' 5.  Calculates a global `rv_coefficient` measuring the overall similarity between the new matrix and the group compromise matrix.
#' 6.  Calculates the `distance_to_compromise`, the Frobenius distance between the new matrix and the full compromise matrix.
#' 
#' @return A list containing the following elements:
#'   \item{subject_scores}{A matrix (`n_subjects` × `ncomp`) of subject-level coordinates in the compromise space.}
#'   \item{subject_cosines}{A matrix (`n_subjects` × `ncomp`) of cosines, indicating alignment with each dimension.}
#'   \item{scalar_summaries}{A `data.frame` with one row per subject, containing the `rv_coefficient` and `distance_to_compromise` (Frobenius distance to the full compromise).}
#'   \item{roi_scores}{A list of matrices, where each element contains the ROI-level scores for a subject.}
#'
#' @export
#' @md
project_subjects.covstatis <- function(x, new_data, subject_ids = NULL, ...) {
  if (!is.list(new_data)) {
    new_data <- list(new_data)
  }
  
  if (is.null(subject_ids)) {
    subject_ids <- names(new_data) %||% paste0("Subject_", seq_along(new_data))
  }
  
  # pre-compute compromise and projector
  compromise <- reconstruct(x)
  denom_compromise <- sum(compromise^2)
  U <- coefficients(x)
  P <- U %*% t(U)
  if (is.null(x$table_matrix)) {
    stop(
      "This covstatis object does not store the preprocessed training tables needed ",
      "for subject projection. Refit the model with the current muscal version.",
      call. = FALSE
    )
  }
  n_axes <- ncol(U)
  rv_vectors <- x$eigC$vectors[, seq_len(n_axes), drop = FALSE]
  rv_values <- x$eigC$values[seq_len(n_axes)]
  rv_scale <- numeric(n_axes)
  keep_scale <- abs(rv_values) > 1e-12
  rv_scale[keep_scale] <- x$sdev[keep_scale] / rv_values[keep_scale]
  
  res_list <- lapply(seq_along(new_data), function(i) {
    S_raw <- new_data[[i]]
    
    # 1. Pre-process
    S_new <- .pre_process_new_cov(x, S_raw)
    
    # 2. ROI coordinates ("Partial-F")
    F_new <- S_new %*% x$projmat
    
    # 3. Subject coordinates in RV/table space
    rv_profile <- drop(crossprod(x$table_matrix, c(S_new)))
    subject_scores <- drop(rv_profile %*% rv_vectors) * rv_scale
    
    # 4. Cosine alignment
    norm_g <- sqrt(sum(subject_scores^2))
    cosines <- if (norm_g > 1e-9) subject_scores / norm_g else subject_scores * 0
    
    # 5A. RV coefficient (against compromise)
    numer <- sum(S_new * compromise)
    denom_new <- sum(S_new^2)
    rv <- numer / sqrt(denom_new * denom_compromise)
    
    # 5B. Euclidean distance to compromise subspace
    S_proj <- P %*% S_new %*% P
    Residual <- S_new - S_proj
    dist <- sqrt(sum(Residual^2))
    
    list(
      subject_scores = subject_scores,
      cosines = cosines,
      rv_coefficient = rv,
      distance_to_compromise = dist,
      roi_scores = F_new
    )
  })
  
  # Consolidate results
  g_scores_mat <- do.call(rbind, lapply(res_list, `[[`, "subject_scores"))
  rownames(g_scores_mat) <- subject_ids
  colnames(g_scores_mat) <- paste0("Dim", 1:ncol(g_scores_mat))
  
  cosines_mat <- do.call(rbind, lapply(res_list, `[[`, "cosines"))
  rownames(cosines_mat) <- subject_ids
  colnames(cosines_mat) <- paste0("Dim", 1:ncol(cosines_mat))
  
  scalar_stats <- data.frame(
    subject_id = subject_ids,
    rv_coefficient = sapply(res_list, `[[`, "rv_coefficient"),
    distance_to_compromise = sapply(res_list, `[[`, "distance_to_compromise"),
    row.names = subject_ids
  )
  
  roi_scores_list <- lapply(res_list, `[[`, "roi_scores")
  names(roi_scores_list) <- subject_ids
  
  list(
    subject_scores = g_scores_mat,
    subject_cosines = cosines_mat,
    scalar_summaries = scalar_stats,
    roi_scores = roi_scores_list
  )
}

#' Retrieve subject scores from a covstatis object
#'
#' This internal helper returns the subject-level coordinates used for plotting
#' and supplementary covariate projection. The historical "G-score" name is
#' retained for compatibility, but the coordinates live in the RV/table space.
#' Barycentric means of partial ROI scores collapse to zero after double-centering,
#' so they are not used here.
#'
#' @param x A fitted `covstatis` object, which contains `partial_scores`.
#' @return A numeric matrix of table scores with dimensions `n_tables` x `n_components`.
#' @noRd
#' @keywords internal
.get_table_scores <- function(x) {
  rv_subject_scores(x)
}

#' @rdname project_covariate
#' @description
#' This function projects a table-level or feature-level covariate into the
#' space of a fitted `covstatis` model, treating it as a supplementary variable
#' without re-fitting the model.
#' @param side Which side the covariate belongs to: `"table"` for one value per
#'   input table/block, or `"feature"` for one value per matrix feature.
#' @param what For `side = "table"`, `"dimension"` (default) returns cosine / β
#'   per compromise dimension and `"feature"` returns a feature-wise pattern.
#'   `"observation"` is accepted as a backwards-compatible alias for `"feature"`.
#'   For `side = "feature"`, only `"dimension"` is currently supported.
#' @param scale  `"cosine"` (−1…1) or `"beta"` (regression coeff.).
#' @return If `what = "dimension"`, a named numeric vector of length `ncomp(x)`.
#'         If `side = "table"` and `what = "feature"`, a `p × ncomp` matrix.
#' @details
#' The interpretation of the output depends on `side`:
#' - **Table-side covariate**: `y` has length equal to the number of input
#'   tables/blocks. Dimension-wise output quantifies how strongly the covariate
#'   aligns with the RV/table-space axes. Feature-wise output is a weighted sum
#'   of partial factor scores and gives a feature pattern associated with the
#'   covariate.
#' - **Feature-side covariate**: `y` has length equal to the matrix dimension
#'   (`p`). Dimension-wise output quantifies how strongly the covariate aligns
#'   with the compromise feature-space axes.
#' 
#' @examples
#' # Create a list of correlation matrices
#' Xs <- lapply(1:5, function(i) matrix(rnorm(10*10), 10, 10))
#' Xs <- lapply(Xs, cor)
#' 
#' # Apply COVSTATIS
#' res <- covstatis(Xs, ncomp=3)
#' 
#' # Table-side covariate: one value per table/block
#' y <- rnorm(length(Xs))
#' 
#' # Project the covariate to get dimension-wise coordinates
#' dim_cos <- project_table_covariate(res, y, what = "dimension", scale = "cosine")
#' 
#' # Project the covariate to get a feature-wise pattern
#' feature_beta <- project_table_covariate(res, y, what = "feature", scale = "beta")
#' 
#' # Feature-side covariate: one value per feature/row-column entity
#' z <- rnorm(nrow(Xs[[1]]))
#' feature_dim_cos <- project_feature_covariate(res, z, scale = "cosine")
#'
#' @export
project_covariate.covstatis <- function(x, y, side = c("table", "feature"),
                              what = NULL, scale = c("cosine", "beta"), ...) {
  side <- match.arg(side)
  scale <- match.arg(scale)
  
  if (side == "table") {
    if (is.null(what)) {
      what <- "dimension"
    }
    what <- match.arg(what, c("dimension", "feature", "observation"))
    if (what == "observation") {
      what <- "feature"
    }
    
    ns <- length(x$partial_scores)
    chk::chk_vector(y)
    chk::chk_equal(length(y), ns)
    
    G <- .get_table_scores(x)        # ns × D
    if (what == "dimension") {
      if (scale == "cosine") {
        numer <- as.numeric(t(G) %*% y)            # length D
        denom <- sqrt(colSums(G^2)) * sqrt(sum(y^2))
        denom[denom < 1e-12] <- Inf
        out   <- numer / denom
      } else { # beta
        denom <- colSums(G^2)
        denom[denom < 1e-12] <- Inf
        out <- as.numeric(t(G) %*% y / denom)
      }
      out[!is.finite(out)] <- 0
      names(out) <- paste0("Dim", seq_along(out))
      return(out)
    }
    
    ## --- feature-wise pattern -----------------------------------------
    R  <- nrow(x$partial_scores[[1]])
    D  <- ncol(x$partial_scores[[1]])
    out <- matrix(0, R, D)
    for (k in seq_len(ns)) out <- out + y[k] * x$partial_scores[[k]]
    if (scale == "cosine") {
      norm_out <- sqrt(sum(out^2))
      if (norm_out > 1e-12) out <- out / norm_out
    } else {                     # beta-style scaling
      denom <- colSums(G^2)      # D-vector
      denom[denom < 1e-12] <- Inf
      out   <- out / matrix(denom, R, D, byrow = TRUE)
      out[!is.finite(out)] <- 0
    }
    dimnames(out) <- list(x$labels, paste0("Dim", 1:D))
    return(out)
  }
  
  if (is.null(what)) {
    what <- "dimension"
  }
  what <- match.arg(what, c("dimension"))
  feature_scores <- x$roi_scores
  chk::chk_vector(y)
  chk::chk_equal(length(y), nrow(feature_scores))
  
  if (scale == "cosine") {
    numer <- as.numeric(t(feature_scores) %*% y)
    denom <- sqrt(colSums(feature_scores^2)) * sqrt(sum(y^2))
    denom[denom < 1e-12] <- Inf
    out <- numer / denom
  } else {
    denom <- colSums(feature_scores^2)
    denom[denom < 1e-12] <- Inf
    out <- as.numeric(t(feature_scores) %*% y / denom)
  }
  out[!is.finite(out)] <- 0
  names(out) <- paste0("Dim", seq_along(out))
  out
}

# Compatibility shim for historical internal name.
.get_G_scores <- function(x) {
  .get_table_scores(x)
}

#' Subject coordinates in the RV/table space
#'
#' Returns a K × D matrix whose columns are the subject (table) scores
#' associated with each compromise dimension.
#'
#' @param x A fitted `covstatis` object.
#' @return A K x D matrix of subject scores in the RV space.
#' @keywords internal
#' @noRd
rv_subject_scores <- function(x) {
  eigC <- x$eigC
  V    <- eigC$vectors[, seq_len(ncol(coefficients(x))), drop = FALSE]
  sigma <- x$sdev[seq_len(ncol(coefficients(x)))]
  sweep(V, 2, sigma, `*`)            # K × D  (same Σ as in ROI space)
}

#' Verify duality between ROI- and RV-space coordinates
#' 
#' This function checks whether the squared norms of the ROI-space coordinates 
#' (F-scores) are equal to the squared norms of the RV-space coordinates 
#' (T-scores) for each dimension, which is a fundamental property of COVSTATIS.
#'
#' @param x A fitted `covstatis` object.
#' @param tol The tolerance for the comparison.
#' @return `TRUE` if the duality property holds within the given tolerance.
#' @export
#' @examples
#' # Create a list of correlation matrices
#' Xs <- lapply(1:5, function(i) matrix(rnorm(10*10), 10, 10))
#' Xs <- lapply(Xs, cor)
#' # Apply COVSTATIS
#' res <- covstatis(Xs, ncomp=3)
#' # Check duality
#' check_duality(res)
check_duality <- function(x, tol = 1e-10) {
  F_scores  <- x$roi_scores                         # ROI × D  (U Σ)
  T_scores  <- rv_subject_scores(x)        # subj × D (V Σ)
  isTRUE(all.equal(colSums(F_scores^2),
                   colSums(T_scores^2),
                   tolerance = tol))
}

#' Extract ROI map for a compromise dimension
#' 
#' @param x A fitted `covstatis` object.
#' @param d The dimension to extract.
#' @return A numeric vector representing the ROI map for the specified dimension.
#' @export
roi_map <- function(x, d) {
  x$roi_scores[, d]
}

#' Plot subject scores in RV space
#'
#' @param x A fitted `covstatis` object.
#' @param dim A numeric vector of length 2 specifying the dimensions to plot.
#' @importFrom graphics par plot
#' @export
plot_subjects <- function(x, dim = c(1,2)) {
  T_scores <- rv_subject_scores(x)
  plot(T_scores[,dim], main = "Subject scores (RV space)",
       xlab = paste("Dim", dim[1]), ylab = paste("Dim", dim[2]))
}

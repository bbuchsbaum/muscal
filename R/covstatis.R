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
  n <- nrow(x)
  H <- diag(n) - matrix(1/n, n, n)          # Ξ with uniform masses
  H %*% x %*% H                     # (1/2) Ξ R Ξ, 1/2 factor removed for covariance
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
  S / sqrt(c(crossprod(c(S))))     # ~15% quicker than sum(S^2)
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
#' 4. Determines optimal weights via first eigenvector of the RV matrix
#' 5. Creates a compromise matrix as weighted sum of normalized matrices
#' 6. Performs eigendecomposition on the compromise matrix
#' 
#' @param labels Optional character vector of labels for the rows/columns of the
#'               covariance matrices. If NULL, tries to use row names from the 
#'               first matrix, or generates sequential labels.
#' @param norm_method The normalization method to apply to each matrix. One of
#'               `"frobenius"` (default), `"mfa"`, or `"none"`. Frobenius norm
#'               scales each matrix to have a total sum-of-squares of 1. MFA
#'               (multiple factor analysis) scales each matrix so its first
#'               eigenvalue is 1. `"none"` applies no normalization.
#' @param dcenter Logical; if TRUE (default), each matrix is double-centered before
#'               analysis. This removes constant modes before cross-table comparison
#'               (standard in STATIS tradition).
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
#' @export
covstatis.list <- function(data, ncomp=2, 
                           norm_method=c("frobenius", "mfa", "none"), 
                           dcenter=TRUE, labels=NULL) {
  
  norm_method <- match.arg(norm_method)
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
    mfa_norm <- function(mat) mat / eigen(mat, symmetric = TRUE, only.values = TRUE)$values[1]
    S <- lapply(S, mfa_norm)
  }
  
  # RV-matrix and weights
  C <- compute_prodmat(S, normalize = (norm_method == "frobenius"))
  alpha <- abs(eigen(C)$vectors[,1])
  alpha <- alpha/(sum(alpha))
  
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
                                 s=scores,
                                 projmat=projmat,
                                 sdev=sqrt(fit$values[keep]),
                                 norm_method=norm_method, 
                                 dcenter=dcenter, 
                                 alpha=alpha, 
                                 C=C,
                                 partial_scores=partial_scores,
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
    if (abs(eig1) < 1e-10) {
      stop("First eigenvalue of a new matrix is too small for MFA normalization.")
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
project_cov.covstatis <- function(x, new_data) {
  processed_data <- .pre_process_new_cov(x, new_data)
  processed_data %*% x$projmat
}

#' Low-rank reconstruction of the compromise matrix
#'
#' Recreates the compromise matrix using only the selected components.
#'
#' @param x A fitted `covstatis` model
#' @param comp Integer vector of component indices to retain
#' @return A matrix of the same size as each input matrix, representing 
#'         a low-rank approximation of the compromise matrix in the preprocessed space
#'         (after double-centering and normalization, if they were applied)
#' @export
reconstruct.covstatis <- function(x, comp = 1:multivarious::ncomp(x)) {
  chk::chk_numeric(comp)
  chk::chk_true(max(comp) <= multivarious::ncomp(x))

  V <- coefficients(x)[, comp, drop = FALSE]  # eigen-vectors
  lambda <- as.numeric(x$sdev[comp]^2)        # eigen-values
  V %*% (lambda * t(V))                       # fast diag trick
}

#' Project and Summarize New Data in a COVSTATIS Space
#' 
#' @description 
#' This function provides a comprehensive analysis for one or more new subject matrices 
#' by projecting them into the compromise space of a fitted `covstatis` model. It
#' computes several key metrics for each new matrix, following the projection logic
#' of `DISTATIS`.
#' 
#' @param x A fitted `covstatis` object.
#' @param new_data A single matrix or a list of matrices to project. Each matrix must
#'   have the same dimensions as the data used to train the model.
#' @param subject_ids Optional character vector of identifiers for the new subjects.
#'   If not provided, names will be taken from `new_data` or generated automatically.
#' @param ... other arguments (not used).
#' 
#' @details
#' The function performs the following steps for each new matrix:
#' 1.  Applies the same pre-processing (double-centering, normalization) as the original model.
#' 2.  Calculates ROI-level factor scores (`roi_scores`), representing the coordinates of each ROI in the compromise space.
#' 3.  Calculates subject-level scores (`subject_scores` or "g-scores") as the barycentric mean of the ROI scores.
#' 4.  Calculates `subject_cosines` indicating the alignment of the subject's scores with each compromise dimension.
#' 5.  Calculates a global `rv_coefficient` measuring the overall similarity between the new matrix and the group compromise matrix.
#' 6.  Calculates the `distance_to_compromise`, the Frobenius distance between the new matrix and its projection onto the compromise subspace.
#' 
#' @return A list containing the following elements:
#'   \item{subject_scores}{A matrix (`n_subjects` × `ncomp`) of subject-level coordinates in the compromise space.}
#'   \item{subject_cosines}{A matrix (`n_subjects` × `ncomp`) of cosines, indicating alignment with each dimension.}
#'   \item{scalar_summaries}{A `data.frame` with one row per subject, containing the `rv_coefficient` and `distance_to_compromise`.}
#'   \item{roi_scores}{A list of matrices, where each element contains the ROI-level scores for a subject.}
#'
#' @export
#' @md
project_subjects <- function(x, ...) {
  UseMethod("project_subjects")
}

#' @rdname project_subjects
#' @export
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
  
  res_list <- lapply(seq_along(new_data), function(i) {
    S_raw <- new_data[[i]]
    
    # 1. Pre-process
    S_new <- .pre_process_new_cov(x, S_raw)
    
    # 2. ROI coordinates ("Partial-F")
    F_new <- S_new %*% x$projmat
    
    # 3. Subject coordinate ("G-scores")
    g_scores <- colMeans(F_new)
    
    # 4. Cosine alignment
    norm_g <- sqrt(sum(g_scores^2))
    cosines <- if (norm_g > 1e-9) g_scores / norm_g else g_scores * 0
    
    # 5A. RV coefficient
    numer <- sum(S_new * compromise)
    denom_new <- sum(S_new^2)
    rv <- numer / sqrt(denom_new * denom_compromise)
    
    # 5B. Euclidean distance
    S_proj <- P %*% S_new %*% P
    Residual <- S_new - S_proj
    dist <- sqrt(sum(Residual^2))
    
    list(
      g_scores = g_scores,
      cosines = cosines,
      rv_coefficient = rv,
      distance_to_compromise = dist,
      roi_scores = F_new
    )
  })
  
  # Consolidate results
  g_scores_mat <- do.call(rbind, lapply(res_list, `[[`, "g_scores"))
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


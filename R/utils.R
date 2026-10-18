#' Prepare Block Preprocessors and Apply to Data
#' 
#' This utility function handles the setup and application of preprocessing steps 
#' for list-based multi-block data, accommodating NULL, a single preprocessor 
#' definition, or a list of definitions.
#' 
#' @param data_list A list of data matrices (blocks).
#' @param preproc_arg The user-provided `preproc` argument. Can be NULL, a single
#'   `prepper`/`pre_processor` object, or a list of them.
#' @param check_consistent_ncol Logical (default: TRUE). If TRUE, checks if all blocks 
#'   have the same number of columns after preprocessing and issues a warning if not.
#' 
#' @return A list containing:
#'   * `proclist`: A list of fitted `pre_processor` objects (one per block), or NULL 
#'       if `preproc_arg` was NULL.
#'   * `Xp`: The list of preprocessed data blocks (identical to `data_list` if 
#'       `proclist` is NULL).
#'   * `p_post`: The number of columns in the first block after preprocessing. If 
#'       `check_consistent_ncol` is TRUE, this assumes all blocks have this dimension.
#' 
#' @importFrom multivarious fit fresh
#' @keywords internal
#' @noRd
prepare_block_preprocessors <- function(data_list, preproc_arg, check_consistent_ncol = TRUE) {
  
  S <- length(data_list)
  Xp <- data_list # Default: use original data if no preproc
  proclist <- NULL
  p_post <- if (S > 0) ncol(data_list[[1]]) else 0 # Initial dimension
  block_names <- names(data_list)
  if (is.null(block_names)) block_names <- paste0("Block_", 1:S)
  
  if (!is.null(preproc_arg)) {
    if (!requireNamespace("multivarious", quietly = TRUE)) {
      stop("Package 'multivarious' needed for preprocessing. Please install it.", call. = FALSE)
    }
    
    # Check if preproc_arg is a single definition or a list
    if (!is.list(preproc_arg) ||
        inherits(preproc_arg, "prepper") ||
        inherits(preproc_arg, "pre_processor")) {
      # Case 1: Single preproc definition - replicate for each block
      message("Applying the same preprocessor definition independently to each block.")
      single_preproc_def <- preproc_arg
      proclist <- vector("list", S)
      for (i in seq_len(S)) {
          proclist[[i]] <- multivarious::fresh(single_preproc_def) %>%
                               multivarious::fit(data_list[[i]])
      }

    } else if (is.list(preproc_arg)) {
      # Case 2: List of preproc definitions
      if (length(preproc_arg) != S) {
          stop(sprintf("If 'preproc' is a list, its length (%d) must match the number of data blocks (%d).",
                       length(preproc_arg), S), call. = FALSE)
      }
      message("Applying preprocessor list: one definition per block.")
      proclist <- vector("list", S)
      for (i in seq_len(S)) {
          if (!inherits(preproc_arg[[i]], "prepper") &&
              !inherits(preproc_arg[[i]], "pre_processor")) {
              stop(sprintf("Element %d of 'preproc' list is not a valid prepper/pre_processor object.", i), call.=FALSE)
          }
          proclist[[i]] <- multivarious::fresh(preproc_arg[[i]]) %>%
                               multivarious::fit(data_list[[i]])
      }

    } else {
        stop("'preproc' must be NULL, a single prepper/pre_processor object, or a list of them.", call.=FALSE)
    }

    names(proclist) <- block_names

    # Apply the fitted preprocessors
    Xp <- vector("list", S)
    names(Xp) <- block_names
    for (i in seq_along(data_list)) {
      p_i_fitted <- proclist[[i]]
      Xp[[i]] <- multivarious::transform(p_i_fitted, data_list[[i]])
    }
    
    # Update post-preprocessing dimension
    p_post <- if (S > 0) ncol(Xp[[1]]) else 0

    # Check for consistent dimensions after preprocessing if requested
    if (check_consistent_ncol) {
        dims_post_preproc <- sapply(Xp, ncol)
        if (S > 0 && !all(dims_post_preproc == p_post)) {
            warning("Preprocessing resulted in blocks with different numbers of columns. Subsequent steps might require consistent dimensions.", call.=FALSE)
        }
    }

  } else {
    # No preprocessing: Ensure Xp still has names 
    names(Xp) <- block_names
    # Check original dimensions if requested
    if (check_consistent_ncol) {
        dims_orig <- sapply(data_list, ncol)
        if (S > 0 && !all(dims_orig == p_post)) {
             warning("Input blocks have different numbers of columns and no preprocessing was applied. Subsequent steps might require consistent dimensions.", call.=FALSE)
        }
    }
  }
  
  return(list(proclist = proclist, Xp = Xp, p_post = p_post))
}

#' Select an SVD backend for exact vs. truncated decompositions
#'
#' @param mat Numeric matrix-like object.
#' @param ncomp Requested number of components.
#' @param exact_dim_threshold Minimum dimension at or below which base SVD is used.
#' @param truncation_fraction Maximum fraction of the smaller matrix dimension to
#'   request from a truncated solver.
#'
#' @return `"base"` or `"irlba"`.
#' @keywords internal
#' @noRd
.muscal_svd_method <- function(mat,
                               ncomp = 1L,
                               exact_dim_threshold = 2L,
                               truncation_fraction = 0.5) {
  mdim <- min(dim(mat))
  if (mdim <= exact_dim_threshold) {
    return("base")
  }

  ncomp <- min(as.integer(ncomp), mdim)
  if (ncomp >= truncation_fraction * mdim) {
    "base"
  } else {
    "irlba"
  }
}

#' Materialize Fitted Identity Preprocessors When None Were Requested
#'
#' @param data_list A list of raw data blocks.
#' @param proclist A fitted preprocessor list, or `NULL`.
#'
#' @return A named list of fitted `pre_processor` objects, one per block.
#' @keywords internal
#' @noRd
.muscal_materialize_block_preprocessors <- function(data_list, proclist) {
  if (!is.null(proclist)) {
    return(proclist)
  }

  block_names <- names(data_list)
  if (is.null(block_names)) {
    block_names <- paste0("Block_", seq_along(data_list))
  }

  fitted_pass <- lapply(data_list, function(x) {
    multivarious::fit(multivarious::pass(), x)
  })
  names(fitted_pass) <- block_names
  fitted_pass
}

#' Identify Significant Components via RMT and Inter-block Coherence Tests
#'
#' Assess declared supplied-noise spectral edges and optional common-feature
#' loading coherence. Arbitrary multiblock loadings are not assumed to share a
#' feature space, and an asymptotic support edge is not a calibrated p-value.
#'
#' @param fit A fitted multi-block model object containing at minimum a `V_list`
#'   attribute (list of loading matrices) and optionally `sdev` (singular values).
#' @param n Integer; the number of observations (rows) used to fit the model.
#' @param k_vec Optional integer vector of block dimensions (number of columns
#'   per original data block). If NULL, inferred from loading-matrix row counts.
#'   Supplied counts and block-index metadata must agree with those fitted feature
#'   counts, including any preprocessing that changes the feature dimension.
#' @param alpha Finite probability strictly between 0 and 1; the threshold for
#'   explicitly enabled experimental ICC (default 0.05). Explicit alpha without
#'   ICC warns that it is inapplicable; an MP support edge is not alpha-level.
#' @param check_rmt Logical; if TRUE (default), performs the Marchenko-Pastur
#'   edge test to check if eigenvalues exceed the noise threshold.
#' @param tail_frac Retained for compatibility. Automatic noise estimation from
#'   truncated fitted spectra is unavailable; this argument does not estimate noise.
#' @param noise_variance Positive finite noise variance, supplied independently
#'   in covariance-eigenvalue units. No automatic estimate is assumed.
#' @param spectrum Spectral units: `"covariance"` means `sdev^2` are covariance
#'   eigenvalues; `"singular"` means raw singular values, divided by `n_eff`
#'   after squaring. The default `"unknown"` leaves RMT unavailable.
#' @param noise_model `"isotropic"` explicitly declares the IID isotropic noise
#'   model. The default `"unknown"` leaves RMT unavailable, including for
#'   arbitrary weighted MFA compromise spectra.
#' @param n_eff Positive effective degrees of freedom of the covariance
#'   normalization, at most `n`; defaults to `n` for an uncentered model.
#'   For a centered unweighted sample covariance, supply `n - 1`.
#' @param check_icc Logical; explicitly request the experimental common-feature
#'   coherence criterion (default FALSE). It is independent of RMT applicability.
#' @param common_features Unique character feature identities declaring a common
#'   coordinate domain for ICC. Every loading matrix must have exactly these row
#'   identities; rows are matched by identity before comparison. Equal widths
#'   alone do not establish a common space.
#'
#' @return A list with the following elements:
#'   \describe{
#'     \item{keep}{Component indices passing all enabled criteria. A requested
#'       unavailable criterion cannot retain components. With both disabled,
#'       no components are declared retained.}
#'     \item{rmt_pass}{Logical edge decisions, or NA if RMT is unavailable/disabled.}
#'     \item{icc_pass}{Experimental ICC decisions, or NA when disabled.}
#'     \item{icc}{Numeric vector of inter-block coherence values per component.}
#'     \item{icc_pvalue}{Numeric vector of p-values for the ICC test.}
#'     \item{mp_edge}{The Marchenko-Pastur edge threshold (NA if RMT skipped).}
#'     \item{sigma2_est}{Supplied noise variance (legacy field name; NA if unavailable).}
#'     \item{lambda}{Eigenvalues (squared singular values).}
#'     \item{n, k_vec, alpha}{Input parameters echoed back.}
#'     \item{rmt_status, rmt_reason}{Whether the requested edge was evaluated,
#'       unavailable, or explicitly disabled, and why. Unavailable is never a pass.}
#'     \item{gamma, n_eff, spectrum, noise_model}{Resolved aspect ratio and declared units/model.}
#'     \item{icc_status}{`"disabled"` or `"evaluated_experimental"`.}
#'   }
#'
#' @details
#' The supplied-noise RMT criterion uses the asymptotic Marchenko-Pastur upper
#' support edge. This is not an alpha-level finite-sample significance test.
#' Its isotropic covariance representation must be declared explicitly; no
#' noise variance is inferred from a few retained fitted components.
#'
#' The ICC test measures the squared cosine similarity of loading vectors
#' across pairs of explicitly matched common-feature blocks. Its Gaussian null
#' approximation with Bonferroni correction remains experimental; feature-domain
#' validation does not establish statistical calibration.
#'
#' @importFrom stats median quantile pnorm p.adjust
#' @importFrom multivarious block_indices
#' @export
significant_components <- function(fit, n, k_vec = NULL,
                                   alpha = 0.05, 
                                   check_rmt = TRUE, 
                                   tail_frac = 0.3,
                                   noise_variance = NULL,
                                   spectrum = c("unknown", "covariance", "singular"),
                                   noise_model = c("unknown", "isotropic"),
                                   n_eff = n,
                                   check_icc = FALSE,
                                   common_features = NULL) {
  
  spectrum <- match.arg(spectrum)
  noise_model <- match.arg(noise_model)
  chk::chk_flag(check_rmt)
  chk::chk_flag(check_icc)
  if (!is.numeric(alpha) || length(alpha) != 1L || !is.finite(alpha) ||
      alpha <= 0 || alpha >= 1) {
    stop("`alpha` must be a finite scalar strictly between 0 and 1.", call. = FALSE)
  }
  if (!missing(alpha) && !check_icc) {
    warning("`alpha` applies only to explicitly enabled ICC; a support edge is not an alpha-level test.",
            call. = FALSE)
  }

  # ---- 0. Input Extraction and Validation -----------------------------
  
  # Extract V_list (essential for ICC)
  V_list <- NULL
  if (!is.null(fit$V_list)) {
      V_list <- fit$V_list
  } else if (!is.null(attr(fit, "V_list"))) {
      V_list <- attr(fit, "V_list")
  }
  if (is.null(V_list) || !is.list(V_list) || length(V_list) == 0) {
      stop("Could not extract a valid V_list (list of loading matrices) from the fit object.")
  }
  S <- length(V_list)
  if (any(!vapply(V_list, function(V) is.matrix(V) && is.numeric(V) &&
                  nrow(V) > 0L && ncol(V) > 0L && all(is.finite(V)), logical(1)))) {
    stop("V_list must contain nonempty finite numeric loading matrices.", call. = FALSE)
  }
  
  feature_counts <- vapply(V_list, nrow, integer(1))
  if (!is.null(fit$block_indices) &&
      (!is.list(fit$block_indices) ||
       !identical(unname(lengths(fit$block_indices)), unname(feature_counts)))) {
    stop("Stored block indices disagree with loading feature counts.", call. = FALSE)
  }
  if (is.null(k_vec)) k_vec <- feature_counts
  if (length(k_vec) != S) {
    stop("Length of k_vec does not match the number of blocks in V_list.", call. = FALSE)
  }
  if (!is.numeric(k_vec) || any(!is.finite(k_vec)) || any(k_vec <= 0) ||
      any(k_vec != feature_counts)) {
    stop("k_vec must match the positive loading feature counts.", call. = FALSE)
  }

  # Extract sdev (singular values, optional for RMT)
  sdev <- NULL
  if (check_rmt) {
      if (!is.null(fit$sdev)) {
          sdev <- fit$sdev
      } else if (!is.null(attr(fit, "sdev"))) {
          sdev <- attr(fit, "sdev")
      }
      if (is.null(sdev)) {
           warning("Singular values (sdev) not found in fit object. Skipping RMT test.", call.=FALSE)
           # Keep the requested criterion unavailable; missing data is not a pass.
      }
  }
  
  # Get number of components (r)
  r <- ncol(V_list[[1]]) # Assume all V in list have same ncol
  if (any(sapply(V_list, ncol) != r)) {
      stop("Inconsistent number of components found across matrices in V_list.")
  }
  
  # Validate n
  chk::chk_number(n)
  chk::chk_gt(n, 1)
  if (!is.numeric(n_eff) || length(n_eff) != 1L || !is.finite(n_eff) ||
      n_eff <= 0 || n_eff > n) {
    stop("`n_eff` must be positive, finite and no greater than n.", call. = FALSE)
  }

  # A supplied-noise edge is evaluated only for declared spectral units/model.
  # Automatic noise estimation from truncated fitted roots is not qualified.
  pass_rmt <- rep(NA, r)
  mp_edge <- sigma2_est <- NA_real_
  lambda <- numeric(0)
  gamma <- sum(k_vec) / n_eff
  rmt_status <- if (check_rmt) "unavailable" else "disabled"
  rmt_reason <- if (check_rmt) {
    "RMT requires supplied noise variance, declared spectral units and isotropic noise."
  } else "RMT was explicitly disabled."
  if (check_rmt && !is.null(sdev)) {
    if (!is.numeric(sdev) || length(sdev) != r || any(!is.finite(sdev)) || any(sdev < 0)) {
      stop("RMT requires one finite nonnegative sdev per component.", call. = FALSE)
    }
    if (!is.null(noise_variance) &&
        (!is.numeric(noise_variance) || length(noise_variance) != 1L ||
         !is.finite(noise_variance) || noise_variance <= 0)) {
      stop("`noise_variance` must be a positive finite scalar.", call. = FALSE)
    }
    lambda <- sdev^2
    if (spectrum == "singular") lambda <- lambda / n_eff
    if (any(!is.finite(lambda))) {
      stop("Squared spectral values must be finite; rescale the input.", call. = FALSE)
    }
    if (!is.null(noise_variance) && spectrum != "unknown" && noise_model == "isotropic") {
      mp_edge <- noise_variance * (1 + sqrt(gamma))^2
      if (!is.finite(mp_edge)) stop("MP edge is not finite; rescale the input.", call. = FALSE)
      sigma2_est <- noise_variance
      pass_rmt <- lambda > mp_edge
      rmt_status <- "evaluated"
      rmt_reason <- "Supplied-noise asymptotic support edge; not a finite-sample p-value."
    }
  } else if (!check_rmt) {
    message("RMT check skipped.")
  }

  icc <- pval <- rep(NA_real_, r)
  pass_icc <- rep(NA, r)
  icc_status <- "disabled"
  if (check_icc) {
    if (!is.character(common_features) || length(common_features) < 2L ||
        anyNA(common_features) || any(common_features == "") || anyDuplicated(common_features)) {
      stop("ICC requires unique nonempty `common_features` identities for at least two features.", call. = FALSE)
    }
    if (S < 2L) stop("ICC requires at least two blocks in a common feature domain.", call. = FALSE)
    common <- lapply(V_list, function(V) {
      ids <- rownames(V)
      if (is.null(ids) || anyNA(ids) || anyDuplicated(ids) ||
          length(ids) != length(common_features) || !setequal(ids, common_features)) {
        stop("ICC blocks must declare the same common feature domain using unique row identities.", call. = FALSE)
      }
      V <- V[match(common_features, ids), , drop = FALSE]
      # Normalize without squaring very large/small original magnitudes.
      amplitude <- apply(abs(V), 2L, max)
      if (any(amplitude == 0)) stop("ICC cannot compare zero-norm loading directions.", call. = FALSE)
      V <- sweep(V, 2L, amplitude, "/")
      sweep(V, 2L, sqrt(colSums(V^2)), "/")
    })
    pairs <- utils::combn(seq_len(S), 2L, simplify = FALSE)
    squared <- lapply(pairs, function(pair) {
      pmin(1, colSums(common[[pair[[1L]]]] * common[[pair[[2L]]]])^2)
    })
    icc <- Reduce(`+`, squared) / length(pairs)
    p <- length(common_features)
    z <- (icc - 1 / p) / sqrt(2 / (length(pairs) * p^2))
    pval <- stats::pnorm(z, lower.tail = FALSE)
    pass_icc <- stats::p.adjust(pval, method = "bonferroni") < alpha
    icc_status <- "evaluated_experimental"
  } else if (!is.null(common_features)) {
    stop("`common_features` requires `check_icc = TRUE`.", call. = FALSE)
  }

  # ---- 3. Combined Decision ------------------------------------------
  rmt_filter <- if (rmt_status == "disabled") rep(TRUE, r) else pass_rmt
  icc_filter <- if (check_icc) pass_icc else rep(TRUE, r)
  keep <- if (!check_rmt && !check_icc) integer(0) else which(rmt_filter & icc_filter)
  
  # ---- 4. Return Results --------------------------------------------
  return(list(keep = keep,
              rmt_pass = pass_rmt,
              icc_pass = pass_icc,
              icc = icc,
              icc_pvalue = pval,
              mp_edge = mp_edge,
              sigma2_est = sigma2_est,
              lambda = lambda,
              n = n,
              k_vec = k_vec,
              alpha = alpha,
              rmt_status = rmt_status,
              rmt_reason = rmt_reason,
              gamma = gamma,
              n_eff = n_eff,
              spectrum = spectrum,
              noise_model = noise_model,
              icc_status = icc_status))
}

#' @noRd
`%||%` <- function(a, b) {
  if (is.null(a)) b else a
}

#' Calculate Within-Class Scatter Matrix
#'
#' Computes the within-class scatter matrix (Sw), which is the sum of scatter
#' matrices for each class.
#'
#' @param X A numeric matrix where rows are observations and columns are features.
#' @param class_vector A factor or vector specifying the class for each row of X.
#' @return The within-class scatter matrix (p x p).
#' @noRd
#' @keywords internal
within_class_scatter <- function(X, class_vector) {
  p <- ncol(X)
  # Using as.data.frame is necessary for split to work on a matrix by rows
  X_split_by_class <- split(as.data.frame(X), class_vector)
  
  # Calculate scatter matrix for each class and sum them up
  scatter_matrices <- lapply(X_split_by_class, function(class_df) {
    # scale() is efficient for centering
    crossprod(scale(as.matrix(class_df), center = TRUE, scale = FALSE))
  })
  
  # Sum the scatter matrices from all classes
  Reduce("+", scatter_matrices)
}


#' Calculate Between-Class Scatter Matrix
#'
#' Computes the between-class scatter matrix (Sb), which measures the scatter
#' of class means around the grand mean.
#'
#' @param X A numeric matrix where rows are observations and columns are features.
#' @param class_vector A factor or vector specifying the class for each row of X.
#' @param grand_mean A numeric vector representing the grand mean of X.
#' @return The between-class scatter matrix (p x p).
#' @importFrom stats aggregate
#' @noRd
#' @keywords internal
between_class_scatter <- function(X, class_vector, grand_mean) {
  p <- ncol(X)
  Sb <- matrix(0, p, p)
  
  # Get class counts and means efficiently
  class_counts <- table(class_vector)
  # aggregate is a good way to get class means
  class_means <- aggregate(X ~ class_vector, FUN = mean)
  
  # Iterate over classes to compute Sb
  for (i in seq_along(class_counts)) {
    class_name <- names(class_counts)[i]
    n_c <- class_counts[class_name]
    # class_means has 'class_vector' as first column
    mu_c <- as.numeric(class_means[i, -1])
    
    diff_vec <- mu_c - grand_mean
    # Add this class's contribution to Sb
    Sb <- Sb + n_c * tcrossprod(diff_vec)
  }
  
  Sb
}

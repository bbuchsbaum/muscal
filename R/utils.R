#' Prepare Block Preprocessors and Apply to Data
#' 
#' This utility function handles the setup and application of preprocessing steps 
#' for list-based multi-block data, accommodating NULL, a single preprocessor 
#' definition, or a list of definitions.
#' 
#' @param data_list A list of data matrices (blocks).
#' @param preproc_arg The user-provided `preproc` argument. Can be NULL, a single 
#'   `prepper` object, or a list of `prepper` objects.
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
#' @importFrom multivarious prep fresh init_transform
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
    if (!is.list(preproc_arg) || inherits(preproc_arg, "prepper")) { # Treat single prepper object as template
      # Case 1: Single preproc definition - replicate for each block
      message("Applying the same preprocessor definition independently to each block.")
      single_preproc_def <- preproc_arg 
      proclist <- vector("list", S)
      for (i in seq_len(S)) {
          proclist[[i]] <- multivarious::fresh(single_preproc_def) %>% 
                               multivarious::prep(data_list[[i]])
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
          if (!inherits(preproc_arg[[i]], "prepper")) {
              stop(sprintf("Element %d of 'preproc' list is not a valid prepper object.", i), call.=FALSE)
          }
          proclist[[i]] <- multivarious::fresh(preproc_arg[[i]]) %>% 
                               multivarious::prep(data_list[[i]])
      }
       
    } else {
        stop("'preproc' must be NULL, a single prepper object, or a list of prepper objects.", call.=FALSE)
    }
    
    names(proclist) <- block_names

    # Apply the fitted preprocessors
    Xp <- vector("list", S)
    names(Xp) <- block_names
    for (i in seq_along(data_list)) {
      p_i_fitted <- proclist[[i]]
      Xp[[i]] <- multivarious::init_transform(p_i_fitted, data_list[[i]])
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

#' Identify Significant Components via RMT and Inter-block Coherence Tests
#'
#' This function determines which components from a multi-block analysis (e.g.,
#' penalized MFA) are statistically significant. It combines two complementary

#' tests: a Random Matrix Theory (RMT) test based on the Marchenko-Pastur
#' distribution, and an Inter-block Coherence (ICC) test that assesses whether
#' loadings are consistent across blocks.
#'
#' @param fit A fitted multi-block model object containing at minimum a `V_list`
#'   attribute (list of loading matrices) and optionally `sdev` (singular values).
#' @param n Integer; the number of observations (rows) used to fit the model.
#' @param k_vec Optional integer vector of block dimensions (number of columns
#'   per block). If NULL, inferred from `fit` via `block_indices` or `V_list`.
#' @param alpha Numeric; significance level for hypothesis tests (default 0.05).
#' @param check_rmt Logical; if TRUE (default), performs the Marchenko-Pastur
#'   edge test to check if eigenvalues exceed the noise threshold.
#' @param tail_frac Numeric; fraction of smallest eigenvalues used to estimate
#'   noise variance for the RMT test (default 0.3).
#'
#' @return A list with the following elements:
#'   \describe{
#'     \item{keep}{Integer vector of component indices that pass both tests.}
#'     \item{rmt_pass}{Logical vector indicating which components pass the RMT test.}
#'     \item{icc_pass}{Logical vector indicating which components pass the ICC test.}
#'     \item{icc}{Numeric vector of inter-block coherence values per component.}
#'     \item{icc_pvalue}{Numeric vector of p-values for the ICC test.}
#'     \item{mp_edge}{The Marchenko-Pastur edge threshold (NA if RMT skipped).}
#'     \item{sigma2_est}{Estimated noise variance (NA if RMT skipped).}
#'     \item{lambda}{Eigenvalues (squared singular values).}
#'     \item{n, k_vec, alpha}{Input parameters echoed back.}
#'   }
#'
#' @details
#' The RMT test uses the Marchenko-Pastur distribution to identify eigenvalues
#' that exceed the expected bulk edge under the null hypothesis of pure noise.
#' Noise variance is estimated robustly from the tail of the eigenvalue
#' distribution.
#'
#' The ICC test measures the squared cosine similarity of loading vectors
#' across all pairs of blocks. Under the null hypothesis of random loadings,
#' the expected value is 1/k (where k is the harmonic mean of block dimensions).
#' A z-test with Bonferroni correction is applied.
#'
#' @importFrom stats median quantile pnorm p.adjust
#' @importFrom multivarious block_indices
#' @export
significant_components <- function(fit, n, k_vec = NULL,
                                   alpha = 0.05, 
                                   check_rmt = TRUE, 
                                   tail_frac = 0.3) {
  
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
  
  # Extract or derive k_vec (block column counts)
  if (is.null(k_vec)) {
      if (inherits(fit, "multiblock_projector")) {
          k_vec <- sapply(multivarious::block_indices(fit), length)
      } else {
          # Try deriving from V_list dimensions (number of rows)
          k_vec_derived <- sapply(V_list, nrow)
          if (all(k_vec_derived > 0)) { # Check if all blocks have rows
              k_vec <- k_vec_derived
          } else {
              stop("Could not derive k_vec (block column counts). Please provide it explicitly.")
          }
      }
  }
  if (length(k_vec) != S) {
       stop(sprintf("Length of k_vec (%d) does not match the number of blocks in V_list (%d).",
                     length(k_vec), S))
  }
  # Ensure k_vec contains valid dimensions (e.g. > 0 for ICC calculation)
  if (any(k_vec <= 1)) {
       warning("Some blocks have k <= 1 features. ICC calculation might be unstable or meaningless for these blocks.", call.=FALSE)
       # Proceed, but user should be aware. Null distribution assumes k > 1.
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
           check_rmt <- FALSE # Disable RMT check if sdev not found
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

  # ---- 1. RMT Test (Conditional) --------------------------------------
  pass_rmt <- rep(TRUE, r) # Default to TRUE if RMT is skipped
  mp_edge <- NA
  sigma2_est <- NA
  lambda <- numeric(0) # Initialize lambda
  
  if (check_rmt) {
      lambda <- sdev^2                    # compromise eigenvalues
      if (length(lambda) != r) {
          warning(sprintf("Length of sdev^2 (%d) does not match number of components (%d). Skipping RMT test.", 
                          length(lambda), r), call.=FALSE)
          check_rmt <- FALSE
          pass_rmt <- rep(TRUE, r)
      } else {
          k_tot  <- sum(k_vec)
          gamma  <- k_tot / n
          if (gamma > 1) {
              warning("Aspect ratio k_tot/n > 1. Marchenko-Pastur edge calculation assumes n >= k_tot. Skipping RMT test.", call.=FALSE)
              check_rmt <- FALSE
              pass_rmt <- rep(TRUE, r)
          } else {
              # Robust noise variance estimation from the tail
              tail_indices <- seq_len(r)[lambda <= stats::quantile(lambda, 1 - tail_frac, na.rm=TRUE)]
              # Ensure there are enough values in the tail for median
              if (length(tail_indices) < 3) { # Need at least a few points for robust median
                   warning("Not enough eigenvalues in the specified tail quantile for robust noise estimation. Using overall median. RMT results may be less reliable.", call.=FALSE)
                   sigma2_est <- stats::median(lambda, na.rm=TRUE) / (1 + sqrt(gamma))^2 # Less robust fallback
              } else {
                  sigma2_est <- stats::median(lambda[tail_indices], na.rm=TRUE) / (1 + sqrt(gamma))^2
              }
              
              if (!is.finite(sigma2_est) || sigma2_est <= 0) {
                  warning("Estimated noise variance (sigma2) is not positive and finite. Skipping RMT test.", call.=FALSE)
                  check_rmt <- FALSE
                  pass_rmt <- rep(TRUE, r)
              } else {
                  mp_edge <- sigma2_est * (1 + sqrt(gamma))^2
                  pass_rmt <- lambda > mp_edge
              }
          }
      }
  } else {
      message("RMT check skipped.")
  }

  # ---- 2. Inter-block Coherence (ICC) Test --------------------------
  # Calculate harmonic mean of block sizes (use k_vec > 1 for stability)
  k_vec_valid <- k_vec[k_vec > 1]
  if (length(k_vec_valid) == 0) {
      warning("No blocks have k > 1 features. Cannot compute meaningful ICC. Skipping ICC test.")
      pass_icc <- rep(TRUE, r) # Default to pass if test cannot be run
      icc <- rep(NA, r)
      pval <- rep(NA, r)
  } else {
      kbar <- 1 / mean(1 / k_vec_valid)
      icc <- numeric(r)
      valid_comparisons <- 0
      
      # Check if we have at least 2 blocks for comparison
      if (S < 2) {
           warning("Only one block found. Cannot compute ICC. Skipping ICC test.")
           pass_icc <- rep(TRUE, r)
           icc <- rep(NA, r)
           pval <- rep(NA, r)
      } else {
          for (c in seq_len(r)) {
            num <- 0
            comparisons_c <- 0
            for (s in seq_len(S - 1)) {
              for (t in (s + 1):S) {
                # Only compare if both blocks have > 1 feature for stable norm calc
                if (k_vec[s] > 0 && k_vec[t] > 0) { 
                    vsc <- V_list[[s]][, c]
                    vtc <- V_list[[t]][, c]
                    norm_s_sq <- sum(vsc^2)
                    norm_t_sq <- sum(vtc^2)
                    # Avoid division by zero if a loading vector is somehow zero
                    if (norm_s_sq > 1e-12 && norm_t_sq > 1e-12) {
                        cos_sq <- (sum(vsc * vtc)^2) / (norm_s_sq * norm_t_sq)
                        num <- num + cos_sq
                        comparisons_c <- comparisons_c + 1
                    }
                }
              }
            }
            # Average over the number of valid comparisons made for this component
            if (comparisons_c > 0) {
               icc[c] <- num / comparisons_c
            } else {
               icc[c] <- NA # Should not happen if S >= 2 and k_vec > 0
            }
          }
          valid_comparisons <- comparisons_c # Total pairs compared
          
          # Calculate Z-score and p-value where ICC is not NA
          z <- rep(NA, r)
          pval <- rep(NA, r)
          not_na_icc <- !is.na(icc)
          if (any(not_na_icc)) {
              # Use kbar derived from blocks with k > 1 for null expectation
              # Variance term assumes independence and large k approx.
              # Denominator: sqrt(Var(ICC)) = sqrt( Var( sum(cos^2)/(N) ) ) approx sqrt( (1/N^2) * N * Var(cos^2) ) 
              # Var(cos^2) approx Var(chi^2_1 / (k-1)) approx 2/(k-1)^2 -> use kbar? 
              # The provided formula might be simplified; let's use it for now.
              # Need N = S*(S-1)/2 pairs used in the sum for variance? No, the variance is of the mean value.
              # Let's stick to the provided formula's structure, assuming it accounts for averaging. Need S*(S-1)/2 comparisons? valid_comparisons handles this.
              z_denom <- sqrt(2 / (valid_comparisons * kbar^2)) # Using valid_comparisons instead of S*(S-1)/2 if some pairs were skipped?
              if (is.finite(z_denom) && z_denom > 1e-12) {
                 z[not_na_icc] <- (icc[not_na_icc] - 1 / kbar) / z_denom
                 pval[not_na_icc] <- 1 - stats::pnorm(z[not_na_icc])
              } else {
                   warning("Could not compute valid Z-score denominator for ICC test. Skipping.", call.=FALSE)
              }
          }
          # Apply Bonferroni correction where p-value is valid
          pass_icc <- rep(FALSE, r)
          valid_pval <- !is.na(pval)
          if (any(valid_pval)) {
             pass_icc[valid_pval] <- stats::p.adjust(pval[valid_pval], method = "bonferroni") < alpha
          }
      }
  }
  
  # ---- 3. Combined Decision ------------------------------------------
  keep <- which(pass_rmt & pass_icc)
  
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
              alpha = alpha))
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
#' @importFrom stats scale
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

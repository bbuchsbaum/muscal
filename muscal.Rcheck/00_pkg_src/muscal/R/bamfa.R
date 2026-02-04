#' @useDynLib muscal, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @import RcppArmadillo
#' @importFrom stats rnorm sd
#' @importFrom chk chk_list chk_integer chk_numeric chk_gte chk_flag chk_matrix chk_not_empty
#' @importFrom chk chk_true
#' @importFrom Matrix tcrossprod Diagonal sparseMatrix rankMatrix
#' @importFrom MASS ginv
#' @importFrom crayon bold magenta green blue yellow
#' @importFrom purrr map map2
#' @importFrom multivarious fit fit_transform transform inverse_transform fresh center concat_pre_processors pass
#' @importFrom multidesign multiblock
#' @importFrom RSpectra svds
#' @importFrom rlang enquo quo_is_missing quo_get_expr
#' @importFrom dplyr select pull
#' @importFrom multivarious multiblock_projector
#' @importFrom RSpectra svds
#' @importFrom stats rnorm
#' @import Rcpp
#' @importFrom rlang .data
NULL

#' Internal ridge least-squares fit
#'
#' Solves min ||X - ZB||_F^2 + lambda*||B||_F^2 for B.
#' Regularizes the unpenalized case with a tiny ridge penalty for stability.
#'
#' @param Z Predictor matrix (n x k)
#' @param X Response matrix (n x p)
#' @param lambda Ridge penalty (scalar, non-negative)
#' @return Coefficient matrix B (k x p)
#' @noRd
#' @keywords internal
ls_ridge <- function(Z, X, lambda = 0) {
  k <- ncol(Z)
  if (k == 0) { # Handle case with zero predictors
    return(matrix(0.0, 0, ncol(X)))
  }

  # Add tiny ridge for stability in unpenalized case
  effective_lambda <- if (lambda == 0) 1e-8 else lambda
  
  # Call the fast Rcpp implementation
  coefs <- tryCatch(
    ridge_solve(Z, X, effective_lambda),
    error = function(e) {
      warning("Rcpp ridge_solve failed: ", e$message, ". Falling back to R.", call. = FALSE)
      # Fallback R implementation
      ZtZ <- crossprod(Z)
      ZtX <- crossprod(Z, X)
      I_k <- diag(k)
      solve(ZtZ + effective_lambda * I_k, ZtX)
    }
  )

  if (!is.matrix(coefs)) {
      coefs <- matrix(coefs, nrow = k, ncol = ncol(X))
  }
  return(coefs) # Returns B (k x p)
}


#' @keywords internal
prep_bamfa_blocks <- function(data, preproc) {
  proclist_fitted <- lapply(seq_along(data), function(i) {
    multivarious::fresh(preproc) %>% multivarious::fit(data[[i]])
  })

  X_p <- lapply(seq_along(data), function(i) {
    multivarious::transform(proclist_fitted[[i]], data[[i]])
  })
  
  p_vec <- vapply(X_p, ncol, integer(1))
  
  block_indices <- list()
  current_start <- 1
  for (i in seq_along(p_vec)) {
    block_indices[[i]] <- current_start:(current_start + p_vec[i] - 1)
    current_start <- current_start + p_vec[i]
  }
  names(block_indices) <- names(data)
  
  list(
    X_p = X_p,
    proclist_fitted = proclist_fitted,
    block_indices = block_indices,
    p_vec = p_vec
  )
}

#' Barycentric Multiple Factor Analysis (BaMFA)
#'
#' @md
#' @description
#' Performs Barycentric Multiple Factor Analysis (BaMFA) using an alternating
#' optimization approach based on a two-level factor model.
#' It decomposes multi-block data (e.g., multiple subjects) into a shared global
#' subspace (`G`) and block-specific subspaces (`B_i`).
#'
#' @details
#' The algorithm models each data block \\(X_i\\) as:
#' \\\[ X_i = S_i G^T + U_i B_i^T + E_i \\\]
#' where \\(G\\) represents shared global loadings, \\(B_i\\) represents block-specific
#' local loadings, \\(S_i\\) and \\(U_i\\) are the corresponding scores, and \\(E_i\\) is noise.
#' Loadings are constrained to be orthonormal (\\(G^T G = I, B_i^T B_i = I\\)).
#'
#' The algorithm aims to minimize the total reconstruction error:
#' \\\[ \\sum_{i=1}^{m} ||X_i - S_i G^T - U_i B_i^T ||_F^2 \\\]
#' using an iterative alternating optimization strategy (similar to Expectation-Maximization):
#' 1. **Preprocessing:** Each block is preprocessed using the provided `preproc` pipeline.
#' 2. **Initialization:** Initialize global loadings `G` (via SVD on the mean block).
#' 3. **Iterate (E-step like):** For each block `i`, holding `G` fixed, update scores `S_i` (global projection), calculate residuals, find the local basis `B_i` from residuals (via SVD), and estimate local scores `U_i` (via projection, potentially regularized by `lambda_l`).
#' 4. **Iterate (M-step like):** Holding scores `S_i` fixed, update the global loadings `G` by performing SVD on an aggregated cross-product matrix `sum(X_i^T S_i)`.
#' 5. **Repeat steps 3-4** for `niter` iterations or until convergence based on `tol`.
#'
#' @param data A `list` of matrices (each a block, e.g., subject), a `multiblock` object,
#'   or a `multidesign` object. If a list or multiblock, all blocks must have the same
#'   number of rows (observations, e.g., conditions). If a multidesign, the `subject` parameter
#'   must be specified to indicate how to split the data.
#' @param k_g Integer: Number of global components to extract (shared across blocks).
#' @param k_l Integer: Number of local components to extract (specific to each block).
#' @param niter Integer: Maximum number of iterations (default: 10).
#' @param preproc A preprocessing pipeline object from `multivarious` (default: `multivarious::center()`).
#' @param lambda_l Numeric: Non-negative ridge penalty applied when estimating the local scores `U_i`
#'   using the internal `ls_ridge` function. Default: 0 (no penalty).
#' @param tol Numeric: Convergence tolerance. Stops if the relative change in the objective function
#'   (Mean Squared Error) and/or the global loadings `G` is less than `tol`. Default: 1e-5.
#' @param subject Optional: A variable name identifying the blocking/subject variable when using a
#'   multidesign object. Required only for the multidesign method.
#' @param ... Additional arguments (currently unused).
#'
#' @return A `multiblock_projector` object with class `"bamfa"`.
#'   The base projector stores the global loading matrix in `v`, the
#'   concatenated preprocessor in `preproc`, and block mappings in
#'   `block_indices`. Additional named elements passed via `...` include:
#'   * `B_list` -- block-specific local loading matrices.
#'   * `S_list` -- block-specific global score matrices.
#'   * `U_list` -- block-specific local score matrices.
#'   * `k_g` -- number of global components in the final model.
#'   * `k_l` -- requested number of local components.
#'   * `lambda_l` -- regularization parameter used for local scores `U_i`.
#'   * `niter_actual` -- actual number of iterations performed.
#'   * `objective_trace` -- objective function value at each iteration.
#'   * `g_change_trace` -- relative change in global loadings at each iteration.
#'   * `data_names` -- names of the input blocks/subjects.
#'   * `proclist_fitted` -- list of fitted preprocessors for each block.
#'
#' @section Caveats and Limitations:
#' * **Model Choice:** Assumes a linear factor model with orthogonal global and local components.
#' * **Initialization:** Results can be sensitive to the initialization of `G`.
#' * **Local Minima:** The alternating optimization algorithm may converge to a local minimum.
#' * **Interpretation:** The separation into global and local components depends on the model fit and ranks (`k_g`, `k_l`).
#' * **Regularization (`lambda_l`):** Penalizes the squared Frobenius norm of the local scores `U_i` via ridge regression.
#' * **Inference:** The method does not provide p-values or confidence intervals.
#'
#' @examples
#' # Generate example multi-block data (e.g., 3 subjects, 10 conditions, 50 features)
#' set.seed(123)
#' n_obs <- 10
#' n_features <- 50
#' n_subjects <- 3
#' data_list <- lapply(1:n_subjects, function(i) {
#'   matrix(rnorm(n_obs * n_features), n_obs, n_features) +
#'   matrix(rnorm(n_obs * 1, mean=i), n_obs, n_features) # Add subject offset
#' })
#' names(data_list) <- paste0("Subject_", 1:n_subjects)
#'
#' # Run BaMFA with k_g=3 global, k_l=2 local components
#' result <- bamfa(data_list, k_g = 3, k_l = 2, niter = 10)
#' print(result)
#'
#' @seealso \code{\link[multivarious]{pca}}, \code{\link[multidesign]{multiblock}}
#' @export
#' @rdname bamfa
bamfa <- function(data, k_g = 2, k_l = 2, niter = 10,
                  preproc = multivarious::center(), lambda_l = 0, tol = 1e-5, 
                  subject = NULL, ...) {
  UseMethod("bamfa")
}

#' @md
#' @rdname bamfa
#' @export
bamfa.default <- function(data, k_g = 2, k_l = 2, niter = 10,
                          preproc = multivarious::center(), lambda_l = 0, tol = 1e-5, ...) {
  # -----------------------------------------------------------------------
  # 0. Basic checks and setup
  # -----------------------------------------------------------------------
  # Input validation
  chk::chk_not_empty(data)
  chk::chk_list(data)
  k_g <- as.integer(k_g)
  k_l <- as.integer(k_l)
  niter <- as.integer(niter)
  chk::chk_integer(k_g)
  chk::chk_integer(k_l)
  chk::chk_integer(niter)
  chk::chk_numeric(lambda_l)
  chk::chk_gte(lambda_l, 0)
  chk::chk_numeric(tol)
  chk::chk_gte(tol, 0)
  
  # Check for consistent number of rows
  row_counts <- vapply(data, nrow, integer(1))
  if (length(unique(row_counts)) > 1) {
    stop("All blocks must have the same number of rows (observations).")
  }
  n_obs <- row_counts[1]
  
  # Preprocessing
  prep_res <- prep_bamfa_blocks(data, preproc)
  X_p <- prep_res$X_p
  proclist_fitted <- prep_res$proclist_fitted
  block_indices <- prep_res$block_indices
  p_vec <- prep_res$p_vec
  p_tot <- sum(p_vec)
  
  # Check if k_g is too large
  k_g <- min(k_g, n_obs, p_tot)
  
  # Check if k_l is too large for any block
  k_l <- min(k_l, n_obs, min(p_vec))
  
  # -----------------------------------------------------------------------
  # 1. Initialization
  # -----------------------------------------------------------------------
  # Initialize G from a random orthogonal matrix
  G <- matrix(rnorm(p_tot * k_g), p_tot, k_g)
  G <- qr.Q(qr(G))
  
  # Lists to store block-specific components
  B_list <- vector("list", length(data)) # Local loadings
  S_list <- vector("list", length(data)) # Global scores
  U_list <- vector("list", length(data)) # Local scores
  
  obj_trace <- numeric(niter) # To track convergence
  g_change_trace <- numeric(niter) # To track G convergence
  
  # -----------------------------------------------------------------------
  # 2. Iterative Optimization
  # -----------------------------------------------------------------------
  for (iter in 1:niter) {
    
    # E-step: Update scores (S_i, U_i) and local bases (B_i) for each block
    obj_sum <- 0
    G_old <- G
    
    for (i in seq_along(X_p)) {
      G_i <- G[block_indices[[i]], , drop = FALSE]
      S_list[[i]] <- X_p[[i]] %*% G_i # n_obs x k_g
      R_i <- X_p[[i]] - S_list[[i]] %*% t(G_i)
      
      # Use RSpectra for faster SVD when k_l is smaller than both dims, with fallback
      if (k_l > 0 && k_l < min(dim(R_i))) {
        svd_R <- tryCatch(
          RSpectra::svds(R_i, k = k_l, nu = k_l, nv = k_l),
          error = function(e) svd(R_i, nu = k_l, nv = k_l)
        )
      } else {
        svd_R <- svd(R_i, nu = k_l, nv = k_l)
      }
      B_list[[i]] <- svd_R$v # p_i x k_l
      
      # Estimate local scores with optional ridge shrinkage
      if (k_l > 0) {
        U_proj <- R_i %*% B_list[[i]]
        if (lambda_l > 0) {
          U_list[[i]] <- U_proj / (1 + lambda_l)  # B_i is orthonormal
        } else {
          U_list[[i]] <- U_proj
        }
      } else {
        U_list[[i]] <- matrix(0, n_obs, 0)
      }
      
      recon_i <- S_list[[i]] %*% t(G_i) + U_list[[i]] %*% t(B_list[[i]])
      obj_sum <- obj_sum + sum((X_p[[i]] - recon_i)^2)
    }
    
    obj_trace[iter] <- obj_sum / (n_obs * p_tot) # Mean Squared Error
    
    # M-step: Update global basis G
    C <- matrix(0, p_tot, k_g)
    for (i in seq_along(X_p)) {
      C[block_indices[[i]], ] <- crossprod(X_p[[i]], S_list[[i]])
    }
    
    svd_C <- svd(C, nu = k_g)
    G <- svd_C$u # Update G with orthogonal vectors
    
    # Convergence check
    denom_g <- max(norm(G_old, "F"), 1e-12)
    g_change_trace[iter] <- norm(G - G_old, "F") / denom_g
    if (iter > 1) {
      rel_change_obj <- abs(obj_trace[iter] - obj_trace[iter - 1]) / obj_trace[iter - 1]
      if (rel_change_obj < tol && g_change_trace[iter] < tol) {
        obj_trace <- obj_trace[1:iter]
        g_change_trace <- g_change_trace[1:iter]
        break
      }
    }
  }
  
  # -----------------------------------------------------------------------
  # 3. Construct result object
  # -----------------------------------------------------------------------
  preproc_concat <- multivarious::concat_pre_processors(proclist_fitted, block_indices)
  
  multivarious::multiblock_projector(
    v = G,
    preproc = preproc_concat,
    B_list = B_list,
    S_list = S_list,
    U_list = U_list,
    k_g = k_g,
    k_l = k_l,
    lambda_l = lambda_l,
    niter_actual = iter,
    objective_trace = obj_trace,
    g_change_trace = g_change_trace,
    data_names = names(data),
    proclist_fitted = proclist_fitted,
    block_indices = block_indices,
    classes = "bamfa"
  )
}


#' @md
#' @rdname bamfa
#' @export
bamfa.list <- function(data, k_g = 2, k_l = 2, niter = 10,
                       preproc = multivarious::center(), lambda_l = 0, tol = 1e-5, ...) {
  if (is.null(names(data))) {
    names(data) <- paste0("Block_", seq_along(data))
  } else if (any(names(data) == "")) {
     names(data)[names(data) == ""] <- paste0("Block_", which(names(data) == ""))
  }
  
  bamfa.default(
    data = data, k_g = k_g, k_l = k_l, niter = niter,
    preproc = preproc, lambda_l = lambda_l, tol = tol, ...
  )
}

#' @md
#' @rdname bamfa
#' @export
bamfa.multiblock <- function(data, k_g = 2, k_l = 2, niter = 10,
                             preproc = multivarious::center(), lambda_l = 0, tol = 1e-5, ...) {
  
  data_list <- as.list(data)
  
  bamfa.list(
    data = data_list, k_g = k_g, k_l = k_l, niter = niter,
    preproc = preproc, lambda_l = lambda_l, tol = tol, ...
  )
}

#' @md
#' @rdname bamfa
#' @importFrom rlang enquo quo_is_missing quo_get_expr
#' @importFrom dplyr select pull
#' @export
bamfa.multidesign <- function(data, k_g = 2, k_l = 2, niter = 10,
                              preproc = multivarious::center(), lambda_l = 0, 
                              tol = 1e-5, subject, ...) {
  subject_quo <- rlang::enquo(subject)
  
  if (rlang::quo_is_missing(subject_quo)) {
    stop("The 'subject' parameter is required for bamfa.multidesign(). Please specify the subject variable to split the data by.")
  }
  
  subjects <- factor(data$design %>% dplyr::pull(!!subject_quo))
  sdat <- split(data, subjects, drop = TRUE)
  
  # Since multidesign objects handle preprocessing internally, we pass preprocessed data
  # and a 'pass-through' preprocessor to the next method.
  
  result <- bamfa.list(
    data     = sdat,
    k_g      = k_g,
    k_l      = k_l,
    niter    = niter,
    preproc  = preproc,
    lambda_l = lambda_l,
    tol      = tol,
    ...
  )
  
  result$subject_variable <- as.character(rlang::quo_get_expr(subject_quo))
  result$subject_levels <- levels(subjects)
  
  return(result)
}


#' Predict method for BaMFA objects
#' 
#' Reconstructs data blocks from a fitted BaMFA model.
#' 
#' @param object A fitted `bamfa` object.
#' @param new_data Optional new data to reconstruct. If NULL (default), reconstructs the training data.
#' @param ... Additional arguments (unused).
#' @return A list of reconstructed data matrices.
#' @export
predict.bamfa <- function(object, new_data = NULL, ...) {
  
  if (is.null(new_data)) {
    # Reconstruct training data
    G <- object$v
    recon_list <- lapply(seq_along(object$data_names), function(i) {
      G_i <- G[object$block_indices[[i]], , drop = FALSE]
      S_i <- object$S_list[[i]]
      U_i <- object$U_list[[i]]
      B_i <- object$B_list[[i]]
      
      recon <- S_i %*% t(G_i) + U_i %*% t(B_i)
      
      # Inverse transform to original data space
      multivarious::inverse_transform(object$proclist_fitted[[i]], recon)
    })
    names(recon_list) <- object$data_names
    return(recon_list)
  }
  
  # Reconstruct new data
  if (!is.list(new_data)) {
    new_data <- list(new_data)
  }
  
  # Preprocess new data using fitted preprocessors
  X_p <- lapply(seq_along(new_data), function(i) {
    multivarious::transform(object$proclist_fitted[[i]], new_data[[i]])
  })
  
  G <- object$v
  k_l <- object$k_l
  lambda_l <- object$lambda_l
  
  recon_list <- lapply(seq_along(X_p), function(i) {
    G_i <- G[object$block_indices[[i]], , drop = FALSE]
    B_i <- object$B_list[[i]]
    if (!is.null(B_i) && ncol(B_i) != k_l) {
      stop("Stored local loadings B_list[[", i, "]] do not match k_l.")
    }
    S_new <- X_p[[i]] %*% G_i
    R_new <- X_p[[i]] - S_new %*% t(G_i)
    
    if (k_l > 0) {
      if (is.null(B_i)) {
        stop("Stored local loadings missing; cannot project new data without refitting.")
      }
      U_proj <- R_new %*% B_i
      if (lambda_l > 0) {
        U_new <- U_proj / (1 + lambda_l)
      } else {
        U_new <- U_proj
      }
    } else {
      B_i <- matrix(0, ncol(R_new), 0)
      U_new <- matrix(0, nrow(R_new), 0)
    }
    
    recon <- S_new %*% t(G_i) + U_new %*% t(B_i)
    
    # Inverse transform to original data space
    multivarious::inverse_transform(object$proclist_fitted[[i]], recon)
  })
  names(recon_list) <- names(new_data) %||% paste0("Block_", seq_along(new_data))
  
  recon_list
}


#' Print Method for BaMFA Objects
#'
#' @md
#' @description
#' Prints a concise summary of a `bamfa` object fitted with the EM-like algorithm,
#' highlighting the model structure and fit.
#'
#' @param x An object of class `bamfa`.
#' @param ... Additional parameters (unused).
#'
#' @return Invisibly returns the input object.
#'
#' @export
print.bamfa <- function(x, ...) {
  header <- crayon::bold(crayon::blue("Barycentric Multiple Factor Analysis (BaMFA)"))
  cat(header, "\n\n")

  k_g_val <- x$k_g
  k_l_val <- x$k_l
  niter_actual_val <- x$niter_actual
  lambda_l_val <- x$lambda_l
  data_names_val <- x$data_names
  objective_trace_val <- x$objective_trace
  n_blocks <- length(x$block_indices)
  
  # Model structure
  cat(crayon::green("Model Structure:"), "\n")
  cat("  Global components (k_g):", crayon::bold(k_g_val), "\n")
  cat("  Local components (k_l requested):", crayon::bold(k_l_val), "\n")
  cat("  Convergence after", crayon::bold(niter_actual_val), "iterations\n")
  cat("  Local score regularization (lambda_l):", crayon::bold(lambda_l_val), "\n\n")
  
  # Block information
  cat(crayon::green("Block Information:"), "\n")
  for (i in seq_len(n_blocks)) {
    cat("  Block", crayon::bold(i), "(", crayon::blue(data_names_val[i]), "):", "\n")
    num_features <- length(x$block_indices[[i]])
    num_local_comp <- if (!is.null(x$B_list[[i]])) ncol(x$B_list[[i]]) else NA
    cat("    Features:", crayon::yellow(num_features), "\n")
    cat("    Local components retained:", crayon::yellow(num_local_comp), "\n")
  }

  # Final objective
  final_obj_idx <- length(objective_trace_val)
  cat("\nFinal reconstruction error (per feature):", 
      crayon::bold(format(objective_trace_val[final_obj_idx], digits=6)), "\n")
  
  invisible(x)
}

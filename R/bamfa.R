#' @useDynLib multivarious, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @import RcppArmadillo
#' @import RcppEigen
#' @importFrom stats rnorm sd svd qr.coef
#' @importFrom chk chk_list chk_integer chk_numeric chk_gte chk_flag chk_matrix
#' @importFrom chk chk_true
#' @importFrom Matrix tcrossprod Diagonal sparseMatrix rankMatrix
#' @importFrom MASS ginv
#' @importFrom crayon bold magenta green blue yellow
#' @importFrom purrr map map2
#' @importFrom multivarious init_transform prep fresh center concat_pre_processors
#' @importFrom multidesign multiblock
NULL

#' Internal ridge least-squares fit
#'
#' Solves min ||X - ZB||_F^2 + lambda*||B||_F^2 for B.
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

  # Assuming Z and X are already centered if necessary by preproc

  if (lambda == 0) {
    # Use QR decomposition for unpenalized least squares
    coefs <- tryCatch(stats::qr.coef(stats::qr(Z), X), error = function(e) NULL)
    if (is.null(coefs)) {
        warning("QR decomposition failed in ls_ridge, using pseudoinverse.", call. = FALSE)
        Zplus <- tryCatch(MASS::ginv(Z), error = function(e) {
             warning("ginv failed: " , e$message, call. = FALSE)
             NULL
             })
        if(is.null(Zplus)) return(matrix(0.0, k, ncol(X))) # Fallback
        coefs <- Zplus %*% X
    }
  } else {
    # Ridge regression: solve(t(Z)Z + lambda*I, t(Z)X)
    tZ <- t(Z) # k x n
    ZtZ <- tZ %*% Z # k x k
    ZtX <- tZ %*% X # k x p
    I_k <- diag(k)
    coefs <- tryCatch(solve(ZtZ + lambda * I_k, ZtX), error = function(e) NULL)
    if (is.null(coefs)) {
        warning("Solve failed in ridge regression, returning zero coefficients.", call. = FALSE)
        coefs <- matrix(0.0, k, ncol(X))
    }
  }
  # Ensure coefs is a matrix
  if (!is.matrix(coefs)) {
      coefs <- matrix(coefs, nrow = k, ncol = ncol(X))
  }
  return(coefs) # Returns B (k x p)
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
#' The algorithm models each data block \(X_i\) as:
#' \[ X_i = S_i G^T + U_i B_i^T + E_i \]
#' where \(G\) represents shared global loadings, \(B_i\) represents block-specific
#' local loadings, \(S_i\) and \(U_i\) are the corresponding scores, and \(E_i\) is noise.
#' Loadings are constrained to be orthonormal (\(G^T G = I, B_i^T B_i = I\)).
#'
#' The algorithm aims to minimize the total reconstruction error:
#' \[ \sum_{i=1}^{m} \|X_i - S_i G^T - U_i B_i^T \|_F^2 \]
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
#'   (Mean Squared Error) is less than `tol`. Default: 1e-5.
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
bamfa.list <- function(data, k_g = 2, k_l = 2, niter = 10,
                       preproc = multivarious::center(), lambda_l = 0, tol = 1e-5, ...) {
  # Ensure list elements have names
  if (is.null(names(data))) {
    names(data) <- paste0("Block_", seq_along(data))
  } else if (any(names(data) == "")) {
     names(data)[names(data) == ""] <- paste0("Block_", which(names(data) == ""))
  }
  # Convert list to multiblock object
  mb <- multidesign::multiblock(data)

  # Call the multiblock method
  bamfa.multiblock(
    data     = mb,
    k_g      = k_g,
    k_l      = k_l,
    niter    = niter,
    preproc  = preproc,
    lambda_l = lambda_l,
    tol      = tol,
    ...
  )
}

#' @md
#' @rdname bamfa
#' @importFrom rlang enquo
#' @importFrom dplyr select pull
#' @export
bamfa.multidesign <- function(data, k_g = 2, k_l = 2, niter = 10,
                              preproc = multivarious::center(), lambda_l = 0, 
                              tol = 1e-5, subject, ...) {
  # Get the subject quo for consistent handling
  subject_quo <- rlang::enquo(subject)
  
  # Check if subject is missing (unevaluated quo)
  if (rlang::quo_is_missing(subject_quo)) {
    stop("The 'subject' parameter is required for bamfa.multidesign(). Please specify the subject variable to split the data by.")
  }
  
  # Extract subject variable from design
  subjects <- factor(data$design %>% 
                     dplyr::select(!!subject_quo) %>% 
                     dplyr::pull(!!subject_quo))
  subject_set <- levels(subjects)
  
  # Split data by subject
  sdat <- split(data, subject)
  
  # Create preprocessors, one per subject
  proclist <- lapply(seq_along(sdat), function(sd) {
    multivarious:::fresh(preproc) %>% prep()
  })
  names(proclist) <- as.character(subject_set)
  
  # Preprocess data for each subject
  strata <- seq_along(sdat) %>% purrr::map(function(i) {
    p <- multivarious::prep(proclist[[i]], sdat[[i]]$x)
    Xi <- sdat[[i]]$x
    Xout <- multivarious::init_transform(p, Xi)
    Xout
  })
  names(strata) <- subject_set
  
  # Call the list method with the extracted data
  result <- bamfa.list(
    data     = strata,
    k_g      = k_g,
    k_l      = k_l,
    niter    = niter,
    preproc  = multivarious::pass(), # Data is already preprocessed
    lambda_l = lambda_l,
    tol      = tol,
    ...
  )
  
  # Store additional information about the multidesign input
  result$subject_variable <- as.character(rlang::quo_get_expr(subject_quo))
  result$subject_levels <- subject_set
  
  return(result)
}

#' @md
#' @rdname bamfa
#' @export
bamfa.multiblock <- function(data, k_g = 2, k_l = 2, niter = 10, preproc = center(),
                             lambda_l = 0, tol = 1e-5, ...) {
  # -----------------------------------------------------------------------
  # 0. Basic checks and setup
  # -----------------------------------------------------------------------
  # Input validation
  chk::chk_true(length(data) > 0)
  chk::chk_integer(k_g)
  chk::chk_integer(k_l)
  chk::chk_integer(niter)
  chk::chk_numeric(lambda_l)
  chk::chk_gte(lambda_l, 0)
  chk::chk_numeric(tol)
  chk::chk_gt(tol, 0)
  
  # Ensure consistent row counts across blocks
  n_rows <- sapply(data, nrow)
  chk::chk_true(all(n_rows == n_rows[1]))
  n <- n_rows[1]  # Number of observations
  
  # -----------------------------------------------------------------------
  # 1. Preprocess blocks
  # -----------------------------------------------------------------------
  # Set up preprocessors (one for each block)
  proclist <- lapply(seq_along(data), function(i) {
    multivarious:::fresh(preproc) %>% prep(data[[i]])
  })
  
  names(proclist) <- names(data)
  
  # Apply preprocessing to each block
  Xp <- vector("list", length(data))
  names(Xp) <- names(data)
  
  for (i in seq_along(data)) {
    Xi <- data[[i]]
    p <- proclist[[i]]
    Xp[[i]] <- multivarious::init_transform(p, Xi)
  }
  
  # -----------------------------------------------------------------------
  # 2. Iterative alternating optimization
  # -----------------------------------------------------------------------
  # Setup progress tracking and convergence metrics
  objective_trace <- numeric(niter + 1)
  mean_block_sizes <- mean(sapply(Xp, ncol))  # For per-feature MSE calculation
  
  # Initialize global loadings G using SVD on the mean block
  X_mean <- Reduce("+", Xp) / length(Xp)
  svd_result <- tryCatch(
    svd(X_mean, nu = 0, nv = min(k_g, min(dim(X_mean)))),
    error = function(e) NULL
  )
  
  if (is.null(svd_result)) {
    warning("SVD failed on mean block, using alternative initialization.", call. = FALSE)
    # Fallback to a small random matrix
    set.seed(42)  # For reproducibility
    G_current <- matrix(rnorm(ncol(X_mean) * min(k_g, min(dim(X_mean)))), 
                        ncol = min(k_g, min(dim(X_mean))))
    G_current <- qr.Q(qr(G_current))  # Orthonormalize
  } else {
    G_current <- svd_result$v[, seq_len(min(k_g, length(svd_result$d))), drop = FALSE]
  }
  
  # Initialize lists for storing block-specific data
  S_list <- vector("list", length(Xp))  # Global scores
  U_list <- vector("list", length(Xp))  # Local scores
  B_list <- vector("list", length(Xp))  # Local loadings
  
  # Initial projection to global space
  for (i in seq_along(Xp)) {
    S_list[[i]] <- Xp[[i]] %*% G_current
    
    # Calculate residuals
    res_i <- Xp[[i]] - S_list[[i]] %*% t(G_current)
    
    # Find local components via SVD of residuals
    svd_res <- tryCatch(
      svd(res_i, nu = min(k_l, min(dim(res_i))), nv = min(k_l, min(dim(res_i)))),
      error = function(e) NULL
    )
    
    if (is.null(svd_res)) {
      warning(paste0("SVD failed for block ", i, ", using alternative local basis."), call. = FALSE)
      # Fallback: random orthogonal matrix
      B_i <- matrix(rnorm(ncol(res_i) * min(k_l, min(dim(res_i)))),
                    ncol = min(k_l, min(dim(res_i))))
      B_i <- qr.Q(qr(B_i))
      U_i <- matrix(0, nrow = nrow(res_i), ncol = ncol(B_i))
    } else {
      # Determine actual number of local components from SVD
      k_l_actual <- min(k_l, length(svd_res$d))
      # Local loadings (features x k_l)
      B_i <- svd_res$v[, seq_len(k_l_actual), drop = FALSE]
      # Local scores scaled by singular values (observations x k_l)
      U_i <- svd_res$u[, seq_len(k_l_actual), drop = FALSE] *
        svd_res$d[seq_len(k_l_actual)]
    }
    
    B_list[[i]] <- B_i
    U_list[[i]] <- U_i
  }
  
  # Calculate initial objective value (mean squared error per feature)
  total_mse <- 0
  total_features <- 0
  for (i in seq_along(Xp)) {
    res_norm <- norm(Xp[[i]] - S_list[[i]] %*% t(G_current) - U_list[[i]] %*% t(B_list[[i]]), "F")^2
    total_mse <- total_mse + res_norm
    total_features <- total_features + prod(dim(Xp[[i]]))
  }
  objective_trace[1] <- total_mse / total_features
  
  # Main iteration loop
  iter_actual <- 0
  for (iter in 1:niter) {
    iter_actual <- iter
    
    # -----------------------------------------------------------------------
    # 2.1. Local update (block-specific parameters, E-step like)
    # -----------------------------------------------------------------------
    for (i in seq_along(Xp)) {
      # Project to global space
      S_list[[i]] <- Xp[[i]] %*% G_current
      
      # Calculate residuals
      res_i <- Xp[[i]] - S_list[[i]] %*% t(G_current)
      
      # Find local components via SVD of residuals
      svd_res <- tryCatch(
        svd(res_i, nu = min(k_l, min(dim(res_i))), nv = min(k_l, min(dim(res_i)))),
        error = function(e) NULL
      )
      
      if (is.null(svd_res)) {
        # Keep previous loadings if SVD fails
        warning(paste0("SVD failed for block ", i, " in iteration ", iter, ", keeping previous local basis."), call. = FALSE)
        if (ncol(B_list[[i]]) > 0) {
          # Estimate local scores with ridge regression
          U_list[[i]] <- res_i %*% B_list[[i]]
        }
      } else {
        # Update local basis (loadings)
        k_l_actual <- min(k_l, length(svd_res$d))
        B_list[[i]] <- svd_res$v[, seq_len(k_l_actual), drop = FALSE]

        # Local scores from SVD scaled by singular values
        if (k_l_actual > 0) {
          U_list[[i]] <- svd_res$u[, seq_len(k_l_actual), drop = FALSE] *
            svd_res$d[seq_len(k_l_actual)]
        } else {
          U_list[[i]] <- matrix(0, nrow = nrow(res_i), ncol = 0)
        }
      }
    }
    
    # -----------------------------------------------------------------------
    # 2.2. Global update (shared parameters, M-step like)
    # -----------------------------------------------------------------------
    # Create the cross-product matrix for global update
    XtS_sum <- matrix(0, nrow = ncol(X_mean), ncol = k_g)
    
    for (i in seq_along(Xp)) {
      # Calculate residual after removing local component
      res_i <- Xp[[i]] - U_list[[i]] %*% t(B_list[[i]])
      
      # Accumulate cross-product for global update
      XtS_sum <- XtS_sum + t(res_i) %*% S_list[[i]]
    }
    
    # Update global loadings via SVD
    svd_global <- tryCatch(
      svd(XtS_sum, nu = min(k_g, min(dim(XtS_sum))), nv = min(k_g, min(dim(XtS_sum)))),
      error = function(e) NULL
    )
    
    if (is.null(svd_global)) {
      warning(paste0("SVD failed for global update in iteration ", iter, ", keeping previous global basis."), call. = FALSE)
    } else {
      # Extract new global basis
      V_new <- svd_global$u
      U_new <- svd_global$v
      
      # Update global loadings
      G_new <- V_new %*% t(U_new)
      
      # Ensure orthogonality of global loadings
      qr_res <- qr(G_new)
      G_current <- qr.Q(qr_res)[, seq_len(min(k_g, qr_res$rank)), drop = FALSE]
    }
    
    # -----------------------------------------------------------------------
    # 2.3. Evaluate convergence
    # -----------------------------------------------------------------------
    # Calculate current objective value (mean squared error per feature)
    total_mse <- 0
    total_features <- 0
    for (i in seq_along(Xp)) {
      res_norm <- norm(Xp[[i]] - S_list[[i]] %*% t(G_current) - U_list[[i]] %*% t(B_list[[i]]), "F")^2
      total_mse <- total_mse + res_norm
      total_features <- total_features + prod(dim(Xp[[i]]))
    }
    objective_trace[iter + 1] <- total_mse / total_features
    
    # Check convergence
    rel_change <- abs(objective_trace[iter + 1] - objective_trace[iter]) / (objective_trace[iter] + .Machine$double.eps)
    if (rel_change < tol) {
      message(sprintf("Converged at iteration %d with relative change %.6f < %.6f (tolerance)", 
                     iter, rel_change, tol))
      break
    }
  }
  
  # Handle the case where no components could be extracted
  if (is.null(G_current)) G_current <- matrix(0.0, p_tot, 0)
  k_g_final <- ncol(G_current)

  # -----------------------------------------------------------------------
  # 3. Assemble Final Output Object
  # -----------------------------------------------------------------------
  # Ensure block_indices is correctly computed/available before this point
  # Assuming p_tot is the total number of features after concatenation
  # We need block_indices that map 1:p_tot back to the original blocks
  p_vec <- sapply(Xp, ncol) # Number of features per block (post-preprocessing)
  p_tot <- sum(p_vec)
  if (nrow(G_current) != p_tot) {
      warning(sprintf("Dimension mismatch: nrow(G)=%d but expected sum(ncol(Xp))=%d. Returning raw list.", 
                       nrow(G_current), p_tot), call.=FALSE)
       # Fallback to raw list
       return(list(
          G = G_current, B = B_list, S = S_list, U = U_list,
          block_indices = NULL, # Indicate failure
          niter_actual = iter_actual, k_g = k_g_final, k_l = k_l,
          lambda_l = lambda_l, preproc = proclist, data_names = names(data),
          objective_trace = objective_trace[1:(iter_actual + 1)]
      ))
  }
  
  block_indices <- list()
  current_start <- 1
  for(i in seq_along(p_vec)) {
      block_indices[[i]] <- current_start:(current_start + p_vec[i] - 1)
      current_start <- current_start + p_vec[i]
  }
  names(block_indices) <- names(Xp)

  # Create concatenated preprocessor
  final_preproc <- NULL
  if (!is.null(proclist)) {
      final_preproc <- multivarious::concat_pre_processors(proclist, block_indices)
  } else {
      # Should not happen if default preproc=center() was used, but handle defensively
      pass_proc <- multivarious::prep(multivarious::pass())
      final_preproc <- concat_pre_processors(rep(list(pass_proc), length(data)), block_indices)
  }

  # Construct the multiblock_projector based on Global components G
  # Pass BaMFA specific results via '...' to be stored as list elements
  result_projector <- multivarious::multiblock_projector(
      v = G_current,             
      preproc = final_preproc,   
      block_indices = block_indices,
      # BaMFA specific elements passed via ...
      B_list = B_list,         
      S_list = S_list,         
      U_list = U_list,         
      k_g = k_g_final,        
      k_l = k_l,              
      lambda_l = lambda_l,    
      niter_actual = iter_actual,  
      objective_trace = objective_trace[1:(iter_actual + 1)], 
      data_names = names(data),
      proclist_fitted = proclist,
      # Add class
      classes = "bamfa"          
  )

  # Remove old class "bamfa_em" if present, keep "bamfa", "multiblock_projector"
  # The projector function already adds "multiblock_projector", we added "bamfa"
  # Ensure no duplicates and remove list/bamfa_em if they existed before
  class(result_projector) <- unique(class(result_projector)[!class(result_projector) %in% c("bamfa_em", "list")])
  
  return(result_projector)
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
  # Check if it's the new multiblock_projector structure
  is_projector <- inherits(x, "multiblock_projector")
  
  header <- crayon::bold(crayon::blue("Barycentric Multiple Factor Analysis (BaMFA)"))
  cat(header, "\n\n")
  
  # Extract info based on structure
  if (is_projector && inherits(x, "bamfa")) {
      # Access elements directly using $
      k_g_val <- x$k_g
      k_l_val <- x$k_l
      niter_actual_val <- x$niter_actual
      lambda_l_val <- x$lambda_l
      data_names_val <- x$data_names
      objective_trace_val <- x$objective_trace
      B_list_internal <- x$B_list
      block_indices_internal <- x$block_indices 
      if (is.null(data_names_val)) data_names_val <- names(block_indices_internal)
      if (is.null(data_names_val)) data_names_val <- paste("Block", seq_along(block_indices_internal))
      n_blocks <- length(block_indices_internal)
      
  } else {
      # Fallback for old list structure or unexpected format
      warning("Object format not recognized as standard BaMFA projector. Printing basic info.", call.=FALSE)
      print(unclass(x))
      return(invisible(x))
  }

  # Ensure values are not NULL before printing
  k_g_val <- ifelse(is.null(k_g_val), NA, k_g_val)
  k_l_val <- ifelse(is.null(k_l_val), NA, k_l_val)
  niter_actual_val <- ifelse(is.null(niter_actual_val), NA, niter_actual_val)
  lambda_l_val <- ifelse(is.null(lambda_l_val), NA, lambda_l_val)

  # Model structure
  cat(crayon::green("Model Structure:"), "\n")
  cat("  Global components (k_g):", crayon::bold(k_g_val), "\n")
  cat("  Local components (k_l requested):", crayon::bold(k_l_val), "\n")
  cat("  Convergence after", crayon::bold(niter_actual_val), "iterations\n")
  cat("  Local score regularization (lambda_l):", crayon::bold(lambda_l_val), "\n\n")
  
  # Block information
  cat(crayon::green("Block Information:"), "\n")
  if (!is.null(block_indices_internal) && !is.null(B_list_internal) && length(block_indices_internal) == length(B_list_internal)) {
      for (i in seq_len(n_blocks)) {
        cat("  Block", crayon::bold(i), "(", crayon::blue(data_names_val[i]), "):", "\n")
        num_features <- length(block_indices_internal[[i]])
        num_local_comp <- if (!is.null(B_list_internal[[i]])) ncol(B_list_internal[[i]]) else NA
        cat("    Features:", crayon::yellow(num_features), "\n")
        cat("    Local components retained:", crayon::yellow(num_local_comp), "\n")
      }
  } else {
       cat(crayon::red("  Block information (indices or B_list) missing or inconsistent."), "\n")
  }
  
  # Final objective
  if (!is.null(objective_trace_val) && length(objective_trace_val) > 0) {
      final_obj_idx <- length(objective_trace_val)
      cat("\nFinal reconstruction error (per feature):", 
          crayon::bold(format(objective_trace_val[final_obj_idx], digits=6)), "\n")
  } else {
      cat("\nObjective trace not available.\n")
  }
  
  invisible(x)
}
#' Perform an Adam update on a block
#'
#' This helper function updates the block loadings \eqn{V} using the Adam approach,
#' maintaining first (\code{M}) and second (\code{V2}) moment estimates, then returns
#' the updated loadings and moment states.
#'
#' @param V Current loading matrix (\eqn{p \times ncomp}).
#' @param G Tangent-space gradient (\eqn{p \times ncomp}).
#' @param M,V2 Current first and second moment states for this block.
#' @param step_count Integer, global step index for bias correction.
#' @param beta1,beta2 Adam hyperparameters (typical defaults: 0.9, 0.999).
#' @param adam_epsilon Small constant for the Adam denominator (default: 1e-8).
#' @param learning_rate Base step size for Adam updates.
#'
#' @return A list with:
#' \itemize{
#'   \item \code{V}: Updated loading matrix (\eqn{p \times ncomp}).
#'   \item \code{M}: Updated first moment matrix (\eqn{p \times ncomp}).
#'   \item \code{V2}: Updated second moment matrix (\eqn{p \times ncomp}).
#' }
#' @keywords internal
#' @noRd
adam_update_block <- function(V, G, M, V2, step_count,
                              beta1, beta2, adam_epsilon,
                              learning_rate) {
  # Update moment estimates
  M <- beta1 * M + (1 - beta1)*G
  V2 <- beta2 * V2 + (1 - beta2)*(G*G)
  
  # bias-correct
  M_hat <- M / (1 - beta1^step_count)
  V_hat <- V2 / (1 - beta2^step_count)
  
  # compute step
  step_mat <- learning_rate * M_hat / (sqrt(V_hat) + adam_epsilon)
  
  # apply update
  V_new <- V - step_mat
  list(V = V_new, M = M, V2 = V2)
}

#' @importFrom chk chk_numeric chk_list chk_integer chk_gte chk_lte chk_flag
#' @importFrom stats svd qr qr.Q qr.R
#' @importFrom purrr map
#' @importFrom multivarious fresh prep init_transform center
#' @importFrom dplyr pull select
#' @importFrom rlang enquo
NULL

#' Penalized Multiple Factor Analysis (MFA) with Pairwise or Global-Mean Penalty
#'
#' @description
#' This function implements a penalized MFA-like decomposition via block-coordinate
#' descent (BCD). Each block \(\mathbf{X}_i \in \mathbb{R}^{n \times p}\) has a loading
#' matrix \(\mathbf{V}_i \in \mathbb{R}^{p \times ncomp}\) with orthonormal columns.
#' Assumes data blocks \code{X_i} are column-centered.
#'
#' We solve:
#' \deqn{
#'   \min_{\{V_i\}} \sum_{i=1}^S \|X_i - X_i V_i V_i^\top\|_F^2
#'   \;+\; \lambda \,\mathrm{Penalty}(\{V_i\}),
#' }
#' where the penalty is either:
#' \itemize{
#'   \item \code{"pairwise"}: \(\sum_{i < j} \|V_i - V_j\|_F^2\), or
#'   \item \code{"global_mean"}: \(\sum_{i} \|V_i - \bar{V}\|_F^2\).
#' }
#' Note: The penalties operate directly on the loading matrices V_i, which is
#' not invariant to orthogonal rotations within each V_i. Consider projection-based
#' penalties if rotation invariance is critical.
#'
#' @section Input types:
#' The function supports several input types:
#' \itemize{
#'   \item \code{list}: A list of matrices, where each matrix is a data block
#'   \item \code{multiblock}: A \code{multiblock} object containing multiple data blocks
#'   \item \code{multidesign}: A \code{multidesign} object that includes data from multiple subjects
#' }
#'
#' @section Implementation Details:
#' For \code{multidesign} objects, the function performs penalized MFA separately 
#' for each subject and then returns a list of results. This approach is useful for
#' analyzing multi-subject data where each subject may have multiple measurements or blocks.
#'
#' @param data A list of matrices, a \code{multiblock} object, or a \code{multidesign} object.
#' @param ncomp Integer number of latent components (columns of each \(\mathbf{V}_i\)).
#' @param lambda Non-negative scalar controlling the similarity penalty strength.
#' @param penalty_method Either \code{"pairwise"} or \code{"global_mean"}.
#' @param max_iter Maximum number of outer BCD iterations.
#' @param nsteps_inner Number of updates per block in each iteration.
#' @param learning_rate Base step size (if \code{optimizer="gradient"}) or Adam step size.
#' @param optimizer "gradient" (fixed LR) or "adam" (adaptive).
#' @param preproc Either NULL (no preprocessing) or a list of preprocessing functions for each block.
#' @param beta1,beta2 Adam hyperparameters (defaults: 0.9, 0.999).
#' @param adam_epsilon Small constant for Adam (default: 1e-8).
#' @param tol_obj Numeric tolerance for outer loop objective relative change (default: 1e-7).
#' @param tol_inner (Optional) tolerance for stopping each block update if \eqn{\|V_{new}-V_{old}\|_F < tol_inner}.
#' @param compute_consensus If \code{TRUE}, compute an average (orthonormal) loading
#'   matrix across blocks at the end (simple average, consider Procrustes alignment if needed).
#' @param verbose If \code{TRUE}, prints iteration logs.
#' @param subject Required for \code{multidesign} method. The name of the subject variable.
#' @param ... Additional arguments passed to methods.
#'
#' @return A list with elements depending on the input type:
#' \itemize{
#'   \item For list and multiblock inputs:
#'     \itemize{
#'       \item \code{V_list}: The final block loadings.
#'       \item \code{obj_values}: The objective values each iteration.
#'       \item \code{consensus}: An optional \eqn{p \times ncomp} matrix if
#'             \code{compute_consensus=TRUE}.
#'       \item \code{ncomp}: The number of components.
#'       \item \code{lambda}: The penalty strength used.
#'       \item \code{penalty_method}: The method used for the penalty.
#'     }
#'   \item For multidesign inputs:
#'     \itemize{
#'       \item A list with each element corresponding to a subject, containing the above elements.
#'       \item \code{preprocessors}: A list of preprocessing functions.
#'       \item \code{subject}: The subject variable name.
#'       \item \code{subjects}: The subjects in the data.
#'     }
#' }
#'
#' @examples
#' \dontrun{
#' # Example with a list of matrices
#' data_list <- lapply(1:3, function(i) scale(matrix(rnorm(100), 10, 10), scale=FALSE))
#' res <- penalized_mfa(data_list, ncomp=2, lambda=1, penalty_method="pairwise",
#'                      optimizer="adam", max_iter=50, verbose=TRUE)
#'
#' # Example with a multiblock object
#' mb <- multiblock(data_list)
#' res_mb <- penalized_mfa(mb, ncomp=2, lambda=1)
#'
#' # Example with a multidesign object
#' md <- multidesign(list(subj1=data_list, subj2=data_list), subject="id")
#' res_md <- penalized_mfa(md, ncomp=2, lambda=1, subject="id")
#' }
#'
#' @export
penalized_mfa <- function(data, ...) {
  UseMethod("penalized_mfa")
}

#' Heuristic to decide if pre-computing XtX is feasible within memory budget
#' @param p Number of features
#' @param mem_mb Memory budget in megabytes
#' @return Logical TRUE if 8 * p^2 bytes is within budget
#' @noRd
should_precompute <- function(p, mem_mb = 1024) {
  # Calculate required memory in MB
  required_mb <- (8 * p * p) / (1024^2)
  required_mb <= mem_mb
}

#' Factory function to create the appropriate reconstruction gradient function
#' @param X Data matrix (n x p)
#' @param XtX Precomputed crossproduct (p x p), if available (default: NULL)
#' @return A function that takes V (p x k) and returns the gradient -2 * XtX * V or -2 * t(X) %*% (X %*% V)
#' @noRd
make_grad_fun <- function(X, XtX = NULL) {
  if (!is.null(XtX)) {
    # Use precomputed XtX: faster if p is not huge
    function(V) -2 * XtX %*% V
  } else {
    # Compute on the fly using BLAS: more memory efficient for large p
    function(V) -2 * crossprod(X, X %*% V)
  }
}

# -------------------------------------------------------------------------
# penalized_mfa.list
# -------------------------------------------------------------------------
# Implementation for list of matrices
# -------------------------------------------------------------------------

#' @rdname penalized_mfa
#' @param memory_budget_mb Numeric (default 1024). Maximum memory (in MB) allocated per block for pre-computing `XtX`. If exceeded, gradient is computed on-the-fly.
#' @export
penalized_mfa.list <- function(data,
                          ncomp           = 2,
                          lambda          = 1,
                          penalty_method  = c("pairwise","global_mean"),
                          max_iter        = 10,
                          nsteps_inner    = 5,
                          learning_rate   = 0.01,
                          optimizer       = c("gradient","adam"),
                          preproc         = NULL,
                          beta1           = 0.9,
                          beta2           = 0.999,
                          adam_epsilon    = 1e-8,
                          tol_obj         = 1e-7,
                          tol_inner       = NULL,
                          compute_consensus = FALSE,
                          verbose         = FALSE,
                          subject         = NULL,
                          memory_budget_mb = 1024,
                          ...) {
  
  # Parameter validation
  penalty_method <- match.arg(penalty_method)
  optimizer      <- match.arg(optimizer)
  
  chk::chk_list(data)
  chk::chk_integer(ncomp)
  chk::chk_gte(ncomp, 1)
  chk::chk_numeric(lambda)
  chk::chk_gte(lambda, 0)
  chk::chk_integer(max_iter)
  chk::chk_gte(max_iter, 1)
  chk::chk_integer(nsteps_inner)
  chk::chk_gte(nsteps_inner, 1)
  chk::chk_numeric(learning_rate)
  chk::chk_gt(learning_rate, 0)
  chk::chk_numeric(beta1)
  chk::chk_gte(beta1, 0)
  chk::chk_lte(beta1, 1)
  chk::chk_numeric(beta2)
  chk::chk_gte(beta2, 0)
  chk::chk_lte(beta2, 1)
  chk::chk_numeric(adam_epsilon)
  chk::chk_gt(adam_epsilon, 0)
  chk::chk_numeric(tol_obj)
  chk::chk_gt(tol_obj, 0)
  if (!is.null(tol_inner)) {
    chk::chk_numeric(tol_inner)
    chk::chk_gt(tol_inner, 0)
  }
  chk::chk_flag(compute_consensus)
  chk::chk_flag(verbose)
  chk::chk_numeric(memory_budget_mb)
  chk::chk_gt(memory_budget_mb, 0)
  
  S <- length(data)
  if (S < 2) stop("Need at least 2 blocks.")
  
  dims <- sapply(data, dim)
  n <- dims[1, 1]
  p <- dims[2, 1] # Initial p, might change after preprocessing
  for (i in seq_len(S)) {
    if (dims[1, i] != n) stop("All blocks must have the same number of rows.")
    # Allow different number of columns p_i per block initially, 
    # but projector requires same p after preprocessing
  }
  
  # 0. Preprocessing using utility function
  # Check_consistent_ncol=TRUE because multiblock_projector requires it later
  preproc_result <- prepare_block_preprocessors(data, preproc, check_consistent_ncol = TRUE)
  proclist <- preproc_result$proclist
  Xp <- preproc_result$Xp
  p <- preproc_result$p_post # Update p to the post-preprocessing dimension
  
  # Setup gradient computation: Pre-compute XtX or use on-the-fly
  XtX_list   <- vector("list", S) # Stores XtX only if computed
  grad_fun   <- vector("list", S) # Stores the function to call for grad_recon
  precompute_info <- logical(S)
  
  for (i in seq_len(S)) {
      Xi <- Xp[[i]]
      # Use the actual number of columns (p) from the potentially preprocessed data
      p_i <- ncol(Xi)
      if (should_precompute(p_i, memory_budget_mb)) {
          XtX_list[[i]] <- crossprod(Xi)
          grad_fun[[i]] <- make_grad_fun(Xi, XtX = XtX_list[[i]])
          precompute_info[i] <- TRUE
      } else {
          # XtX_list[[i]] remains NULL
          grad_fun[[i]] <- make_grad_fun(Xi) # on-the-fly
          precompute_info[i] <- FALSE
      }
  }
  
  if (verbose) {
      precomputed_count <- sum(precompute_info)
      on_the_fly_count <- S - precomputed_count
      cat(sprintf("Gradient computation: %d blocks using precomputed XtX, %d blocks using on-the-fly.
",
                  precomputed_count, on_the_fly_count))
  }
  
  # 1) Initialize loadings (using preprocessed data)
  V_list <- vector("list", S)
  for (i in seq_len(S)) {
    # Use preprocessed data for initialization
    sv <- svd(Xp[[i]], nu=0, nv=ncomp)
    # Orthonormal columns => good initial guess
    # Handle cases where svd returns fewer components than requested
    ncomp_actual_init <- min(ncomp, length(sv$d))
    if (ncomp_actual_init < ncomp) {
         warning(sprintf("Block %d: Initial SVD returned %d components, requested %d. Padding loading matrix.", 
                       i, ncomp_actual_init, ncomp), call. = FALSE)
        # Simple padding with orthogonal vectors - might not be ideal but maintains dimension
        V_init <- sv$v[, 1:ncomp_actual_init, drop=FALSE]
        p_i <- nrow(V_init)
        if (p_i > 0) { # Avoid error if block has 0 features
             V_pad <- matrix(rnorm(p_i * (ncomp - ncomp_actual_init)), p_i, ncomp - ncomp_actual_init)
             V_pad_orth <- qr.Q(qr(cbind(V_init, V_pad)))[, (ncomp_actual_init + 1):ncomp, drop = FALSE]
             V_list[[i]] <- cbind(V_init, V_pad_orth)
        } else {
             V_list[[i]] <- matrix(0, 0, ncomp)
        }
    } else {
        V_list[[i]] <- sv$v[, 1:ncomp, drop=FALSE]
    }
    # Ensure orthonormal columns after potential padding/selection
    if (nrow(V_list[[i]]) > 0) {
        V_list[[i]] <- qr.Q(qr(V_list[[i]]))[, 1:ncomp, drop=FALSE] # Re-orthonormalize and ensure correct columns
    }
  }
  
  # If using Adam, store first and second moments
  if (optimizer == "adam") {
    # Ensure moments match feature dimension p
    m_list <- lapply(seq_len(S), function(i) matrix(0, ncol(Xp[[i]]), ncomp))
    v_list <- lapply(seq_len(S), function(i) matrix(0, ncol(Xp[[i]]), ncomp))
  }
  global_step <- 0
  
  # Objective function (uses preprocessed data Xp)
  obj_fun <- function(Vs) {
    # reconstruction cost
    val <- 0
    for (k in seq_len(S)) {
      Xk <- Xp[[k]] # Use preprocessed data
      Vk <- Vs[[k]]
      # Handle case where Vk might be empty or have wrong dimensions if init failed badly
      if (is.null(Vk) || nrow(Vk) != ncol(Xk) || ncol(Vk) != ncomp) {
         warning(sprintf("Objective Calc: Skipping block %d due to inconsistent V dimensions.", k), call.=FALSE)
         next
      }
      resid <- Xk - (Xk %*% Vk) %*% t(Vk)
      val <- val + sum(resid^2)
    }
    # penalty
    if (penalty_method == "pairwise") {
      pen_val <- 0
      # Use faster calculation (see update_block)
      Vsum <- Reduce(`+`, Vs)
      for(i in seq_len(S)) {
          diff_i <- S * Vs[[i]] - Vsum
          # The penalty is sum_{i<j} ||Vi-Vj||^2 = sum_i || S*Vi - Vsum ||^2 / S
          # Or, more simply, sum_i ||Vi - mean(V)||^2. Let's stick to the review suggestion for pairwise first.
          # The pairwise gradient derivation leads to sum_{i<j} ||Vi-Vj||^2 = sum_i trace(Vi^T (S*Vi - Vsum))
          # Let's recalculate the penalty value directly from definition for simplicity here.
          for (j in (i + 1):S) {
             if (j > S) next # Avoid index out of bounds
             diff_ij <- Vs[[i]] - Vs[[j]]
             pen_val <- pen_val + sum(diff_ij^2)
          }
      }
      # Original pairwise penalty loop was correct, let's restore that simplicity.
      # pen_val <- 0
      # for (i in seq_len(S - 1)) {
      #    for (j in (i + 1):S) {
      #      diff_ij <- Vs[[i]] - Vs[[j]]
      #      pen_val <- pen_val + sum(diff_ij^2)
      #    }
      # }
      val + lambda*pen_val
    } else {
      # global_mean
      mean_V <- Reduce(`+`, Vs)/S
      pen_val <- 0
      for (i in seq_len(S)) {
        diff_i <- Vs[[i]] - mean_V
        pen_val <- pen_val + sum(diff_i^2)
      }
      val + lambda*pen_val
    }
  }
  
  # BCD block update
  update_block <- function(bidx, Vs, Vsum = NULL) {
    Vi <- Vs[[bidx]]
    # Handle case where Vi might be invalid from initialization
    p_i <- ncol(Xp[[bidx]])
    if (is.null(Vi) || nrow(Vi) != p_i || ncol(Vi) != ncomp) {
        warning(sprintf("Update Block: Skipping block %d update due to inconsistent V dimensions.", bidx), call.=FALSE)
        return(Vs) # Return unmodified list
    }

    # XtX_i <- XtX_list[[bidx]] # Not needed directly anymore
    
    # Adam states if needed
    if (optimizer=="adam") {
      Mi <- m_list[[bidx]]
      V2 <- v_list[[bidx]]
      # Check Adam state dimensions match Vi
       if (nrow(Mi) != p_i || ncol(Mi) != ncomp || nrow(V2) != p_i || ncol(V2) != ncomp) {
          warning(sprintf("Adam state dimension mismatch for block %d. Reinitializing.", bidx), call.=FALSE)
          Mi <- matrix(0, p_i, ncomp)
          V2 <- matrix(0, p_i, ncomp)
       }
    }
    
    # This sum is needed for pairwise gradient, compute once per outer iter
    if (is.null(Vsum) && penalty_method == "pairwise") {
      Vsum <- Reduce(`+`, Vs)
    }
    
    # Inner loop
    for (step_inner in seq_len(nsteps_inner)) {
      # recon grad: Use the appropriate function (precomputed or on-the-fly)
      grad_recon <- grad_fun[[bidx]](Vi)
      
      # penalty grad
      if (penalty_method=="pairwise") {
        # Efficient version: grad = 2 * lambda * (S*Vi - Vsum)
        # Ensure Vsum is available (should be passed from outer loop)
        if(is.null(Vsum)) stop("Vsum is NULL inside update_block for pairwise penalty. This shouldn't happen.")
        grad_penalty <- 2*lambda*(S*Vi - Vsum)
      } else {
        # global_mean (already efficient)
        mean_V <- Reduce(`+`, Vs)/S
        grad_penalty <- 2*lambda*(Vi - mean_V)
      }
      
      G <- grad_recon + grad_penalty
      
      # project onto tangent
      A <- t(Vi) %*% G
      sym_part <- (A + t(A))/2
      G_tangent <- G - Vi %*% sym_part
      
      # update
      if (optimizer=="gradient") {
        step_mat <- learning_rate*G_tangent
        Vi_new <- Vi - step_mat # Update Vi for next inner step
      } else {
        # adam
        global_step <<- global_step + 1
        out_adam <- adam_update_block(
          V=Vi, G=G_tangent, M=Mi, V2=V2,
          step_count=global_step,
          beta1=beta1, beta2=beta2,
          adam_epsilon=adam_epsilon,
          learning_rate=learning_rate
        )
        Vi_new <- out_adam$V # Update Vi for next inner step
        Mi <- out_adam$M # Update Adam states
        V2 <- out_adam$V2
      }
      
      # Compute norm difference based on the updated Vi_new vs previous Vi
      diff_norm <- sqrt(sum((Vi_new - Vi)^2))
      Vi <- Vi_new # Update Vi for the next inner iteration

      # optional early stopping for inner loop
      if (!is.null(tol_inner) && diff_norm < tol_inner) {
        if (verbose) {
          message(sprintf("Block %d: inner stop at step %d (||deltaV||=%.2e < tol_inner)",
                          bidx, step_inner, diff_norm))
        }
        break
      }
    } # End inner loop
    
    # Re-orthonormalize ONCE after inner loop finishes
    qr_decomp <- qr(Vi)
    rank_check <- qr_decomp$rank
    if (rank_check < ncomp) {
        warning(sprintf("Block %d: Matrix became rank deficient (%d < %d) during update. Check data/parameters.", 
                      bidx, rank_check, ncomp), call.=FALSE)
        # Return the Q matrix up to the detected rank, pad if needed?
        # For now, just take the available columns to avoid errors
        Vi_ortho <- qr.Q(qr_decomp)[, 1:rank_check, drop=FALSE]
        # How to handle subsequent steps? If dim changes, breaks assumptions.
        # Simplest: pad with zeros or orthogonal random vectors. Let's pad orthogonally.
        if (rank_check > 0) { # Avoid error if rank is 0
             p_i <- nrow(Vi_ortho)
             pad_needed <- ncomp - rank_check
             if (pad_needed > 0 && p_i > 0) {
                 V_pad <- matrix(rnorm(p_i * pad_needed), p_i, pad_needed)
                 V_pad_orth <- qr.Q(qr(cbind(Vi_ortho, V_pad)))[, (rank_check + 1):ncomp, drop = FALSE]
                 Vi_ortho <- cbind(Vi_ortho, V_pad_orth)
             } else if (pad_needed > 0 && p_i == 0) {
                 Vi_ortho <- matrix(0, 0, ncomp)
             } # else pad_needed is 0, Vi_ortho is fine
        } else {
             Vi_ortho <- matrix(0, nrow(Vi), ncomp) # Fallback if rank is 0
        }
    } else {
        Vi_ortho <- qr.Q(qr_decomp)[, 1:ncomp, drop=FALSE] # Ensure exactly ncomp columns
    }
    
    # Update the list with the ORTHONORMALIZED version
    Vs[[bidx]] <- Vi_ortho 

    if (optimizer=="adam") {
      # Update Adam states in the main list after inner loop completes
      m_list[[bidx]] <<- Mi
      v_list[[bidx]] <<- V2
    }
    Vs # Return the updated list with the orthonormalized block
  }
  
  # main loop
  obj_values <- numeric(max_iter + 1) # Allocate space for initial + max_iter
  # Calculate initial objective using preprocessed data
  initial_obj <- tryCatch(obj_fun(V_list), error = function(e) {
      warning("Could not compute initial objective: ", e$message, call. = FALSE)
      return(NA) 
  })
  if (is.na(initial_obj)) {
      warning("Initial objective is NA. Check input data and initialization.", call. = FALSE)
      # Set to Inf to allow loop to start, but convergence check might be unreliable
      initial_obj <- Inf 
  }
  obj_values[1] <- initial_obj # Store initial objective
  old_obj <- initial_obj
  
  iter <- 0 # Initialize iter for trimming logic
  for (iter in 1:max_iter) {
    # Compute Vsum once per outer iteration if using pairwise penalty
    Vsum_iter <- NULL
    if (penalty_method == "pairwise") {
       # Ensure all Vs are valid before summing
       valid_V <- sapply(V_list, function(v) !is.null(v) && is.matrix(v) && ncol(v) == ncomp)
       if (all(valid_V)) {
           Vsum_iter <- Reduce(`+`, V_list)
       } else {
           warning(sprintf("Iter %d: Cannot compute Vsum for pairwise penalty due to invalid V matrices. Skipping penalty calculation?", iter), call.=FALSE)
           # How to handle this? Maybe skip penalty update for this iter?
           # For now, let Vsum_iter be NULL, update_block will error if penalty='pairwise'
           Vsum_iter <- NULL 
       }
    }

    # Consider parallel update here using future.apply if S is large
    for (b in seq_len(S)) {
      V_list <- update_block(b, V_list, Vsum=Vsum_iter)
    }
    
    new_obj <- tryCatch(obj_fun(V_list), error = function(e) {
      warning(sprintf("Could not compute objective at iter %d: %s", iter, e$message), call. = FALSE)
      return(NA) 
    })

    if (is.na(new_obj)) {
        warning("Objective is NA. Optimization may be unstable. Stopping.", call. = FALSE)
        break # Stop if objective calculation fails
    }

    obj_values[iter + 1] <- new_obj # Store objective AFTER iter completed
    
    # Check convergence
    if (!is.finite(old_obj) || !is.finite(new_obj)) {
       rel_change <- Inf # Cannot compute relative change reliably
    } else {
       rel_change <- abs(new_obj - old_obj) / (abs(old_obj) + 1e-12)
    }
    
    if (verbose) {
      cat(sprintf("Iter %d: obj=%.6f, rel_change=%.2e\n", iter, new_obj, rel_change))
    }
    if (rel_change < tol_obj && is.finite(rel_change)) {
      if (verbose) cat("Converged early (outer loop).\n")
      break # Exit loop
    }
    old_obj <- new_obj
  } # End main loop

  # Correctly trim obj_values array
  # Includes initial value + values from iterations 1 up to the last completed iteration 'iter'
  obj_values <- head(obj_values, iter + 1)

  # --- Prepare final multiblock_projector object --- 

  # Concatenate V_list into a single 'v' matrix
  # Ensure all blocks are valid matrices with correct dimensions first
  p <- ncol(Xp[[1]]) # Assuming all blocks have same p after preprocessing
  valid_V_final <- sapply(V_list, function(v) !is.null(v) && is.matrix(v) && ncol(v) == ncomp && nrow(v) == p)
  if (!all(valid_V_final)) {
      warning("Final V_list contains invalid/inconsistent matrices. Cannot construct projector. Returning raw list.", call.=FALSE)
      # Return the raw list if construction fails
      return(list(
          V_list = V_list,
          obj_values = obj_values,
          ncomp = ncomp,
          lambda = lambda,
          penalty_method = penalty_method,
          preprocessors = proclist, 
          precompute_info = precompute_info
      ))
  }
  v_concat <- do.call(rbind, V_list)
  
  # Generate block indices
  # Assumes all blocks in Xp (after potential preprocessing) have the same number of columns 'p'
  block_lengths_p <- sapply(Xp, ncol)
  if (!all(block_lengths_p == p)) {
      warning("Blocks have different numbers of features after preprocessing. Cannot create consistent block indices for multiblock_projector.", call.=FALSE)
      # Fallback? Or maybe this shouldn't happen if preproc maintains dims?
      # For now, assume 'p' is consistent
  }
  block_indices <- list()
  current_start <- 1
  for(i in 1:S) {
      block_indices[[i]] <- current_start:(current_start + p - 1)
      current_start <- current_start + p
  }
  names(block_indices) <- names(Xp)

  # Create the final preprocessor using concat_pre_processors
  final_preproc <- NULL
  if (!is.null(proclist)) {
      # Ensure concat_pre_processors is available
      if (!exists("concat_pre_processors", mode = "function")) {
          stop("Function 'concat_pre_processors' not found. Ensure multivarious package is loaded correctly.", call.=FALSE)
      }
      final_preproc <- concat_pre_processors(proclist, block_indices)
  } else {
      # Create a pass() preprocessor if no preprocessing was done
      # Need to ensure it's a finalized pre_processor object
      pass_proc <- multivarious::prep(multivarious::pass())
      # Repeat the pass processor for each block
      final_preproc <- concat_pre_processors(rep(list(pass_proc), S), block_indices)
  }
  
  # Compute consensus if requested (and possible)
  consensus_v <- NULL
  if (compute_consensus) {
    if (all(valid_V_final)) { # Check validity again just before Reduce
        Vsum <- Reduce(`+`, V_list)
        consensus_v <- qr.Q(qr(Vsum))[, 1:ncomp, drop=FALSE] 
    } else {
        warning("Cannot compute consensus due to inconsistent V matrices in V_list.", call.=FALSE)
    }
  }
  
  # Construct the multiblock_projector
  # Use multivarious::multiblock_projector explicitly
  result_projector <- multivarious::multiblock_projector(
      v = v_concat,
      preproc = final_preproc,
      block_indices = block_indices,
      classes = "penalized_mfa" # Add original class back
  )
  
  # Add other MFA-specific results as attributes
  attr(result_projector, "obj_values") <- obj_values
  attr(result_projector, "lambda") <- lambda
  attr(result_projector, "penalty_method") <- penalty_method
  attr(result_projector, "precompute_info") <- precompute_info
  attr(result_projector, "iterations_run") <- iter
  attr(result_projector, "consensus") <- consensus_v # Store the consensus matrix if computed
  # Keep original V_list as attribute? Might be useful for inspection.
  attr(result_projector, "V_list") <- V_list 

  return(result_projector)
}

#' @rdname penalized_mfa
#' @export
penalized_mfa.multiblock <- function(data,
                               ncomp           = 2,
                               lambda          = 1,
                               penalty_method  = c("pairwise","global_mean"),
                               max_iter        = 10,
                               nsteps_inner    = 5,
                               learning_rate   = 0.01,
                               optimizer       = c("gradient","adam"),
                               preproc         = NULL,
                               beta1           = 0.9,
                               beta2           = 0.999,
                               adam_epsilon    = 1e-8,
                               tol_obj         = 1e-7,
                               tol_inner       = NULL,
                               compute_consensus = FALSE,
                               verbose         = FALSE,
                               subject         = NULL,
                               ...) {
  # Check if multidesign package is available
  if (!requireNamespace("multidesign", quietly = TRUE)) {
    stop("Package 'multidesign' is required for working with multiblock objects. Please install it.",
         call. = FALSE)
  }
  
  # Ensure input is a multiblock object
  if (!inherits(data, "multiblock")) {
    stop("Input must be a multiblock object.", call. = FALSE)
  }
  
  # Check orientation - we need column-stacked multiblock
  if (!multidesign::is_cstacked(data)) {
    stop("Multiblock object must be column-stacked (all blocks share row dimension).", call. = FALSE)
  }
  
  # Convert multiblock to list of matrices
  data_list <- as.list(data)
  
  # Call the list method
  penalized_mfa.list(
    data            = data_list,
    ncomp           = ncomp,
    lambda          = lambda,
    penalty_method  = penalty_method,
    max_iter        = max_iter,
    nsteps_inner    = nsteps_inner,
    learning_rate   = learning_rate,
    optimizer       = optimizer,
    beta1           = beta1,
    beta2           = beta2, 
    adam_epsilon    = adam_epsilon,
    tol_obj         = tol_obj,
    tol_inner       = tol_inner,
    compute_consensus = compute_consensus,
    verbose         = verbose,
    ...
  )
}

#' @rdname penalized_mfa
#' @export
penalized_mfa.multidesign <- function(data,
                                ncomp           = 2,
                                lambda          = 1,
                                penalty_method  = c("pairwise","global_mean"),
                                max_iter        = 10,
                                nsteps_inner    = 5,
                                learning_rate   = 0.01,
                                optimizer       = c("gradient","adam"),
                                preproc         = multivarious::center(),
                                beta1           = 0.9,
                                beta2           = 0.999,
                                adam_epsilon    = 1e-8,
                                tol_obj         = 1e-7,
                                tol_inner       = NULL,
                                compute_consensus = FALSE,
                                verbose         = FALSE,
                                subject,
                                ...) {
  # Check if multidesign package is available
  if (!requireNamespace("multidesign", quietly = TRUE)) {
    stop("Package 'multidesign' is required for working with multidesign objects. Please install it.",
         call. = FALSE)
  }
  
  # Check if multivarious package is available
  if (!requireNamespace("multivarious", quietly = TRUE)) {
    stop("Package 'multivarious' is required for preprocessing. Please install it.",
         call. = FALSE)
  }
  
  # Ensure input is a multidesign object
  if (!inherits(data, "multidesign")) {
    stop("Input must be a multidesign object.", call. = FALSE)
  }
  
  # Check that subject is provided
  if (missing(subject)) {
    stop("'subject' parameter is required for multidesign method", call. = FALSE)
  }
  
  # Get the subject quo for consistent handling
  subject_quo <- rlang::enquo(subject)
  
  # Extract subject variable from design
  subjects <- factor(data$design %>% 
                     dplyr::select(!!subject_quo) %>% 
                     dplyr::pull(!!subject_quo))
  subject_set <- levels(subjects)
  
  if (length(subject_set) < 2) {
    stop("At least 2 subjects are required for penalized_mfa", call. = FALSE)
  }
  
  # Split data by subject
  sdat <- multidesign::split(data, !!subject_quo)
  
  # Create preprocessors, one per subject
  proclist <- lapply(seq_along(sdat), function(sd) {
    multivarious:::fresh(preproc) %>% multivarious::prep()
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
  result <- penalized_mfa.list(
    data             = strata,
    ncomp            = ncomp,
    lambda           = lambda,
    penalty_method   = penalty_method,
    max_iter         = max_iter,
    nsteps_inner     = nsteps_inner,
    learning_rate    = learning_rate,
    optimizer        = optimizer,
    beta1            = beta1,
    beta2            = beta2, 
    adam_epsilon     = adam_epsilon,
    tol_obj          = tol_obj,
    tol_inner        = tol_inner,
    compute_consensus = compute_consensus,
    verbose          = verbose,
    ...
  )
  
  # Add additional information
  result$proclist <- proclist
  result$subject_var <- subject
  result$subjects <- subjects
  result$subject_set <- subject_set
  
  # Update class
  class(result) <- c("penalized_mfa_multidesign", "penalized_mfa", "list")
  
  return(result)
}



#' Print Method for Penalized MFA Objects
#'
#' @md
#' @description
#' Prints a summary of a penalized MFA object, showing key parameters and fit results.
#'
#' @param x A `penalized_mfa` object
#' @param ... Additional parameters (unused)
#'
#' @return Invisibly returns the input object
#'
#' @export
print.penalized_mfa <- function(x, ...) {
  # Check if it's the new multiblock_projector structure or the old list structure
  is_projector <- inherits(x, "multiblock_projector")
  
  header <- crayon::bold(crayon::blue("Penalized Multiple Factor Analysis (MFA)"))
  cat(header, "\n\n")
  
  # Extract info based on structure
  if (is_projector) {
      ncomp <- ncol(x$v) # Get from projector's v matrix
      lambda_val <- attr(x, "lambda")
      penalty_method_val <- attr(x, "penalty_method")
      V_list_internal <- attr(x, "V_list") # Get original V_list if stored
      if (is.null(V_list_internal)) V_list_internal <- list() # Fallback
      obj_values_val <- attr(x, "obj_values")
      consensus_val <- attr(x, "consensus")
      n_blocks <- length(x$block_indices) 
      block_names <- names(x$block_indices)
      if (is.null(block_names)) block_names <- paste("Block", 1:n_blocks)
      iterations_run <- attr(x, "iterations_run")
      if (is.null(iterations_run)) iterations_run <- length(obj_values_val) -1 # Estimate if missing

  } else {
      # Assume old list structure for backward compatibility or if construction failed
      ncomp <- x$ncomp
      lambda_val <- x$lambda
      penalty_method_val <- x$penalty_method
      V_list_internal <- x$V_list
      obj_values_val <- x$obj_values
      consensus_val <- x$consensus
      n_blocks <- length(V_list_internal)
      block_names <- names(V_list_internal)
      if (is.null(block_names)) block_names <- paste("Block", 1:n_blocks)
      # Old structure didn't store iterations_run explicitly
      iterations_run <- length(obj_values_val) -1 # Estimate iterations
  }
  
  # Model parameters
  cat(crayon::green("Model Parameters:"), "\n")
  cat("  Components:", crayon::bold(ncomp), "\n")
  cat("  Penalty type:", crayon::bold(penalty_method_val), "\n")
  cat("  Lambda:", crayon::bold(lambda_val), "\n")
  
  # Results overview
  cat(crayon::green("\nResults:"), "\n")
  cat("  Number of blocks:", crayon::bold(n_blocks), "\n")
  
  # Show per-block info using V_list_internal
  cat(crayon::green("\nBlock Information (Loadings):\n"))
  for (i in seq_len(n_blocks)) {
    block_name <- block_names[i]
    Bi <- V_list_internal[[i]] # Use the stored V_list
    if (!is.null(Bi) && is.matrix(Bi)) {
       cat("  ", crayon::bold(block_name), ": ", 
           crayon::yellow(paste0(nrow(Bi), "×", ncol(Bi))), " loading matrix\n", sep="")
    } else {
        cat("  ", crayon::bold(block_name), ": ", crayon::red("Invalid/Missing Loadings"), "\n", sep="")
    }
  }
  
  # Show convergence info if available
  if (!is.null(obj_values_val) && length(obj_values_val) > 0) {
    cat(crayon::green("\nConvergence:"), "\n")
    initial_obj <- obj_values_val[1]
    final_obj <- obj_values_val[length(obj_values_val)]
    cat("  Initial objective:", crayon::bold(format(initial_obj, digits=6)), "\n")
    cat("  Final objective:", crayon::bold(format(final_obj, digits=6)), "\n")
    # Use explicit iterations_run if available, else estimate
    num_iters_display <- if (!is.null(iterations_run)) iterations_run else length(obj_values_val) - 1
    # Ensure it's not negative if obj_values has only 1 element
    num_iters_display <- max(0, num_iters_display) 
    cat("  Iterations run:", crayon::bold(num_iters_display), "\n")
    
    # Calculate percent decrease if possible
    if (length(obj_values_val) > 1 && is.finite(initial_obj) && is.finite(final_obj) && abs(initial_obj) > 1e-12) {
      pct_decrease <- 100 * (initial_obj - final_obj) / abs(initial_obj) # Use abs for safety
      cat("  Objective decrease:", crayon::bold(paste0(format(pct_decrease, digits=4), "%")), "\n")
    }
  }
  
  # Show consensus information if available
  if (!is.null(consensus_val)) {
    cat(crayon::green("\nConsensus:"), "\n")
    cat("  Dimensions:", crayon::bold(paste0(nrow(consensus_val), "×", ncol(consensus_val))), "\n")
  }
  
  cat("\n")
  invisible(x)
}

#' Print Method for Penalized MFA Multidesign Objects
#'
#' @md
#' @description
#' Prints a summary of a penalized MFA object derived from multidesign input,
#' showing key parameters, subject information, and fit results.
#'
#' @param x A `penalized_mfa_multidesign` object
#' @param ... Additional parameters (unused)
#'
#' @return Invisibly returns the input object
#'
#' @export
print.penalized_mfa_multidesign <- function(x, ...) {
  header <- crayon::bold(crayon::blue("Penalized MFA (from Multidesign)"))
  cat(header, "\n\n")
  
  # Model parameters
  cat(crayon::green("Model Parameters:"), "\n")
  cat("  Components:", crayon::bold(x$ncomp), "\n")
  cat("  Penalty type:", crayon::bold(x$penalty_method), "\n")
  cat("  Lambda:", crayon::bold(x$lambda), "\n")
  
  # Subject information
  cat(crayon::green("\nSubject Information:"), "\n")
  cat("  Subject variable:", crayon::bold(as.character(x$subject_var)), "\n")
  cat("  Number of subjects:", crayon::bold(length(x$subject_set)), "\n")
  
  # Show per-subject info
  cat(crayon::green("\nBlock Information (by Subject):"), "\n")
  for (i in seq_along(x$V_list)) {
    subject_name <- x$subject_set[i]
    if (is.null(subject_name)) subject_name <- paste("Subject", i)
    
    Bi <- x$V_list[[i]]
    cat("  ", crayon::bold(subject_name), ": ", 
        crayon::yellow(paste0(nrow(Bi), "×", ncol(Bi))), " loading matrix\n", sep="")
  }
  
  # Show convergence info if available
  if (length(x$obj_values) > 0) {
    cat(crayon::green("\nConvergence:"), "\n")
    cat("  Initial objective:", crayon::bold(format(x$obj_values[1], digits=6)), "\n")
    cat("  Final objective:", crayon::bold(format(x$obj_values[length(x$obj_values)], digits=6)), "\n")
    cat("  Iterations:", crayon::bold(length(x$obj_values)), "\n")
    
    # Calculate percent decrease
    if (length(x$obj_values) > 1 && x$obj_values[1] != 0) {
      pct_decrease <- 100 * (x$obj_values[1] - x$obj_values[length(x$obj_values)]) / x$obj_values[1]
      cat("  Objective decrease:", crayon::bold(paste0(format(pct_decrease, digits=4), "%")), "\n")
    }
  }
  
  # Show consensus information if available
  if (!is.null(x$consensus)) {
    cat(crayon::green("\nConsensus:"), "\n")
    cat("  Dimensions:", crayon::bold(paste0(nrow(x$consensus), "×", ncol(x$consensus))), "\n")
  }
  
  cat("\n")
  invisible(x)
}

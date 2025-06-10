#' Penalized MFA with Clusterwise Spatial Smoothness Constraints
#'
#' @description
#' A penalized MFA-like method for multi-subject cluster data, where each subject
#' \eqn{X_s \in \mathbb{R}^{n_s \times k_s}} has \code{k_s} clusters (columns).
#' Cluster coordinates (\code{coords_list}) are used to build a graph structure,
#' and the penalty encourages smoothness of loadings \eqn{V_s} across connected
#' clusters using the graph Laplacian. Assumes data blocks \code{X_s} are column-centered.
#'
#' We solve:
#' \deqn{
#'   \min_{\{V_s\}} \sum_{s=1}^S \|X_s - X_s V_s V_s^\top\|_F^2
#'   + \lambda \,\mathrm{trace}(V^\top L\,V),
#' }
#' where \(\mathbf{V}\) is the row-wise concatenation of \(\{\mathbf{V}_s\}\) and
#' \eqn{L = D - A} is the graph Laplacian derived from the adjacency matrix \eqn{A}
#' constructed by \code{spatial_constraints}. Minimizing this term encourages
#' loadings \eqn{v_i} and \eqn{v_j} to be similar if clusters \eqn{i} and \eqn{j}
#' are connected (\eqn{A_{ij}=1}).
#'
#' @param data_list A list of length \eqn{S}, each \eqn{n_s \times k_s}. Assumed column-centered.
#' @param coords_list A list of length \eqn{S}, each \eqn{k_s \times 3} of cluster coords.
#' @param ncomp Number of latent components per block.
#' @param lambda Nonnegative penalty weight for the spatial smoothness constraint.
#' @param adjacency_opts A list passed to \code{\link{spatial_constraints}} for adjacency construction.
#' @param max_iter Outer loop (BCD) steps.
#' @param nsteps_inner Inner steps (gradient or Adam) per block update.
#' @param learning_rate Base step size.
#' @param optimizer "gradient" or "adam".
#' @param beta1,beta2 Adam hyperparameters (defaults: 0.9, 0.999).
#' @param adam_epsilon Small constant for Adam denominator (default: 1e-8).
#' @param tol_obj Numeric tolerance for outer loop objective relative change (default: 1e-7).
#' @param tol_inner Tolerance for stopping block updates if \eqn{\|V_{new}-V_{old}\|_F < tol_inner}.
#' @param verbose If \code{TRUE}, prints iteration logs.
#' @param preproc Optional preprocessing for each block. Can be NULL (default), a single `prepper` object (applied independently to each block), or a list of `prepper` objects (one per block).
#' @param memory_budget_mb Numeric (default 1024). Maximum memory (in MB) allocated per block for pre-computing `XtX`. If exceeded, gradient is computed on-the-fly, and objective uses a different formula.
#'
#' @return A list with:
#' \itemize{
#'   \item \code{V_list}: A list of \eqn{S} final loadings \(\mathbf{V}_s\).
#'   \item \code{Sadj}: The graph Laplacian matrix \eqn{L=D-A} used for the penalty.
#'   \item \code{obj_values}: Outer iteration objective values.
#' }
#' @importFrom Matrix Diagonal rowSums
#' @importFrom pracma orth
#' @importFrom crayon bold blue green yellow red cyan
#' @export
penalized_mfa_clusterwise <- function(data_list,
                                      coords_list,
                                      ncomp         = 2,
                                      lambda        = 1,
                                      adjacency_opts= list(),
                                      max_iter      = 10,
                                      nsteps_inner  = 5,
                                      learning_rate = 0.01,
                                      optimizer     = c("gradient","adam"),
                                      beta1         = 0.9,
                                      beta2         = 0.999,
                                      adam_epsilon  = 1e-8,
                                      tol_obj       = 1e-7,
                                      tol_inner     = NULL,
                                      verbose       = FALSE,
                                      preproc       = NULL,
                                      memory_budget_mb = 1024) {
  
  # Check for Matrix package needed for sparse matrices
  if (!requireNamespace("Matrix", quietly = TRUE)) {
      stop("Package 'Matrix' needed for this function to work with sparse Laplacians. Please install it.", call. = FALSE)
  }

  optimizer <- match.arg(optimizer)
  # Parameter validation
  chk::chk_list(data_list)
  chk::chk_list(coords_list)
  chk::chk_integer(ncomp)
  chk::chk_gte(ncomp, 1)
  chk::chk_numeric(lambda)
  chk::chk_gte(lambda, 0)
  chk::chk_list(adjacency_opts)
  chk::chk_integer(max_iter)
  chk::chk_gte(max_iter, 1)
  chk::chk_integer(nsteps_inner)
  chk::chk_gte(nsteps_inner, 1)
  chk::chk_numeric(learning_rate)
  chk::chk_gt(learning_rate, 0)
  chk::chk_numeric(beta1); chk::chk_range(beta1, c(0, 1))
  chk::chk_numeric(beta2); chk::chk_range(beta2, c(0, 1))
  chk::chk_numeric(adam_epsilon); chk::chk_gt(adam_epsilon, 0)
  chk::chk_numeric(tol_obj); chk::chk_gt(tol_obj, 0)
  if (!is.null(tol_inner)) {
      chk::chk_numeric(tol_inner); chk::chk_gt(tol_inner, 0)
  }
  chk::chk_flag(verbose)
  chk::chk_numeric(memory_budget_mb); chk::chk_gt(memory_budget_mb, 0)


  S <- length(data_list)
  if (S < 2) stop("Need >=2 subjects.")
  if (length(coords_list) != S) stop("coords_list must match data_list length.")
  
  # Check data dimensions early using number of columns only
  k_s_vec <- sapply(data_list, ncol) # Number of columns (clusters) per subject
  if(!all(sapply(coords_list, nrow) == k_s_vec)) {
      stop("Number of rows in each coords_list element must match number of columns in corresponding data_list element.")
  }
  
  # Allow different number of columns k_s_vec per block initially
  # Preprocessing check below ensures consistency if needed by later steps

  # 0. Preprocessing using utility function
  # Check_consistent_ncol=FALSE initially, as clusterwise doesn't strictly require it,
  # but spatial constraints assume alignment based on original k_s_vec
  preproc_result <- prepare_block_preprocessors(data_list, preproc, check_consistent_ncol = FALSE)
  proclist <- preproc_result$proclist
  Xp <- preproc_result$Xp
  # Note: k_s_vec still refers to original column counts for spatial indexing
  # Preprocessing must NOT change the number of columns for spatial logic to work
  k_s_post_preproc <- sapply(Xp, ncol)
  if (!all(k_s_post_preproc == k_s_vec)) {
      stop("Preprocessing changed the number of columns (clusters) for some blocks. This invalidates the spatial constraints. Please use preprocessing steps that preserve the number of columns.", call.=FALSE)
  }
  
  # Build adjacency matrix A (using original k_s_vec and coords_list)
  # Ensure spatial_constraints is available (assuming it's in the same package or loaded)
  if (!exists("spatial_constraints", mode = "function")) {
      stop("Function 'spatial_constraints' not found. Ensure it is loaded.", call. = FALSE)
  }
  adj_opts <- modifyList(list(nblocks=S), adjacency_opts)
  Aadj <- do.call(spatial_constraints, c(list(coords_list), adj_opts))
  
  # Calculate Graph Laplacian L = D - A
  deg <- Matrix::rowSums(Aadj)
  Sadj <- Matrix::Diagonal(x = deg) - Aadj # Sadj now holds the Laplacian L
  bigK <- nrow(Sadj)
  
  # Offsets for indexing into Vbig and L
  offset <- cumsum(c(0, k_s_vec))
  if (offset[S+1] != bigK) {
      stop(sprintf("Total number of clusters (%d) does not match Laplacian dimension (%d).", offset[S+1], bigK))
  }
  # Create row index list for easy slicing of Vbig/LV
  row_indices_list <- lapply(seq_len(S), function(s) (offset[s]+1):offset[s+1])

  # Precompute XtX or NormX2 based on memory budget
  XtX_list <- vector("list", S)
  normX2_list <- numeric(S)
  grad_fun_list <- vector("list", S)
  precompute_info <- logical(S)

  for (s in seq_len(S)) {
      Xs <- Xp[[s]]
      k_s <- k_s_vec[s]
      # Skip if block is empty
      if (k_s == 0) {
          normX2_list[s] <- 0
          grad_fun_list[[s]] <- function(V) matrix(0, 0, ncomp) # Placeholder grad fun
          precompute_info[s] <- FALSE # Explicitly false for empty block
          next
      }

      if (should_precompute(k_s, memory_budget_mb)) {
          XtX_s <- crossprod(Xs)
          XtX_list[[s]] <- XtX_s
          normX2_list[s] <- sum(diag(XtX_s)) # Faster than sum(Xs^2)
          grad_fun_list[[s]] <- make_grad_fun(Xs, XtX = XtX_s)
          precompute_info[s] <- TRUE
      } else {
          # XtX_list[[s]] remains NULL
          normX2_list[s] <- sum(Xs^2)
          grad_fun_list[[s]] <- make_grad_fun(Xs) # on-the-fly
          precompute_info[s] <- FALSE
      }
  }

  if (verbose) {
      precomputed_count <- sum(precompute_info)
      on_the_fly_count <- S - precomputed_count
      cat(sprintf("Gradient/Objective: %d blocks using precomputed XtX, %d blocks using on-the-fly.
",
                  precomputed_count, on_the_fly_count))
  }
  
  # init loadings (handle empty blocks and padding)
  V_list <- vector("list", S)
  for (s in seq_len(S)) {
    Xs <- Xp[[s]]
    k_s <- k_s_vec[s]
    if(k_s == 0) {
        V_list[[s]] <- matrix(0, 0, ncomp) # Placeholder for empty block
        next
    }
    # Perform SVD safely
    sv <- tryCatch(svd(Xs, nu=0, nv=ncomp), error = function(e) {
        warning(sprintf("SVD failed during initialization for block %d: %s", s, e$message), call.=FALSE)
        NULL
    })
    
    if (is.null(sv) || length(sv$d) == 0) {
        warning(sprintf("Block %d: SVD init failed or returned no singular values. Initializing with random orthogonal matrix.", s), call.=FALSE)
        Vtemp <- pracma::orth(matrix(rnorm(k_s * ncomp), k_s, ncomp))
    } else {
        Vtemp <- sv$v
    }

    # Pad if rank < ncomp
    ncomp_svd <- ncol(Vtemp)
    if (ncomp_svd < ncomp) {
      extras <- ncomp - ncomp_svd
      Vextra <- matrix(rnorm(k_s * extras), k_s, extras)
      # Orthogonalize padded part against existing Vtemp and itself
      # Using QR for simplicity and robustness vs pracma::gramSchmidt
      Q_all <- qr.Q(qr(cbind(Vtemp, Vextra)))
      # Ensure correct dimensions even if original rank was < ncomp_svd
      Vtemp <- Q_all[, 1:ncomp, drop = FALSE]
    } else if (ncomp_svd > ncomp) {
        # Should not happen with svd(nv=ncomp), but check
        Vtemp <- Vtemp[, 1:ncomp, drop=FALSE]
    }
    
    # Final orthonormalization for the initial V_list[s]
    V_list[[s]] <- qr.Q(qr(Vtemp))[, 1:ncomp, drop=FALSE]
  }
  
  # Adam storage (initialize based on actual k_s)
  if (optimizer=="adam") {
    m_list <- lapply(seq_len(S), function(s) matrix(0, k_s_vec[s], ncomp))
    v_list <- lapply(seq_len(S), function(s) matrix(0, k_s_vec[s], ncomp))
  }
  global_step <- 0
  
  # Helper to combine list of V_s into one big V matrix
  # This version assumes Vl elements correspond to blocks 1:S and handles empty ones
  combine_loadings <- function(Vl) {
    # Check for consistency before combining
    if (length(Vl) != S) stop("Internal error: Vl length mismatch in combine_loadings")
    k_s_current <- sapply(Vl, nrow)
    if (!all(k_s_current == k_s_vec)) {
         warning("Row counts in Vl differ from expected k_s_vec. Using expected.", call.=FALSE)
         # This might indicate an issue, but proceed by creating Vbig based on k_s_vec
    }
    
    # Create Vbig using original offsets, filling with zeros for empty blocks
    out <- matrix(0, bigK, ncomp)
    for (s in seq_len(S)) {
        idx_out <- row_indices_list[[s]]
        # Ensure the block Vl[[s]] has the correct dimensions before assigning
        if(length(idx_out) > 0 && !is.null(Vl[[s]]) && nrow(Vl[[s]]) == length(idx_out) && ncol(Vl[[s]]) == ncomp) {
           out[idx_out, ] <- Vl[[s]]
        } else if (length(idx_out) > 0) {
            # If dimensions mismatch or Vl[[s]] is NULL, leave as zeros but warn
             warning(sprintf("Dimension mismatch or NULL V matrix for block %d in combine_loadings. Filling with zeros.", s), call.=FALSE)
        }
    }
    out
  }
  
  # Initial objective calculation requires a first pass of combine+LV
  Vbig_current <- combine_loadings(V_list)
  LV_current <- Sadj %*% Vbig_current # Initial LV for first iter spatial grad
  
  # Objective function using cached components
  calculate_objective <- function(Vl, Vbig, LV) {
      recon_total <- 0
      for(s in seq_len(S)) {
          Vs <- Vl[[s]]
          k_s <- k_s_vec[s]
          if (k_s == 0) next # Skip empty blocks
          
          # Check Vs validity before calculation
          if (is.null(Vs) || nrow(Vs) != k_s || ncol(Vs) != ncomp) {
              warning(sprintf("Objective: Skipping block %d due to invalid V dimensions.", s), call.=FALSE)
              recon_total <- recon_total + normX2_list[s] # Penalize fully if V is bad
              next
          }
          
          if (precompute_info[s]) {
              # Use XtX if available: ||X||^2 - tr(V' XtX V)
              XtXs <- XtX_list[[s]]
              recon_term <- normX2_list[s] - sum(diag(crossprod(Vs, XtXs %*% Vs)))
          } else {
              # Use on-the-fly: ||X||^2 - ||XV||^2
              Xs <- Xp[[s]]
              recon_term <- normX2_list[s] - sum((Xs %*% Vs)^2)
          }
          # Ensure non-negative reconstruction cost (numerical precision issues)
          recon_total <- recon_total + max(0, recon_term)
      }
      
      # Spatial penalty: tr(V' L V) = sum (LV * V)
      penalty_term <- lambda * sum(LV * Vbig) 
      
      return(recon_total + penalty_term)
  }

  # Project gradient G onto tangent space of Stiefel manifold at Vs
  stiefel_project <- function(Vs, G) {
    # Handle empty block
    if (nrow(Vs) == 0) return(matrix(0, 0, ncomp))
    A <- t(Vs) %*% G
    sym_part <- (A + t(A)) / 2
    G - Vs %*% sym_part
  }
  
  # Main BCD loop
  obj_values <- numeric(max_iter + 1)
  # Calculate initial objective value
  initial_obj <- tryCatch(calculate_objective(V_list, Vbig_current, LV_current),
                          error = function(e) {
                              warning("Could not compute initial objective: ", e$message, call. = FALSE)
                              return(NA)
                          })
  if (!is.finite(initial_obj)) { # Check for NA or Inf
      warning("Initial objective is not finite. Check input data and initialization.", call. = FALSE)
      initial_obj <- Inf # Allow loop to start
  }
  obj_values[1] <- initial_obj
  old_obj <- initial_obj

  iter <- 0 # Initialize iter for trimming logic
  for (iter in 1:max_iter) {
    
    # Capture V_list state at the BEGINNING of the iteration
    # This represents V_old for the off-diagonal spatial gradient part
    V_list_begin_iter <- V_list
    Vbig_begin_iter <- combine_loadings(V_list_begin_iter)

    # ---------------------------------------------
    # 1. Block update loop (without QR retraction)
    # ---------------------------------------------
    for (s in seq_len(S)) {
      k_s <- k_s_vec[s]
      # Skip update for empty blocks
      if (k_s == 0) next 

      Vs <- V_list[[s]] # Current (potentially un-retracted) V for block s
      idx_s <- row_indices_list[[s]] # Indices for this block in Vbig/LV
      
      # Skip if Vs is somehow ill-defined
      if (is.null(Vs) || nrow(Vs) != k_s || ncol(Vs) != ncomp) {
          warning(sprintf("Iter %d, Block %d: Skipping update due to inconsistent V dimensions.", iter, s), call.=FALSE)
          next
      }

      # Adam states if needed
      if (optimizer=="adam") {
        Mi <- m_list[[s]]
        V2 <- v_list[[s]]
        # Ensure Adam states match Vs dimensions (already checked in init, maybe redundant?)
        if(nrow(Mi) != k_s || ncol(Mi) != ncomp || nrow(V2) != k_s || ncol(V2) != ncomp) {
             warning(sprintf("Adam state dimension mismatch for block %d. Reinitializing.", s), call.=FALSE)
             Mi <- matrix(0, k_s, ncomp)
             V2 <- matrix(0, k_s, ncomp)
        }
      }

      # Inner loop for updating block s (without QR)
      for (step_inner in seq_len(nsteps_inner)) {
        # Reconstruction gradient (using appropriate function)
        G_recon <- grad_fun_list[[s]](Vs)
        
        # Spatial gradient (using LV from *previous* outer iteration)
        # Gradient of tr(V' L V) = 2 L V
        # Corrected spatial gradient: uses current Vs for diagonal block L_ss
        # and V from start of iteration for off-diagonal L_{s, not s}
        G_spatial_fresh_diag <- Sadj[idx_s, idx_s] %*% Vs
        # Need V_big representing state at start of *this* iteration
        G_spatial_off_diag <- Sadj[idx_s, -idx_s] %*% Vbig_begin_iter[-idx_s, , drop = FALSE]
        G_spatial <- 2 * (G_spatial_fresh_diag + G_spatial_off_diag)
        
        # Combine gradients
        G <- G_recon + lambda * G_spatial
        
        # Project onto tangent space
        G_t <- stiefel_project(Vs, G)
        
        # Update step (Adam or Gradient)
        if (optimizer=="gradient") {
          Vs_new <- Vs - learning_rate * G_t
        } else {
          # adam
          global_step <- global_step + 1 # Increment global step for Adam bias correction
          out_adam <- adam_update_block(
            V=Vs, G=G_t, M=Mi, V2=V2,
            step_count=global_step,
            beta1=beta1, beta2=beta2,
            adam_epsilon=adam_epsilon,
            learning_rate=learning_rate
          )
          Vs_new <- out_adam$V
          Mi <- out_adam$M # Update Adam states for next inner step
          V2 <- out_adam$V2
        }
        
        # Check for convergence in inner loop (optional)
        # Compare Vs_new with Vs BEFORE updating Vs
        if (!is.null(tol_inner)) {
             diff_norm <- sqrt(sum((Vs_new - Vs)^2))
             if (diff_norm < tol_inner) {
                 if (verbose) {
                     message(sprintf("Iter %d, Block %d: inner stop at step %d (||deltaV||=%.2e < tol_inner)",
                                     iter, s, step_inner, diff_norm))
                 }
                 Vs <- Vs_new # Update Vs before breaking
                 break
             }
        }
        
        # Update Vs for next inner iteration
        Vs <- Vs_new

      } # End inner loop

      # Store the potentially un-retracted Vs back into the list
      V_list[[s]] <- Vs 
      
      # Update Adam states in the main list after inner loop completes
      if (optimizer=="adam") {
        m_list[[s]] <- Mi
        v_list[[s]] <- V2
      }
    } # End loop over blocks s for updates

    # ---------------------------------------------
    # 2. QR Retraction Pass
    # ---------------------------------------------
    for (s in seq_len(S)) {
        k_s <- k_s_vec[s]
        if (k_s == 0) next # Skip empty blocks
        
        Vs_unretracted <- V_list[[s]]
        # Check if matrix is valid before QR
        if (is.null(Vs_unretracted) || !is.matrix(Vs_unretracted) || nrow(Vs_unretracted) != k_s || ncol(Vs_unretracted) != ncomp) {
             warning(sprintf("Iter %d, Block %d: Skipping QR retraction due to invalid matrix.", iter, s), call.=FALSE)
             # What to do? Keep the invalid matrix? Revert to old? Set to zero? Let's keep it for now.
             next
        }

        qr_decomp <- qr(Vs_unretracted)
        rank_check <- qr_decomp$rank
        
        if (rank_check < ncomp) {
            warning(sprintf("Iter %d, Block %d: Matrix became rank deficient (%d < %d) during update. Padding after QR.", 
                          iter, s, rank_check, ncomp), call.=FALSE)
            # Take available columns and pad orthogonally
            Vs_ortho_partial <- qr.Q(qr_decomp)[, 1:rank_check, drop=FALSE]
            if (rank_check > 0) { 
                 p_i <- nrow(Vs_ortho_partial)
                 pad_needed <- ncomp - rank_check
                 if (pad_needed > 0 && p_i > 0) {
                     V_pad <- matrix(rnorm(p_i * pad_needed), p_i, pad_needed)
                     V_pad_orth <- qr.Q(qr(cbind(Vs_ortho_partial, V_pad)))[, (rank_check + 1):ncomp, drop = FALSE]
                     Vs_ortho <- cbind(Vs_ortho_partial, V_pad_orth)
                 } else { # pad_needed is 0 or p_i is 0
                     Vs_ortho <- Vs_ortho_partial # Should have ncomp columns if pad_needed=0
                     if (p_i == 0) Vs_ortho <- matrix(0, 0, ncomp)
                 } 
            } else {
                 # Rank is 0
                 Vs_ortho <- matrix(0, k_s, ncomp) # Fallback if rank is 0
            }
        } else {
            # Full rank, ensure exactly ncomp columns
            Vs_ortho <- qr.Q(qr_decomp)[, 1:ncomp, drop=FALSE] 
        }
        # Store the ORTHONORMALIZED version back into V_list
        V_list[[s]] <- Vs_ortho
    }

    # -----------------------------------------------------------
    # 3. Compute Vbig and LV (for next iter grad & current obj)
    # -----------------------------------------------------------
    Vbig_current <- combine_loadings(V_list) # Use the now orthonormalized V_list
    LV_current <- Sadj %*% Vbig_current 

    # -----------------------------------------------------------
    # 4. Calculate Objective and Check Convergence
    # -----------------------------------------------------------
    new_obj <- tryCatch(calculate_objective(V_list, Vbig_current, LV_current),
                         error = function(e) {
                           warning(sprintf("Could not compute objective at iter %d: %s", iter, e$message), call. = FALSE)
                           return(NA) 
                         })

    if (is.na(new_obj)) {
        warning("Objective is NA. Optimization may be unstable. Stopping.", call. = FALSE)
        # Trim obj_values up to the previous iteration
        obj_values <- obj_values[1:iter] # iter because obj_values[1] is initial
        break # Stop if objective calculation fails
    }

    obj_values[iter + 1] <- new_obj # Store objective AFTER iter completed
    
    # Check convergence
    if (!is.finite(old_obj) || !is.finite(new_obj)) {
       rel_change <- Inf # Cannot compute relative change reliably
    } else {
       # Use relative change in objective function
       rel_change <- abs(new_obj - old_obj) / (abs(old_obj) + 1e-12)
    }
    
    if (verbose) {
      cat(sprintf("Iter %d: obj=%.6f, rel_change=%.2e\n", iter, new_obj, rel_change))
    }
    
    # Stop if converged
    if (rel_change < tol_obj && is.finite(rel_change)) {
      if (verbose) cat("Converged early (outer loop objective tolerance reached).\n")
      break # Exit loop
    }
    
    # Optional: Add check for gradient norm convergence? (More complex to implement efficiently here)
    
    old_obj <- new_obj

  } # End main loop (iter)

  # Final trim of obj_values array
  obj_values <- head(obj_values, iter + 1)

  # Prepare results list - converting to multiblock_projector
  
  # Concatenate V_list into a single 'v' matrix
  # Ensure validity first
  valid_V_final <- sapply(V_list, function(v) !is.null(v) && is.matrix(v) && ncol(v) == ncomp)
  # Check if number of rows match k_s_vec
  valid_rows <- mapply(function(v, k_s) { if (is.null(v)) FALSE else nrow(v) == k_s }, V_list, k_s_vec)
  
  if (!all(valid_V_final) || !all(valid_rows)) {
      warning("Final V_list contains invalid/inconsistent matrices. Cannot construct projector. Returning raw list.", call.=FALSE)
      # Return the raw list if construction fails
      return(list(
          V_list = V_list, 
          LV = LV_current, 
          Sadj = Sadj, 
          obj_values = obj_values, 
          ncomp = ncomp,
          lambda = lambda,
          precompute_info = precompute_info,
          iterations_run = iter
      ))
  }
  v_concat <- do.call(rbind, V_list)
  
  # Block indices are already computed in row_indices_list
  # Ensure names are set if data_list had names
  if (!is.null(names(data_list))) {
      names(row_indices_list) <- names(data_list)
  } else {
      names(row_indices_list) <- paste0("Block_", seq_along(row_indices_list))
  }

  # Create a pass() preprocessor, as this function assumes pre-processed input
  final_preproc <- multivarious::prep(multivarious::pass())

  # Construct the multiblock_projector
  result_projector <- multivarious::multiblock_projector(
      v = v_concat,
      preproc = final_preproc,
      block_indices = row_indices_list,
      classes = "penalized_mfa_clusterwise" # Add original class back
  )

  # Add other specific results as attributes
  attr(result_projector, "Sadj") <- Sadj
  attr(result_projector, "LV") <- LV_current
  attr(result_projector, "obj_values") <- obj_values
  attr(result_projector, "lambda") <- lambda
  attr(result_projector, "precompute_info") <- precompute_info
  attr(result_projector, "iterations_run") <- iter
  # Store original V_list for potential inspection
  attr(result_projector, "V_list") <- V_list
  
  return(result_projector)
}

#' Print Method for Penalized MFA Clusterwise Objects
#'
#' @md
#' @description
#' Prints a summary of a `penalized_mfa_clusterwise` object, showing key parameters 
#' (lambda, components), convergence info, and block details.
#'
#' @param x A `penalized_mfa_clusterwise` object (which inherits from `multiblock_projector`)
#' @param ... Additional parameters (unused)
#'
#' @return Invisibly returns the input object
#'
#' @export
print.penalized_mfa_clusterwise <- function(x, ...) {
  # Verify it's the expected structure
  if (!inherits(x, "multiblock_projector") || !inherits(x, "penalized_mfa_clusterwise")) {
    warning("Input object does not seem to be a valid penalized_mfa_clusterwise result.")
    print(unclass(x)) # Fallback print
    return(invisible(x))
  }
  
  header <- crayon::bold(crayon::blue("Penalized MFA with Clusterwise Smoothness"))
  cat(header, "\n\n")
  
  # Extract info from object and attributes
  ncomp <- ncol(x$v) 
  lambda_val <- attr(x, "lambda")
  V_list_internal <- attr(x, "V_list") 
  obj_values_val <- attr(x, "obj_values")
  n_blocks <- length(x$block_indices) 
  block_names <- names(x$block_indices)
  if (is.null(block_names)) block_names <- paste("Block", 1:n_blocks)
  iterations_run <- attr(x, "iterations_run")
  if (is.null(iterations_run)) iterations_run <- length(obj_values_val) - 1 # Estimate if missing
  precompute_info <- attr(x, "precompute_info")
  if (is.null(precompute_info)) precompute_info <- rep(NA, n_blocks) # Handle if missing
  
  # Model parameters
  cat(crayon::green("Model Parameters:"), "\n")
  cat("  Components (ncomp):", crayon::bold(ncomp), "\n")
  cat("  Smoothness Lambda:", crayon::bold(lambda_val), "\n")
  
  # Results overview
  cat(crayon::green("\nData Structure:"), "\n")
  cat("  Number of blocks/subjects:", crayon::bold(n_blocks), "\n")
  
  # Show per-block info using V_list_internal
  cat(crayon::green("\nBlock Information (Loadings & Gradient Method):\n"))
  for (i in seq_len(n_blocks)) {
    block_name <- block_names[i]
    Bi <- V_list_internal[[i]] # Use the stored V_list
    grad_method <- if (is.na(precompute_info[i])) "Unknown" else if (precompute_info[i]) "Precomputed XtX" else "On-the-fly"
    
    if (!is.null(Bi) && is.matrix(Bi)) {
       cat("  ", crayon::bold(block_name), ": ", 
           crayon::yellow(paste0(nrow(Bi), "×", ncol(Bi))), " loading matrix (Grad: ", crayon::cyan(grad_method), ")\n", sep="")
    } else {
        cat("  ", crayon::bold(block_name), ": ", crayon::red("Invalid/Missing Loadings"), 
            " (Grad: ", crayon::cyan(grad_method), ")\n", sep="")
    }
  }
  
  # Show convergence info if available
  if (!is.null(obj_values_val) && length(obj_values_val) > 0) {
    cat(crayon::green("\nConvergence:"), "\n")
    initial_obj <- obj_values_val[1]
    final_obj <- obj_values_val[length(obj_values_val)]
    # Use explicit iterations_run if available, else estimate
    num_iters_display <- if (!is.null(iterations_run)) iterations_run else length(obj_values_val) - 1
    num_iters_display <- max(0, num_iters_display) 
    
    cat("  Iterations run:", crayon::bold(num_iters_display), "\n")
    cat("  Initial objective:", crayon::bold(format(initial_obj, digits=6)), "\n")
    cat("  Final objective:", crayon::bold(format(final_obj, digits=6)), "\n")
    
    # Calculate percent decrease if possible
    if (length(obj_values_val) > 1 && is.finite(initial_obj) && is.finite(final_obj) && abs(initial_obj) > 1e-12) {
      pct_decrease <- 100 * (initial_obj - final_obj) / abs(initial_obj) 
      cat("  Objective decrease:", crayon::bold(paste0(format(pct_decrease, digits=4), "%")), "\n")
    }
  }
  
  return(invisible(x))
}

#' Compute a similarity matrix from blocks of data
#'
#' Creates a symmetric matrix where each element [i,j] is the similarity between
#' blocks i and j, calculated using the supplied function.
#'
#' @param blocks A list of numeric matrices or data frames
#' @param FUN Function to compute similarity between two blocks
#' @param ... Additional arguments passed to FUN
#' @return A symmetric similarity matrix with dimensions length(blocks) × length(blocks)
#' @importFrom utils combn
#' @noRd
#' @keywords internal
compute_sim_mat <- function(blocks, FUN, ...) {
  pairs <- combn(length(blocks),2)
  M <- matrix(0, length(blocks), length(blocks))
  
  for (i in 1:ncol(pairs)) {
    p1 <- pairs[1,i]
    p2 <- pairs[2,i]
    sim <- FUN(blocks[[p1]], blocks[[p2]], ...)
    M[p1,p2] <- sim
    M[p2,p1] <- sim
  }
  
  M
}

#' Calculate normalization factors for blocks in MFA
#'
#' Determines weighting factors for each block depending on the selected normalization method.
#'
#' @param blocks A list of preprocessed data blocks
#' @param type The normalization method to use: "MFA", "RV", "RV2", "None", or "Frob"
#' @return A numeric vector of normalization factors, one per block
#' @noRd
#' @keywords internal
normalization_factors <- function(blocks, type=c("MFA", "RV", "RV2", "None", "Frob")) {
  type <- match.arg(type)
  message("normalization type:", type)
  
  first_sv <- function(mat, method_pref = "irlba") {
    method <- .muscal_svd_method(mat, ncomp = 1)
    if (!identical(method_pref, "irlba") && !identical(method, "base")) {
      method <- method_pref
    }
    tryCatch(
      multivarious::svd_wrapper(mat, ncomp = 1, method = method)$sdev[1],
      error = function(e) multivarious::svd_wrapper(mat, ncomp = 1, method = "base")$sdev[1]
    )
  }
  
  alpha <- if (type == "MFA") {
    unlist(lapply(blocks, function(X) 1 / (first_sv(X)^2)))
  } else if (type == "RV") {
    smat <- compute_sim_mat(blocks, function(x1,x2) MatrixCorrelation::RV(x1,x2))
    diag(smat) <- 1
    u1 <- multivarious::svd_wrapper(smat, ncomp=1, method = .muscal_svd_method(smat, ncomp = 1))$u[,1]
    abs(u1) / sum(abs(u1))
  } else if (type == "RV2") {
    smat <- compute_sim_mat(blocks, function(x1,x2) MatrixCorrelation::RV2(x1,x2))
    diag(smat) <- 1
    u1 <- multivarious::svd_wrapper(smat, ncomp=1, method = .muscal_svd_method(smat, ncomp = 1))$u[,1]
    abs(u1) / sum(abs(u1))
  } else if (type == "Frob") {
    unlist(lapply(as.list(blocks), function(X) sum(X^2)))
  } else {
    rep(1, length(blocks))
  }
}


#' @md
#' @rdname mfa
#' @details
#' The `mfa.list` method applies the MFA algorithm to a list of data matrices or data frames.
#' This method first converts the list to a multiblock object and then calls `mfa.multiblock`.
#'
#' @examples
#' \donttest{
#' # Apply MFA to a list of matrices
#' X <- replicate(5, { matrix(rnorm(10*10), 10, 10) }, simplify=FALSE)
#' res <- mfa(X, ncomp=3, normalization="MFA")
#' }
#' @export
mfa.list <- function(data, preproc=center(), ncomp=2,
                     normalization=c("MFA", "RV", "None", "Frob", "custom"),
                     A=NULL, M=NULL, ...) {
  data <- multiblock(data)
  mfa.multiblock(data, preproc, ncomp, normalization, A, M,...)
}


#' @md
#' @rdname mfa
#' @details
#' The `mfa.multiblock` method implements Multiple Factor Analysis for a collection of 
#' data blocks. This method handles data preprocessing, block normalization, and integration 
#' of multiple data tables that share the same observations.
#' 
#' Normalization options include:
#' * `MFA`: Scales each block by its first singular value (default)
#' * `RV`: Normalizes blocks based on RV matrix correlation
#' * `None`: No scaling applied
#' * `Frob`: Uses Frobenius norm for scaling
#' * `custom`: Uses custom weight matrices provided via A and M parameters
#'
#' @section Missing data:
#' With `missing = "regularized"` or `"em"`, `NA`/`NaN` cells are imputed by
#' iterative MFA (Josse & Husson, 2012; the algorithm of
#' `missMDA::imputeMFA()`), and the returned model is an ordinary MFA fit to
#' the completed blocks. Missing cells start at their column means; each
#' iteration then refits the preprocessing, block weights, and an
#' `ncp_impute`-component MFA to the current completed data and replaces only
#' the missing cells with the low-rank reconstruction. `"em"` uses the plain
#' rank-`ncp_impute` reconstruction; `"regularized"` shrinks each retained
#' singular value `d_k` to `(d_k^2 - sigma2) / d_k`, where `sigma2` is the
#' residual inertia per residual degree of freedom, which guards against
#' overfitting when data are noisy or many cells are missing.
#'
#' Rows may be missing entire blocks, but every column needs at least one
#' observed value, blocks must share the same rows (see [linked_mfa()] for
#' blocks observed on different row sets), and infinite values are rejected.
#' `ncp_impute` is capped at `min(n - 2, p - 1)` with a warning, and a warning
#' is issued if the imputation does not converge within `missing_maxiter`
#' iterations. Convergence is the root-mean-square change of the imputed cells,
#' in units of each variable's observed standard deviation, falling below
#' `missing_tol`. Imputed cells influence the final scores and loadings, so
#' choose `ncp_impute` with care (it defaults to `ncomp`). The result records
#' `missing = list(method, ncp_impute, iterations, converged, delta, mask)`.
#'
#' @param missing Missing-data handling mode. `"error"` (default) rejects `NA`
#'   values. `"regularized"` and `"em"` run iterative MFA imputation before the
#'   final complete-data MFA fit (see the *Missing data* section). `"nipals"`
#'   is reserved for a future direct available-data backend.
#' @param ncp_impute Integer; number of MFA components used inside the iterative
#'   imputation loop. Values above `min(n - 2, p - 1)` are capped.
#' @param missing_tol Numeric convergence tolerance for iterative imputation.
#' @param missing_maxiter Integer maximum number of imputation iterations.
#' @param return_imputed Logical; if `TRUE`, attach completed blocks as
#'   `imputed_data`.
#'
#' @examples
#' \donttest{
#' # Create 5 random matrices of the same size
#' X <- replicate(5, { matrix(rnorm(10*10), 10, 10) }, simplify=FALSE)
#'
#' # Apply MFA with MFA normalization
#' res <- mfa(X, ncomp=3, normalization="MFA")
#'
#' # Project a block onto the model
#' p <- multivarious::project_block(res, X[[1]], 1)
#'
#' # Verify number of components
#' stopifnot(ncol(multivarious::scores(res)) == 3)
#'
#' # Create a classifier
#' labs <- letters[1:10]
#' cfier <- multivarious::classifier(res, new_data=do.call(cbind, X), labels=labs)
#' pred <- predict(cfier, do.call(cbind, X)[1:2,])
#'
#' # Create a classifier using a specific block
#' cfier2 <- multivarious::classifier(res, new_data=X[[2]], labels=labs,
#'                                   colind=res$block_indices[[2]])
#' pred2 <- predict(cfier2, X[[2]][1:2,])
#' }
#' @export
mfa.multiblock <- function(data, preproc=center(), ncomp=2,
                normalization=c("MFA", "RV", "None", "Frob", "custom"),
                A=NULL, M=NULL,
                missing = c("error", "regularized", "em", "nipals"),
                ncp_impute = ncomp,
                missing_tol = 1e-6,
                missing_maxiter = 100,
                return_imputed = FALSE,
                ...) {
  fit_call <- match.call(expand.dots = FALSE)
  fit_dots <- list(...)
  missing <- match.arg(missing)

  has_na <- any(vapply(data, anyNA, logical(1L)))
  if (!has_na) {
    return(.mfa_complete(
      data = data,
      preproc = preproc,
      ncomp = ncomp,
      normalization = normalization,
      A = A,
      M = M,
      fit_call = fit_call,
      fit_dots = fit_dots,
      ...
    ))
  }

  if (identical(missing, "error")) {
    stop(
      "mfa() does not accept NA values unless `missing` is set to ",
      "'regularized', 'em', or 'nipals'.",
      call. = FALSE
    )
  }

  if (identical(missing, "nipals")) {
    stop(
      "missing = 'nipals' is not implemented yet; use missing = 'regularized' or 'em'.",
      call. = FALSE
    )
  }

  missing_controls <- .mfa_validate_missing_controls(
    ncp_impute = ncp_impute,
    missing_tol = missing_tol,
    missing_maxiter = missing_maxiter,
    return_imputed = return_imputed
  )
  ncp_impute <- missing_controls$ncp_impute
  missing_tol <- missing_controls$missing_tol
  missing_maxiter <- missing_controls$missing_maxiter
  return_imputed <- missing_controls$return_imputed

  imp <- .mfa_impute_iterative(
    data = data,
    preproc = preproc,
    ncp = ncp_impute,
    normalization = normalization,
    A = A,
    M = M,
    method = missing,
    tol = missing_tol,
    maxiter = missing_maxiter,
    ...
  )

  fit <- .mfa_complete(
    data = imp$data,
    preproc = preproc,
    ncomp = ncomp,
    normalization = normalization,
    A = A,
    M = M,
    fit_call = fit_call,
    fit_dots = fit_dots,
    ...
  )
  ncp_impute <- imp$ncp
  fit$missing <- list(
    method = missing,
    ncp_impute = ncp_impute,
    iterations = imp$iterations,
    converged = imp$converged,
    delta = imp$delta,
    mask = imp$mask
  )
  fit$fit_spec$refit <- .mfa_missing_refit_spec(
    data = data,
    preproc = preproc,
    ncomp = ncomp,
    normalization = normalization,
    A = A,
    M = M,
    missing = missing,
    ncp_impute = ncp_impute,
    missing_tol = missing_tol,
    missing_maxiter = missing_maxiter,
    return_imputed = return_imputed,
    fit_dots = fit_dots
  )
  fit$fit_spec$refit_supported <- TRUE
  if (isTRUE(return_imputed)) {
    fit$imputed_data <- imp$data
  }
  fit
}

.mfa_complete <- function(data, preproc=center(), ncomp=2,
                normalization=c("MFA", "RV", "None", "Frob", "custom"),
                A=NULL, M=NULL,
                build_diagnostics = TRUE,
                fit_call = NULL,
                fit_dots = list(),
                ...) {
  if (is.null(fit_call)) {
    fit_call <- match.call(expand.dots = FALSE)
  }
  if (is.null(fit_dots)) {
    fit_dots <- list(...)
  }
  A_input <- A
  M_input <- M
  
  
  chk::chk_true(length(data) > 1)
  for (i in 1:length(data)) {
    chk::chkor_vld(chk::vld_matrix(data[[i]]), chk::vld_s4_class(data[[i]], "Matrix"))
  }
  
  nrs <- sapply(data, nrow)
  chk::chk_true(all(nrs == nrs[1]))
  nr <- nrs[1]
  
  normalization <- match.arg(normalization)
  
  if (normalization == "custom" && is.null(A) && is.null(M)) {
    stop('At least one of A or M must be provided when normalization = "custom"')
  }
  
  S <- length(data)
  if (is.null(names(data))) {
    names(data) <- paste0("B", 1:S)
  }
  data_refit <- lapply(data, function(x) as.matrix(x))
  
  # Preprocessing using utility function
  # Set check_consistent_ncol=FALSE as MFA handles potential concatenation later
  preproc_result <- prepare_block_preprocessors(data, preproc, check_consistent_ncol = FALSE)
  proclist <- preproc_result$proclist
  fitted_proclist <- .muscal_materialize_block_preprocessors(data, proclist)
  strata <- preproc_result$Xp # Renamed from Xp for consistency with original MFA code
  
  ## calculate block normalization factors
  if (normalization != "custom") {
    alpha <- normalization_factors(strata, type=normalization)
    A <- rep(alpha, sapply(strata, ncol))
  } else {
    alpha <- rep(1, length(strata))
  }
  
  ## compute block indicees
  block_indices <- list()
  ind <- 1
  for (i in 1:length(strata)) {
    block_indices[[i]] <- seq(ind, ind+ncol(strata[[i]])-1)
    ind <- ind + ncol(strata[[i]])
  }
  
  proc <- multivarious::concat_pre_processors(fitted_proclist, block_indices)

    ## fit genpca
  if (!requireNamespace("genpca", quietly = TRUE)) {
    stop("Package 'genpca' needed for MFA analysis. Please install it.", call. = FALSE)
  }
  Xp <- do.call(cbind, strata)
  fit <- genpca::genpca(Xp, 
            preproc=multivarious::pass(),
            A=A, 
            M=M,
            ncomp=ncomp,
            ...)
  
  fit[["block_indices"]] <- block_indices
  fit[["alpha"]] <- alpha
  fit[["normalization"]] <- normalization
  fit[["names"]] <- names(data)
  
  ## this is awkward...
  ## instead, we need a "delegation" mechanism, where a multiblock projector simply wraps a projector
  ## here, we rely on the fact that we use "pass()" pre-processing for inner genpca fit
  fit[["preproc"]] <- proc

  # -------------------------------------------------------------------------
  # Derived quantities for plotting / diagnostics
  # -------------------------------------------------------------------------
  if (isTRUE(build_diagnostics)) {
    partial_scores <- lapply(seq_along(strata), function(i) {
      idx <- block_indices[[i]]
      strata[[i]] %*% fit$v[idx, , drop = FALSE]
    })
    names(partial_scores) <- names(data)

    cor_loadings <- tryCatch(
      {
        C <- stats::cor(Xp, fit$s)
        C[!is.finite(C)] <- 0
        C
      },
      error = function(e) NULL
    )

    rv <- tryCatch(
      {
        M <- compute_sim_mat(strata, function(x1, x2) MatrixCorrelation::RV(x1, x2))
        diag(M) <- 1
        rownames(M) <- colnames(M) <- names(data)
        M
      },
      error = function(e) NULL
    )

    rv2 <- tryCatch(
      {
        M <- compute_sim_mat(strata, function(x1, x2) MatrixCorrelation::RV2(x1, x2))
        diag(M) <- 1
        rownames(M) <- colnames(M) <- names(data)
        M
      },
      error = function(e) NULL
    )
  } else {
    partial_scores <- NULL
    cor_loadings <- NULL
    rv <- NULL
    rv2 <- NULL
  }

  # Construct the final multiblock_biprojector
  mfa_result <- multivarious::multiblock_biprojector(
      v = fit$v,              # Loadings from genpca (concatenated space)
      s = fit$s,              # Scores from genpca
      sdev = fit$sdev,        # Singular values from genpca
      preproc = proc,         # Use the concatenated block-aware preprocessor
      block_indices = block_indices,
      # Pass MFA specific info via ...
      alpha = alpha,
      normalization = normalization,
      names = names(data),
      partial_scores = partial_scores,
      cor_loadings = cor_loadings,
      rv = rv,
      rv2 = rv2,
      # Pass relevant genpca info via ... as well
      ou = fit$ou,
      ov = fit$ov,
      A_genpca = fit$A, # Renamed to avoid conflict if A was input
      M_genpca = fit$M,
      # Add the class
      classes = "mfa"
  )
  
  mfa_result <- .muscal_attach_fit_contract(
    mfa_result,
    method = "mfa",
    task = "reconstruction",
    oos_types = c("scores", "reconstruction"),
    fit_call = fit_call,
    refit_supported = TRUE,
    prediction_target = "blocks",
    refit = .muscal_make_refit_spec(
      data = data_refit,
      fit_fn = function(data) {
        do.call(
          mfa,
          c(
            list(
              data = data,
              preproc = preproc,
              ncomp = ncomp,
              normalization = normalization,
              A = A_input,
              M = M_input
            ),
            fit_dots
          )
        )
      },
      bootstrap_fn = function(data) {
        n <- nrow(data[[1L]])
        idx <- sample.int(n, size = n, replace = TRUE)
        lapply(data, function(block) block[idx, , drop = FALSE])
      },
      permutation_fn = function(data) {
        lapply(data, function(block) block[sample.int(nrow(block)), , drop = FALSE])
      },
      resample_unit = "rows"
    )
  )

  return(mfa_result)
}

.mfa_impute_iterative <- function(data, preproc, ncp, normalization, A, M,
                                  method = c("regularized", "em"),
                                  tol = 1e-6, maxiter = 100, ...) {
  method <- match.arg(method)
  X <- lapply(data, function(x) as.matrix(x))
  if (is.null(names(X))) {
    names(X) <- paste0("B", seq_along(X))
  }
  nrs <- vapply(X, nrow, integer(1L))
  if (any(nrs != nrs[[1L]])) {
    stop("All blocks must have the same number of rows.", call. = FALSE)
  }
  mask <- lapply(X, is.na)

  for (b in seq_along(X)) {
    if (any(is.infinite(X[[b]]))) {
      stop(
        sprintf(
          "Block '%s' contains infinite values; only NA/NaN cells are treated as missing.",
          names(X)[[b]]
        ),
        call. = FALSE
      )
    }
    bad_cols <- which(colSums(!mask[[b]]) == 0L)
    if (length(bad_cols) > 0L) {
      stop(
        sprintf(
          "Cannot impute variable/column %d with all values missing in block '%s'.",
          bad_cols[[1L]], names(X)[[b]]
        ),
        call. = FALSE
      )
    }
  }

  empty_rows <- which(Reduce(`&`, lapply(mask, function(m) rowSums(!m) == 0L)))
  if (length(empty_rows) > 0L) {
    warning(
      sprintf(
        paste0(
          "%d row(s) have no observed values in any block (first: row %d); ",
          "their imputed values and scores are uninformative."
        ),
        length(empty_rows), empty_rows[[1L]]
      ),
      call. = FALSE
    )
  }

  # The residual-variance estimate behind the regularized step needs positive
  # degrees of freedom, (n - 1 - ncp) * (p - ncp) > 0, and an (ncp + 1)-th
  # component; this is the same bound missMDA::imputeMFA enforces.
  n <- nrs[[1L]]
  p <- sum(vapply(X, ncol, integer(1L)))
  max_ncp <- min(n - 2L, p - 1L)
  if (max_ncp < 1L) {
    stop("Iterative imputation needs at least 3 rows and 2 columns.", call. = FALSE)
  }
  if (ncp > max_ncp) {
    warning(
      sprintf(
        "`ncp_impute` = %d is too large for %d rows and %d columns; using %d.",
        ncp, n, p, max_ncp
      ),
      call. = FALSE
    )
    ncp <- max_ncp
  }

  # Convergence is measured on the imputed cells in units of each variable's
  # observed standard deviation, so it does not depend on location or scale.
  miss_scale <- unlist(Map(function(block, m) {
    s <- apply(block, 2L, stats::sd, na.rm = TRUE)
    s[!is.finite(s) | s <= 0] <- 1
    s[col(m)[m]]
  }, X, mask), use.names = FALSE)

  Ximp <- lapply(X, function(block) {
    out <- block
    cm <- colMeans(out, na.rm = TRUE)
    miss <- is.na(out)
    if (any(miss)) {
      out[miss] <- cm[col(out)[miss]]
    }
    out
  })

  # The regularized step needs the (ncp + 1)-th singular value.
  fit_ncomp <- if (identical(method, "regularized")) ncp + 1L else ncp
  converged <- FALSE
  delta <- Inf
  iter <- 0L
  for (iter in seq_len(maxiter)) {
    fit <- suppressMessages(.mfa_complete(
      data = Ximp,
      preproc = preproc,
      ncomp = fit_ncomp,
      normalization = normalization,
      A = A,
      M = M,
      build_diagnostics = FALSE,
      ...
    ))
    Xhat <- .mfa_reconstruct_blocks(fit, Ximp, ncomp = ncp, method = method)

    old_missing <- unlist(Map(function(x, m) x[m], Ximp, mask), use.names = FALSE)
    for (b in seq_along(Ximp)) {
      Ximp[[b]][mask[[b]]] <- Xhat[[b]][mask[[b]]]
    }
    new_missing <- unlist(Map(function(x, m) x[m], Ximp, mask), use.names = FALSE)
    delta <- sqrt(mean(((new_missing - old_missing) / miss_scale)^2))
    if (is.finite(delta) && delta < tol) {
      converged <- TRUE
      break
    }
  }

  if (!converged) {
    warning(
      sprintf(
        paste0(
          "MFA missing-data imputation did not converge in %d iterations ",
          "(delta = %.3g, missing_tol = %.3g); consider increasing ",
          "`missing_maxiter` or lowering `ncp_impute`."
        ),
        maxiter, delta, tol
      ),
      call. = FALSE
    )
  }

  list(
    data = Ximp,
    mask = mask,
    ncp = ncp,
    iterations = iter,
    converged = converged,
    delta = delta
  )
}

.mfa_validate_missing_controls <- function(ncp_impute, missing_tol,
                                           missing_maxiter, return_imputed) {
  if (!is.numeric(ncp_impute) ||
      length(ncp_impute) != 1L ||
      is.na(ncp_impute) ||
      !is.finite(ncp_impute) ||
      ncp_impute < 1 ||
      abs(ncp_impute - round(ncp_impute)) > 1e-8) {
    stop("`ncp_impute` must be a single integer-like value >= 1.", call. = FALSE)
  }

  if (!is.numeric(missing_tol) ||
      length(missing_tol) != 1L ||
      is.na(missing_tol) ||
      !is.finite(missing_tol) ||
      missing_tol <= 0) {
    stop("`missing_tol` must be a single finite value > 0.", call. = FALSE)
  }

  if (!is.numeric(missing_maxiter) ||
      length(missing_maxiter) != 1L ||
      is.na(missing_maxiter) ||
      !is.finite(missing_maxiter) ||
      missing_maxiter < 1 ||
      abs(missing_maxiter - round(missing_maxiter)) > 1e-8) {
    stop("`missing_maxiter` must be a single integer-like value >= 1.", call. = FALSE)
  }

  if (!is.logical(return_imputed) || length(return_imputed) != 1L || is.na(return_imputed)) {
    stop("`return_imputed` must be `TRUE` or `FALSE`.", call. = FALSE)
  }

  list(
    ncp_impute = as.integer(round(ncp_impute)),
    missing_tol = as.numeric(missing_tol),
    missing_maxiter = as.integer(round(missing_maxiter)),
    return_imputed = isTRUE(return_imputed)
  )
}

.mfa_reconstruct_blocks <- function(fit, data, ncomp, method = c("regularized", "em")) {
  method <- match.arg(method)
  d <- fit$sdev
  comp <- seq_len(min(as.integer(ncomp), length(d)))

  shrink <- if (identical(method, "em")) {
    rep(1, length(comp))
  } else {
    .mfa_regularized_shrinkage(
      d = d,
      ncp = length(comp),
      n = nrow(data[[1L]]),
      p = sum(vapply(data, ncol, integer(1L))),
      total_inertia = .mfa_total_inertia(fit, data)
    )
  }

  # `fit` was trained on `data`, so its scores are the in-sample `ou D`.
  scores <- multivarious::scores(fit)[, comp, drop = FALSE]
  ov <- .mfa_data_axes(fit)[, comp, drop = FALSE]
  Xhat_p <- sweep(scores, 2L, shrink, `*`) %*% t(ov)
  Xhat <- multivarious::inverse_transform(fit$preproc, Xhat_p)

  lapply(fit$block_indices, function(idx) Xhat[, idx, drop = FALSE])
}

# Regularized iterative MFA shrinkage (Josse & Husson, 2012; missMDA::imputeMFA):
# each retained singular value d_k becomes (d_k^2 - sigma2) / d_k, where sigma2
# is the residual inertia per residual degree of freedom, capped at d_{ncp+1}^2.
# Returned as multiplicative factors on the scores `ou D`.
.mfa_regularized_shrinkage <- function(d, ncp, n, p, total_inertia) {
  d_k <- d[seq_len(ncp)]
  resid <- max(total_inertia - sum(d_k^2), 0)
  df <- (n - 1 - ncp) * (p - ncp)
  sigma2 <- if (df > 0) n * p / min(p, n - 1) * resid / df else 0
  d_next <- if (length(d) > ncp) d[[ncp + 1L]] else 0
  sigma2 <- min(sigma2, d_next^2)
  shrink <- ifelse(d_k > 0, (d_k^2 - sigma2) / d_k^2, 0)
  pmax(shrink, 0)
}

# Total generalized inertia trace(X' M X A) of the preprocessed data, i.e. the
# sum of all squared generalized singular values of the fitted MFA.
.mfa_total_inertia <- function(fit, data) {
  Xp <- as.matrix(multivarious::reprocess(fit, do.call(cbind, data)))
  XA <- if (is.null(fit$A_genpca)) Xp else as.matrix(Xp %*% fit$A_genpca)
  MX <- if (is.null(fit$M_genpca)) Xp else as.matrix(fit$M_genpca %*% Xp)
  sum(MX * XA)
}

# genpca stores loadings in the column metric, `components(x) = A %*% ov`, so
# `scores %*% t(components)` reconstructs `Xhat %*% A`. Back-projecting with
# the A-orthonormal axes `ov` gives the generalized SVD approximation `ou D ov'`.
.mfa_data_axes <- function(x) {
  ov <- x$ov
  if (is.null(ov) || nrow(ov) != nrow(multivarious::components(x))) {
    return(NULL)
  }
  as.matrix(ov)
}

#' Reconstruct Data from an MFA Fit
#'
#' MFA loadings are expressed in the block-weighted column metric
#' (`components(x) = A %*% ov`), so the generic bi-projector reconstruction
#' `scores %*% t(components)` would return the data multiplied by the MFA block
#' weights. These methods back-project with the metric-orthonormal axes `ov`
#' instead, giving the rank-`length(comp)` generalized SVD approximation of the
#' preprocessed data, and then undo the preprocessing.
#'
#' @param x A fitted `mfa` object.
#' @param new_data Matrix of new observations: all concatenated variables, or
#'   only the variables in `colind`.
#' @param comp Integer vector of components to use.
#' @param rowind Training rows to reconstruct.
#' @param colind Concatenated variables (columns) to reconstruct.
#' @param ... Unused.
#'
#' @return A numeric matrix on the original data scale.
#' @importFrom multivarious reconstruct reconstruct_new
#' @export
reconstruct.mfa <- function(x, comp = 1:multivarious::ncomp(x),
                            rowind = 1:nrow(multivarious::scores(x)),
                            colind = 1:nrow(multivarious::components(x)), ...) {
  ov <- .mfa_data_axes(x)
  if (is.null(ov)) {
    return(NextMethod())
  }
  chk::chk_subset(comp, seq_len(multivarious::ncomp(x)))
  chk::chk_subset(rowind, seq_len(nrow(multivarious::scores(x))))
  chk::chk_subset(colind, seq_len(nrow(ov)))

  rec <- multivarious::scores(x)[rowind, comp, drop = FALSE] %*%
    t(ov[colind, comp, drop = FALSE])
  multivarious::inverse_transform(x$preproc, rec, colind = colind)
}

#' @rdname reconstruct.mfa
#' @export
reconstruct_new.mfa <- function(x, new_data, comp = 1:multivarious::ncomp(x),
                                colind = 1:nrow(multivarious::components(x)), ...) {
  ov <- .mfa_data_axes(x)
  if (is.null(ov)) {
    return(NextMethod())
  }
  if (is.vector(new_data)) {
    new_data <- matrix(new_data, nrow = 1)
  }
  new_data <- as.matrix(new_data)
  chk::chk_subset(comp, seq_len(multivarious::ncomp(x)))
  chk::chk_subset(colind, seq_len(nrow(ov)))

  if (length(colind) == nrow(ov) && all(colind == seq_len(nrow(ov)))) {
    chk::chk_equal(ncol(new_data), nrow(ov))
    scores_new <- multivarious::project(x, new_data)[, comp, drop = FALSE]
    rec <- scores_new %*% t(ov[, comp, drop = FALSE])
    multivarious::inverse_transform(x$preproc, rec)
  } else {
    chk::chk_equal(ncol(new_data), length(colind))
    scores_new <- multivarious::partial_project(x, new_data, colind)[, comp, drop = FALSE]
    rec <- scores_new %*% t(ov[colind, comp, drop = FALSE])
    multivarious::inverse_transform(x$preproc, rec, colind = colind)
  }
}

.mfa_missing_refit_spec <- function(data, preproc, ncomp, normalization, A, M,
                                    missing, ncp_impute, missing_tol,
                                    missing_maxiter, return_imputed, fit_dots) {
  data_refit <- lapply(data, function(x) as.matrix(x))
  .muscal_make_refit_spec(
    data = data_refit,
    fit_fn = function(data) {
      do.call(
        mfa,
        c(
          list(
            data = data,
            preproc = preproc,
            ncomp = ncomp,
            normalization = normalization,
            A = A,
            M = M,
            missing = missing,
            ncp_impute = ncp_impute,
            missing_tol = missing_tol,
            missing_maxiter = missing_maxiter,
            return_imputed = return_imputed
          ),
          fit_dots
        )
      )
    },
    bootstrap_fn = function(data) {
      n <- nrow(data[[1L]])
      idx <- sample.int(n, size = n, replace = TRUE)
      lapply(data, function(block) block[idx, , drop = FALSE])
    },
    permutation_fn = function(data) {
      lapply(data, function(block) block[sample.int(nrow(block)), , drop = FALSE])
    },
    resample_unit = "rows"
  )
}

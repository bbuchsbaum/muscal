#' @import abind
#' @importFrom stats quantile sd
#' @noRd
# A helper function to summarize a matrix
summarize_matrix <- function(matrices, .f) {
  stacked_mat <- abind::abind(matrices, along = 3)
  unname(apply(stacked_mat, c(1, 2), .f))
}


#' Summarize Bootstrap Resampling
#'
#' This function summarizes the results of bootstrap resampling by computing the mean, standard deviation,
#' and percentiles for each list of matrices in the input tibble.
#'
#' @param boot_ret A tibble containing the results of bootstrap resampling. It should contain two columns:
#'        - boot_i: A list column where each list element is a matrix.
#'        - boot_j: A list column where each list element is a matrix.
#' @param alpha The percentile level for the computation of the lower and upper percentiles (default is 0.05).
#' @return A list containing the mean, standard deviation, upper and lower percentiles for each set of matrices.
#' @importFrom stats quantile sd
#' @noRd
summarize_boot <- function(boot_ret, alpha=.05){
  # stack once per type to avoid repeated abind allocations
  boot_i_stack <- abind::abind(boot_ret$boot_i, along = 3)
  boot_j_stack <- abind::abind(boot_ret$boot_j, along = 3)
  
  boot_i_mean  <- unname(apply(boot_i_stack, c(1, 2), mean))
  boot_i_sd    <- unname(apply(boot_i_stack, c(1, 2), sd))
  boot_i_upper <- unname(apply(boot_i_stack, c(1, 2), quantile, probs = 1 - alpha, type = 1))
  boot_i_lower <- unname(apply(boot_i_stack, c(1, 2), quantile, probs = alpha, type = 1))
  
  boot_j_mean  <- unname(apply(boot_j_stack, c(1, 2), mean))
  boot_j_sd    <- unname(apply(boot_j_stack, c(1, 2), sd))
  boot_j_upper <- unname(apply(boot_j_stack, c(1, 2), quantile, probs = 1 - alpha, type = 1))
  boot_j_lower <- unname(apply(boot_j_stack, c(1, 2), quantile, probs = alpha, type = 1))
  
  list(
    boot_scores_mean  = boot_i_mean,
    boot_scores_sd    = boot_i_sd,
    boot_scores_upper = boot_i_upper,
    boot_scores_lower = boot_i_lower,
    boot_lds_mean     = boot_j_mean,
    boot_lds_sd       = boot_j_sd,
    boot_lds_upper    = boot_j_upper,
    boot_lds_lower    = boot_j_lower
  )
}


#' Bootstrap Resampling for bada Multivariate Models
#'
#' Perform bootstrap resampling on a multivariate bada model to estimate the variability
#' of components and scores.
#'
#' @param x A fitted bada model object that has been fit to a training dataset.
#' @param nboot An integer specifying the number of bootstrap resamples to perform (default is 500).
#' @param ... Additional arguments:
#'   \describe{
#'     \item{data}{The dataset on which the bootstrap resampling should be performed (required).}
#'     \item{alpha}{The percentile level for the computation of the lower and upper percentiles (default is 0.05).}
#'     \item{verbose}{Logical; print draw indices (default FALSE).}
#'     \item{seed}{Optional integer seed passed to furrr.}
#'     \item{resample_indices}{Optional list of nboot subject-index vectors, each of length the number of subjects.}
#'     \item{keep_replicates}{Retain aligned loading/score matrices and draw metadata (default FALSE).}
#'   }
#' @details Whole subjects are sampled with replacement, with unique replicate subject
#' identities for duplicate draws. Each draw refits the point estimator and freshly
#' fits the original subject preprocessing recipe. This assumes a deterministic,
#' row-preserving preprocessing recipe and supports `rescomp = 0` only.
#' Permutations and signs are aligned to the fitted loading basis. Tied retained
#' reference roots have NA coordinate summaries and subspace diagnostics instead.
#' This experimental percentile summary does not establish coverage or subject
#' exchangeability. Failed/rank-deficient draws are reported and excluded; selective
#' failures may bias inference. At least two valid draws are required.
#' `alpha` is a one-sided tail probability: the default 0.05 gives central 90%
#' intervals; use 0.025 for central 95% intervals. Additional fields report attempted
#' and valid counts, failures, resampling indices, and alignment diagnostics.
#' The returned list contains eight summary matrices:
#'    * `boot_scores_mean`: A matrix which is the mean of all bootstrapped scores matrices.
#'    * `boot_scores_sd`: A matrix which is the standard deviation of all bootstrapped scores matrices.
#'    * `boot_scores_upper`: A matrix which is the upper alpha percentile of all bootstrapped scores matrices.
#'    * `boot_scores_lower`: A matrix which is the lower alpha percentile of all bootstrapped scores matrices.
#'    * `boot_lds_mean`: A matrix which is the mean of all bootstrapped loadings matrices.
#'    * `boot_lds_sd`: A matrix which is the standard deviation of all bootstrapped loadings matrices.
#'    * `boot_lds_upper`: A matrix which is the upper alpha percentile of all bootstrapped loadings matrices.
#'    * `boot_lds_lower`: A matrix which is the lower alpha percentile of all bootstrapped loadings matrices.
#' The dimensions of each matrix in the list correspond to the dimensions of the respective matrices in the input data.
#' @section Bootstrap controls:
#' `bootstrap(x, nboot, data=...)` requires the original fitted multidesign data.
#' Optional `seed` controls reproducibility, `verbose=TRUE` prints draw indices,
#' `keep_replicates=TRUE` retains aligned loading/score matrices, and
#' `resample_indices` supplies a list of nboot subject-index vectors, each of
#' length the original subject count. `alpha` is the one-sided percentile tail.
#' @importFrom furrr future_map furrr_options
#' @importFrom dplyr tibble bind_rows
#' @importFrom purrr map
#' @importFrom multivarious bootstrap
#' @importFrom rlang %||%
#' @importFrom stats sd
#' @rdname bada
#' @export
#' @method bootstrap bada
bootstrap.bada <- function(x, nboot = 500, ...) {
  dots <- list(...)
  data <- dots$data
  alpha <- dots$alpha %||% 0.05
  if (is.null(data)) stop("'data' argument is required for bootstrap.bada")
  if (!is.numeric(nboot) || length(nboot) != 1L || !is.finite(nboot) ||
      nboot < 2 || nboot != floor(nboot)) stop("`nboot` must be an integer >= 2.", call. = FALSE)
  if (!is.numeric(alpha) || length(alpha) != 1L || !is.finite(alpha) ||
      alpha <= 0 || alpha >= 0.5) stop("`alpha` must be strictly between 0 and 0.5.", call. = FALSE)
  if ((x$rescomp %||% 0L) != 0L) stop("BaDA bootstrap supports `rescomp = 0` only.", call. = FALSE)
  y_name <- rlang::as_name(rlang::quo_get_expr(x$y_var))
  subject_name <- rlang::as_name(rlang::quo_get_expr(x$subject_var))
  if (!inherits(data, "multidesign") || !is.numeric(data$x) || any(!is.finite(data$x)) ||
      !identical(as.character(data$design[[y_name]]), as.character(x$labels)) ||
      !identical(as.character(data$design[[subject_name]]), as.character(x$subjects)) ||
      ncol(data$x) != nrow(x$v) ||
      (!is.null(rownames(x$v)) && !identical(colnames(data$x), rownames(x$v)))) {
    stop("`data` must match the fitted BaDA row identities and feature schema and be finite.", call. = FALSE)
  }
  blocks <- .bada_subject_barycenter_blocks(x, data)
  reference <- .bada_average_barycenter_blocks(blocks)
  if (!isTRUE(all.equal(reference, x$barycenters, check.attributes = FALSE))) {
    stop("`data` does not reproduce the fitted BaDA barycenters.", call. = FALSE)
  }
  K <- ncol(x$v)
  reference_d <- svd(reference, nu = 0L, nv = 0L)$d
  full_rank <- function(B) {
    d <- svd(B, nu = 0L, nv = 0L)$d
    length(d) >= K && d[1L] > reference_d[1L] * 1e-10 && d[K] > d[1L] * 1e-8
  }
  if (!full_rank(reference)) stop("Fitted BaDA barycenters lack the requested rank.", call. = FALSE)
  sdat <- split(data, x$subjects)
  S <- length(sdat)
  schedule <- dots$resample_indices
  if (!is.null(schedule) && (!is.list(schedule) || length(schedule) != nboot ||
      any(!vapply(schedule, function(idx) is.numeric(idx) && length(idx) == S &&
        all(is.finite(idx)) && all(idx == floor(idx)) && all(idx >= 1 & idx <= S), logical(1))))) {
    stop("`resample_indices` must contain nboot vectors of S valid subject indices.", call. = FALSE)
  }
  draws <- furrr::future_map(seq_len(nboot), function(i) {
    idx <- if (is.null(schedule)) sample.int(S, S, replace = TRUE) else schedule[[i]]
    if (isTRUE(dots$verbose)) message(i)
    tryCatch({
      # Reject cancelled/rank-deficient means before the point estimator's jitter fallback.
      if (!full_rank(.bada_average_barycenter_blocks(blocks[idx]))) {
        stop("Resampled class barycenters have insufficient rank.")
      }
      selected <- sdat[idx]
      designs <- lapply(seq_along(selected), function(j) {
        design <- selected[[j]]$design
        design[[subject_name]] <- sprintf(".boot_%08d", j)
        design
      })
      boot_data <- multidesign::multidesign(do.call(rbind, lapply(selected, `[[`, "x")),
                                           do.call(rbind, designs))
      fit <- rlang::inject(bada(boot_data, y = !!rlang::sym(y_name),
        subject = !!rlang::sym(subject_name), preproc = x$proclist[[1L]],
        ncomp = K, rescomp = 0))
      if (!full_rank(fit$barycenters)) stop("Refitted class barycenters have insufficient rank.")
      if (!identical(dim(fit$v), dim(x$v)) || any(!is.finite(fit$v)) ||
          !identical(rownames(fit$fscores), rownames(x$fscores)) ||
          any(!is.finite(fit$fscores))) stop("Refitted BaDA target schema is incompatible.")
      list(index = i, indices = idx, loadings = fit$v, scores = fit$fscores,
           n_subjects = length(unique(fit$subjects)))
    }, error = function(e) list(index = i, indices = idx, stage = "refit", message = conditionMessage(e)))
  }, .options = furrr::furrr_options(seed = dots$seed %||% TRUE))
  failed <- vapply(draws, function(draw) !is.null(draw$message), logical(1))
  valid <- draws[!failed]
  if (length(valid) < 2L) stop(sprintf("Only %d valid BaDA bootstrap refits; %d failed. At least two are required.",
                                     length(valid), sum(failed)), call. = FALSE)
  reliability <- loading_reliability(list(sdev = reference_d[seq_len(K)]),
    V_ref = x$v, boot_loadings = lapply(valid, `[[`, "loadings"), alpha = 2 * alpha)
  aligned <- lapply(seq_along(valid), function(i) {
    draw <- valid[[i]]
    p <- reliability$alignment$permutations[[i]]
    signs <- reliability$alignment$signs[[i]]
    draw$loadings <- sweep(draw$loadings[, p, drop = FALSE], 2L, signs, "*")
    draw$scores <- sweep(draw$scores[, p, drop = FALSE], 2L, signs, "*")
    draw
  })
  ret <- summarize_boot(list(boot_i = lapply(aligned, `[[`, "scores"),
                             boot_j = lapply(aligned, `[[`, "loadings")), alpha)
  for (name in names(ret)) ret[[name]][, reliability$alignment$ambiguous_components] <- NA_real_
  ret$nboot <- nboot
  ret$n_success <- length(valid)
  ret$failures <- draws[failed]
  ret$resample_indices <- lapply(draws, `[[`, "indices")
  ret$alignment <- reliability$alignment
  ret$target <- "Whole-subject BaDA refits with fresh subject preprocessing; conditional on retained reference alignment."
  if (isTRUE(dots$keep_replicates)) ret$replicates <- aligned
  class(ret) <- c("bootstrap_bada_result", "list")
  ret
}

#' @import chk
#' @param resdim pca dimensionality for residual analysis (only relevant if `rescomp` > 0)
#' @param rescomp number of final residual components (default = 0, no residual aanalysis)
#' @param ... Additional arguments passed to methods (currently unused).
#' @rdname bada
#' @importFrom multidesign summarize_by multidesign
#' @importFrom multivarious pca sdev pass fit fit_transform transform discriminant_projector scores shape
#' @importFrom stats coef
#' @importFrom rlang enquo quo_get_expr
#' @importFrom dplyr select pull %>%
#' @importFrom chk chk_true
#' @export
bada.multidesign <- function(data, y, subject, preproc=multivarious::center(), ncomp=2,
                             resdim=20, rescomp=0, ...) {
  y_quo <- rlang::enquo(y)
  subject_quo <- rlang::enquo(subject)
  
  labels <- factor(data$design %>% dplyr::select(!!y_quo) %>% dplyr::pull(!!y_quo))
  label_set <- levels(labels)
  
  subjects <- factor(data$design %>% dplyr::select(!!subject_quo) %>% dplyr::pull(!!subject_quo))
  subject_set <- levels(subjects)
  

  
  ## data split by subject
  sdat <- base::split(data, !!subject_quo)


  ## pre-processors, one per subject - fit and transform in one pass
  proclist <- vector("list", length(sdat))
  strata <- seq_along(sdat) %>% purrr::map(function(i) {
    p <- multivarious::fresh(preproc) %>% multivarious::fit(sdat[[i]]$x)
    proclist[[i]] <<- p
    Xout <- multivarious::transform(p, sdat[[i]]$x)
    multidesign::multidesign(Xout, sdat[[i]]$design)
  })

  names(proclist) <- as.character(subject_set)
  
  block_indices <- list()
  ind <- 1
  for (i in 1:length(strata)) {
    block_indices[[i]] <- seq(ind, ind+ncol(strata[[i]]$x)-1)
    ind <- ind + ncol(strata[[i]]$x)
  }
  
  names(block_indices) <- as.character(subject_set)
  
  
  ## group designs
  Dc <- strata %>% purrr::map(function(s) {
    summarize_by(s, !!y_quo)
  })
  
  ## group barycenters
  Xc <-Reduce("+", lapply(Dc, "[[", "x"))/length(Dc)
  Xc[!is.finite(Xc)] <- 0
  if (max(abs(Xc)) < 1e-12) {
    Xc <- Xc + diag(1e-8, nrow(Xc))
  }
  
  ## requires consistent ordering. Better to extract from design
  row.names(Xc) <- label_set
  ncomp <- min(ncomp, nrow(Xc))

  ## group pca
  pca_group <- multivarious::pca(Xc, ncomp=ncomp, preproc=multivarious::pass(), method = "base")
  
  ## residual analysis
  if (rescomp > 0) {
    chk::chk_true(resdim > 0)
    chk::chk_true(rescomp <= resdim)
    residual_strata <- strata %>% purrr::map(function(s) {
      levs <- s$design %>% dplyr::pull(!!y_quo)
      s$x <- s$x - Xc[levs,,drop=FALSE]
      s
    })
  
    Xresid <- do.call(rbind, residual_strata %>% purrr::map( ~ .x$x))
  
    pca_resid <- pca(Xresid, ncomp=resdim, method=.muscal_svd_method(Xresid, ncomp = resdim))
    Xpca_resid <- scores(pca_resid)
  
    Sw <- within_class_scatter(Xpca_resid, interaction(subjects, labels))
    Sb <- between_class_scatter(Xpca_resid, interaction(subjects, labels), 
                              colMeans(Xpca_resid))
  
    eigout <- eigen(solve(Sw, Sb))
    resid_vecs <- eigout$vectors[, seq_len(rescomp), drop = FALSE]
    resid_v <- pca_resid$v %*% resid_vecs
    v_pre <- cbind(pca_group$v, resid_v)
    v <- qr.Q(qr(v_pre))
  } else {
    resdim <- 0
    rescomp <- 0
    v <- pca_group$v
  }
  
  ## compute projections, one subject at a time.
  s <- do.call(rbind, strata %>% purrr::map(function(s) {
    s$x %*% v
  }))
  
  proc <- concat_pre_processors(proclist, block_indices)
  
  discriminant_projector(v=v,
                         s=s,
                         fscores=Xc %*% v,
                         sdev=apply(s, 2, stats::sd),
                         preproc = proc,
                         proclist = proclist,
                         labels=labels,
                         label_set=label_set,
                         resdim=resdim,
                         rescomp=rescomp,
                         subjects=subjects,
                         barycenters=Xc,
                         block_indices=block_indices,
                         subject_var=subject_quo,
                         y_var=y_quo,
                         classes=c("bada", "projector"))
 
}


#' Reprocess New Data for a Barycentric Discriminant Analysis (BaDA) Model
#'
#' This function transforms new data using the preprocessing pipeline from a fitted BaDA model.
#' It handles different preprocessing scenarios based on the provided parameters.
#'
#' @param x A fitted BaDA model object.
#' @param new_data A numeric matrix of new data to be processed.
#' @param colind An optional integer vector specifying column indices to use within blocks.
#'   If NULL and block is also NULL, all blocks are used. If NULL but block is provided,
#'   all columns in the specified block are used.
#' @param ... Additional arguments:
#'   \describe{
#'     \item{block}{An optional character string specifying which block's preprocessing to apply.
#'       If NULL and colind is also NULL, preprocessing is averaged across all blocks.
#'       If provided, only the specified block's preprocessing is applied.}
#'   }
#'
#' @details
#' The function handles three scenarios:
#' 1. When both colind and block are NULL: The function applies preprocessing for each block
#'    and averages the results.
#' 2. When block is provided: The function applies preprocessing specific to the named block.
#'    If colind is also provided, it acts as a relative subset within the block.
#' 3. When only colind is provided: The function applies preprocessing for each block
#'    using the specified column indices and averages the results.
#'
#' @return A preprocessed numeric matrix with the same number of rows as the input data.
#'
#' @rdname bada
#' @export
reprocess.bada <- function(x, new_data, colind = NULL, ...) {
  dots <- list(...)
  block <- dots$block

  if (is.null(colind) && is.null(block)) {
    ## how to pre-process when you don't know the subject?
    ## we pre-process every way and average.
    chk::chk_equal(ncol(new_data), shape(x)[1])
    
    ## pre-process every which way...
    Reduce("+", lapply(seq_along(x$block_indices), function(i) {
      multivarious::transform(x$preproc, new_data, colind=x$block_indices[[i]])
    }))/length(x$block_indices)
    
  } else if (!is.null(block)) {
    ## pre-process along one block
    chk::chk_character(block)
    chk::chk_equal(ncol(new_data), shape(x)[1])
    sind <- x$block_indices[[block]]
    if (!is.null(colind)) {
      ## relative subset using colind
      chk::chk_true(max(colind) <= length(sind))
      sind <- sind[colind]
    }
    multivarious::transform(x$preproc, new_data, colind=sind)
  } else {
    ## colind not null. pre-process every which way using colind per block
    Reduce("+", lapply(seq_along(x$block_indices), function(i) {
      sind <- x$block_indices[[i]]
      chk::chk_true(max(colind) <= length(sind))
      multivarious::transform(x$preproc, new_data, colind=sind[colind])
    }))/length(x$block_indices)
  }
  
}

#' Project New Data onto a Barycentric Discriminant Analysis (BaDA) Model
#'
#' This function projects new data onto a previously fitted BaDA model,
#' returning the scores of the new data in the space of the original model.
#'
#' @param x A fitted BaDA model object.
#' @param new_data A numeric matrix of new data to be projected.
#' @param ... Additional arguments:
#'   \describe{
#'     \item{block}{An optional character string specifying which block's preprocessing
#'       to apply before projection. If not provided, the generic method is used.}
#'   }
#'
#' @details
#' When a specific block is provided, the function first reprocesses the data using
#' that block's preprocessing pipeline, then projects it onto the model space.
#' If no block is specified, it delegates to the default project method.
#'
#' @return A numeric matrix of projected scores.
#'
#' @rdname bada
#' @export
project.bada <- function(x, new_data, ...) {
  dots <- list(...)
  block <- dots$block

  if (is.null(block)) {
    NextMethod()
  } else {
    #Xp <- multivarious::transform(preproc, new_data)
    reprocess(x, new_data, block = block) %*% coef(x)
  }
}

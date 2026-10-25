# ---------------------------------------------------------------------------
# Block-source (out-of-core) entry points built on the shared Gram engine in
# R/gram_engine.R. Each adapter states which options keep the method exact on
# the Gram path, how to call the in-memory fit on compressed blocks, and which
# result fields are feature-indexed.
# ---------------------------------------------------------------------------

# Run a method's Gram path: compress, fit in memory on the compressed blocks,
# map feature-indexed fields back, and point refits/resampling at the source.
#   inner_fn(comp)   in-memory fit on comp$R (plus any in-memory blocks)
#   spec             field spec for .muscal_gram_finish()
#   method_fn        the user-facing method, called again by refits
#   args             evaluated arguments of the user-facing call
#   field            argument holding the block source ("X", or "data")
# The refit data are the inner fit's own refit data with the compressed blocks
# replaced by the source, so their names are the method's argument names.
#   post(fit, inner, comp) optional method-specific fix-ups
.muscal_gram_fit <- function(src, kinds, gram_rank, gram_tol, inner_fn, spec,
                             method_fn, args, field, fit_call,
                             fit_dots = list(), post = NULL) {
  comp <- .muscal_gram_compress(src, kinds, gram_rank = gram_rank,
                                gram_tol = gram_tol %||% 1e-10)
  inner <- inner_fn(comp)
  fit <- .muscal_gram_finish(inner, comp, spec)
  if (!is.null(post)) fit <- post(fit, inner, comp)
  fit$fit_spec$fit_call <- fit_call
  is_data <- identical(field, "data")
  .muscal_gram_refit(
    fit, inner$fit_spec$refit, src,
    data = if (is_data) NULL else inner$fit_spec$refit$data,
    field = if (is_data) NULL else field,
    fit_fn = function(data) {
      call_args <- args
      if (is_data) call_args$data <- data else call_args[names(data)] <- data
      do.call(method_fn, c(call_args, fit_dots), quote = TRUE)
    }
  )
}

.muscal_gram_reject <- function(cond, what) {
  if (isTRUE(cond)) {
    stop(
      sprintf("%s is not supported for block sources (it acts on individual features).", what),
      call. = FALSE
    )
  }
}

# Scores correlated with block `b`'s features: its own score rows when the fit
# keeps them (Z_list), else the reference rows it maps to, else all scores.
.muscal_gram_cor_rows <- function(fit, b) {
  if (!is.null(fit$Z_list[[b]])) return(fit$Z_list[[b]])
  idx <- fit$row_index[[b]]
  if (is.null(idx)) multivarious::scores(fit) else multivarious::scores(fit)[idx, , drop = FALSE]
}

# Fields shared by the mcca family.
.mcca_gram_spec <- function() {
  list(
    concat_lin = "v",
    list_lin = "canonical_weights",
    concat_cor = c("cor_loadings", "scaled_loadings"),
    list_cor = "scaled_loadings_by_block",
    cor_scores = .muscal_gram_cor_rows
  )
}

#' @rdname mfa
#' @details
#' The `mfa.muscal_block_source` method fits MFA out of core from a
#' [block_source()]: each block is read once, chunk by chunk, compressed to its
#' row Gram matrix, fitted in memory, and read once more to map loadings and
#' correlation loadings back to the original features. Results match the
#' in-memory fit (up to the sign of each component). Only column-local
#' preprocessing is supported, missing values are not, and a custom
#' per-feature column metric `A` is rejected because it acts on individual
#' features.
#' @param gram_rank Block sources only: optional maximum number of Gram
#'   eigenpairs kept per block (one value or one per block); `NULL` keeps all,
#'   which is exact. The fraction of each block's variance left out is returned
#'   as `gram_discarded`.
#' @param gram_tol Block sources only: relative eigenvalue tolerance for
#'   dropping numerically null directions of each Gram matrix.
#' @export
mfa.muscal_block_source <- function(data, preproc = center(), ncomp = 2,
                                    normalization = c("MFA", "RV", "None", "Frob", "custom"),
                                    A = NULL, M = NULL,
                                    gram_rank = NULL, gram_tol = 1e-10, ...) {
  fit_call <- match.call(expand.dots = FALSE)
  args <- .muscal_gram_args()
  normalization <- match.arg(normalization)
  .muscal_gram_reject(!is.null(A), "A custom column metric `A`")
  src <- data
  .muscal_gram_fit(
    src = src,
    kinds = .muscal_gram_kinds(preproc, length(src$names)),
    gram_rank = gram_rank, gram_tol = gram_tol,
    inner_fn = function(comp) {
      mfa(comp$R, preproc = multivarious::pass(), ncomp = ncomp,
          normalization = normalization, M = M, ...)
    },
    spec = list(
      concat_lin = c("v", "ov"),
      concat_cor = "cor_loadings",
      cor_scores = function(fit, b) multivarious::scores(fit)
    ),
    method_fn = mfa, args = args, field = "data",
    fit_call = fit_call, fit_dots = list(...),
    post = function(fit, inner, comp) {
      # The column metric is constant within each block.
      a_inner <- Matrix::diag(inner$A_genpca)
      a_block <- vapply(inner$block_indices, function(idx) a_inner[[idx[[1L]]]], numeric(1L),
                        USE.NAMES = FALSE)
      fit$A_genpca <- Matrix::Diagonal(x = rep(a_block, lengths(fit$block_indices)))
      fit
    }
  )
}

#' @rdname mcca
#' @details
#' The `mcca.muscal_block_source` method fits MCCA out of core from a
#' [block_source()]. MCCA already works through each block's row Gram matrix,
#' so the result matches the in-memory fit (up to component signs); each block
#' is read twice. Only column-local preprocessing is supported.
#' @param gram_rank,gram_tol Block sources only; see [mfa()].
#' @export
mcca.muscal_block_source <- function(data, preproc = multivarious::center(), ncomp = 2,
                                     ridge = 1e-6, block_weights = NULL, use_future = FALSE,
                                     gram_rank = NULL, gram_tol = 1e-10, ...) {
  fit_call <- match.call(expand.dots = FALSE)
  args <- .muscal_gram_args()
  src <- data
  .muscal_gram_fit(
    src = src,
    kinds = .muscal_gram_kinds(preproc, length(src$names)),
    gram_rank = gram_rank, gram_tol = gram_tol,
    inner_fn = function(comp) {
      mcca(comp$R, preproc = multivarious::pass(), ncomp = ncomp, ridge = ridge,
           block_weights = block_weights, use_future = use_future, ...)
    },
    spec = .mcca_gram_spec(),
    method_fn = mcca, args = args, field = "data",
    fit_call = fit_call, fit_dots = list(...)
  )
}

.amcca_fit_gram <- function(args, fit_call, fit_dots) {
  src <- args$X
  .muscal_gram_fit(
    src = src,
    kinds = .muscal_gram_kinds(.muscal_gram_effective(args, "preproc"), length(src$names)),
    gram_rank = args$gram_rank, gram_tol = args$gram_tol,
    inner_fn = function(comp) {
      inner_args <- args
      inner_args$X <- comp$R
      inner_args$preproc <- multivarious::pass()
      inner_args$gram_rank <- NULL
      inner_args$gram_tol <- NULL
      do.call(aligned_mcca, c(inner_args, fit_dots), quote = TRUE)
    },
    spec = .mcca_gram_spec(),
    method_fn = aligned_mcca, args = args,
    field = "X", fit_call = fit_call, fit_dots = fit_dots
  )
}

# Methods with an in-memory reference block `Y` followed by source blocks `X`
# (anchored_mcca, anchored_mfa/linked_mfa, graph_anchored_mfa and the coupled
# variant). `spec` adds method-specific fields to the finish spec.
.muscal_anchored_gram <- function(method_fn, args, fit_call, fit_dots, spec) {
  src <- args$X
  n_src <- length(src$names)
  Y <- as.matrix(args$Y)
  preproc <- .muscal_gram_effective(args, "preproc")
  inner_preproc <- .muscal_gram_inner_preproc(preproc, n_src, skip = 1L)
  y_proc <- multivarious::fit(multivarious::fresh(inner_preproc[[1L]]), Y)
  .muscal_gram_fit(
    src = src,
    kinds = .muscal_gram_kinds(preproc, n_src, skip = 1L),
    gram_rank = args$gram_rank, gram_tol = args$gram_tol,
    inner_fn = function(comp) {
      inner_args <- args
      inner_args$X <- comp$R
      inner_args$preproc <- inner_preproc
      inner_args$gram_rank <- NULL
      inner_args$gram_tol <- NULL
      do.call(method_fn, c(inner_args, fit_dots), quote = TRUE)
    },
    spec = c(spec, list(cor_scores = .muscal_gram_cor_rows, keep_preproc = list(Y = y_proc))),
    method_fn = method_fn, args = args, field = "X",
    fit_call = fit_call, fit_dots = fit_dots
  )
}

.anmcca_fit_gram <- function(args, fit_call, fit_dots) {
  .muscal_anchored_gram(anchored_mcca, args, fit_call, fit_dots, .mcca_gram_spec())
}

# Fields shared by the ridge-ALS MFA variants (V_list loadings per block).
.muscal_als_gram_spec <- function() {
  list(
    concat_lin = "v",
    list_lin = "V_list",
    concat_cor = "cor_loadings",
    block_preproc = "block_preproc"
  )
}

.anmfa_fit_gram <- function(args, fit_call, fit_dots) {
  .muscal_gram_reject(!is.null(args$feature_groups) && isTRUE(args$feature_lambda > 0),
                      "The feature prior (feature_groups with feature_lambda > 0)")
  .muscal_anchored_gram(anchored_mfa, args, fit_call, fit_dots, .muscal_als_gram_spec())
}

.gamfa_fit_gram <- function(method_fn, args, fit_call, fit_dots) {
  .muscal_gram_reject(!is.null(args$feature_graph) && isTRUE(args$graph_lambda > 0),
                      "A feature graph (feature_graph with graph_lambda > 0)")
  .muscal_anchored_gram(method_fn, args, fit_call, fit_dots, .muscal_als_gram_spec())
}

.ramfa_fit_gram <- function(args, fit_call, fit_dots) {
  .muscal_gram_reject(!is.null(args$feature_groups) && isTRUE(args$feature_lambda > 0),
                      "The feature prior (feature_groups with feature_lambda > 0)")
  .muscal_gram_reject(!is.null(args$feature_graph) && isTRUE(args$graph_lambda > 0),
                      "A feature graph (feature_graph with graph_lambda > 0)")
  src <- args$X
  .muscal_gram_fit(
    src = src,
    kinds = .muscal_gram_kinds(.muscal_gram_effective(args, "preproc"), length(src$names)),
    gram_rank = args$gram_rank, gram_tol = args$gram_tol,
    inner_fn = function(comp) {
      inner_args <- args
      inner_args$X <- comp$R
      inner_args$preproc <- multivarious::pass()
      inner_args$gram_rank <- NULL
      inner_args$gram_tol <- NULL
      do.call(response_aligned_mfa, c(inner_args, fit_dots), quote = TRUE)
    },
    spec = c(.muscal_als_gram_spec(), list(cor_scores = .muscal_gram_cor_rows)),
    method_fn = response_aligned_mfa, args = args, field = "X",
    fit_call = fit_call, fit_dots = fit_dots
  )
}

# aligned_interbattery: blocks in a domain share one loading matrix, so the
# compression unit is the row-stacked domain. Each domain becomes one block of
# a virtual stacked source named "<side>::<loading id>" (the keys the fit uses
# for its loadings); its compressed rows are split back into the member
# blocks, and the inner fit is told the same grouping through block_info.
.aib_fit_gram <- function(args, fit_call, fit_dots) {
  if (!is.null(args$gram_rank)) {
    stop("`gram_rank` is not supported for aligned_interbattery(): with shared domain ",
         "loadings, truncation would make per-block fit statistics inexact.", call. = FALSE)
  }
  units <- list()
  plans <- list()
  for (side in c("x", "y")) {
    data <- args[[toupper(side)]]
    if (!inherits(data, "muscal_block_source")) next
    info <- .aib_validate_flat_block_info(args[[paste0(side, "_block_info")]], data$names, side = side)
    dummy <- lapply(seq_along(data$names), function(k) {
      matrix(0, 0L, data$ncol[[k]], dimnames = list(NULL, data$colnames[[k]]))
    })
    names(dummy) <- data$names
    loading_id <- .aib_make_loading_ids(dummy, info, side)
    for (id in unique(loading_id)) {
      members <- match(names(loading_id)[loading_id == id], data$names)
      units[[paste(side, id, sep = "::")]] <- list(side = side, id = id, src = data, members = members)
    }
    plans[[side]] <- list(src = data, info = info, loading_id = loading_id)
  }
  vsrc <- block_source(
    loader = function(g, cols) {
      u <- units[[g]]
      do.call(rbind, lapply(u$members, function(k) .muscal_source_read(u$src, k, cols)))
    },
    nrow = vapply(units, function(u) sum(u$src$nrow[u$members]), integer(1L)),
    ncol = vapply(units, function(u) u$src$ncol[[u$members[[1L]]]], integer(1L)),
    names = names(units),
    colnames = lapply(units, function(u) u$src$colnames[[u$members[[1L]]]]),
    chunk_size = min(vapply(units, function(u) u$src$chunk_size, integer(1L)))
  )
  effective <- attr(args, "effective")
  kinds <- vapply(units, function(u) {
    .muscal_stream_preproc_kinds(effective[[paste0(u$side, "_preproc")]], 1L)
  }, character(1L))
  comp <- .muscal_gram_compress(vsrc, kinds, gram_tol = args$gram_tol %||% 1e-10)

  inner_args <- args
  inner_args$gram_rank <- NULL
  inner_args$gram_tol <- NULL
  for (side in names(plans)) {
    pl <- plans[[side]]
    blocks <- vector("list", length(pl$src$names))
    names(blocks) <- pl$src$names
    for (key in names(units)[vapply(units, `[[`, character(1L), "side") == side]) {
      u <- units[[key]]
      ends <- cumsum(u$src$nrow[u$members])
      starts <- ends - u$src$nrow[u$members] + 1L
      for (j in seq_along(u$members)) {
        blocks[[u$members[[j]]]] <- comp$R[[key]][starts[[j]]:ends[[j]], , drop = FALSE]
      }
    }
    info <- pl$info
    info$domain <- unname(pl$loading_id[info$block])
    inner_args[[toupper(side)]] <- blocks
    inner_args[[paste0(side, "_block_info")]] <- info
    inner_args[paste0(side, "_preproc")] <- list(NULL)
  }
  inner <- do.call(aligned_interbattery, c(inner_args, fit_dots), quote = TRUE)

  fit <- .muscal_gram_finish(inner, comp, list(concat_lin = "v", rebuild_preproc = FALSE))
  procs <- .muscal_gram_fitted_preprocs(comp)
  for (key in names(units)) {
    u <- units[[key]]
    field <- paste0(u$side, "_loadings")
    fit[[field]][[u$id]] <- fit$v[fit$block_indices[[key]], , drop = FALSE]
    fit[[paste0(u$side, "_loading_preproc")]][[u$id]] <- procs[[key]]
  }
  for (side in names(plans)) {
    pl <- plans[[side]]
    info_field <- paste0(side, "_block_info")
    fit[[info_field]]$domain <- pl$info$domain[match(fit[[info_field]]$block, pl$info$block)]
    rows <- which(fit$block_fit$side == side)
    blk <- fit$block_fit$block[rows]
    fit$block_fit$p[rows] <- pl$src$ncol[match(blk, pl$src$names)]
    fit$block_fit$domain[rows] <- pl$info$domain[match(blk, pl$info$block)]
  }
  fit$preproc <- multivarious::fit(
    multivarious::fresh(multivarious::pass()),
    matrix(0, nrow = 1L, ncol = nrow(fit$v))
  )
  fit$fit_spec$fit_call <- fit_call
  fit
}

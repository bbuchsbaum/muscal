# ---------------------------------------------------------------------------
# Shared Gram-path engine for out-of-core fits
#
# For methods whose loadings lie in each block's row space (ridge/least-squares
# loadings with rotation-invariant penalties, or kernel methods built on
# X_k X_k'), the fit depends on a preprocessed block X_k only through its row
# Gram matrix. Writing the thin SVD X_k = U_k D_k W_k', the engine
#
#   1. reads every block once, chunk by chunk, and compresses it to
#      R_k = U_k D_k (n_k x r_k) with the same rows as X_k;
#   2. lets the method fit in memory on the R_k, which yields identical
#      row-space quantities (scores, partial scores, X_k V_k);
#   3. reads every block once more and maps each feature-indexed output back:
#      a loading matrix V_R in R_k coordinates becomes
#      V_X = W_k V_R = X_k' U_k D_k^{-1} V_R, and feature-score correlations are
#      recomputed from X_k' (T_k - mean(T_k)).
#
# Methods supply a small `spec` naming their feature-indexed fields; see
# .muscal_gram_finish().
# ---------------------------------------------------------------------------

# Pass 1: column-local preprocessing statistics and the compressed blocks.
# Only one n_k x n_k Gram matrix exists at a time.
.muscal_gram_compress <- function(src, kinds, gram_rank = NULL, gram_tol = 1e-10) {
  n_blocks <- length(src$names)
  chk::chk_number(gram_tol)
  chk::chk_gte(gram_tol, 0)
  if (!is.null(gram_rank)) {
    gram_rank <- as.integer(gram_rank)
    if (length(gram_rank) == 1L) gram_rank <- rep(gram_rank, n_blocks)
    if (length(gram_rank) != n_blocks || anyNA(gram_rank) || any(gram_rank < 1L)) {
      stop("`gram_rank` must be NULL or positive integers (one, or one per block).",
           call. = FALSE)
    }
  }
  roots <- vector("list", n_blocks)
  gram <- .muscal_stream_gram(src, kinds, keep_gram = FALSE, on_gram = function(k, Gk) {
    roots[[k]] <<- .muscal_gram_root(Gk, tol = gram_tol, max_rank = gram_rank[k])
  })
  names(roots) <- src$names
  list(
    src = src,
    kinds = kinds,
    prep = gram$prep,
    R = lapply(roots, `[[`, "R"),
    lift = lapply(roots, `[[`, "lift"),
    lead = vapply(roots, `[[`, numeric(1L), "lead"),
    rank = vapply(roots, function(r) ncol(r$R), integer(1L)),
    discarded = vapply(roots, `[[`, numeric(1L), "discarded"),
    dropped_ss = vapply(roots, `[[`, numeric(1L), "dropped_ss")
  )
}

# Pass 2: for each source block k, X_k' (lift_k lin_k) for loadings given in
# R_k coordinates, and X_k' (T_k - mean(T_k)) for feature-score correlations.
# `lin` and `cor` are lists over source blocks (NULL entries allowed).
.muscal_gram_lift <- function(comp, lin, cor) {
  src <- comp$src
  n_blocks <- length(src$names)
  rhs <- lapply(seq_len(n_blocks), function(k) {
    parts <- list()
    if (!is.null(lin[[k]])) parts$lin <- comp$lift[[k]] %*% lin[[k]]
    if (!is.null(cor[[k]])) parts$cor <- sweep(cor[[k]], 2L, colMeans(cor[[k]]))
    parts
  })
  out <- lapply(seq_len(n_blocks), function(k) {
    m <- sum(vapply(rhs[[k]], ncol, integer(1L)))
    matrix(0, src$ncol[[k]], m)
  })
  M_list <- lapply(rhs, function(parts) if (length(parts)) do.call(cbind, parts) else NULL)
  todo <- which(!vapply(M_list, is.null, logical(1L)))
  if (length(todo) > 0L) {
    .muscal_stream_pass(src, comp$kinds, comp$prep, function(k, cols, Xp) {
      out[[k]][cols, ] <<- crossprod(Xp, M_list[[k]])
    }, blocks = todo)
  }
  lapply(seq_len(n_blocks), function(k) {
    res <- list(lin = NULL, cor = NULL)
    at <- 0L
    if (!is.null(rhs[[k]]$lin)) {
      m <- ncol(rhs[[k]]$lin)
      res$lin <- out[[k]][, at + seq_len(m), drop = FALSE]
      at <- at + m
    }
    if (!is.null(rhs[[k]]$cor)) {
      Tc <- rhs[[k]]$cor
      C <- out[[k]][, at + seq_len(ncol(Tc)), drop = FALSE] /
        outer(sqrt(comp$prep[[k]]$col_ss), sqrt(colSums(Tc^2)))
      C[!is.finite(C)] <- 0
      res$cor <- C
    }
    res
  })
}

# Map a fit made on compressed blocks back to the original feature space.
#
# spec fields (all optional):
#   concat_lin  names of matrices whose rows follow fit$block_indices and that
#               are linear in the row space (loadings, e.g. "v")
#   list_lin    names of per-block lists of such matrices, keyed by block name
#   concat_cor  names of concatenated feature-score correlation matrices
#   list_cor    names of per-block lists of correlation matrices
#   cor_scores  function(fit, block) giving the n_k x m scores correlated with
#               block `block`'s features (required if any *_cor field is set)
#   keep_preproc named list of fitted preprocessors for fit blocks that are not
#               source blocks (e.g. an in-memory reference block "Y")
#   block_preproc names of per-block lists of fitted preprocessors to replace
# Blocks are matched by name: fit blocks whose names are not source blocks are
# left untouched.
.muscal_gram_finish <- function(fit, comp, spec) {
  src <- comp$src
  indices_named <- !is.null(names(fit$block_indices))
  fit_blocks <- names(fit$block_indices) %||% fit$names
  if (is.null(fit_blocks) && length(fit$block_indices) == length(src$names)) {
    fit_blocks <- src$names
  }
  names(fit$block_indices) <- fit_blocks
  in_src <- fit_blocks %in% src$names
  src_of <- match(fit_blocks, src$names)

  # Collect R-coordinate matrices per source block, column-bound across fields.
  lin_fields <- c(spec$concat_lin, spec$list_lin)
  cor_fields <- c(spec$concat_cor, spec$list_cor)
  piece <- function(field, b) {
    if (field %in% spec$concat_lin || field %in% spec$concat_cor) {
      M <- fit[[field]]
      if (is.null(M)) return(NULL)
      as.matrix(M)[fit$block_indices[[b]], , drop = FALSE]
    } else {
      L <- fit[[field]]
      if (is.null(L) || is.null(L[[b]])) return(NULL)
      as.matrix(L[[b]])
    }
  }
  n_src <- length(src$names)
  lin <- vector("list", n_src)
  cor <- vector("list", n_src)
  lin_cols <- cor_cols <- vector("list", n_src)
  for (b in fit_blocks[in_src]) {
    k <- match(b, src$names)
    mats <- lapply(lin_fields, piece, b = b)
    lin_cols[[k]] <- vapply(mats, function(M) if (is.null(M)) 0L else ncol(M), integer(1L))
    keep <- !vapply(mats, is.null, logical(1L))
    if (any(keep)) lin[[k]] <- do.call(cbind, mats[keep])
    if (length(cor_fields)) {
      cor[[k]] <- as.matrix(spec$cor_scores(fit, b))
    }
  }
  lifted <- .muscal_gram_lift(comp, lin, cor)

  # Feature-space dimensions per fit block.
  dims <- vapply(seq_along(fit_blocks), function(i) {
    if (in_src[[i]]) src$ncol[[src_of[[i]]]] else length(fit$block_indices[[i]])
  }, integer(1L))
  new_indices <- split(seq_len(sum(dims)), rep(seq_along(dims), dims))
  names(new_indices) <- fit_blocks
  lifted_field <- function(field, b, kind) {
    k <- match(b, src$names)
    if (identical(kind, "cor")) {
      C <- lifted[[k]]$cor
      rownames(C) <- src$colnames[[k]]
      return(C)
    }
    widths <- lin_cols[[k]]
    j <- match(field, lin_fields)
    if (widths[[j]] == 0L) return(NULL)
    start <- sum(widths[seq_len(j - 1L)])
    M <- lifted[[k]]$lin[, start + seq_len(widths[[j]]), drop = FALSE]
    rownames(M) <- src$colnames[[k]]
    M
  }
  rebuild_concat <- function(field, kind) {
    old <- fit[[field]]
    if (is.null(old)) return(NULL)
    old <- as.matrix(old)
    parts <- lapply(seq_along(fit_blocks), function(i) {
      b <- fit_blocks[[i]]
      if (in_src[[i]]) {
        M <- lifted_field(field, b, kind)
      } else {
        M <- old[fit$block_indices[[i]], , drop = FALSE]
      }
      colnames(M) <- colnames(old)
      M
    })
    do.call(rbind, parts)
  }
  rebuild_list <- function(field, kind) {
    old <- fit[[field]]
    if (is.null(old)) return(NULL)
    for (b in intersect(names(old), src$names)) {
      M <- lifted_field(field, b, kind)
      if (!is.null(M)) {
        colnames(M) <- colnames(old[[b]])
        old[[b]] <- M
      }
    }
    old
  }
  for (f in spec$concat_lin) fit[[f]] <- rebuild_concat(f, "lin")
  for (f in spec$list_lin) fit[[f]] <- rebuild_list(f, "lin")
  for (f in spec$concat_cor) fit[[f]] <- rebuild_concat(f, "cor")
  for (f in spec$list_cor) fit[[f]] <- rebuild_list(f, "cor")

  # Preprocessors: streamed ones for source blocks, kept ones otherwise
  # (spec$rebuild_preproc = FALSE leaves fit$preproc to the caller).
  src_procs <- .muscal_gram_fitted_preprocs(comp)
  if (!isFALSE(spec$rebuild_preproc)) {
    procs <- lapply(seq_along(fit_blocks), function(i) {
      if (in_src[[i]]) src_procs[[src_of[[i]]]] else spec$keep_preproc[[fit_blocks[[i]]]]
    })
    if (any(vapply(procs, is.null, logical(1L)))) {
      stop("Internal error: missing fitted preprocessor for a non-source block.", call. = FALSE)
    }
    fit$preproc <- multivarious::concat_pre_processors(procs, unname(new_indices))
  }
  fit$block_indices <- if (indices_named) new_indices else unname(new_indices)
  for (f in spec$block_preproc) {
    if (!is.null(fit[[f]])) {
      for (b in intersect(names(fit[[f]]), src$names)) fit[[f]][[b]] <- src_procs[[b]]
    }
  }

  # Goodness-of-fit tables: compressed blocks have the same rows, and under
  # truncation the residual gains exactly the variance that was dropped.
  if (is.data.frame(fit$block_fit) && "block" %in% names(fit$block_fit)) {
    bf <- fit$block_fit
    rows <- which(bf$block %in% src$names)
    for (r in rows) {
      k <- match(bf$block[[r]], src$names)
      if ("p" %in% names(bf)) bf$p[[r]] <- src$ncol[[k]]
      if (all(c("sse", "tss") %in% names(bf))) {
        bf$sse[[r]] <- bf$sse[[r]] + comp$dropped_ss[[k]]
        bf$tss[[r]] <- comp$prep[[k]]$tss
        if ("r2" %in% names(bf)) {
          bf$r2[[r]] <- if (bf$tss[[r]] > 0) 1 - bf$sse[[r]] / bf$tss[[r]] else NA_real_
        }
      }
    }
    fit$block_fit <- bf
  }

  fit$strategy <- "gram"
  fit$gram_rank <- comp$rank
  fit$gram_discarded <- comp$discarded
  fit
}

# Fitted preprocessors for every source block, from the streamed statistics.
.muscal_gram_fitted_preprocs <- function(comp) {
  src <- comp$src
  out <- lapply(seq_along(src$names), function(k) {
    .muscal_stream_fitted_preproc(comp$kinds[[k]], comp$prep[[k]], src$ncol[[k]])
  })
  names(out) <- src$names
  out
}

# Wrap an in-memory resampler (bootstrap or permutation) so it works on data
# whose blocks live in a block source. The resampler runs on stand-in blocks
# holding row ids, and the rows it picked become a row-subset view of the
# source. `field` names the element of `data` holding the blocks, or NULL when
# `data` is itself the list of blocks.
.muscal_source_resampler <- function(resample_fn, field = "X") {
  force(resample_fn)
  function(data) {
    src <- if (is.null(field)) data else data[[field]]
    proxies <- lapply(seq_along(src$names), function(k) {
      matrix(seq_len(src$nrow[[k]]), ncol = 1L)
    })
    names(proxies) <- src$names
    data_proxy <- data
    if (is.null(field)) data_proxy <- proxies else data_proxy[[field]] <- proxies
    out <- resample_fn(data_proxy)
    picked <- if (is.null(field)) out else out[[field]]
    rows <- lapply(src$names, function(b) as.integer(picked[[b]][, 1L]))
    if (any(lengths(rows) == 0L)) {
      stop("A resampling draw left a block with no rows.", call. = FALSE)
    }
    view <- .muscal_source_subset_rows(src, rows)
    if (is.null(field)) {
      view
    } else {
      out[[field]] <- view
      out
    }
  }
}

# Replace a fit's refit metadata so that refits and resampling read the block
# source instead of the compressed blocks the inner fit saw.
.muscal_gram_refit <- function(fit, inner_spec, src, data, fit_fn, field = "X") {
  if (!is.null(field)) data[[field]] <- src else data <- src
  fit$fit_spec$refit <- .muscal_make_refit_spec(
    data = data,
    fit_fn = fit_fn,
    bootstrap_fn = .muscal_source_resampler(inner_spec$bootstrap_fn, field),
    permutation_fn = .muscal_source_resampler(inner_spec$permutation_fn, field),
    resample_unit = inner_spec$resample_unit
  )
  fit$fit_spec$refit_supported <- TRUE
  fit
}

# Evaluated values of the arguments the caller was actually given (without
# `...`), used to call a method again on compressed blocks or on the source.
# Unsupplied arguments are left out so their defaults (and any missing()
# logic) apply as in the original call.
.muscal_gram_args <- function(env = parent.frame(), fun = sys.function(sys.parent()),
                              call = sys.call(sys.parent())) {
  mc <- match.call(fun, call, expand.dots = FALSE)
  nms <- setdiff(names(mc)[-1L], c("", "..."))
  out <- mget(nms, envir = env)
  # All arguments with their effective (default-filled) values, for adapters
  # that need e.g. the preprocessing actually in force.
  attr(out, "effective") <- mget(setdiff(names(formals(fun)), "..."), envir = env)
  out
}

.muscal_gram_effective <- function(args, name) {
  attr(args, "effective")[[name]]
}

# Preprocessing definitions for source blocks, from a single definition or a
# list (optionally skipping leading in-memory blocks such as "Y").
.muscal_gram_kinds <- function(preproc, n_src, skip = 0L) {
  defs <- if (is.list(preproc) && !inherits(preproc, "pre_processor") &&
              !inherits(preproc, "prepper")) {
    if (skip > 0L) preproc[-seq_len(skip)] else preproc
  } else {
    preproc
  }
  .muscal_stream_preproc_kinds(defs, n_src)
}

# Inner preprocessing for a fit on compressed blocks: compressed blocks are
# already preprocessed, so they get pass(); leading in-memory blocks keep theirs.
.muscal_gram_inner_preproc <- function(preproc, n_src, skip = 0L) {
  passes <- rep(list(multivarious::pass()), n_src)
  if (skip == 0L) return(passes)
  lead <- if (is.list(preproc) && !inherits(preproc, "pre_processor") &&
              !inherits(preproc, "prepper")) {
    preproc[seq_len(skip)]
  } else {
    rep(list(preproc), skip)
  }
  c(lead, passes)
}

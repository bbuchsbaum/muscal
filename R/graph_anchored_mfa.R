#' Graph-Anchored Multiple Factor Analysis
#'
#' @description
#' `graph_anchored_mfa()` extends [anchored_mfa()] to settings where auxiliary
#' blocks do not share aligned columns and borrowing across blocks is induced by
#' a sparse feature graph. Rows of each auxiliary block are linked to a common
#' anchor matrix `Y` through `row_index`, while a graph Laplacian penalty on the
#' concatenated auxiliary loadings encourages similar features across blocks to
#' have similar latent representations.
#'
#' Missing domains are handled by omission: each observed subject-domain pair is
#' treated as one auxiliary block. This supports subjects with only `D1`,
#' subjects with `D1` and `D2`, and mixtures of observed-domain patterns.
#'
#' @details
#' ## Model
#' The fitted model has the form:
#' \deqn{Y \approx S B^\top}
#' \deqn{X_k \approx S[\mathrm{idx}_k,] V_k^\top}
#' where `S` is `N × ncomp`, `B` is `q × ncomp`, and each `V_k` is
#' `p_k × ncomp`.
#'
#' Let `V` denote the row-wise concatenation of the auxiliary loading matrices
#' `V_k`, and let `L` be a feature-graph Laplacian over all auxiliary features.
#' The estimator minimizes the anchored-MFA reconstruction loss plus the graph
#' smoothness term
#' \deqn{\lambda_G \mathrm{tr}(V^\top L V).}
#'
#' When `graph_lambda = 0` (or `feature_graph = NULL`), the method reduces to
#' [anchored_mfa()] up to numerical tolerance.
#'
#' ## Input organization
#' `X` and `row_index` may be supplied either as:
#' * flat lists of observed blocks, or
#' * nested subject/domain lists, e.g. `X[[subject]][[domain]]`.
#'
#' Nested input is flattened internally into one block per observed
#' subject-domain pair. The resulting mapping is recorded in `block_info`.
#'
#' ## Feature graph
#' `feature_graph` may be:
#' * `NULL` (no graph penalty),
#' * `"colnames"` to connect identical auxiliary column names across blocks,
#' * a data frame with columns `block1`, `feature1`, `block2`, `feature2` and
#'   optional `weight`, or
#' * a square sparse/dense matrix interpreted according to `graph_form`.
#'
#' @param Y Numeric matrix/data.frame (`N × q`) serving as the anchored target
#'   space.
#' @param X Auxiliary blocks. Either a flat named list of matrices/data frames,
#'   or a nested list `X[[subject]][[domain]]`.
#' @param row_index A structure mirroring `X`. Each vector maps rows of the
#'   corresponding auxiliary block to rows of `Y`.
#' @param block_info Optional data frame describing flattened blocks. If
#'   supplied, it must have one row per flattened block. Recommended columns are
#'   `block`, `subject`, and `domain`.
#' @param preproc A `multivarious` preprocessing pipeline (a `pre_processor` or
#'   `prepper`) or a list of them. If a list, it must have length
#'   `1 + length(flattened_X)` and will be applied to `c(list(Y), flattened_X)`.
#' @param ncomp Integer number of latent components.
#' @param normalization Block weighting scheme. `"MFA"` uses inverse squared
#'   first singular values; `"None"` uses uniform weights; `"custom"` uses
#'   `alpha`.
#' @param alpha Optional numeric vector of per-block weights. When
#'   `normalization = "custom"`, it must have length `1 + length(flattened_X)`,
#'   with the first weight corresponding to `Y`.
#' @param feature_graph Feature-graph specification; see Details.
#' @param graph_lambda Non-negative scalar controlling graph-penalty strength.
#' @param graph_form Interpretation of `feature_graph` when it is matrix-like, or
#'   the Laplacian construction used for edge-based inputs.
#' @param max_iter Maximum ALS iterations.
#' @param tol Relative convergence tolerance on the objective.
#' @param ridge Non-negative ridge stabilization applied to loading and score
#'   updates.
#' @param verbose Logical; if `TRUE`, prints iteration diagnostics.
#' @param ... Unused (reserved for future extensions).
#'
#' @return An object inheriting from `multivarious::multiblock_biprojector` with
#'   additional classes `graph_anchored_mfa`, `anchored_mfa`, and `linked_mfa`.
#'   The object contains anchored scores in `s`, auxiliary loadings in `V_list`,
#'   anchor loadings in `B`, block metadata in `block_info`, and graph metadata
#'   in `graph_laplacian` and `graph_lambda`.
#'
#' @export
graph_anchored_mfa <- function(Y,
                               X,
                               row_index,
                               block_info = NULL,
                               preproc = multivarious::center(),
                               ncomp = 2,
                               normalization = c("MFA", "None", "custom"),
                               alpha = NULL,
                               feature_graph = NULL,
                               graph_lambda = 0,
                               graph_form = c("laplacian", "adjacency", "normalized_laplacian"),
                               max_iter = 50,
                               tol = 1e-6,
                               ridge = 1e-8,
                               verbose = FALSE,
                               ...) {
  fit_call <- match.call(expand.dots = FALSE)
  fit_dots <- list(...)
  normalization <- match.arg(normalization)
  graph_form <- match.arg(graph_form)

  chk::chk_true(is.matrix(Y) || is.data.frame(Y))
  Y <- as.matrix(Y)
  chk::chk_list(X)
  chk::chk_list(row_index)
  ncomp <- as.integer(ncomp)
  chk::chk_integer(ncomp)
  chk::chk_gte(ncomp, 1)
  chk::chk_numeric(graph_lambda)
  chk::chk_gte(graph_lambda, 0)
  max_iter <- as.integer(max_iter)
  chk::chk_integer(max_iter)
  chk::chk_gte(max_iter, 1)
  chk::chk_numeric(tol)
  chk::chk_gte(tol, 0)
  chk::chk_numeric(ridge)
  chk::chk_gte(ridge, 0)
  chk::chk_flag(verbose)

  flat <- .gamfa_flatten_blocks(X, row_index, block_info)
  X <- flat$X
  row_index <- flat$row_index
  block_info <- flat$block_info
  data_refit <- list(
    Y = as.matrix(Y),
    X = lapply(X, function(x) as.matrix(x)),
    row_index = lapply(row_index, as.integer),
    block_info = block_info
  )

  if (length(X) < 1L) {
    stop("At least one auxiliary block is required.", call. = FALSE)
  }

  N <- nrow(Y)
  for (k in seq_along(X)) {
    idx <- row_index[[k]]
    chk::chk_true(is.numeric(idx) || is.integer(idx))
    idx <- as.integer(idx)
    chk::chk_equal(length(idx), nrow(X[[k]]))
    if (anyNA(idx)) stop("row_index cannot contain NA.", call. = FALSE)
    if (any(idx < 1L) || any(idx > N)) {
      stop("row_index values must be in 1..nrow(Y).", call. = FALSE)
    }
    row_index[[k]] <- idx
  }

  blocks <- c(list(Y = Y), X)
  if (is.list(preproc) && !inherits(preproc, "pre_processor") && !inherits(preproc, "prepper")) {
    chk::chk_equal(length(preproc), length(blocks))
  }
  prep_res <- prepare_block_preprocessors(blocks, preproc, check_consistent_ncol = FALSE)
  Yp <- prep_res$Xp[[1]]
  Xp <- prep_res$Xp[-1]
  proclist <- prep_res$proclist
  N <- nrow(Yp)

  alpha_blocks <- if (normalization == "custom") {
    if (is.null(alpha)) {
      stop("When normalization = 'custom', 'alpha' must be provided.", call. = FALSE)
    }
    alpha <- as.numeric(alpha)
    chk::chk_equal(length(alpha), length(blocks))
    if (any(!is.finite(alpha)) || any(alpha < 0)) {
      stop("'alpha' must be finite and non-negative.", call. = FALSE)
    }
    alpha
  } else if (normalization == "None") {
    rep(1, length(blocks))
  } else {
    first_sv <- function(mat) {
      mdim <- min(dim(mat))
      method <- if (mdim < 3) "base" else "svds"
      tryCatch(
        multivarious::svd_wrapper(mat, ncomp = 1, method = method)$sdev[1],
        error = function(e) multivarious::svd_wrapper(mat, ncomp = 1, method = "base")$sdev[1]
      )
    }
    sapply(prep_res$Xp, function(M) 1 / (first_sv(M)^2 + ridge))
  }
  names(alpha_blocks) <- c("Y", names(X))

  min_p <- min(c(ncol(Yp), vapply(Xp, ncol, integer(1))))
  K <- min(as.integer(ncomp), nrow(Yp), min_p)
  if (K < 1L) stop("ncomp must be at least 1 and <= min(nrow(Y), ncol(Y), ncol(X_k)).", call. = FALSE)

  graph <- .gamfa_parse_feature_graph(
    feature_graph = feature_graph,
    X_list = Xp,
    graph_lambda = graph_lambda,
    graph_form = graph_form
  )

  init <- .lmfa_init_from_Y(Yp, K, ridge = ridge)
  S <- init$S
  B <- init$B

  V_list <- lapply(seq_along(Xp), function(k) {
    Sk <- S[row_index[[k]], , drop = FALSE]
    .lmfa_update_loadings(Sk, Xp[[k]], ridge = ridge)$V
  })
  names(V_list) <- names(Xp)

  objective_trace <- numeric(0)
  prev_obj <- Inf

  for (iter in seq_len(max_iter)) {
    B <- .lmfa_update_loadings(S, Yp, ridge = ridge)$V

    if (isTRUE(graph$enabled) && isTRUE(graph_lambda > 0)) {
      S_list <- lapply(seq_along(Xp), function(k) S[row_index[[k]], , drop = FALSE])
      V_list <- .gamfa_update_V_graph(
        S_list = S_list,
        X_list = Xp,
        alpha_blocks = alpha_blocks[-1L],
        graph_laplacian = graph$L,
        graph_lambda = graph_lambda,
        ridge = ridge
      )
      names(V_list) <- names(Xp)
    } else {
      for (k in seq_along(Xp)) {
        Sk <- S[row_index[[k]], , drop = FALSE]
        V_list[[k]] <- .lmfa_update_V_block(
          Sk = Sk,
          Xk = Xp[[k]],
          alpha_block = alpha_blocks[k + 1L],
          fg_block = rep.int(NA_character_, ncol(Xp[[k]])),
          weights_block = rep.int(1, ncol(Xp[[k]])),
          centers = NULL,
          feature_lambda = 0,
          ridge = ridge
        )
      }
      names(V_list) <- names(Xp)
    }

    S <- .lmfa_update_scores(
      Y = Yp,
      B = B,
      X_list = Xp,
      V_list = V_list,
      row_index = row_index,
      alpha_y = alpha_blocks[1],
      alpha_blocks = alpha_blocks[-1],
      ridge = ridge
    )

    ortho <- .lmfa_orthonormalize_scores(S)
    S <- ortho$S
    rot <- ortho$R
    B <- B %*% t(rot)
    V_list <- lapply(V_list, function(V) V %*% t(rot))

    obj <- .gamfa_objective(
      Y = Yp,
      S = S,
      B = B,
      X_list = Xp,
      V_list = V_list,
      row_index = row_index,
      alpha_y = alpha_blocks[1],
      alpha_blocks = alpha_blocks[-1],
      graph_laplacian = graph$L,
      graph_lambda = graph_lambda
    )
    objective_trace <- c(objective_trace, obj)

    rel_change <- abs(obj - prev_obj) / (abs(prev_obj) + ridge)
    if (isTRUE(verbose)) {
      message(
        sprintf(
          "graph_anchored_mfa iter %d: obj=%.6g, rel_change=%.3g",
          iter, obj, rel_change
        )
      )
    }
    if (is.finite(prev_obj) && rel_change < tol) break
    prev_obj <- obj
  }

  block_indices <- list()
  current <- 1L
  block_indices[[1]] <- current:(current + ncol(Yp) - 1L)
  current <- current + ncol(Yp)
  for (k in seq_along(Xp)) {
    block_indices[[k + 1L]] <- current:(current + ncol(Xp[[k]]) - 1L)
    current <- current + ncol(Xp[[k]])
  }
  names(block_indices) <- c("Y", names(Xp))

  proc <- multivarious::concat_pre_processors(proclist, block_indices)
  v_concat <- do.call(rbind, c(list(B), V_list))
  sdev <- apply(S, 2, stats::sd)

  partial_scores <- lapply(seq_along(Xp), function(k) {
    Vk <- V_list[[k]]
    Xt <- Xp[[k]]
    kdim <- ncol(Vk)
    XtV <- Xt %*% Vk
    G <- crossprod(Vk) + ridge * diag(kdim)
    XtV %*% solve(G)
  })
  names(partial_scores) <- names(Xp)

  cor_y <- tryCatch(stats::cor(Yp, S), error = function(e) NULL)
  if (!is.null(cor_y)) cor_y[!is.finite(cor_y)] <- 0

  cor_x_list <- lapply(seq_along(Xp), function(k) {
    idx <- row_index[[k]]
    Ck <- tryCatch(stats::cor(Xp[[k]], S[idx, , drop = FALSE]), error = function(e) NULL)
    if (!is.null(Ck)) Ck[!is.finite(Ck)] <- 0
    Ck
  })
  cor_x_list <- cor_x_list[!vapply(cor_x_list, is.null, logical(1))]

  cor_loadings <- if (!is.null(cor_y) && length(cor_x_list) > 0) {
    rbind(cor_y, do.call(rbind, cor_x_list))
  } else if (!is.null(cor_y)) {
    cor_y
  } else {
    NULL
  }

  block_fit <- {
    Y_hat <- S %*% t(B)
    sse_y <- sum((Yp - Y_hat)^2)
    tss_y <- sum(Yp^2)
    r2_y <- if (tss_y > 0) 1 - sse_y / tss_y else NA_real_

    df <- data.frame(
      block = "Y",
      n = nrow(Yp),
      p = ncol(Yp),
      sse = sse_y,
      tss = tss_y,
      r2 = r2_y
    )

    for (k in seq_along(Xp)) {
      idx <- row_index[[k]]
      X_hat <- S[idx, , drop = FALSE] %*% t(V_list[[k]])
      sse <- sum((Xp[[k]] - X_hat)^2)
      tss <- sum(Xp[[k]]^2)
      r2 <- if (tss > 0) 1 - sse / tss else NA_real_
      df <- rbind(
        df,
        data.frame(
          block = names(Xp)[k],
          n = nrow(Xp[[k]]),
          p = ncol(Xp[[k]]),
          sse = sse,
          tss = tss,
          r2 = r2
        )
      )
    }
    df
  }

  fit <- multivarious::multiblock_biprojector(
    v = v_concat,
    s = S,
    sdev = sdev,
    preproc = proc,
    block_indices = block_indices,
    B = B,
    V_list = V_list,
    row_index = row_index,
    alpha_blocks = alpha_blocks,
    normalization = normalization,
    objective_trace = objective_trace,
    partial_scores = partial_scores,
    cor_loadings = cor_loadings,
    block_fit = block_fit,
    block_info = block_info,
    graph_laplacian = graph$L,
    graph_adjacency = graph$A,
    graph_lambda = graph_lambda,
    graph_form = graph_form,
    block_preproc = setNames(proclist[-1L], names(Xp)),
    anchor_preproc = proclist[[1L]],
    ridge = ridge,
    classes = c("graph_anchored_mfa", "anchored_mfa", "linked_mfa")
  )

  .muscal_attach_fit_contract(
    fit,
    method = "graph_anchored_mfa",
    task = "response_prediction",
    oos_types = c("response", "scores", "reconstruction"),
    fit_call = fit_call,
    refit_supported = TRUE,
    prediction_target = "Y",
    refit = .muscal_make_refit_spec(
      data = data_refit,
      fit_fn = function(data) {
        do.call(
          graph_anchored_mfa,
          c(
            list(
              Y = data$Y,
              X = data$X,
              row_index = data$row_index,
              block_info = data$block_info,
              preproc = preproc,
              ncomp = ncomp,
              normalization = normalization,
              alpha = alpha,
              feature_graph = feature_graph,
              graph_lambda = graph_lambda,
              graph_form = graph_form,
              max_iter = max_iter,
              tol = tol,
              ridge = ridge,
              verbose = FALSE
            ),
            fit_dots
          )
        )
      },
      bootstrap_fn = .muscal_bootstrap_anchored_data,
      permutation_fn = .muscal_permute_anchored_data,
      resample_unit = "anchor_rows"
    )
  )
}

#' Project New Rows into a Graph-Anchored MFA Space
#'
#' @param x A fitted `graph_anchored_mfa` object.
#' @param new_data Numeric matrix/data.frame of new rows for a known auxiliary
#'   block.
#' @param block Character or integer identifying the auxiliary block whose
#'   loading matrix should be used.
#' @param preprocess Logical; if `TRUE` (default), applies the fitted block
#'   preprocessing pipeline before projection.
#' @param ... Unused.
#'
#' @return A numeric matrix of latent scores.
#' @export
project.graph_anchored_mfa <- function(x, new_data, block, preprocess = TRUE, ...) {
  chk::chk_true(inherits(x, "graph_anchored_mfa"))
  .muscal_project_row_linked_scores(x, new_data, block = block, preprocess = preprocess)
}

#' Predict from a Graph-Anchored MFA Fit
#'
#' @param object A fitted `graph_anchored_mfa` object.
#' @param new_data Numeric matrix/data.frame of new rows from a known auxiliary
#'   block.
#' @param block Character or integer identifying the auxiliary block.
#' @param type One of `"response"`, `"scores"`, or `"reconstruction"`.
#' @param preprocess Logical; if `TRUE` (default), applies the fitted block
#'   preprocessing pipeline before projection.
#' @param ... Unused.
#'
#' @return A numeric matrix of predicted `Y` rows when `type = "response"`,
#'   reconstructed block rows when `type = "reconstruction"`, or latent scores
#'   when `type = "scores"`.
#' @export
predict.graph_anchored_mfa <- function(object,
                                       new_data,
                                       block,
                                       type = c("response", "scores", "reconstruction"),
                                       preprocess = TRUE,
                                       ...) {
  type <- match.arg(type)
  if (type == "scores") {
    return(project.graph_anchored_mfa(object, new_data, block = block, preprocess = preprocess, ...))
  }
  if (type == "reconstruction") {
    return(.muscal_predict_block_reconstruction(object, new_data, block = block, preprocess = preprocess))
  }

  .muscal_predict_anchor_response(object, new_data, block = block, preprocess = preprocess)
}

.gamfa_is_matrix_like <- function(x) {
  is.matrix(x) || is.data.frame(x)
}

.gamfa_flatten_blocks <- function(X, row_index, block_info = NULL) {
  if (length(X) == 0L) {
    return(list(X = list(), row_index = list(), block_info = data.frame()))
  }

  is_flat <- all(vapply(X, .gamfa_is_matrix_like, logical(1)))
  if (is_flat) {
    chk::chk_equal(length(X), length(row_index))
    X_flat <- lapply(X, as.matrix)
    idx_flat <- lapply(row_index, as.integer)
    if (is.null(names(X_flat))) names(X_flat) <- paste0("X", seq_along(X_flat))
    if (is.null(names(idx_flat))) names(idx_flat) <- names(X_flat)

    info <- .gamfa_validate_block_info(
      block_info = block_info,
      block_names = names(X_flat),
      subject_names = names(X_flat),
      domain_names = rep(NA_character_, length(X_flat))
    )
    return(list(X = X_flat, row_index = idx_flat, block_info = info))
  }

  chk::chk_equal(length(X), length(row_index))
  subj_names <- names(X)
  if (is.null(subj_names)) subj_names <- paste0("S", seq_along(X))

  X_flat <- list()
  idx_flat <- list()
  info_rows <- vector("list", 0L)

  for (si in seq_along(X)) {
    Xi <- X[[si]]
    ri <- row_index[[si]]
    if (.gamfa_is_matrix_like(Xi)) {
      Xi <- list(D1 = Xi)
      ri <- list(D1 = ri)
    }
    if (!is.list(Xi) || !is.list(ri)) {
      stop("Nested `X` and `row_index` must be lists of subject-level lists.", call. = FALSE)
    }

    dom_names <- names(Xi)
    if (is.null(dom_names)) dom_names <- paste0("D", seq_along(Xi))
    if (is.null(names(ri))) names(ri) <- dom_names

    for (di in seq_along(Xi)) {
      dom <- dom_names[di]
      if (!.gamfa_is_matrix_like(Xi[[di]])) {
        stop("Every observed nested block must be a matrix or data.frame.", call. = FALSE)
      }
      if (!(dom %in% names(ri))) {
        stop(sprintf("Missing row_index entry for subject '%s', domain '%s'.", subj_names[si], dom), call. = FALSE)
      }
      block_name <- paste(subj_names[si], dom, sep = "__")
      X_flat[[block_name]] <- as.matrix(Xi[[di]])
      idx_flat[[block_name]] <- as.integer(ri[[dom]])
      info_rows[[length(info_rows) + 1L]] <- data.frame(
        block = block_name,
        subject = subj_names[si],
        domain = dom,
        stringsAsFactors = FALSE
      )
    }
  }

  info_auto <- do.call(rbind, info_rows)
  if (is.null(block_info)) {
    info <- info_auto
  } else {
    info <- .gamfa_validate_block_info(
      block_info = block_info,
      block_names = names(X_flat),
      subject_names = info_auto$subject,
      domain_names = info_auto$domain
    )
  }

  list(X = X_flat, row_index = idx_flat, block_info = info)
}

.gamfa_validate_block_info <- function(block_info, block_names, subject_names, domain_names) {
  if (is.null(block_info)) {
    return(data.frame(
      block = block_names,
      subject = subject_names,
      domain = domain_names,
      stringsAsFactors = FALSE
    ))
  }

  if (!is.data.frame(block_info)) {
    stop("`block_info` must be NULL or a data.frame.", call. = FALSE)
  }
  if (nrow(block_info) != length(block_names)) {
    stop("`block_info` must have one row per flattened auxiliary block.", call. = FALSE)
  }
  info <- block_info
  if (!("block" %in% names(info))) info$block <- block_names
  if (anyDuplicated(info$block)) stop("`block_info$block` must be unique.", call. = FALSE)
  if (!setequal(info$block, block_names)) {
    stop("`block_info$block` must match the flattened block names.", call. = FALSE)
  }
  info <- info[match(block_names, info$block), , drop = FALSE]
  rownames(info) <- NULL
  info
}

.gamfa_parse_feature_graph <- function(feature_graph,
                                       X_list,
                                       graph_lambda,
                                       graph_form) {
  p_total <- sum(vapply(X_list, ncol, integer(1)))
  zero_graph <- list(
    enabled = FALSE,
    L = Matrix::Diagonal(p_total, x = 0),
    A = Matrix::Diagonal(p_total, x = 0)
  )

  if (isTRUE(graph_lambda == 0) || is.null(feature_graph)) {
    return(zero_graph)
  }

  if (is.character(feature_graph) && length(feature_graph) == 1L && identical(feature_graph, "colnames")) {
    A <- .gamfa_graph_from_colnames(X_list)
    L <- .gamfa_laplacian_from_adjacency(A, normalized = identical(graph_form, "normalized_laplacian"))
    return(list(enabled = TRUE, L = L, A = A))
  }

  if (is.data.frame(feature_graph)) {
    A <- .gamfa_graph_from_edges(feature_graph, X_list)
    L <- .gamfa_laplacian_from_adjacency(A, normalized = identical(graph_form, "normalized_laplacian"))
    return(list(enabled = TRUE, L = L, A = A))
  }

  if (inherits(feature_graph, "Matrix") || is.matrix(feature_graph)) {
    G <- Matrix::Matrix(feature_graph, sparse = TRUE)
    if (nrow(G) != p_total || ncol(G) != p_total) {
      stop(
        sprintf("feature_graph has dimension %d x %d, expected %d x %d.", nrow(G), ncol(G), p_total, p_total),
        call. = FALSE
      )
    }
    G <- (G + Matrix::t(G)) / 2
    if (graph_form == "adjacency") {
      A <- G
      diag(A) <- 0
      L <- .gamfa_laplacian_from_adjacency(A, normalized = FALSE)
      return(list(enabled = TRUE, L = L, A = A))
    }
    L <- G
    A <- NULL
    return(list(enabled = TRUE, L = L, A = A))
  }

  stop(
    "feature_graph must be NULL, 'colnames', a data.frame, or a square matrix.",
    call. = FALSE
  )
}

.gamfa_graph_from_colnames <- function(X_list) {
  cn_all <- lapply(X_list, colnames)
  if (any(vapply(cn_all, is.null, logical(1)))) {
    stop("feature_graph='colnames' requires non-NULL colnames for every auxiliary block.", call. = FALSE)
  }

  offsets <- .gamfa_feature_offsets(X_list)
  name_map <- new.env(parent = emptyenv())
  for (k in seq_along(X_list)) {
    for (j in seq_len(ncol(X_list[[k]]))) {
      nm <- cn_all[[k]][j]
      gi <- offsets$starts[k] + j - 1L
      cur <- if (exists(nm, envir = name_map, inherits = FALSE)) {
        get(nm, envir = name_map, inherits = FALSE)
      } else {
        integer(0)
      }
      assign(nm, c(cur, gi), envir = name_map)
    }
  }

  i_idx <- integer(0)
  j_idx <- integer(0)
  x_val <- numeric(0)
  for (nm in ls(envir = name_map)) {
    ids <- get(nm, envir = name_map, inherits = FALSE)
    if (length(ids) < 2L) next
    cmb <- utils::combn(ids, 2L)
    i_idx <- c(i_idx, cmb[1, ], cmb[2, ])
    j_idx <- c(j_idx, cmb[2, ], cmb[1, ])
    x_val <- c(x_val, rep(1, ncol(cmb) * 2L))
  }

  Matrix::sparseMatrix(
    i = i_idx,
    j = j_idx,
    x = x_val,
    dims = c(offsets$total, offsets$total)
  )
}

.gamfa_graph_from_edges <- function(feature_graph, X_list) {
  req <- c("block1", "feature1", "block2", "feature2")
  if (!all(req %in% names(feature_graph))) {
    stop("feature_graph data.frame must contain columns: block1, feature1, block2, feature2.", call. = FALSE)
  }
  if (!("weight" %in% names(feature_graph))) feature_graph$weight <- 1

  x_names <- names(X_list)
  offsets <- .gamfa_feature_offsets(X_list)
  i_idx <- integer(nrow(feature_graph) * 2L)
  j_idx <- integer(nrow(feature_graph) * 2L)
  x_val <- numeric(nrow(feature_graph) * 2L)

  for (i in seq_len(nrow(feature_graph))) {
    k1 <- .gamfa_resolve_block(feature_graph$block1[[i]], x_names)
    k2 <- .gamfa_resolve_block(feature_graph$block2[[i]], x_names)
    f1 <- .gamfa_resolve_feature(feature_graph$feature1[[i]], X_list[[k1]], k1)
    f2 <- .gamfa_resolve_feature(feature_graph$feature2[[i]], X_list[[k2]], k2)
    w <- as.numeric(feature_graph$weight[[i]])
    if (!is.finite(w) || w < 0) {
      stop("feature_graph$weight must be finite and non-negative.", call. = FALSE)
    }
    g1 <- offsets$starts[k1] + f1 - 1L
    g2 <- offsets$starts[k2] + f2 - 1L
    pos <- 2L * i - 1L
    i_idx[pos] <- g1
    j_idx[pos] <- g2
    x_val[pos] <- w
    i_idx[pos + 1L] <- g2
    j_idx[pos + 1L] <- g1
    x_val[pos + 1L] <- w
  }

  Matrix::sparseMatrix(
    i = i_idx,
    j = j_idx,
    x = x_val,
    dims = c(offsets$total, offsets$total)
  )
}

.gamfa_feature_offsets <- function(X_list) {
  p_vec <- vapply(X_list, ncol, integer(1))
  starts <- cumsum(c(1L, p_vec))[seq_along(p_vec)]
  ends <- starts + p_vec - 1L
  list(p_vec = p_vec, starts = starts, ends = ends, total = sum(p_vec))
}

.gamfa_resolve_block <- function(block, block_names) {
  if (is.numeric(block) || is.integer(block)) {
    bi <- as.integer(block)
    if (bi < 1L || bi > length(block_names)) {
      stop("Block index out of range in feature_graph.", call. = FALSE)
    }
    return(bi)
  }
  bi <- match(as.character(block), block_names)
  if (is.na(bi)) {
    stop("Block name not found in feature_graph.", call. = FALSE)
  }
  bi
}

.gamfa_resolve_feature <- function(feature, X_block, block_index) {
  if (is.numeric(feature) || is.integer(feature)) {
    fi <- as.integer(feature)
    if (fi < 1L || fi > ncol(X_block)) {
      stop(sprintf("Feature index out of range in block %d.", block_index), call. = FALSE)
    }
    return(fi)
  }
  cn <- colnames(X_block)
  if (is.null(cn)) {
    stop("Feature names were used in feature_graph, but some block has NULL colnames.", call. = FALSE)
  }
  fi <- match(as.character(feature), cn)
  if (is.na(fi)) {
    stop(sprintf("Feature name '%s' not found in block %d.", as.character(feature), block_index), call. = FALSE)
  }
  fi
}

.gamfa_laplacian_from_adjacency <- function(A, normalized = FALSE) {
  A <- Matrix::Matrix(A, sparse = TRUE)
  A <- (A + Matrix::t(A)) / 2
  diag(A) <- 0
  deg <- Matrix::rowSums(A)
  L <- Matrix::Diagonal(x = as.numeric(deg)) - A
  if (!normalized) return(L)

  inv_sqrt <- ifelse(deg > 0, 1 / sqrt(deg), 0)
  D_half <- Matrix::Diagonal(x = as.numeric(inv_sqrt))
  D_half %*% L %*% D_half
}

.gamfa_update_V_graph <- function(S_list,
                                  X_list,
                                  alpha_blocks,
                                  graph_laplacian,
                                  graph_lambda,
                                  ridge = 1e-8) {
  K <- ncol(S_list[[1]])
  offsets <- .gamfa_feature_offsets(X_list)
  P_total <- offsets$total

  rhs_mat <- matrix(0, nrow = P_total, ncol = K)
  A <- graph_lambda * Matrix::kronecker(Matrix::Diagonal(K), graph_laplacian) +
    ridge * Matrix::Diagonal(P_total * K)

  for (k in seq_along(X_list)) {
    rows_k <- offsets$starts[k]:offsets$ends[k]
    Sk <- S_list[[k]]
    Xk <- X_list[[k]]
    a_k <- alpha_blocks[k]
    Gk <- a_k * crossprod(Sk)
    mask <- rep.int(0, P_total)
    mask[rows_k] <- 1
    Pk <- Matrix::Diagonal(x = mask)
    A <- A + Matrix::kronecker(Gk, Pk)
    rhs_mat[rows_k, ] <- a_k * crossprod(Xk, Sk)
  }

  sol <- Matrix::solve(A, as.vector(rhs_mat))
  V_big <- matrix(as.numeric(sol), nrow = P_total, ncol = K)

  out <- lapply(seq_along(X_list), function(k) {
    Vk <- V_big[offsets$starts[k]:offsets$ends[k], , drop = FALSE]
    rownames(Vk) <- colnames(X_list[[k]])
    Vk
  })
  names(out) <- names(X_list)
  out
}

.gamfa_objective <- function(Y,
                             S,
                             B,
                             X_list,
                             V_list,
                             row_index,
                             alpha_y,
                             alpha_blocks,
                             graph_laplacian,
                             graph_lambda) {
  err_y <- alpha_y * sum((Y - S %*% t(B))^2)
  err_x <- 0
  for (k in seq_along(X_list)) {
    Sk <- S[row_index[[k]], , drop = FALSE]
    resid <- X_list[[k]] - Sk %*% t(V_list[[k]])
    err_x <- err_x + alpha_blocks[k] * sum(resid^2)
  }

  pen <- 0
  if (isTRUE(graph_lambda > 0) && !is.null(graph_laplacian)) {
    V_big <- do.call(rbind, V_list)
    LV <- graph_laplacian %*% V_big
    pen <- graph_lambda * sum(LV * V_big)
  }

  err_y + err_x + pen
}

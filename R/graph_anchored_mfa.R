#' Graph-Anchored Multiple Factor Analysis
#'
#' @description
#' `graph_anchored_mfa()` extends [anchored_mfa()] to settings where auxiliary
#' blocks do not share aligned columns and borrowing across blocks is induced by
#' sparse graphs over both auxiliary features and anchor rows. Rows of each
#' auxiliary block are linked to a common anchor matrix `Y` through
#' `row_index`, while graph Laplacian penalties can encourage both similar
#' auxiliary features across blocks and similar anchor rows in `Y` to share
#' latent structure.
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
#' `V_k`, let `L_G` be a feature-graph Laplacian over all auxiliary features,
#' and let `L_S` be a score-graph Laplacian over anchor rows of `Y`. The
#' estimator minimizes the anchored-MFA reconstruction loss plus the graph
#' smoothness terms
#' \deqn{\lambda_G \mathrm{tr}(V^\top L_G V) + \lambda_S \mathrm{tr}(S^\top L_S S)}
#' and ridge penalties used to identify and stabilize the fitted factors,
#' \deqn{\mathrm{ridge} \left(\|S\|_F^2 + \|B\|_F^2 + \sum_k \|V_k\|_F^2\right).}
#' The anchored score matrix can be fit either with the historical
#' unconstrained/QR update (`score_constraint = "none"`) or with an explicit
#' orthonormal constraint (`score_constraint = "orthonormal"`).
#' The score-graph term is equivalent to
#' \deqn{\frac{\lambda_S}{2} \sum_{i,j} w_{ij} \|S_{i\cdot} - S_{j\cdot}\|_2^2}
#' for adjacency weights `w_{ij}`, so nearby rows in `Y` are encouraged to have
#' nearby latent scores.
#'
#' When both `graph_lambda = 0` (or `feature_graph = NULL`) and
#' `score_graph_lambda = 0` (or `score_graph = NULL`), the method reduces to
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
#' ## Score graph
#' `score_graph` may be:
#' * `NULL` (no score penalty),
#' * `"knn"` to construct a symmetric k-nearest-neighbor graph on preprocessed
#'   rows of `Y`,
#' * a data frame with columns `row1`, `row2`, and optional `weight`, or
#' * a square sparse/dense matrix interpreted according to `score_graph_form`.
#'
#' For more control over row-graph construction, an external graph builder such
#' as the `adjoin` package can be used to create a weighted kNN adjacency or
#' Laplacian from `Y`, then supplied here as `score_graph`.
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
#' @param score_constraint Identification strategy for the anchored score
#'   matrix. `"none"` uses the historical unconstrained update followed by QR
#'   normalization inside each ALS iteration. `"orthonormal"` enforces
#'   `S^\top S = I` with a constrained majorization/polar update.
#' @param feature_graph Feature-graph specification; see Details.
#' @param graph_lambda Non-negative scalar controlling graph-penalty strength.
#' @param graph_form Interpretation of `feature_graph` when it is matrix-like, or
#'   the Laplacian construction used for edge-based inputs.
#' @param score_graph Optional score-graph specification; see Details.
#' @param score_graph_lambda Non-negative scalar controlling row/score-graph
#'   smoothing strength.
#' @param score_graph_form Interpretation of `score_graph` when it is
#'   matrix-like, or the Laplacian construction used for edge-based inputs.
#' @param score_graph_k Integer number of neighbors used when
#'   `score_graph = "knn"`.
#' @param score_graph_weight_mode Weighting applied when `score_graph = "knn"`.
#'   `"heat"` uses a Gaussian similarity kernel on neighbor distances;
#'   `"binary"` assigns weight 1 to every retained neighbor edge.
#' @param score_graph_sigma Optional positive bandwidth used when
#'   `score_graph = "knn"` and `score_graph_weight_mode = "heat"`. If `NULL`, a
#'   robust value is estimated from the k-th neighbor distances.
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
#'   in `graph_laplacian`, `graph_lambda`, `score_graph_laplacian`, and
#'   `score_graph_lambda`.
#'
#' @examples
#' set.seed(1)
#' N <- 30
#' Y <- matrix(rnorm(N * 3), N, 3)
#' X1 <- matrix(rnorm(15 * 5), 15, 5)
#' X2 <- matrix(rnorm(15 * 4), 15, 4)
#' idx <- list(X1 = sample.int(N, 15), X2 = sample.int(N, 15))
#'
#' fit <- graph_anchored_mfa(
#'   Y = Y,
#'   X = list(X1 = X1, X2 = X2),
#'   row_index = idx,
#'   ncomp = 2,
#'   score_graph = "knn",
#'   score_graph_k = 5,
#'   score_graph_weight_mode = "heat",
#'   score_graph_lambda = 1
#' )
#'
#' if (requireNamespace("adjoin", quietly = TRUE)) {
#'   A <- adjoin::graph_weights(
#'     Y,
#'     neighbor_mode = "knn",
#'     k = 5,
#'     weight_mode = "heat"
#'   )
#'
#'   fit_adjoin <- graph_anchored_mfa(
#'     Y = Y,
#'     X = list(X1 = X1, X2 = X2),
#'     row_index = idx,
#'     ncomp = 2,
#'     score_graph = A,
#'     score_graph_form = "adjacency",
#'     score_graph_lambda = 1
#'   )
#' }
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
                               score_constraint = c("none", "orthonormal"),
                               feature_graph = NULL,
                               graph_lambda = 0,
                               graph_form = c("laplacian", "adjacency", "normalized_laplacian"),
                               score_graph = NULL,
                               score_graph_lambda = 0,
                               score_graph_form = c("laplacian", "adjacency", "normalized_laplacian"),
                               score_graph_k = 10,
                               score_graph_weight_mode = c("heat", "binary"),
                               score_graph_sigma = NULL,
                               max_iter = 50,
                               tol = 1e-6,
                               ridge = 1e-8,
                               verbose = FALSE,
                               use_future = FALSE,
                               ...) {
  fit_call <- match.call(expand.dots = FALSE)
  fit_dots <- list(...)
  normalization <- match.arg(normalization)
  score_constraint <- match.arg(score_constraint)
  graph_form <- match.arg(graph_form)
  score_graph_form <- match.arg(score_graph_form)
  score_graph_weight_mode <- match.arg(score_graph_weight_mode)
  chk::chk_flag(use_future)
  if (isTRUE(use_future) && !requireNamespace("furrr", quietly = TRUE)) {
    stop("use_future = TRUE requires the 'furrr' package.", call. = FALSE)
  }

  Y <- .gamfa_as_numeric_matrix(Y, what = "Y")
  chk::chk_list(X)
  chk::chk_list(row_index)
  ncomp <- .gamfa_as_positive_integer(ncomp, what = "ncomp")
  chk::chk_numeric(graph_lambda)
  chk::chk_gte(graph_lambda, 0)
  chk::chk_numeric(score_graph_lambda)
  chk::chk_gte(score_graph_lambda, 0)
  score_graph_k <- .gamfa_as_positive_integer(score_graph_k, what = "score_graph_k")
  if (!is.null(score_graph_sigma)) {
    chk::chk_numeric(score_graph_sigma)
    chk::chk_gt(score_graph_sigma, 0)
  }
  max_iter <- .gamfa_as_positive_integer(max_iter, what = "max_iter")
  chk::chk_numeric(tol)
  chk::chk_gte(tol, 0)
  chk::chk_numeric(ridge)
  chk::chk_gte(ridge, 0)
  chk::chk_flag(verbose)

  flat <- .gamfa_flatten_blocks(X, row_index, block_info)
  X <- flat$X
  row_index <- flat$row_index
  block_info <- flat$block_info
  X <- lapply(X, function(x) .gamfa_as_numeric_matrix(x, what = "X block"))
  data_refit <- list(
    Y = Y,
    X = X,
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
  fitted_proclist <- .muscal_materialize_block_preprocessors(blocks, proclist)
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
      method <- .muscal_svd_method(mat, ncomp = 1)
      tryCatch(
        multivarious::svd_wrapper(mat, ncomp = 1, method = method)$sdev[1],
        error = function(e) multivarious::svd_wrapper(mat, ncomp = 1, method = "base")$sdev[1]
      )
    }
    sapply(prep_res$Xp, function(M) 1 / (first_sv(M)^2 + ridge))
  }
  names(alpha_blocks) <- c("Y", names(X))

  K <- min(as.integer(ncomp), nrow(Yp), ncol(Yp))
  if (K < 1L) stop("ncomp must be at least 1 and <= min(nrow(Y), ncol(Y)).", call. = FALSE)

  graph <- .gamfa_parse_feature_graph(
    feature_graph = feature_graph,
    X_list = Xp,
    graph_lambda = graph_lambda,
    graph_form = graph_form
  )
  score_graph_spec <- .gamfa_parse_score_graph(
    score_graph = score_graph,
    Y = Yp,
    score_graph_lambda = score_graph_lambda,
    score_graph_form = score_graph_form,
    score_graph_k = score_graph_k,
    score_graph_weight_mode = score_graph_weight_mode,
    score_graph_sigma = score_graph_sigma
  )

  core <- .gamfa_fit_anchor_engine(
    Y = Yp,
    X_list = Xp,
    row_index = row_index,
    ncomp = K,
    alpha_y = alpha_blocks[1],
    alpha_blocks = alpha_blocks[-1L],
    graph = graph,
    graph_lambda = graph_lambda,
    score_graph_spec = score_graph_spec,
    score_graph_lambda = score_graph_lambda,
    score_constraint = score_constraint,
    max_iter = max_iter,
    tol = tol,
    ridge = ridge,
    verbose = verbose,
    mode = "hard"
  )
  S <- core$S
  B <- core$B
  V_list <- core$V_list
  Z_list <- core$Z_list
  objective_trace <- core$objective_trace

  block_indices <- list()
  current <- 1L
  block_indices[[1]] <- current:(current + ncol(Yp) - 1L)
  current <- current + ncol(Yp)
  for (k in seq_along(Xp)) {
    block_indices[[k + 1L]] <- current:(current + ncol(Xp[[k]]) - 1L)
    current <- current + ncol(Xp[[k]])
  }
  names(block_indices) <- c("Y", names(Xp))

  proc <- multivarious::concat_pre_processors(fitted_proclist, block_indices)
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
    sse_y <- .lmfa_reconstruction_sse(Yp, S, B)
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
      sse <- .lmfa_reconstruction_sse(Xp[[k]], S[idx, , drop = FALSE], V_list[[k]])
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

  score_stack <- .ramfa_stack_scores(Z_list)

  fit <- multivarious::multiblock_biprojector(
    v = v_concat,
    s = S,
    sdev = sdev,
    preproc = proc,
    block_indices = block_indices,
    B = B,
    S = S,
    V_list = V_list,
    Z_list = Z_list,
    score_index = score_stack$score_index,
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
    score_graph_laplacian = score_graph_spec$L,
    score_graph_adjacency = score_graph_spec$A,
    score_graph_lambda = score_graph_lambda,
    score_graph_form = score_graph_form,
    score_graph_k = score_graph_k,
    score_graph_weight_mode = score_graph_weight_mode,
    score_graph_sigma = score_graph_sigma,
    block_preproc = setNames(fitted_proclist[-1L], names(Xp)),
    anchor_preproc = fitted_proclist[[1L]],
    ridge = ridge,
    score_representation = "anchor_scores",
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
              score_constraint = score_constraint,
              feature_graph = feature_graph,
              graph_lambda = graph_lambda,
              graph_form = graph_form,
              score_graph = score_graph,
              score_graph_lambda = score_graph_lambda,
              score_graph_form = score_graph_form,
              score_graph_k = score_graph_k,
              score_graph_weight_mode = score_graph_weight_mode,
              score_graph_sigma = score_graph_sigma,
              max_iter = max_iter,
              tol = tol,
              ridge = ridge,
              verbose = FALSE,
              use_future = use_future
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

#' Coupled Graph-Anchored Multiple Factor Analysis
#'
#' @description
#' `coupled_graph_anchored_mfa()` extends [graph_anchored_mfa()] by introducing
#' block-specific row scores `Z_k` that are softly coupled to a shared
#' stimulus-level score matrix `S`. This is useful when each auxiliary block is
#' expected to express the common stimulus structure with subject- or
#' domain-specific deviations.
#'
#' @details
#' The fitted model has the form:
#' \deqn{Y \approx S B^\top}
#' \deqn{X_k \approx Z_k V_k^\top}
#' with a consensus penalty
#' \deqn{\mu \sum_k \|Z_k - S[\mathrm{idx}_k,]\|_F^2.}
#'
#' Feature-side and score-side graph penalties are inherited from
#' [graph_anchored_mfa()]:
#' \deqn{\lambda_G \mathrm{tr}(V^\top L_G V) + \lambda_S \mathrm{tr}(S^\top L_S S).}
#' As in [graph_anchored_mfa()], `score_constraint = "none"` uses the
#' historical unconstrained/QR score update, while
#' `score_constraint = "orthonormal"` enforces `S^\top S = I` directly during
#' fitting.
#'
#' When `coupling_lambda` is large, block-specific scores are strongly tied to
#' the shared anchor score space; when it is small, blocks can deviate more
#' freely while still being linked through the common response target `Y`.
#'
#' @inheritParams graph_anchored_mfa
#' @param coupling_lambda Non-negative scalar controlling the strength of the
#'   consensus penalty tying each block-specific score matrix `Z_k` to the
#'   shared anchor scores `S[row_index[[k]], ]`.
#'
#' @return An object inheriting from `multivarious::multiblock_biprojector` with
#'   additional classes `coupled_graph_anchored_mfa` and `linked_mfa`. The
#'   object contains shared anchor scores in `s`, block-specific row scores in
#'   `Z_list`, auxiliary loadings in `V_list`, and anchor loadings in `B`.
#'
#' @export
coupled_graph_anchored_mfa <- function(Y,
                                       X,
                                       row_index,
                                       block_info = NULL,
                                       preproc = multivarious::center(),
                                       ncomp = 2,
                                       normalization = c("MFA", "None", "custom"),
                                       alpha = NULL,
                                       score_constraint = c("none", "orthonormal"),
                                       feature_graph = NULL,
                                       graph_lambda = 0,
                                       graph_form = c("laplacian", "adjacency", "normalized_laplacian"),
                                       score_graph = NULL,
                                       score_graph_lambda = 0,
                                       score_graph_form = c("laplacian", "adjacency", "normalized_laplacian"),
                                       score_graph_k = 10,
                                       score_graph_weight_mode = c("heat", "binary"),
                                       score_graph_sigma = NULL,
                                       coupling_lambda = 1,
                                       max_iter = 50,
                                       tol = 1e-6,
                                       ridge = 1e-8,
                                       verbose = FALSE,
                                       use_future = FALSE,
                                       ...) {
  fit_call <- match.call(expand.dots = FALSE)
  fit_dots <- list(...)
  normalization <- match.arg(normalization)
  score_constraint <- match.arg(score_constraint)
  graph_form <- match.arg(graph_form)
  score_graph_form <- match.arg(score_graph_form)
  score_graph_weight_mode <- match.arg(score_graph_weight_mode)
  chk::chk_flag(use_future)
  if (isTRUE(use_future) && !requireNamespace("furrr", quietly = TRUE)) {
    stop("use_future = TRUE requires the 'furrr' package.", call. = FALSE)
  }

  Y <- .gamfa_as_numeric_matrix(Y, what = "Y")
  chk::chk_list(X)
  chk::chk_list(row_index)
  ncomp <- .gamfa_as_positive_integer(ncomp, what = "ncomp")
  chk::chk_numeric(graph_lambda)
  chk::chk_gte(graph_lambda, 0)
  chk::chk_numeric(score_graph_lambda)
  chk::chk_gte(score_graph_lambda, 0)
  chk::chk_numeric(coupling_lambda)
  chk::chk_gte(coupling_lambda, 0)
  score_graph_k <- .gamfa_as_positive_integer(score_graph_k, what = "score_graph_k")
  if (!is.null(score_graph_sigma)) {
    chk::chk_numeric(score_graph_sigma)
    chk::chk_gt(score_graph_sigma, 0)
  }
  max_iter <- .gamfa_as_positive_integer(max_iter, what = "max_iter")
  chk::chk_numeric(tol)
  chk::chk_gte(tol, 0)
  chk::chk_numeric(ridge)
  chk::chk_gte(ridge, 0)
  chk::chk_flag(verbose)

  flat <- .gamfa_flatten_blocks(X, row_index, block_info)
  X <- flat$X
  row_index <- flat$row_index
  block_info <- flat$block_info
  X <- lapply(X, function(x) .gamfa_as_numeric_matrix(x, what = "X block"))
  data_refit <- list(
    Y = Y,
    X = X,
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
  fitted_proclist <- .muscal_materialize_block_preprocessors(blocks, proclist)

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
      method <- .muscal_svd_method(mat, ncomp = 1)
      tryCatch(
        multivarious::svd_wrapper(mat, ncomp = 1, method = method)$sdev[1],
        error = function(e) multivarious::svd_wrapper(mat, ncomp = 1, method = "base")$sdev[1]
      )
    }
    sapply(prep_res$Xp, function(M) 1 / (first_sv(M)^2 + ridge))
  }
  names(alpha_blocks) <- c("Y", names(X))

  K <- min(as.integer(ncomp), nrow(Yp), ncol(Yp))
  if (K < 1L) stop("ncomp must be at least 1 and <= min(nrow(Y), ncol(Y)).", call. = FALSE)

  graph <- .gamfa_parse_feature_graph(
    feature_graph = feature_graph,
    X_list = Xp,
    graph_lambda = graph_lambda,
    graph_form = graph_form
  )
  score_graph_spec <- .gamfa_parse_score_graph(
    score_graph = score_graph,
    Y = Yp,
    score_graph_lambda = score_graph_lambda,
    score_graph_form = score_graph_form,
    score_graph_k = score_graph_k,
    score_graph_weight_mode = score_graph_weight_mode,
    score_graph_sigma = score_graph_sigma
  )

  core <- .gamfa_fit_anchor_engine(
    Y = Yp,
    X_list = Xp,
    row_index = row_index,
    ncomp = K,
    alpha_y = alpha_blocks[1],
    alpha_blocks = alpha_blocks[-1L],
    graph = graph,
    graph_lambda = graph_lambda,
    score_graph_spec = score_graph_spec,
    score_graph_lambda = score_graph_lambda,
    score_constraint = score_constraint,
    max_iter = max_iter,
    tol = tol,
    ridge = ridge,
    verbose = verbose,
    mode = "coupled",
    coupling_lambda = coupling_lambda
  )
  S <- core$S
  B <- core$B
  V_list <- core$V_list
  Z_list <- core$Z_list
  objective_trace <- core$objective_trace

  block_indices <- list()
  current <- 1L
  block_indices[[1]] <- current:(current + ncol(Yp) - 1L)
  current <- current + ncol(Yp)
  for (k in seq_along(Xp)) {
    block_indices[[k + 1L]] <- current:(current + ncol(Xp[[k]]) - 1L)
    current <- current + ncol(Xp[[k]])
  }
  names(block_indices) <- c("Y", names(Xp))

  proc <- multivarious::concat_pre_processors(fitted_proclist, block_indices)
  v_concat <- do.call(rbind, c(list(B), V_list))
  sdev <- apply(S, 2, stats::sd)
  partial_scores <- Z_list

  cor_y <- tryCatch(stats::cor(Yp, S), error = function(e) NULL)
  if (!is.null(cor_y)) cor_y[!is.finite(cor_y)] <- 0

  cor_x_list <- lapply(seq_along(Xp), function(k) {
    Ck <- tryCatch(stats::cor(Xp[[k]], Z_list[[k]]), error = function(e) NULL)
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
    sse_y <- .lmfa_reconstruction_sse(Yp, S, B)
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
      sse <- .lmfa_reconstruction_sse(Xp[[k]], Z_list[[k]], V_list[[k]])
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

  score_stack <- .ramfa_stack_scores(Z_list)

  fit <- multivarious::multiblock_biprojector(
    v = v_concat,
    s = S,
    sdev = sdev,
    preproc = proc,
    block_indices = block_indices,
    B = B,
    S = S,
    V_list = V_list,
    Z_list = Z_list,
    score_index = score_stack$score_index,
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
    score_graph_laplacian = score_graph_spec$L,
    score_graph_adjacency = score_graph_spec$A,
    score_graph_lambda = score_graph_lambda,
    score_graph_form = score_graph_form,
    score_graph_k = score_graph_k,
    score_graph_weight_mode = score_graph_weight_mode,
    score_graph_sigma = score_graph_sigma,
    coupling_lambda = coupling_lambda,
    block_preproc = setNames(fitted_proclist[-1L], names(Xp)),
    anchor_preproc = fitted_proclist[[1L]],
    ridge = ridge,
    score_representation = "anchor_scores",
    classes = c("coupled_graph_anchored_mfa", "linked_mfa")
  )

  .muscal_attach_fit_contract(
    fit,
    method = "coupled_graph_anchored_mfa",
    task = "response_prediction",
    oos_types = c("response", "scores", "reconstruction"),
    fit_call = fit_call,
    refit_supported = TRUE,
    prediction_target = "Y",
    refit = .muscal_make_refit_spec(
      data = data_refit,
      fit_fn = function(data) {
        do.call(
          coupled_graph_anchored_mfa,
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
              score_constraint = score_constraint,
              feature_graph = feature_graph,
              graph_lambda = graph_lambda,
              graph_form = graph_form,
              score_graph = score_graph,
              score_graph_lambda = score_graph_lambda,
              score_graph_form = score_graph_form,
              score_graph_k = score_graph_k,
              score_graph_weight_mode = score_graph_weight_mode,
              score_graph_sigma = score_graph_sigma,
              coupling_lambda = coupling_lambda,
              max_iter = max_iter,
              tol = tol,
              ridge = ridge,
              verbose = FALSE,
              use_future = use_future
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

#' Project New Rows into a Coupled Graph-Anchored MFA Space
#'
#' @param x A fitted `coupled_graph_anchored_mfa` object.
#' @param new_data Numeric matrix/data.frame of new rows for a known auxiliary
#'   block.
#' @param block Character or integer identifying the auxiliary block whose
#'   loading matrix should be used.
#' @param preprocess Logical; if `TRUE` (default), applies the fitted block
#'   preprocessing pipeline before projection.
#' @param ... Unused.
#'
#' @return A numeric matrix of block-specific latent scores.
#' @export
project.coupled_graph_anchored_mfa <- function(x, new_data, block, preprocess = TRUE, ...) {
  chk::chk_true(inherits(x, "coupled_graph_anchored_mfa"))
  .muscal_project_row_linked_scores(x, new_data, block = block, preprocess = preprocess)
}

#' Predict from a Coupled Graph-Anchored MFA Fit
#'
#' @param object A fitted `coupled_graph_anchored_mfa` object.
#' @param new_data Numeric matrix/data.frame of new rows from a known auxiliary
#'   block.
#' @param block Character or integer identifying the auxiliary block.
#' @param type One of `"response"`, `"scores"`, or `"reconstruction"`.
#' @param preprocess Logical; if `TRUE` (default), applies the fitted block
#'   preprocessing pipeline before projection.
#' @param ... Unused.
#'
#' @return A numeric matrix of predicted `Y` rows when `type = "response"`,
#'   reconstructed block rows when `type = "reconstruction"`, or block-specific
#'   latent scores when `type = "scores"`.
#' @export
predict.coupled_graph_anchored_mfa <- function(object,
                                               new_data,
                                               block,
                                               type = c("response", "scores", "reconstruction"),
                                               preprocess = TRUE,
                                               ...) {
  type <- match.arg(type)
  if (type == "scores") {
    return(project.coupled_graph_anchored_mfa(object, new_data, block = block, preprocess = preprocess, ...))
  }
  if (type == "reconstruction") {
    return(.muscal_predict_block_reconstruction(object, new_data, block = block, preprocess = preprocess))
  }

  .muscal_predict_anchor_response(object, new_data, block = block, preprocess = preprocess)
}

.gamfa_is_matrix_like <- function(x) {
  is.matrix(x) || is.data.frame(x)
}

.gamfa_as_numeric_matrix <- function(x, what) {
  if (!.gamfa_is_matrix_like(x)) {
    stop(sprintf("`%s` must be a matrix or data.frame.", what), call. = FALSE)
  }
  x <- as.matrix(x)
  if (!is.numeric(x)) {
    stop(sprintf("`%s` must be numeric.", what), call. = FALSE)
  }
  if (anyNA(x) || any(!is.finite(x))) {
    stop(sprintf("`%s` must contain only finite values.", what), call. = FALSE)
  }
  x
}

.gamfa_as_positive_integer <- function(x, what) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x < 1 || x != round(x)) {
    stop(sprintf("`%s` must be a positive integer scalar.", what), call. = FALSE)
  }
  as.integer(x)
}

.gamfa_flatten_blocks <- function(X, row_index, block_info = NULL) {
  if (length(X) == 0L) {
    return(list(X = list(), row_index = list(), block_info = data.frame()))
  }

  is_flat <- all(vapply(X, .gamfa_is_matrix_like, logical(1)))
  if (is_flat) {
    chk::chk_equal(length(X), length(row_index))
    X_flat <- lapply(X, as.matrix)
    if (is.null(names(X_flat))) names(X_flat) <- paste0("X", seq_along(X_flat))
    idx_flat <- .gamfa_align_named_index_list(
      index_list = row_index,
      target_names = names(X_flat),
      what = "flat row_index"
    )
    idx_flat <- lapply(idx_flat, as.integer)

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
  row_index <- .gamfa_align_named_index_list(
    index_list = row_index,
    target_names = subj_names,
    what = "top-level row_index"
  )

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
    ri <- .gamfa_align_named_index_list(
      index_list = ri,
      target_names = dom_names,
      what = sprintf("row_index for subject '%s'", subj_names[si])
    )

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

.gamfa_align_named_index_list <- function(index_list, target_names, what) {
  if (is.null(names(index_list))) {
    names(index_list) <- target_names
    return(index_list)
  }

  if (!setequal(names(index_list), target_names)) {
    stop(sprintf("%s names must match the corresponding X names exactly.", what), call. = FALSE)
  }

  index_list[target_names]
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
    G <- .gamfa_validate_graph_matrix(G, form = graph_form, what = "feature_graph")
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

.gamfa_parse_score_graph <- function(score_graph,
                                     Y,
                                     score_graph_lambda,
                                     score_graph_form,
                                     score_graph_k,
                                     score_graph_weight_mode,
                                     score_graph_sigma) {
  n <- nrow(Y)
  zero_graph <- list(
    enabled = FALSE,
    L = Matrix::Diagonal(n, x = 0),
    A = Matrix::Diagonal(n, x = 0)
  )

  if (isTRUE(score_graph_lambda == 0) || is.null(score_graph) || n < 2L) {
    return(zero_graph)
  }

  if (is.character(score_graph) && length(score_graph) == 1L && identical(score_graph, "knn")) {
    A <- .gamfa_score_graph_from_knn(
      Y,
      k = score_graph_k,
      weight_mode = score_graph_weight_mode,
      sigma = score_graph_sigma
    )
    L <- .gamfa_laplacian_from_adjacency(A, normalized = identical(score_graph_form, "normalized_laplacian"))
    return(list(enabled = TRUE, L = L, A = A))
  }

  if (is.data.frame(score_graph)) {
    A <- .gamfa_score_graph_from_edges(score_graph, n = n, row_names = rownames(Y))
    L <- .gamfa_laplacian_from_adjacency(A, normalized = identical(score_graph_form, "normalized_laplacian"))
    return(list(enabled = TRUE, L = L, A = A))
  }

  if (inherits(score_graph, "Matrix") || is.matrix(score_graph)) {
    G <- Matrix::Matrix(score_graph, sparse = TRUE)
    if (nrow(G) != n || ncol(G) != n) {
      stop(
        sprintf("score_graph has dimension %d x %d, expected %d x %d.", nrow(G), ncol(G), n, n),
        call. = FALSE
      )
    }
    G <- .gamfa_validate_graph_matrix(G, form = score_graph_form, what = "score_graph")
    if (score_graph_form == "adjacency") {
      A <- G
      diag(A) <- 0
      L <- .gamfa_laplacian_from_adjacency(A, normalized = FALSE)
      return(list(enabled = TRUE, L = L, A = A))
    }
    return(list(enabled = TRUE, L = G, A = NULL))
  }

  stop(
    "score_graph must be NULL, 'knn', a data.frame, or a square matrix.",
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
        data.frame(block = integer(0), id = integer(0))
      }
      cur[nrow(cur) + 1L, ] <- list(block = k, id = gi)
      assign(nm, cur, envir = name_map)
    }
  }

  i_idx <- integer(0)
  j_idx <- integer(0)
  x_val <- numeric(0)
  for (nm in ls(envir = name_map)) {
    members <- get(nm, envir = name_map, inherits = FALSE)
    if (nrow(members) < 2L) next
    for (a in seq_len(nrow(members) - 1L)) {
      for (b in (a + 1L):nrow(members)) {
        if (members$block[a] == members$block[b]) next
        i_idx <- c(i_idx, members$id[a], members$id[b])
        j_idx <- c(j_idx, members$id[b], members$id[a])
        x_val <- c(x_val, 1, 1)
      }
    }
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

.gamfa_validate_graph_matrix <- function(G, form, what, tol = 1e-8) {
  G <- Matrix::Matrix(G, sparse = TRUE)
  G <- (G + Matrix::t(G)) / 2

  vals <- G@x
  if (length(vals) > 0L && any(!is.finite(vals))) {
    stop(sprintf("%s must contain only finite values.", what), call. = FALSE)
  }

  diag_vals <- Matrix::diag(G)
  off <- G
  diag(off) <- 0

  if (form == "adjacency") {
    if (length(off@x) > 0L && any(off@x < -tol)) {
      stop(sprintf("%s adjacency weights must be non-negative.", what), call. = FALSE)
    }
    if (max(abs(diag_vals)) > tol) {
      stop(sprintf("%s adjacency matrices must have a zero diagonal.", what), call. = FALSE)
    }
    return(G)
  }

  if (length(off@x) > 0L && any(off@x > tol)) {
    stop(sprintf("%s %s matrices must have non-positive off-diagonal entries.", what, form), call. = FALSE)
  }

  eig_min <- .gamfa_smallest_eigenvalue(G)
  if (!is.na(eig_min) && eig_min < -sqrt(tol)) {
    stop(sprintf("%s %s matrices must be positive semidefinite.", what, form), call. = FALSE)
  }

  if (form == "laplacian") {
    row_sums <- Matrix::rowSums(G)
    if (max(abs(row_sums)) > 1e-6) {
      stop(sprintf("%s Laplacian matrices must have row sums near zero.", what), call. = FALSE)
    }
  } else if (form == "normalized_laplacian") {
    if (any(diag_vals < -tol) || any(diag_vals > 1 + 1e-6)) {
      stop(sprintf("%s normalized Laplacians must have diagonal entries in [0, 1].", what), call. = FALSE)
    }
  }

  G
}

.gamfa_smallest_eigenvalue <- function(G) {
  n <- nrow(G)
  if (n < 1L) return(NA_real_)
  if (n == 1L) return(as.numeric(G[1, 1]))

  out <- tryCatch(
    RSpectra::eigs_sym(G, k = 1L, which = "SA", opts = list(retvec = FALSE))$values[[1L]],
    error = function(e) NA_real_
  )
  if (is.na(out)) {
    out <- min(eigen(as.matrix(G), symmetric = TRUE, only.values = TRUE)$values)
  }
  as.numeric(out)
}

.gamfa_score_graph_from_knn <- function(Y, k, weight_mode = c("heat", "binary"), sigma = NULL) {
  n <- nrow(Y)
  if (n < 2L) {
    return(Matrix::Diagonal(n, x = 0))
  }

  weight_mode <- match.arg(weight_mode)
  k_eff <- min(as.integer(k), n - 1L)
  nn <- RANN::nn2(data = Y, query = Y, k = k_eff + 1L)
  idx <- nn$nn.idx[, -1, drop = FALSE]
  dists <- nn$nn.dists[, -1, drop = FALSE]

  weights <- if (identical(weight_mode, "binary")) {
    matrix(1, nrow = nrow(idx), ncol = ncol(idx))
  } else {
    sigma_use <- if (is.null(sigma)) {
      sig <- stats::median(dists[, ncol(dists)], na.rm = TRUE)
      if (!is.finite(sig) || sig <= 0) {
        sig <- stats::quantile(as.vector(dists), probs = 0.25, na.rm = TRUE, names = FALSE)
      }
      if (!is.finite(sig) || sig <= 0) {
        stop("Unable to estimate score_graph_sigma from kNN distances.", call. = FALSE)
      }
      sig
    } else {
      if (!is.finite(sigma) || sigma <= 0) {
        stop("score_graph_sigma must be a positive scalar.", call. = FALSE)
      }
      sigma
    }
    exp(-(dists^2) / (2 * sigma_use^2))
  }

  A <- Matrix::sparseMatrix(
    i = rep(seq_len(nrow(idx)), times = ncol(idx)),
    j = as.vector(idx),
    x = as.vector(weights),
    dims = c(n, n)
  )
  A <- (A + Matrix::t(A) + abs(A - Matrix::t(A))) / 2
  Matrix::drop0(Matrix::Matrix(A, sparse = TRUE))
}

.gamfa_score_graph_from_edges <- function(score_graph, n, row_names = NULL) {
  req <- c("row1", "row2")
  if (!all(req %in% names(score_graph))) {
    stop("score_graph data.frame must contain columns: row1, row2.", call. = FALSE)
  }
  if (!("weight" %in% names(score_graph))) score_graph$weight <- 1

  i_idx <- integer(nrow(score_graph) * 2L)
  j_idx <- integer(nrow(score_graph) * 2L)
  x_val <- numeric(nrow(score_graph) * 2L)

  for (i in seq_len(nrow(score_graph))) {
    r1 <- .gamfa_resolve_anchor_row(score_graph$row1[[i]], n = n, row_names = row_names)
    r2 <- .gamfa_resolve_anchor_row(score_graph$row2[[i]], n = n, row_names = row_names)
    w <- as.numeric(score_graph$weight[[i]])
    if (!is.finite(w) || w < 0) {
      stop("score_graph$weight must be finite and non-negative.", call. = FALSE)
    }
    pos <- 2L * i - 1L
    i_idx[pos] <- r1
    j_idx[pos] <- r2
    x_val[pos] <- w
    i_idx[pos + 1L] <- r2
    j_idx[pos + 1L] <- r1
    x_val[pos + 1L] <- w
  }

  Matrix::sparseMatrix(
    i = i_idx,
    j = j_idx,
    x = x_val,
    dims = c(n, n)
  )
}

.gamfa_resolve_anchor_row <- function(row, n, row_names = NULL) {
  if (is.numeric(row) || is.integer(row)) {
    ri <- as.integer(row)
    if (ri < 1L || ri > n) {
      stop("Anchor-row index out of range in score_graph.", call. = FALSE)
    }
    return(ri)
  }

  if (is.null(row_names)) {
    stop("Character row identifiers in score_graph require non-NULL rownames(Y).", call. = FALSE)
  }
  ri <- match(as.character(row), row_names)
  if (is.na(ri)) {
    stop(sprintf("Anchor-row name '%s' not found in score_graph.", as.character(row)), call. = FALSE)
  }
  ri
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

.gamfa_fit_anchor_engine <- function(Y,
                                     X_list,
                                     row_index,
                                     ncomp,
                                     alpha_y,
                                     alpha_blocks,
                                     graph,
                                     graph_lambda,
                                     score_graph_spec,
                                     score_graph_lambda,
                                     score_constraint,
                                     max_iter,
                                     tol,
                                     ridge,
                                     verbose,
                                     mode = c("hard", "coupled"),
                                     coupling_lambda = 0) {
  mode <- match.arg(mode)
  .muscal_fit_common_space_engine(
    engine = if (identical(mode, "hard")) "hard_anchor" else "coupled_anchor",
    X_list = X_list,
    ncomp = ncomp,
    alpha_blocks = alpha_blocks,
    ridge = ridge,
    max_iter = max_iter,
    tol = tol,
    verbose = verbose,
    anchor_response = Y,
    anchor_response_alpha = alpha_y,
    anchor_response_weights = rep(1, nrow(Y)),
    row_index = row_index,
    fg = .lmfa_parse_feature_groups(NULL, X_list, 0),
    feature_lambda = 0,
    graph = graph,
    graph_lambda = graph_lambda,
    score_graph_spec = score_graph_spec,
    score_graph_lambda = score_graph_lambda,
    score_constraint = score_constraint,
    coupling_lambda = coupling_lambda
  )
}

.gamfa_update_scores_graph <- function(Y,
                                       B,
                                       X_list,
                                       V_list,
                                       row_index,
                                       alpha_y,
                                       alpha_blocks,
                                       score_graph_laplacian,
                                       score_graph_lambda,
                                       local = NULL,
                                       ridge = 1e-8) {
  K <- ncol(B)
  if (is.null(local)) {
    local <- .gamfa_score_system(
      Y = Y,
      B = B,
      X_list = X_list,
      V_list = V_list,
      row_index = row_index,
      alpha_y = alpha_y,
      alpha_blocks = alpha_blocks
    )
  }

  apply_A <- function(S) {
    .muscal_apply_row_operator(
      S = S,
      A_list = local$A_list,
      graph_laplacian = score_graph_laplacian,
      graph_lambda = score_graph_lambda,
      ridge = ridge
    )
  }
  X0 <- .lmfa_update_scores(
    Y = Y,
    B = B,
    X_list = X_list,
    V_list = V_list,
    row_index = row_index,
    alpha_y = alpha_y,
    alpha_blocks = alpha_blocks,
    ridge = ridge
  )
  .muscal_matrix_cg(
    rhs = local$rhs,
    apply_A = apply_A,
    X0 = X0,
    tol = 1e-8,
    max_iter = max(50L, 5L * K)
  )
}

.gamfa_score_gradient <- function(S,
                                  Y,
                                  B,
                                  X_list,
                                  V_list,
                                  row_index,
                                  alpha_y,
                                  alpha_blocks,
                                  score_graph_laplacian,
                                  score_graph_lambda,
                                  ridge = 0) {
  local <- .gamfa_score_system(
    Y = Y,
    B = B,
    X_list = X_list,
    V_list = V_list,
    row_index = row_index,
    alpha_y = alpha_y,
    alpha_blocks = alpha_blocks
  )
  grad <- matrix(0, nrow = nrow(S), ncol = ncol(S))
  for (i in seq_len(nrow(S))) {
    grad[i, ] <- 2 * (as.matrix(local$A_list[[i]]) %*% S[i, ] - local$rhs[i, ])
  }
  if (isTRUE(score_graph_lambda > 0) && !is.null(score_graph_laplacian)) {
    grad <- grad + 2 * score_graph_lambda * as.matrix(score_graph_laplacian %*% S)
  }
  if (ridge > 0) {
    grad <- grad + 2 * ridge * S
  }
  grad
}

.gamfa_score_system <- function(Y,
                                B,
                                X_list,
                                V_list,
                                row_index,
                                alpha_y,
                                alpha_blocks) {
  N <- nrow(Y)
  K <- ncol(B)

  BtB <- crossprod(B)
  YB <- Y %*% B
  VtV_list <- lapply(V_list, crossprod)
  XV_list <- lapply(seq_along(X_list), function(k) X_list[[k]] %*% V_list[[k]])
  agg_by_i <- lapply(seq_along(X_list), function(k) {
    .muscal_rowsum_counts(XV_list[[k]], row_index[[k]], N)
  })
  sums_by_i <- lapply(agg_by_i, `[[`, "sums")
  counts_by_i <- lapply(agg_by_i, `[[`, "counts")

  A_list <- vector("list", N)
  rhs <- matrix(0, nrow = N, ncol = K)
  for (i in seq_len(N)) {
    A_i <- alpha_y * BtB
    b_i <- alpha_y * YB[i, ]
    for (k in seq_along(X_list)) {
      cnt <- counts_by_i[[k]][i]
      if (cnt == 0L) next
      a_k <- alpha_blocks[k]
      A_i <- A_i + a_k * cnt * VtV_list[[k]]
      b_i <- b_i + a_k * sums_by_i[[k]][i, ]
    }
    A_list[[i]] <- A_i
    rhs[i, ] <- b_i
  }

  list(A_list = A_list, rhs = rhs)
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
                             graph_lambda,
                             score_graph_laplacian,
                             score_graph_lambda,
                             ridge = 0) {
  err_y <- alpha_y * .lmfa_reconstruction_sse(Y, S, B)
  err_x <- 0
  for (k in seq_along(X_list)) {
    Sk <- S[row_index[[k]], , drop = FALSE]
    err_x <- err_x + alpha_blocks[k] * .lmfa_reconstruction_sse(X_list[[k]], Sk, V_list[[k]])
  }

  pen <- 0
  if (isTRUE(graph_lambda > 0) && !is.null(graph_laplacian)) {
    V_big <- do.call(rbind, V_list)
    LV <- graph_laplacian %*% V_big
    pen <- graph_lambda * sum(LV * V_big)
  }

  if (isTRUE(score_graph_lambda > 0) && !is.null(score_graph_laplacian)) {
    LS <- score_graph_laplacian %*% S
    pen <- pen + score_graph_lambda * sum(LS * S)
  }

  ridge_pen <- ridge * (sum(S^2) + sum(B^2) + sum(vapply(V_list, function(V) sum(V^2), numeric(1))))

  err_y + err_x + pen + ridge_pen
}

.cgamfa_update_Z_list <- function(S,
                                  X_list,
                                  V_list,
                                  row_index,
                                  alpha_blocks,
                                  coupling_lambda,
                                  ridge = 1e-8) {
  K <- ncol(S)
  Ik <- diag(K)

  out <- lapply(seq_along(X_list), function(k) {
    Vk <- V_list[[k]]
    Xk <- X_list[[k]]
    idx <- row_index[[k]]
    a_k <- alpha_blocks[k]
    Atot <- a_k * crossprod(Vk) + (coupling_lambda + ridge) * Ik
    rhs <- a_k * Xk %*% Vk + coupling_lambda * S[idx, , drop = FALSE]
    Zk <- rhs %*% solve(Atot)
    rownames(Zk) <- rownames(Xk)
    Zk
  })
  names(out) <- names(X_list)
  out
}

.cgamfa_score_system <- function(Y,
                                 B,
                                 Z_list,
                                 row_index,
                                 alpha_y,
                                 coupling_lambda) {
  N <- nrow(Y)
  K <- ncol(B)
  BtB <- crossprod(B)
  YB <- Y %*% B
  agg_by_i <- lapply(seq_along(Z_list), function(k) {
    .muscal_rowsum_counts(Z_list[[k]], row_index[[k]], N)
  })
  counts_by_i <- lapply(agg_by_i, `[[`, "counts")
  sums_by_i <- lapply(agg_by_i, `[[`, "sums")

  A_list <- vector("list", N)
  rhs <- matrix(0, nrow = N, ncol = K)
  Ik <- diag(K)
  for (i in seq_len(N)) {
    coupling_count <- 0
    coupling_sum <- rep(0, K)
    for (k in seq_along(Z_list)) {
      cnt <- counts_by_i[[k]][i]
      if (cnt == 0L) next
      coupling_count <- coupling_count + cnt
      coupling_sum <- coupling_sum + sums_by_i[[k]][i, ]
    }

    A_i <- alpha_y * BtB + coupling_lambda * coupling_count * Ik
    b_i <- alpha_y * YB[i, ] + coupling_lambda * coupling_sum
    A_list[[i]] <- A_i
    rhs[i, ] <- b_i
  }

  list(A_list = A_list, rhs = rhs)
}

.cgamfa_update_scores <- function(Y,
                                  B,
                                  Z_list,
                                  row_index,
                                  alpha_y,
                                  coupling_lambda,
                                  local = NULL,
                                  ridge = 1e-8) {
  if (is.null(local)) {
    local <- .cgamfa_score_system(
      Y = Y,
      B = B,
      Z_list = Z_list,
      row_index = row_index,
      alpha_y = alpha_y,
      coupling_lambda = coupling_lambda
    )
  }
  N <- nrow(Y)
  K <- ncol(B)
  S <- matrix(0, nrow = N, ncol = K)
  for (i in seq_len(N)) {
    A_i <- as.matrix(local$A_list[[i]]) + ridge * diag(K)
    S[i, ] <- solve(A_i, local$rhs[i, ])
  }
  S
}

.cgamfa_update_scores_graph <- function(Y,
                                        B,
                                        Z_list,
                                        row_index,
                                        alpha_y,
                                        coupling_lambda,
                                        score_graph_laplacian,
                                        score_graph_lambda,
                                        local = NULL,
                                        ridge = 1e-8) {
  K <- ncol(B)
  if (is.null(local)) {
    local <- .cgamfa_score_system(
      Y = Y,
      B = B,
      Z_list = Z_list,
      row_index = row_index,
      alpha_y = alpha_y,
      coupling_lambda = coupling_lambda
    )
  }

  apply_A <- function(S) {
    .muscal_apply_row_operator(
      S = S,
      A_list = local$A_list,
      graph_laplacian = score_graph_laplacian,
      graph_lambda = score_graph_lambda,
      ridge = ridge
    )
  }
  X0 <- .cgamfa_update_scores(
    Y = Y,
    B = B,
    Z_list = Z_list,
    row_index = row_index,
    alpha_y = alpha_y,
    coupling_lambda = coupling_lambda,
    ridge = ridge
  )
  .muscal_matrix_cg(
    rhs = local$rhs,
    apply_A = apply_A,
    X0 = X0,
    tol = 1e-8,
    max_iter = max(50L, 5L * K)
  )
}

.cgamfa_score_gradient <- function(S,
                                   Y,
                                   B,
                                   Z_list,
                                   row_index,
                                   alpha_y,
                                   coupling_lambda,
                                   score_graph_laplacian,
                                   score_graph_lambda,
                                   ridge = 0) {
  local <- .cgamfa_score_system(
    Y = Y,
    B = B,
    Z_list = Z_list,
    row_index = row_index,
    alpha_y = alpha_y,
    coupling_lambda = coupling_lambda
  )
  grad <- matrix(0, nrow = nrow(S), ncol = ncol(S))
  for (i in seq_len(nrow(S))) {
    grad[i, ] <- 2 * (as.matrix(local$A_list[[i]]) %*% S[i, ] - local$rhs[i, ])
  }
  if (isTRUE(score_graph_lambda > 0) && !is.null(score_graph_laplacian)) {
    grad <- grad + 2 * score_graph_lambda * as.matrix(score_graph_laplacian %*% S)
  }
  if (ridge > 0) {
    grad <- grad + 2 * ridge * S
  }
  grad
}

.cgamfa_objective <- function(Y,
                              S,
                              B,
                              X_list,
                              Z_list,
                              V_list,
                              row_index,
                              alpha_y,
                              alpha_blocks,
                              coupling_lambda,
                              graph_laplacian,
                              graph_lambda,
                              score_graph_laplacian,
                              score_graph_lambda,
                              ridge = 0) {
  err_y <- alpha_y * .lmfa_reconstruction_sse(Y, S, B)
  err_x <- 0
  coupling_pen <- 0
  for (k in seq_along(X_list)) {
    err_x <- err_x + alpha_blocks[k] * .lmfa_reconstruction_sse(X_list[[k]], Z_list[[k]], V_list[[k]])
    dev <- Z_list[[k]] - S[row_index[[k]], , drop = FALSE]
    coupling_pen <- coupling_pen + coupling_lambda * sum(dev^2)
  }

  pen <- 0
  if (isTRUE(graph_lambda > 0) && !is.null(graph_laplacian)) {
    V_big <- do.call(rbind, V_list)
    LV <- graph_laplacian %*% V_big
    pen <- graph_lambda * sum(LV * V_big)
  }

  if (isTRUE(score_graph_lambda > 0) && !is.null(score_graph_laplacian)) {
    LS <- score_graph_laplacian %*% S
    pen <- pen + score_graph_lambda * sum(LS * S)
  }

  ridge_pen <- ridge * (
    sum(S^2) +
      sum(B^2) +
      sum(vapply(V_list, function(V) sum(V^2), numeric(1))) +
      sum(vapply(Z_list, function(Z) sum(Z^2), numeric(1)))
  )

  err_y + err_x + coupling_pen + pen + ridge_pen
}

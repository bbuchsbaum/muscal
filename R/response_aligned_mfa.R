#' Response-Aligned Multiple Factor Analysis
#'
#' @description
#' `response_aligned_mfa()` learns a shared response-informed latent space for a
#' set of multiblock predictors `X` paired with blockwise multivariate responses
#' `Y`. Unlike [anchored_mfa()], it does **not** assume that rows of `Y` are
#' shared across blocks. Instead, each block `k` receives its own latent score
#' matrix `Z_k`, the `X_k` blocks reconstruct from `Z_k`, and every `Y_k`
#' projects through a common response loading matrix `B`.
#'
#' @details
#' ## Model
#' The fitted model minimizes
#' \deqn{
#' \sum_k \alpha_k \|X_k - Z_k V_k^\top\|_F^2 +
#' \sum_k \eta_k \|R_k^{1/2}(Y_k - Z_k B^\top)\|_F^2 +
#' \eta_0 \|R_0^{1/2}(Y^{(a)} - S B^\top)\|_F^2 +
#' \sum_k \mu_k \|M_k^{1/2}(Z_k - A_k S)\|_F^2
#' }
#' plus optional feature-side coupling on the rows of the block loading
#' matrices `V_k` and ridge stabilization terms on `Z_k`, `V_k`, `B`, and the
#' optional anchor scores `S`. In v1, feature coupling may be supplied either
#' through `feature_groups`/`feature_lambda` or through
#' `feature_graph`/`graph_lambda`, but not both simultaneously.
#'
#' Rowwise response weights `R_k` are diagonal matrices represented
#' internally by non-negative weight vectors `response_weights[[k]]`.
#'
#' ## Identifiability
#' The latent basis is fixed by a pooled score normalization rather than by
#' constraining `B`. After each iteration, the fitted block scores are rotated so
#' that the weighted pooled score stack has orthonormal columns. This yields a
#' score-side identification that remains available even when the response
#' dimension is smaller than the latent dimension. When anchor structure is
#' active, the anchor score matrix `S` is rotated in the same basis. A
#' deterministic sign convention is then applied componentwise so that refits do
#' not differ only by arbitrary sign flips.
#'
#' ## Preprocessing
#' Predictor blocks `X` are preprocessed blockwise using `preproc`, following the
#' same conventions as [aligned_mfa()]. Response blocks `Y` are preprocessed in a
#' **shared** response space: `response_preproc` is fit once on the row-wise
#' concatenation of all response blocks plus `anchor_response` when present, and
#' then applied to each response source.
#'
#' ## Prediction contract
#' By default, out-of-sample prediction uses only information available at
#' prediction time:
#' `x_new -> z_hat -> y_hat`, optionally refined by `new_anchor_map` when test
#' rows have anchor-state information. If `new_response` is supplied and
#' `conditional = TRUE` is requested in [project()] or [predict()], the latent
#' score solve performs **conditional completion** rather than pure prediction by
#' incorporating observed response information explicitly.
#'
#' @param Y A list of numeric matrices/data.frames. Each element `Y[[k]]` is
#'   `n_k x q` and must have the same number of columns across blocks.
#' @param X A list of numeric matrices/data.frames. Each element `X[[k]]` is
#'   `n_k x p_k`.
#' @param preproc A `multivarious` preprocessing pipeline (a `pre_processor` or
#'   `prepper`) or a list of them for the `X` blocks. If a list, it must have
#'   length `length(X)` and is applied to `X` in block order.
#' @param response_preproc A single `multivarious` preprocessing pipeline used to
#'   define the shared response space. It is fit once on the pooled response
#'   rows and then applied to every block-specific `Y_k`.
#' @param ncomp Integer number of latent components.
#' @param normalization Block weighting scheme for the `X` blocks. `"MFA"` uses
#'   inverse squared first singular values; `"None"` uses uniform weights;
#'   `"custom"` uses `alpha`.
#' @param alpha Optional numeric vector of per-block `X` weights (length
#'   `length(X)`) used when `normalization = "custom"`.
#' @param response_alpha Optional numeric scalar or vector of per-block response
#'   weights. Defaults to `1` for every block.
#' @param response_weights Optional rowwise response weights. One of:
#'   * `NULL` (all rows weighted equally),
#'   * a single non-negative scalar applied to every row of every block,
#'   * a numeric vector of length `length(X)` giving a constant row weight per
#'     block, or
#'   * a list mirroring `X`, where each element is a non-negative numeric vector
#'     of length `nrow(Y[[k]])`.
#' @param anchor_response Optional anchor-level multivariate response table with
#'   one row per anchor state. Requires active `anchor_map` information and must
#'   have the same number of columns as the blockwise responses.
#' @param anchor_response_alpha Optional non-negative scalar controlling the
#'   contribution of `anchor_response` to the shared response loading update.
#' @param anchor_response_weights Optional rowwise weights for
#'   `anchor_response`. May be `NULL`, a scalar, or a non-negative vector of
#'   length `nrow(anchor_response)`.
#' @param anchor_map Optional anchor structure. One of:
#'   * `NULL` (no anchor coupling),
#'   * a list of integer vectors with one element per block, where values in
#'     `1..N_anchor` identify anchor states and `NA` indicates an unanchored row,
#'   * or a list of non-negative matrices `A_k` whose rows encode soft anchor
#'     assignments for each block row.
#' @param anchor_weight Optional rowwise anchor weights. One of:
#'   * `NULL` (defaults to `1` for anchored rows and `0` for unanchored rows),
#'   * a single non-negative scalar,
#'   * a numeric vector of length `length(X)` giving a constant row weight per
#'     block, or
#'   * a list mirroring `X`, where each element is a non-negative numeric vector
#'     of length `nrow(X[[k]])`.
#' @param coupling_lambda Optional non-negative scalar or vector of per-block
#'   anchor-coupling strengths. When all values are zero, the anchor machinery is
#'   dropped.
#' @param feature_groups Feature prior specification for the `X` loadings. One
#'   of:
#'   * `NULL` (no feature prior),
#'   * `"colnames"` to group features with identical column names across blocks,
#'   * a `data.frame` with columns `block`, `feature`, `group` and optional
#'     `weight`.
#' @param feature_graph Optional graph specification for the `X` loadings. One
#'   of:
#'   * `NULL` (no graph penalty),
#'   * `"colnames"` to connect identically named features across blocks,
#'   * a `data.frame` with columns `block1`, `feature1`, `block2`, `feature2`
#'     and optional `weight`,
#'   * or a square sparse/dense matrix interpreted according to `graph_form`.
#'   This is an alternative to `feature_groups`, not an additional penalty.
#' @param feature_lambda Non-negative scalar controlling the strength of the
#'   feature-group shrinkage.
#' @param graph_lambda Non-negative scalar controlling the strength of the
#'   feature-graph penalty.
#' @param graph_form Interpretation of `feature_graph` when it is matrix-like,
#'   or the Laplacian construction used for edge-based inputs.
#' @param max_iter Maximum ALS iterations.
#' @param tol Relative convergence tolerance on the objective.
#' @param ridge Non-negative ridge stabilization applied to score and loading
#'   updates.
#' @param verbose Logical; if `TRUE`, prints iteration diagnostics.
#' @param ... Unused (reserved for future extensions).
#'
#' @return An object inheriting from `multivarious::multiblock_biprojector` with
#'   additional class `"response_aligned_mfa"`. The object stores block-specific
#'   scores in `Z_list`, shared response loadings in `B`, block loadings in
#'   `V_list`, rowwise response weights in `response_weights`, optional anchor
#'   scores in `S`, optional anchor-level responses in `anchor_response`, parsed
#'   anchor maps in `anchor_map`, rowwise anchor weights in `anchor_weight`,
#'   optional `anchor_response_fit` diagnostics, pooled score normalization
#'   weights in `score_weights`, and optional feature-graph metadata in
#'   `graph_laplacian`, `graph_adjacency`, `graph_lambda`, and `graph_form`.
#'   For compatibility with generic score accessors, the `s` slot stores a
#'   concatenated stack of `Z_list`, and `score_index` maps stacked score rows
#'   back to blocks. The `score_representation` field records that `s` is this
#'   stacked block-score view rather than a shared observation-level score table.
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' X1 <- matrix(rnorm(30 * 6), 30, 6)
#' X2 <- matrix(rnorm(28 * 5), 28, 5)
#' Y1 <- matrix(rnorm(30 * 3), 30, 3)
#' Y2 <- matrix(rnorm(28 * 3), 28, 3)
#'
#' fit <- response_aligned_mfa(
#'   Y = list(X1 = Y1, X2 = Y2),
#'   X = list(X1 = X1, X2 = X2),
#'   ncomp = 2
#' )
#'
#' pred <- predict(fit, X1[1:5, , drop = FALSE], block = "X1", type = "response")
#' stopifnot(nrow(pred) == 5)
#' }
#' @export
response_aligned_mfa <- function(Y,
                                 X,
                                 preproc = multivarious::center(),
                                 response_preproc = multivarious::center(),
                                 ncomp = 2,
                                 normalization = c("MFA", "None", "custom"),
                                 alpha = NULL,
                                 response_alpha = 1,
                                 response_weights = NULL,
                                 anchor_response = NULL,
                                 anchor_response_alpha = 1,
                                 anchor_response_weights = NULL,
                                 anchor_map = NULL,
                                 anchor_weight = NULL,
                                 coupling_lambda = 0,
                                 feature_groups = NULL,
                                 feature_graph = NULL,
                                 feature_lambda = 0,
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

  chk::chk_list(X)
  chk::chk_true(length(X) >= 2L)
  chk::chk_list(Y)
  chk::chk_equal(length(X), length(Y))
  ncomp <- as.integer(ncomp)
  chk::chk_integer(ncomp)
  chk::chk_gte(ncomp, 1)
  chk::chk_numeric(feature_lambda)
  chk::chk_gte(feature_lambda, 0)
  chk::chk_numeric(graph_lambda)
  chk::chk_gte(graph_lambda, 0)
  chk::chk_numeric(coupling_lambda)
  chk::chk_gte(coupling_lambda, 0)
  chk::chk_numeric(anchor_response_alpha)
  chk::chk_gte(anchor_response_alpha, 0)
  if (length(anchor_response_alpha) != 1L) {
    stop("`anchor_response_alpha` must be a single non-negative scalar.", call. = FALSE)
  }
  max_iter <- as.integer(max_iter)
  chk::chk_integer(max_iter)
  chk::chk_gte(max_iter, 1)
  chk::chk_numeric(tol)
  chk::chk_gte(tol, 0)
  chk::chk_numeric(ridge)
  chk::chk_gte(ridge, 0)
  chk::chk_flag(verbose)

  X <- lapply(X, function(x) {
    chk::chk_true(is.matrix(x) || is.data.frame(x))
    as.matrix(x)
  })
  Y <- lapply(Y, function(y) {
    chk::chk_true(is.matrix(y) || is.data.frame(y))
    as.matrix(y)
  })
  if (!is.null(anchor_response)) {
    chk::chk_true(is.matrix(anchor_response) || is.data.frame(anchor_response))
    anchor_response <- as.matrix(anchor_response)
  }

  if (is.null(names(X))) names(X) <- paste0("X", seq_along(X))
  if (is.null(names(Y))) {
    names(Y) <- names(X)
  } else {
    Y <- .lmfa_align_named_index_list(Y, names(X), what = "Y")
  }

  for (k in seq_along(X)) {
    chk::chk_equal(nrow(X[[k]]), nrow(Y[[k]]))
  }

  q_dims <- unique(vapply(Y, ncol, integer(1)))
  if (length(q_dims) != 1L) {
    stop("All response blocks in `Y` must have the same number of columns.", call. = FALSE)
  }
  if (!is.null(anchor_response) && ncol(anchor_response) != q_dims[[1]]) {
    stop("`anchor_response` must have the same number of columns as the blockwise responses.", call. = FALSE)
  }

  total_rows <- sum(vapply(X, nrow, integer(1)))
  if (total_rows < ncomp) {
    stop("`ncomp` must be <= total number of observed rows across blocks.", call. = FALSE)
  }

  if (is.list(response_preproc) &&
      !inherits(response_preproc, "pre_processor") &&
      !inherits(response_preproc, "prepper")) {
    stop("`response_preproc` must be a single shared preprocessing pipeline, not a list.", call. = FALSE)
  }

  if (is.list(preproc) && !inherits(preproc, "pre_processor") && !inherits(preproc, "prepper")) {
    chk::chk_equal(length(preproc), length(X))
  }
  prep_res <- prepare_block_preprocessors(X, preproc, check_consistent_ncol = FALSE)
  Xp <- prep_res$Xp
  proclist <- prep_res$proclist

  pooled_response_rows <- if (is.null(anchor_response)) {
    do.call(rbind, Y)
  } else {
    do.call(rbind, c(Y, list(anchor_response)))
  }
  response_preproc_fit <- multivarious::fit(
    multivarious::fresh(response_preproc),
    pooled_response_rows
  )
  Yp <- lapply(Y, function(y) multivarious::transform(response_preproc_fit, y))
  anchor_response_p <- if (is.null(anchor_response)) {
    NULL
  } else {
    multivarious::transform(response_preproc_fit, anchor_response)
  }

  alpha_blocks <- if (normalization == "custom") {
    if (is.null(alpha)) {
      stop("When normalization = 'custom', 'alpha' must be provided.", call. = FALSE)
    }
    alpha <- as.numeric(alpha)
    chk::chk_equal(length(alpha), length(X))
    if (any(!is.finite(alpha)) || any(alpha < 0)) {
      stop("'alpha' must be finite and non-negative.", call. = FALSE)
    }
    alpha
  } else if (normalization == "None") {
    rep(1, length(X))
  } else {
    first_sv <- function(mat) {
      mdim <- min(dim(mat))
      method <- if (mdim < 3) "base" else "svds"
      tryCatch(
        multivarious::svd_wrapper(mat, ncomp = 1, method = method)$sdev[1],
        error = function(e) multivarious::svd_wrapper(mat, ncomp = 1, method = "base")$sdev[1]
      )
    }
    sapply(Xp, function(M) 1 / (first_sv(M)^2 + ridge))
  }
  names(alpha_blocks) <- names(X)

  response_alpha <- .ramfa_parse_block_weights(
    response_alpha,
    lengths = vapply(Yp, nrow, integer(1)),
    block_names = names(X),
    what = "response_alpha"
  )
  response_weights <- .ramfa_parse_response_weights(
    response_weights = response_weights,
    lengths = vapply(Yp, nrow, integer(1)),
    block_names = names(X)
  )
  coupling_lambda <- .ramfa_parse_block_weights(
    coupling_lambda,
    lengths = vapply(Xp, nrow, integer(1)),
    block_names = names(X),
    what = "coupling_lambda"
  )
  anchor_spec <- .ramfa_parse_anchor_map(
    anchor_map = anchor_map,
    lengths = vapply(Xp, nrow, integer(1)),
    block_names = names(X),
    coupling_lambda = coupling_lambda,
    anchor_weight = anchor_weight
  )
  if (!is.null(anchor_response_p) && !isTRUE(anchor_spec$active)) {
    stop("`anchor_response` requires active anchor structure via `anchor_map` and positive `coupling_lambda`.", call. = FALSE)
  }
  anchor_response_weights <- .ramfa_parse_single_row_weights(
    weights = anchor_response_weights,
    n = if (is.null(anchor_response_p)) 0L else nrow(anchor_response_p),
    what = "anchor_response_weights"
  )
  if (!is.null(anchor_response_p) && anchor_spec$n_anchor != nrow(anchor_response_p)) {
    stop("`anchor_response` must have one row per anchor state implied by `anchor_map`.", call. = FALSE)
  }

  has_block_response <- any(response_alpha > 0) && any(vapply(response_weights, sum, numeric(1)) > 0)
  has_anchor_response <- !is.null(anchor_response_p) &&
    as.numeric(anchor_response_alpha) > 0 &&
    sum(anchor_response_weights) > 0
  if (!has_block_response && !has_anchor_response) {
    stop("At least one supervised response source must have positive weight.", call. = FALSE)
  }

  data_refit <- list(
    Y = Y,
    X = X,
    response_weights = response_weights,
    anchor_response = anchor_response,
    anchor_response_weights = anchor_response_weights,
    anchor_map = if (isTRUE(anchor_spec$active)) anchor_spec$A_list else NULL,
    anchor_weight = if (isTRUE(anchor_spec$active)) anchor_spec$weights else NULL
  )

  fg <- .lmfa_parse_feature_groups(feature_groups, Xp, feature_lambda)
  feature_lambda_eff <- fg$feature_lambda
  graph <- .gamfa_parse_feature_graph(
    feature_graph = feature_graph,
    X_list = Xp,
    graph_lambda = graph_lambda,
    graph_form = graph_form
  )
  if (isTRUE(fg$enabled) && isTRUE(feature_lambda_eff > 0) &&
      isTRUE(graph$enabled) && isTRUE(graph_lambda > 0)) {
    stop("`feature_groups` and `feature_graph` cannot be combined in v1; choose one feature-coupling penalty.", call. = FALSE)
  }

  core <- .muscal_fit_common_space_engine(
    engine = "response",
    X_list = Xp,
    ncomp = as.integer(ncomp),
    alpha_blocks = alpha_blocks,
    ridge = ridge,
    max_iter = max_iter,
    tol = tol,
    verbose = verbose,
    Y_list = Yp,
    response_alpha = response_alpha,
    response_weights = response_weights,
    anchor_response = anchor_response_p,
    anchor_response_alpha = anchor_response_alpha,
    anchor_response_weights = anchor_response_weights,
    anchor_map = if (isTRUE(anchor_spec$active)) anchor_spec$A_list else NULL,
    anchor_weight = if (isTRUE(anchor_spec$active)) anchor_spec$weights else NULL,
    coupling_lambda = coupling_lambda,
    n_anchor = anchor_spec$n_anchor,
    fg = fg,
    feature_lambda = feature_lambda_eff,
    graph = graph,
    graph_lambda = graph_lambda,
    score_graph_spec = list(enabled = FALSE, L = NULL, A = NULL),
    score_graph_lambda = 0
  )
  B <- core$B
  V_list <- core$V_list
  Z_list <- core$Z_list
  S <- core$S
  score_weights <- core$score_weights
  objective_trace <- core$objective_trace

  block_indices <- list()
  current <- 1L
  for (k in seq_along(Xp)) {
    block_indices[[k]] <- current:(current + ncol(Xp[[k]]) - 1L)
    current <- current + ncol(Xp[[k]])
  }
  names(block_indices) <- names(Xp)

  score_stack <- .ramfa_stack_scores(Z_list)
  score_index <- score_stack$score_index

  proc <- multivarious::concat_pre_processors(proclist, block_indices)
  v_concat <- do.call(rbind, V_list)
  s_stack <- score_stack$s
  sdev <- apply(s_stack, 2, stats::sd)

  cor_x_list <- lapply(seq_along(Xp), function(k) {
    Ck <- tryCatch(stats::cor(Xp[[k]], Z_list[[k]]), error = function(e) NULL)
    if (!is.null(Ck)) Ck[!is.finite(Ck)] <- 0
    Ck
  })
  cor_x_list <- cor_x_list[!vapply(cor_x_list, is.null, logical(1))]
  cor_loadings <- if (length(cor_x_list) > 0) do.call(rbind, cor_x_list) else NULL

  block_fit <- do.call(rbind, lapply(seq_along(Xp), function(k) {
    X_hat <- Z_list[[k]] %*% t(V_list[[k]])
    sse <- sum((Xp[[k]] - X_hat)^2)
    tss <- sum(Xp[[k]]^2)
    data.frame(
      block = names(Xp)[k],
      n = nrow(Xp[[k]]),
      p = ncol(Xp[[k]]),
      sse = sse,
      tss = tss,
      r2 = if (tss > 0) 1 - sse / tss else NA_real_
    )
  }))

  response_fit <- do.call(rbind, lapply(seq_along(Yp), function(k) {
    Y_hat <- Z_list[[k]] %*% t(B)
    resid <- Yp[[k]] - Y_hat
    w <- response_weights[[k]]
    sse <- sum(w * rowSums(resid^2))
    tss <- sum(w * rowSums(Yp[[k]]^2))
    data.frame(
      block = names(Yp)[k],
      n = nrow(Yp[[k]]),
      q = ncol(Yp[[k]]),
      sse = sse,
      tss = tss,
      r2 = if (tss > 0) 1 - sse / tss else NA_real_,
      mean_weight = mean(w)
    )
  }))
  anchor_response_fit <- if (!is.null(anchor_response_p)) {
    Y_hat <- S %*% t(B)
    resid <- anchor_response_p - Y_hat
    sse <- sum(anchor_response_weights * rowSums(resid^2))
    tss <- sum(anchor_response_weights * rowSums(anchor_response_p^2))
    data.frame(
      block = ".anchor",
      n = nrow(anchor_response_p),
      q = ncol(anchor_response_p),
      sse = sse,
      tss = tss,
      r2 = if (tss > 0) 1 - sse / tss else NA_real_,
      mean_weight = mean(anchor_response_weights)
    )
  } else {
    NULL
  }

  fit <- multivarious::multiblock_biprojector(
    v = v_concat,
    s = s_stack,
    sdev = sdev,
    preproc = proc,
    block_indices = block_indices,
    B = B,
    V_list = V_list,
    Z_list = Z_list,
    S = S,
    score_index = score_index,
    alpha_blocks = alpha_blocks,
    response_alpha = response_alpha,
    response_weights = response_weights,
    anchor_response = anchor_response_p,
    anchor_response_alpha = as.numeric(anchor_response_alpha),
    anchor_response_weights = anchor_response_weights,
    anchor_map = if (isTRUE(anchor_spec$active)) anchor_spec$A_list else NULL,
    anchor_weight = if (isTRUE(anchor_spec$active)) anchor_spec$weights else NULL,
    coupling_lambda = coupling_lambda,
    n_anchor = anchor_spec$n_anchor,
    graph_laplacian = graph$L,
    graph_adjacency = graph$A,
    graph_lambda = graph_lambda,
    graph_form = graph_form,
    score_weights = score_weights,
    normalization = normalization,
    feature_groups = feature_groups,
    feature_lambda = feature_lambda_eff,
    objective_trace = objective_trace,
    partial_scores = Z_list,
    cor_loadings = cor_loadings,
    block_fit = block_fit,
    response_fit = response_fit,
    anchor_response_fit = anchor_response_fit,
    block_preproc = setNames(proclist, names(Xp)),
    response_preproc = response_preproc_fit,
    ridge = ridge,
    score_representation = "stacked_block_scores",
    classes = "response_aligned_mfa"
  )

  .muscal_attach_fit_contract(
    fit,
    method = "response_aligned_mfa",
    task = "response_prediction",
    oos_types = c("response", "scores", "reconstruction"),
    fit_call = fit_call,
    refit_supported = TRUE,
    prediction_target = "Y",
    refit = .muscal_make_refit_spec(
      data = data_refit,
      fit_fn = function(data) {
        do.call(
          response_aligned_mfa,
          c(
            list(
              Y = data$Y,
              X = data$X,
              preproc = preproc,
              response_preproc = response_preproc,
              ncomp = ncomp,
              normalization = normalization,
              alpha = alpha,
              response_alpha = response_alpha,
              response_weights = data$response_weights,
              anchor_response = data$anchor_response,
              anchor_response_alpha = anchor_response_alpha,
              anchor_response_weights = data$anchor_response_weights,
              anchor_map = data$anchor_map,
              anchor_weight = data$anchor_weight,
              coupling_lambda = coupling_lambda,
              feature_groups = feature_groups,
              feature_graph = feature_graph,
              feature_lambda = feature_lambda,
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
      bootstrap_fn = .muscal_bootstrap_response_aligned_data,
      permutation_fn = .muscal_permute_response_aligned_data,
      resample_unit = "block_rows"
    )
  )
}

.ramfa_parse_block_weights <- function(weights, lengths, block_names, what) {
  K <- length(lengths)
  if (length(block_names) != K) {
    stop("Internal error: block name/length mismatch.", call. = FALSE)
  }

  if (length(weights) == 1L) {
    out <- rep(as.numeric(weights), K)
  } else {
    out <- as.numeric(weights)
    chk::chk_equal(length(out), K)
  }

  if (any(!is.finite(out)) || any(out < 0)) {
    stop(sprintf("`%s` must be finite and non-negative.", what), call. = FALSE)
  }
  stats::setNames(out, block_names)
}

.ramfa_parse_response_weights <- function(response_weights, lengths, block_names) {
  K <- length(lengths)
  if (is.null(response_weights)) {
    out <- lapply(lengths, function(n) rep(1, n))
    names(out) <- block_names
    return(out)
  }

  if (is.list(response_weights) &&
      !inherits(response_weights, "pre_processor") &&
      !inherits(response_weights, "prepper")) {
    if (length(response_weights) != K) {
      stop("`response_weights` list must have one element per block.", call. = FALSE)
    }
    if (is.null(names(response_weights))) {
      names(response_weights) <- block_names
    } else {
      response_weights <- .lmfa_align_named_index_list(response_weights, block_names, what = "response_weights")
    }
    out <- lapply(seq_along(response_weights), function(k) {
      w <- as.numeric(response_weights[[k]])
      chk::chk_equal(length(w), lengths[[k]])
      if (any(!is.finite(w)) || any(w < 0)) {
        stop("Each `response_weights[[k]]` must be finite and non-negative.", call. = FALSE)
      }
      w
    })
    names(out) <- block_names
    return(out)
  }

  w <- as.numeric(response_weights)
  if (length(w) == 1L) {
    out <- lapply(lengths, function(n) rep(w, n))
  } else if (length(w) == K) {
    out <- lapply(seq_along(lengths), function(k) rep(w[[k]], lengths[[k]]))
  } else {
    stop("`response_weights` must be NULL, a scalar, a length-K vector, or a blockwise list.", call. = FALSE)
  }

  if (any(!is.finite(unlist(out, use.names = FALSE))) || any(unlist(out, use.names = FALSE) < 0)) {
    stop("`response_weights` must be finite and non-negative.", call. = FALSE)
  }
  names(out) <- block_names
  out
}

.ramfa_parse_single_row_weights <- function(weights, n, what) {
  if (n < 1L) {
    return(NULL)
  }
  if (is.null(weights)) {
    return(rep(1, n))
  }
  weights <- as.numeric(weights)
  if (length(weights) == 1L) {
    weights <- rep(weights, n)
  }
  chk::chk_equal(length(weights), n)
  if (any(!is.finite(weights)) || any(weights < 0)) {
    stop(sprintf("`%s` must be finite and non-negative.", what), call. = FALSE)
  }
  weights
}

.ramfa_parse_anchor_map <- function(anchor_map,
                                    lengths,
                                    block_names,
                                    coupling_lambda,
                                    anchor_weight = NULL) {
  if (is.null(anchor_map) || all(coupling_lambda <= 0)) {
    if (!is.null(anchor_map) && all(coupling_lambda <= 0)) {
      return(list(
        active = FALSE,
        A_list = NULL,
        weights = NULL,
        n_anchor = 0L
      ))
    }
    return(list(
      active = FALSE,
      A_list = NULL,
      weights = NULL,
      n_anchor = 0L
    ))
  }

  K <- length(lengths)
  chk::chk_list(anchor_map)
  chk::chk_equal(length(anchor_map), K)
  if (is.null(names(anchor_map))) {
    names(anchor_map) <- block_names
  } else {
    anchor_map <- .lmfa_align_named_index_list(anchor_map, block_names, what = "anchor_map")
  }

  n_anchor <- 0L
  A_list <- vector("list", K)
  anchored_default <- vector("list", K)

  for (k in seq_along(anchor_map)) {
    Ak <- anchor_map[[k]]
    nk <- lengths[[k]]
    if (is.matrix(Ak) || is.data.frame(Ak)) {
      A_mat <- as.matrix(Ak)
      chk::chk_equal(nrow(A_mat), nk)
      if (any(!is.finite(A_mat)) || any(A_mat < 0)) {
        stop("Matrix `anchor_map` entries must be finite and non-negative.", call. = FALSE)
      }
      row_totals <- rowSums(A_mat)
      if (any(row_totals > 0 & abs(row_totals - 1) > 1e-8)) {
        stop("Rows of matrix `anchor_map` with positive mass must sum to 1.", call. = FALSE)
      }
      n_anchor <- max(n_anchor, ncol(A_mat))
      A_list[[k]] <- A_mat
      anchored_default[[k]] <- as.numeric(row_totals > 0)
    } else if (is.numeric(Ak) || is.integer(Ak)) {
      idx <- as.integer(Ak)
      chk::chk_equal(length(idx), nk)
      bad <- !is.na(idx) & idx < 1L
      if (any(bad)) {
        stop("Integer `anchor_map` entries must be positive or NA.", call. = FALSE)
      }
      if (any(!is.na(idx))) {
        n_anchor <- max(n_anchor, max(idx, na.rm = TRUE))
      }
      A_list[[k]] <- idx
      anchored_default[[k]] <- as.numeric(!is.na(idx))
    } else {
      stop("Each `anchor_map[[k]]` must be an integer vector or numeric matrix.", call. = FALSE)
    }
  }

  if (n_anchor < 1L) {
    stop("`anchor_map` did not define any anchor states.", call. = FALSE)
  }

  A_list <- lapply(A_list, function(Ak) {
    if (is.integer(Ak)) {
      out <- matrix(0, nrow = length(Ak), ncol = n_anchor)
      keep <- which(!is.na(Ak))
      if (length(keep) > 0L) {
        out[cbind(keep, Ak[keep])] <- 1
      }
      out
    } else {
      if (ncol(Ak) == n_anchor) {
        Ak
      } else {
        out <- matrix(0, nrow = nrow(Ak), ncol = n_anchor)
        out[, seq_len(ncol(Ak))] <- Ak
        out
      }
    }
  })
  names(A_list) <- block_names

  weights <- .ramfa_parse_optional_row_weights(
    weights = anchor_weight,
    lengths = lengths,
    block_names = block_names,
    default = anchored_default,
    what = "anchor_weight"
  )

  for (k in seq_along(A_list)) {
    row_totals <- rowSums(A_list[[k]])
    if (any(weights[[k]] > 0 & row_totals <= 0)) {
      stop("Rows with positive `anchor_weight` must have a non-zero `anchor_map` row.", call. = FALSE)
    }
  }

  list(
    active = TRUE,
    A_list = A_list,
    weights = weights,
    n_anchor = n_anchor
  )
}

.ramfa_parse_optional_row_weights <- function(weights,
                                              lengths,
                                              block_names,
                                              default,
                                              what) {
  K <- length(lengths)
  if (is.null(weights)) {
    out <- lapply(default, as.numeric)
    names(out) <- block_names
    return(out)
  }

  if (is.list(weights) &&
      !inherits(weights, "pre_processor") &&
      !inherits(weights, "prepper")) {
    chk::chk_equal(length(weights), K)
    if (is.null(names(weights))) {
      names(weights) <- block_names
    } else {
      weights <- .lmfa_align_named_index_list(weights, block_names, what = what)
    }
    out <- lapply(seq_along(weights), function(k) {
      w <- as.numeric(weights[[k]])
      chk::chk_equal(length(w), lengths[[k]])
      if (any(!is.finite(w)) || any(w < 0)) {
        stop(sprintf("Each `%s[[k]]` must be finite and non-negative.", what), call. = FALSE)
      }
      w
    })
    names(out) <- block_names
    return(out)
  }

  w <- as.numeric(weights)
  if (length(w) == 1L) {
    out <- lapply(lengths, function(n) rep(w, n))
  } else if (length(w) == K) {
    out <- lapply(seq_along(lengths), function(k) rep(w[[k]], lengths[[k]]))
  } else {
    stop(sprintf("`%s` must be NULL, a scalar, a length-K vector, or a blockwise list.", what), call. = FALSE)
  }

  if (any(!is.finite(unlist(out, use.names = FALSE))) || any(unlist(out, use.names = FALSE) < 0)) {
    stop(sprintf("`%s` must be finite and non-negative.", what), call. = FALSE)
  }
  names(out) <- block_names
  out
}

.ramfa_initialize_response_loadings <- function(Y_list, K) {
  q <- ncol(Y_list[[1]])
  Y_pool <- do.call(rbind, Y_list)
  B <- matrix(0, nrow = q, ncol = K)
  k_resp <- min(K, nrow(Y_pool), ncol(Y_pool))
  if (k_resp >= 1L) {
    sv <- multivarious::svd_wrapper(
      Y_pool,
      ncomp = k_resp,
      method = if (min(dim(Y_pool)) < 3L) "base" else "svds"
    )
    B[, seq_len(k_resp)] <- sv$v[, seq_len(k_resp), drop = FALSE]
  }
  if (k_resp < K) {
    B[, seq.int(k_resp + 1L, K)] <- matrix(stats::rnorm(q * (K - k_resp), sd = 0.1), nrow = q)
  }
  B
}

.ramfa_initialize_score_matrix <- function(response_mat,
                                           predictor_mat = NULL,
                                           B,
                                           K) {
  Z <- matrix(0, nrow = nrow(response_mat), ncol = K)
  use_resp <- min(K, ncol(B))
  if (use_resp >= 1L) {
    Z[, seq_len(use_resp)] <- response_mat %*% B[, seq_len(use_resp), drop = FALSE]
  }
  if (use_resp < K) {
    extra <- K - use_resp
    take <- if (is.null(predictor_mat)) 0L else min(extra, min(dim(predictor_mat)))
    if (take >= 1L) {
      sv <- multivarious::svd_wrapper(
        predictor_mat,
        ncomp = take,
        method = if (min(dim(predictor_mat)) < 3L) "base" else "svds"
      )
      Z[, seq.int(use_resp + 1L, length.out = take)] <-
        sv$u[, seq_len(take), drop = FALSE] %*%
        diag(sv$sdev[seq_len(take)], take, take)
    }
    if (take < extra) {
      cols <- seq.int(use_resp + take + 1L, K)
      Z[, cols] <- matrix(stats::rnorm(nrow(Z) * length(cols), sd = 0.1), nrow = nrow(Z))
    }
  }
  Z
}

.ramfa_initialize_scores <- function(X_list, Y_list, B, K) {
  Z_list <- lapply(seq_along(Y_list), function(k) {
    .ramfa_initialize_score_matrix(
      response_mat = Y_list[[k]],
      predictor_mat = X_list[[k]],
      B = B,
      K = K
    )
  })
  names(Z_list) <- names(X_list)
  Z_list
}

.ramfa_initialize_anchor_scores <- function(Z_list,
                                            anchor_map,
                                            anchor_weight,
                                            n_anchor,
                                            ridge = 1e-8) {
  .ramfa_update_anchor_scores(
    Z_list = Z_list,
    anchor_map = anchor_map,
    anchor_weight = anchor_weight,
    coupling_lambda = rep(1, length(Z_list)),
    n_anchor = n_anchor,
    ridge = ridge
  )
}

.ramfa_score_weights <- function(alpha_blocks, response_alpha, response_weights) {
  response_means <- vapply(response_weights, mean, numeric(1))
  weights <- unname(alpha_blocks) + unname(response_alpha) * response_means
  if (all(weights <= 0)) {
    stop("At least one block must contribute positive weight to pooled score normalization.", call. = FALSE)
  }
  stats::setNames(weights, names(alpha_blocks))
}

.ramfa_stack_scores <- function(Z_list) {
  score_index <- vector("list", length(Z_list))
  current <- 1L
  for (k in seq_along(Z_list)) {
    score_index[[k]] <- current:(current + nrow(Z_list[[k]]) - 1L)
    current <- current + nrow(Z_list[[k]])
  }
  names(score_index) <- names(Z_list)
  list(
    s = do.call(rbind, Z_list),
    score_index = score_index
  )
}

.ramfa_orthonormalize_score_list <- function(Z_list, weights, ridge = 1e-8) {
  K <- ncol(Z_list[[1]])
  Zw <- do.call(rbind, lapply(seq_along(Z_list), function(k) {
    sqrt(max(as.numeric(weights[[k]]), 0)) * Z_list[[k]]
  }))
  sv <- svd(Zw, nu = K, nv = K)
  d <- sv$d[seq_len(K)]
  safe_d <- pmax(d, sqrt(ridge))
  transform <- sv$v[, seq_len(K), drop = FALSE] %*% diag(1 / safe_d, K, K)
  loading_rot <- sv$v[, seq_len(K), drop = FALSE] %*% diag(safe_d, K, K)
  list(
    Z_list = lapply(Z_list, function(Z) Z %*% transform),
    transform = transform,
    loading_rot = loading_rot
  )
}

.ramfa_component_signs <- function(B,
                                   V_list = NULL,
                                   S = NULL,
                                   Z_list = NULL,
                                   tol = 1e-12) {
  K <- ncol(B)
  signs <- rep(1, K)

  dominant_sign <- function(x) {
    if (is.null(x)) return(NA_real_)
    x <- as.numeric(x)
    keep <- which(is.finite(x) & abs(x) > tol)
    if (length(keep) == 0L) return(NA_real_)
    s <- sign(x[keep[[which.max(abs(x[keep]))]]])
    if (!is.finite(s) || s == 0) NA_real_ else s
  }

  for (j in seq_len(K)) {
    sign_j <- dominant_sign(B[, j, drop = TRUE])
    if (is.na(sign_j) && !is.null(V_list)) {
      sign_j <- dominant_sign(do.call(c, lapply(V_list, function(V) V[, j, drop = TRUE])))
    }
    if (is.na(sign_j) && !is.null(S)) {
      sign_j <- dominant_sign(S[, j, drop = TRUE])
    }
    if (is.na(sign_j) && !is.null(Z_list)) {
      sign_j <- dominant_sign(do.call(c, lapply(Z_list, function(Z) Z[, j, drop = TRUE])))
    }
    signs[[j]] <- if (is.na(sign_j)) 1 else sign_j
  }

  signs
}

.ramfa_orient_components <- function(B,
                                     V_list = NULL,
                                     Z_list = NULL,
                                     S = NULL) {
  signs <- .ramfa_component_signs(
    B = B,
    V_list = V_list,
    S = S,
    Z_list = Z_list
  )
  D <- diag(signs, nrow = length(signs), ncol = length(signs))
  list(
    B = B %*% D,
    V_list = if (is.null(V_list)) NULL else lapply(V_list, function(V) V %*% D),
    Z_list = if (is.null(Z_list)) NULL else lapply(Z_list, function(Z) Z %*% D),
    S = if (is.null(S)) NULL else S %*% D
  )
}

.lmfa_feature_group_penalty <- function(V_list,
                                        fg,
                                        centers,
                                        feature_lambda) {
  if (!isTRUE(feature_lambda > 0) || !isTRUE(fg$enabled) || is.null(centers) || length(fg$members) == 0) {
    return(0)
  }

  pen <- 0
  groups <- names(fg$members)
  if (is.null(groups)) groups <- seq_along(fg$members)
  for (gi in seq_along(fg$members)) {
    gname <- groups[gi]
    cvec <- centers$centers[match(gname, rownames(centers$centers)), , drop = TRUE]
    for (m in fg$members[[gi]]) {
      v <- V_list[[m$block]][m$feature, , drop = TRUE]
      pen <- pen + m$weight * sum((v - cvec)^2)
    }
  }

  feature_lambda * pen
}

.muscal_fit_common_space_engine <- function(engine = c("response", "hard_anchor", "coupled_anchor"),
                                            X_list,
                                            ncomp,
                                            alpha_blocks,
                                            ridge,
                                            max_iter,
                                            tol,
                                            verbose,
                                            Y_list = NULL,
                                            response_alpha = NULL,
                                            response_weights = NULL,
                                            anchor_response = NULL,
                                            anchor_response_alpha = 0,
                                            anchor_response_weights = NULL,
                                            anchor_map = NULL,
                                            anchor_weight = NULL,
                                            coupling_lambda = 0,
                                            n_anchor = 0L,
                                            row_index = NULL,
                                            fg = NULL,
                                            feature_lambda = 0,
                                            graph = NULL,
                                            graph_lambda = 0,
                                            score_graph_spec = NULL,
                                            score_graph_lambda = 0,
                                            score_constraint = c("none", "orthonormal")) {
  engine <- match.arg(engine)
  score_constraint <- match.arg(score_constraint)
  if (is.null(graph)) {
    graph <- list(enabled = FALSE, L = NULL, A = NULL)
  }
  if (is.null(score_graph_spec)) {
    score_graph_spec <- list(enabled = FALSE, L = NULL, A = NULL)
  }
  if (is.null(fg)) {
    fg <- list(
      enabled = FALSE,
      by_block = lapply(X_list, function(X) rep.int(NA_character_, ncol(X))),
      weights_by_block = lapply(X_list, function(X) rep.int(1, ncol(X))),
      members = list()
    )
  }

  if (identical(engine, "response")) {
    K <- as.integer(ncomp)
    init_response_list <- Y_list
    if (!is.null(anchor_response)) {
      init_response_list <- c(init_response_list, list(.anchor = anchor_response))
    }
    B <- .ramfa_initialize_response_loadings(init_response_list, K = K)
    Z_list <- .ramfa_initialize_scores(X_list = X_list, Y_list = Y_list, B = B, K = K)
    S <- if (!is.null(anchor_map)) {
      S0 <- .ramfa_initialize_anchor_scores(
        Z_list = Z_list,
        anchor_map = anchor_map,
        anchor_weight = anchor_weight,
        n_anchor = n_anchor,
        ridge = ridge
      )
      if (!is.null(anchor_response) && as.numeric(anchor_response_alpha) > 0 && sum(anchor_response_weights) > 0) {
        S_resp <- .ramfa_initialize_score_matrix(
          response_mat = anchor_response,
          predictor_mat = NULL,
          B = B,
          K = K
        )
        S0 <- (S0 + S_resp) / 2
      }
      S0
    } else {
      NULL
    }
    score_weights <- .ramfa_score_weights(
      alpha_blocks = alpha_blocks,
      response_alpha = response_alpha,
      response_weights = response_weights
    )
    ortho <- .ramfa_orthonormalize_score_list(Z_list, weights = score_weights, ridge = ridge)
    Z_list <- ortho$Z_list
    if (!is.null(S)) {
      S <- S %*% ortho$transform
    }
    B <- B %*% ortho$loading_rot
    oriented <- .ramfa_orient_components(B = B, Z_list = Z_list, S = S)
    B <- oriented$B
    Z_list <- oriented$Z_list
    S <- oriented$S

    V_list <- lapply(seq_along(X_list), function(k) {
      .lmfa_update_loadings(Z_list[[k]], X_list[[k]], alpha = alpha_blocks[k], ridge = ridge)$V
    })
    names(V_list) <- names(X_list)
    Z_list <- setNames(Z_list, names(X_list))

    objective_trace <- numeric(0)
    prev_obj <- Inf

    for (iter in seq_len(max_iter)) {
      B <- .ramfa_update_response_loadings(
        Z_list = Z_list,
        Y_list = Y_list,
        response_alpha = response_alpha,
        response_weights = response_weights,
        S_anchor = S,
        anchor_response = anchor_response,
        anchor_response_alpha = anchor_response_alpha,
        anchor_response_weights = anchor_response_weights,
        ridge = ridge
      )

      if (isTRUE(graph$enabled) && isTRUE(graph_lambda > 0)) {
        V_list <- .gamfa_update_V_graph(
          S_list = Z_list,
          X_list = X_list,
          alpha_blocks = alpha_blocks,
          graph_laplacian = graph$L,
          graph_lambda = graph_lambda,
          ridge = ridge
        )
      } else {
        centers <- .lmfa_group_centers(V_list, fg, block_weights = alpha_blocks)
        for (k in seq_along(X_list)) {
          V_list[[k]] <- .lmfa_update_V_block(
            Sk = Z_list[[k]],
            Xk = X_list[[k]],
            alpha_block = alpha_blocks[k],
            fg_block = fg$by_block[[k]],
            weights_block = fg$weights_by_block[[k]],
            centers = centers,
            feature_lambda = feature_lambda,
            ridge = ridge
          )
        }
      }

      Z_list <- lapply(seq_along(X_list), function(k) {
        .ramfa_update_scores_block(
          Xk = X_list[[k]],
          Yk = Y_list[[k]],
          Vk = V_list[[k]],
          B = B,
          alpha_block = alpha_blocks[k],
          response_alpha = response_alpha[k],
          response_weight = response_weights[[k]],
          anchor_target = if (!is.null(S) && !is.null(anchor_map)) anchor_map[[k]] %*% S else NULL,
          coupling_lambda = coupling_lambda[k],
          anchor_weight = if (!is.null(anchor_weight)) anchor_weight[[k]] else NULL,
          ridge = ridge
        )
      })
      names(Z_list) <- names(X_list)

      if (!is.null(S)) {
        S <- .ramfa_update_anchor_scores(
          Z_list = Z_list,
          anchor_map = anchor_map,
          anchor_weight = anchor_weight,
          coupling_lambda = coupling_lambda,
          n_anchor = n_anchor,
          B = B,
          anchor_response = anchor_response,
          anchor_response_alpha = anchor_response_alpha,
          anchor_response_weights = anchor_response_weights,
          score_graph_laplacian = score_graph_spec$L,
          score_graph_lambda = score_graph_lambda,
          ridge = ridge
        )
      }

      ortho <- .ramfa_orthonormalize_score_list(Z_list, weights = score_weights, ridge = ridge)
      Z_list <- ortho$Z_list
      if (!is.null(S)) {
        S <- S %*% ortho$transform
      }
      B <- B %*% ortho$loading_rot
      V_list <- lapply(V_list, function(V) V %*% ortho$loading_rot)
      oriented <- .ramfa_orient_components(
        B = B,
        V_list = V_list,
        Z_list = Z_list,
        S = S
      )
      B <- oriented$B
      V_list <- oriented$V_list
      Z_list <- oriented$Z_list
      S <- oriented$S

      centers <- if (isTRUE(graph$enabled) && isTRUE(graph_lambda > 0)) {
        NULL
      } else {
        .lmfa_group_centers(V_list, fg, block_weights = alpha_blocks)
      }
      obj <- .ramfa_objective(
        Y_list = Y_list,
        Z_list = Z_list,
        B = B,
        X_list = X_list,
        V_list = V_list,
        alpha_blocks = alpha_blocks,
        response_alpha = response_alpha,
        response_weights = response_weights,
        anchor_response = anchor_response,
        anchor_response_alpha = anchor_response_alpha,
        anchor_response_weights = anchor_response_weights,
        S = S,
        anchor_map = anchor_map,
        anchor_weight = anchor_weight,
        coupling_lambda = coupling_lambda,
        graph_laplacian = graph$L,
        graph_lambda = graph_lambda,
        score_graph_laplacian = score_graph_spec$L,
        score_graph_lambda = score_graph_lambda,
        fg = fg,
        centers = centers,
        feature_lambda = feature_lambda,
        ridge = ridge
      )
      objective_trace <- c(objective_trace, obj)

      rel_change <- abs(obj - prev_obj) / (abs(prev_obj) + ridge)
      if (isTRUE(verbose)) {
        message(sprintf("response_aligned_mfa iter %d: obj=%.6g, rel_change=%.3g", iter, obj, rel_change))
      }
      if (is.finite(prev_obj) && rel_change < tol) break
      prev_obj <- obj
    }

    return(list(
      S = S,
      B = B,
      V_list = V_list,
      Z_list = Z_list,
      objective_trace = objective_trace,
      score_weights = score_weights
    ))
  }

  alpha_y <- as.numeric(anchor_response_alpha)
  anchor_weights <- anchor_response_weights
  if (is.null(anchor_weights)) {
    anchor_weights <- rep(1, nrow(anchor_response))
  }
  init <- .lmfa_init_from_Y(anchor_response, min(as.integer(ncomp), nrow(anchor_response), ncol(anchor_response)), ridge = ridge)
  S <- init$S
  B <- init$B
  Z_list <- lapply(seq_along(X_list), function(k) S[row_index[[k]], , drop = FALSE])
  names(Z_list) <- names(X_list)
  V_list <- lapply(seq_along(X_list), function(k) {
    .lmfa_update_loadings(Z_list[[k]], X_list[[k]], alpha = alpha_blocks[k], ridge = ridge)$V
  })
  names(V_list) <- names(X_list)

  objective_trace <- numeric(0)
  prev_obj <- Inf

  for (iter in seq_len(max_iter)) {
    B <- .ramfa_update_shared_response_loadings(
      S_anchor = S,
      anchor_response = anchor_response,
      anchor_response_alpha = alpha_y,
      anchor_response_weights = anchor_weights,
      ridge = ridge
    )

    score_blocks <- if (identical(engine, "coupled_anchor")) {
      Z_list
    } else {
      lapply(seq_along(X_list), function(k) S[row_index[[k]], , drop = FALSE])
    }
    names(score_blocks) <- names(X_list)

    if (isTRUE(graph$enabled) && isTRUE(graph_lambda > 0)) {
      V_list <- .gamfa_update_V_graph(
        S_list = score_blocks,
        X_list = X_list,
        alpha_blocks = alpha_blocks,
        graph_laplacian = graph$L,
        graph_lambda = graph_lambda,
        ridge = ridge
      )
    } else {
      centers <- .lmfa_group_centers(V_list, fg, block_weights = alpha_blocks)
      for (k in seq_along(X_list)) {
        V_list[[k]] <- .lmfa_update_V_block(
          Sk = score_blocks[[k]],
          Xk = X_list[[k]],
          alpha_block = alpha_blocks[k],
          fg_block = fg$by_block[[k]],
          weights_block = fg$weights_by_block[[k]],
          centers = centers,
          feature_lambda = feature_lambda,
          ridge = ridge
        )
      }
    }

    if (identical(engine, "coupled_anchor")) {
      Z_list <- .cgamfa_update_Z_list(
        S = S,
        X_list = X_list,
        V_list = V_list,
        row_index = row_index,
        alpha_blocks = alpha_blocks,
        coupling_lambda = coupling_lambda,
        ridge = ridge
      )
    }

    if (identical(score_constraint, "orthonormal")) {
      local <- if (identical(engine, "coupled_anchor")) {
        .cgamfa_score_system(
          Y = anchor_response,
          B = B,
          Z_list = Z_list,
          row_index = row_index,
          alpha_y = alpha_y,
          coupling_lambda = coupling_lambda
        )
      } else {
        .gamfa_score_system(
          Y = anchor_response,
          B = B,
          X_list = X_list,
          V_list = V_list,
          row_index = row_index,
          alpha_y = alpha_y,
          alpha_blocks = alpha_blocks
        )
      }

      S_warm <- if (isTRUE(score_graph_spec$enabled) && isTRUE(score_graph_lambda > 0)) {
        if (identical(engine, "coupled_anchor")) {
          .cgamfa_update_scores_graph(
            Y = anchor_response,
            B = B,
            Z_list = Z_list,
            row_index = row_index,
            alpha_y = alpha_y,
            coupling_lambda = coupling_lambda,
            score_graph_laplacian = score_graph_spec$L,
            score_graph_lambda = score_graph_lambda,
            local = local,
            ridge = ridge
          )
        } else {
          .gamfa_update_scores_graph(
            Y = anchor_response,
            B = B,
            X_list = X_list,
            V_list = V_list,
            row_index = row_index,
            alpha_y = alpha_y,
            alpha_blocks = alpha_blocks,
            score_graph_laplacian = score_graph_spec$L,
            score_graph_lambda = score_graph_lambda,
            local = local,
            ridge = ridge
          )
        }
      } else {
        if (identical(engine, "coupled_anchor")) {
          .cgamfa_update_scores(
            Y = anchor_response,
            B = B,
            Z_list = Z_list,
            row_index = row_index,
            alpha_y = alpha_y,
            coupling_lambda = coupling_lambda,
            local = local,
            ridge = ridge
          )
        } else {
          .lmfa_update_scores(
            Y = anchor_response,
            B = B,
            X_list = X_list,
            V_list = V_list,
            row_index = row_index,
            alpha_y = alpha_y,
            alpha_blocks = alpha_blocks,
            local = local,
            ridge = ridge
          )
        }
      }

      S_start <- .muscal_stiefel_retract(S)
      S_warm <- .muscal_stiefel_retract(S_warm)
      obj_start <- .muscal_score_objective_from_system(
        S_start,
        local$A_list,
        local$rhs,
        graph_laplacian = score_graph_spec$L,
        graph_lambda = score_graph_lambda,
        ridge = ridge
      )
      obj_warm <- .muscal_score_objective_from_system(
        S_warm,
        local$A_list,
        local$rhs,
        graph_laplacian = score_graph_spec$L,
        graph_lambda = score_graph_lambda,
        ridge = ridge
      )
      if (obj_warm <= obj_start) {
        S_start <- S_warm
      }
      s_opt <- .muscal_stiefel_mm(
        S = S_start,
        A_list = local$A_list,
        rhs = local$rhs,
        graph_laplacian = score_graph_spec$L,
        graph_lambda = score_graph_lambda,
        ridge = ridge
      )
      S <- s_opt$S
    } else {
      S <- if (isTRUE(score_graph_spec$enabled) && isTRUE(score_graph_lambda > 0)) {
        if (identical(engine, "coupled_anchor")) {
          .cgamfa_update_scores_graph(
            Y = anchor_response,
            B = B,
            Z_list = Z_list,
            row_index = row_index,
            alpha_y = alpha_y,
            coupling_lambda = coupling_lambda,
            score_graph_laplacian = score_graph_spec$L,
            score_graph_lambda = score_graph_lambda,
            ridge = ridge
          )
        } else {
          .gamfa_update_scores_graph(
            Y = anchor_response,
            B = B,
            X_list = X_list,
            V_list = V_list,
            row_index = row_index,
            alpha_y = alpha_y,
            alpha_blocks = alpha_blocks,
            score_graph_laplacian = score_graph_spec$L,
            score_graph_lambda = score_graph_lambda,
            ridge = ridge
          )
        }
      } else {
        if (identical(engine, "coupled_anchor")) {
          .cgamfa_update_scores(
            Y = anchor_response,
            B = B,
            Z_list = Z_list,
            row_index = row_index,
            alpha_y = alpha_y,
            coupling_lambda = coupling_lambda,
            ridge = ridge
          )
        } else {
          .lmfa_update_scores(
            Y = anchor_response,
            B = B,
            X_list = X_list,
            V_list = V_list,
            row_index = row_index,
            alpha_y = alpha_y,
            alpha_blocks = alpha_blocks,
            ridge = ridge
          )
        }
      }

      ortho <- .lmfa_orthonormalize_scores(S)
      S <- ortho$S
      rot <- ortho$R
      B <- B %*% t(rot)
      V_list <- lapply(V_list, function(V) V %*% t(rot))
      if (identical(engine, "coupled_anchor")) {
        rot_inv <- solve(rot)
        Z_list <- lapply(Z_list, function(Z) Z %*% rot_inv)
      }
    }

    if (!identical(engine, "coupled_anchor")) {
      Z_list <- lapply(seq_along(X_list), function(k) S[row_index[[k]], , drop = FALSE])
      names(Z_list) <- names(X_list)
    }

    centers <- if (isTRUE(graph$enabled) && isTRUE(graph_lambda > 0)) {
      NULL
    } else {
      .lmfa_group_centers(V_list, fg, block_weights = alpha_blocks)
    }
    group_pen <- .lmfa_feature_group_penalty(V_list, fg = fg, centers = centers, feature_lambda = feature_lambda)
    obj <- if (identical(engine, "coupled_anchor")) {
      .cgamfa_objective(
        Y = anchor_response,
        S = S,
        B = B,
        X_list = X_list,
        Z_list = Z_list,
        V_list = V_list,
        row_index = row_index,
        alpha_y = alpha_y,
        alpha_blocks = alpha_blocks,
        coupling_lambda = coupling_lambda,
        graph_laplacian = graph$L,
        graph_lambda = graph_lambda,
        score_graph_laplacian = score_graph_spec$L,
        score_graph_lambda = score_graph_lambda,
        ridge = ridge
      ) + group_pen
    } else {
      .gamfa_objective(
        Y = anchor_response,
        S = S,
        B = B,
        X_list = X_list,
        V_list = V_list,
        row_index = row_index,
        alpha_y = alpha_y,
        alpha_blocks = alpha_blocks,
        graph_laplacian = graph$L,
        graph_lambda = graph_lambda,
        score_graph_laplacian = score_graph_spec$L,
        score_graph_lambda = score_graph_lambda,
        ridge = ridge
      ) + group_pen
    }
    objective_trace <- c(objective_trace, obj)

    rel_change <- abs(obj - prev_obj) / (abs(prev_obj) + ridge)
    if (isTRUE(verbose)) {
      message(sprintf("%s iter %d: obj=%.6g, rel_change=%.3g", engine, iter, obj, rel_change))
    }
    if (is.finite(prev_obj) && rel_change < tol) break
    prev_obj <- obj
  }

  list(
    S = S,
    B = B,
    V_list = V_list,
    Z_list = Z_list,
    objective_trace = objective_trace,
    score_weights = NULL
  )
}

.ramfa_update_shared_response_loadings <- function(score_list = NULL,
                                                   response_list = NULL,
                                                   response_alpha = NULL,
                                                   response_weights = NULL,
                                                   ridge = 1e-8,
                                                   S_anchor = NULL,
                                                   anchor_response = NULL,
                                                   anchor_response_alpha = 0,
                                                   anchor_response_weights = NULL) {
  score_list <- score_list %||% list()
  response_list <- response_list %||% list()
  response_alpha <- response_alpha %||% numeric(0)
  response_weights <- response_weights %||% list()

  score_ref <- if (!is.null(S_anchor)) S_anchor else if (length(score_list) > 0L) score_list[[1]] else NULL
  response_ref <- if (!is.null(anchor_response)) anchor_response else if (length(response_list) > 0L) response_list[[1]] else NULL
  if (is.null(score_ref) || is.null(response_ref)) {
    stop("At least one score/response source is required to update shared response loadings.", call. = FALSE)
  }

  K <- ncol(score_ref)
  q <- ncol(response_ref)
  Z_pool <- matrix(0, nrow = 0, ncol = K)
  Y_pool <- matrix(0, nrow = 0, ncol = q)

  if (!is.null(S_anchor) && !is.null(anchor_response)) {
    if (is.null(anchor_response_weights)) {
      anchor_response_weights <- rep(1, nrow(anchor_response))
    }
    w_anchor <- sqrt(as.numeric(anchor_response_alpha) * as.numeric(anchor_response_weights))
    keep_anchor <- which(w_anchor > 0)
    if (length(keep_anchor) > 0L) {
      Z_pool <- rbind(Z_pool, S_anchor[keep_anchor, , drop = FALSE] * w_anchor[keep_anchor])
      Y_pool <- rbind(Y_pool, anchor_response[keep_anchor, , drop = FALSE] * w_anchor[keep_anchor])
    }
  }

  for (k in seq_along(score_list)) {
    w <- sqrt(as.numeric(response_alpha[[k]]) * response_weights[[k]])
    keep <- which(w > 0)
    if (length(keep) == 0L) next
    Z_pool <- rbind(Z_pool, score_list[[k]][keep, , drop = FALSE] * w[keep])
    Y_pool <- rbind(Y_pool, response_list[[k]][keep, , drop = FALSE] * w[keep])
  }

  if (nrow(Z_pool) == 0L) {
    return(matrix(0, nrow = q, ncol = K))
  }

  .lmfa_update_loadings(Z_pool, Y_pool, alpha = 1, ridge = ridge)$V
}

.ramfa_update_response_loadings <- function(Z_list,
                                            Y_list,
                                            response_alpha,
                                            response_weights,
                                            S_anchor = NULL,
                                            anchor_response = NULL,
                                            anchor_response_alpha = 0,
                                            anchor_response_weights = NULL,
                                            ridge = 1e-8) {
  .ramfa_update_shared_response_loadings(
    score_list = Z_list,
    response_list = Y_list,
    response_alpha = response_alpha,
    response_weights = response_weights,
    S_anchor = S_anchor,
    anchor_response = anchor_response,
    anchor_response_alpha = anchor_response_alpha,
    anchor_response_weights = anchor_response_weights,
    ridge = ridge
  )
}

.ramfa_update_scores_block <- function(Xk,
                                       Yk,
                                       Vk,
                                       B,
                                       alpha_block,
                                       response_alpha,
                                       response_weight,
                                       anchor_target = NULL,
                                       coupling_lambda = 0,
                                       anchor_weight = NULL,
                                       ridge = 1e-8) {
  K <- ncol(Vk)
  VtV <- crossprod(Vk)
  BtB <- crossprod(B)
  XV <- Xk %*% Vk
  YB <- Yk %*% B
  I_K <- diag(K)

  Z <- matrix(0, nrow = nrow(Xk), ncol = K)
  for (i in seq_len(nrow(Xk))) {
    rw <- response_weight[[i]]
    A <- alpha_block * VtV + response_alpha * rw * BtB + ridge * I_K
    b <- alpha_block * XV[i, ] + response_alpha * rw * YB[i, ]
    if (!is.null(anchor_target) && coupling_lambda > 0) {
      aw <- anchor_weight[[i]]
      A <- A + coupling_lambda * aw * I_K
      b <- b + coupling_lambda * aw * anchor_target[i, ]
    }
    if (all(abs(b) < 1e-14)) {
      Z[i, ] <- 0
    } else {
      Z[i, ] <- solve(A, b)
    }
  }
  Z
}

.ramfa_update_anchor_scores <- function(Z_list,
                                        anchor_map,
                                        anchor_weight,
                                        coupling_lambda,
                                        n_anchor,
                                        B = NULL,
                                        anchor_response = NULL,
                                        anchor_response_alpha = 0,
                                        anchor_response_weights = NULL,
                                        score_graph_laplacian = NULL,
                                        score_graph_lambda = 0,
                                        ridge = 1e-8) {
  K <- ncol(Z_list[[1]])
  lhs <- matrix(0, nrow = n_anchor, ncol = n_anchor)
  rhs <- matrix(0, nrow = n_anchor, ncol = K)

  for (k in seq_along(Z_list)) {
    lam <- as.numeric(coupling_lambda[[k]])
    if (!is.finite(lam) || lam <= 0) next
    Wk <- diag(anchor_weight[[k]], nrow = length(anchor_weight[[k]]))
    Ak <- anchor_map[[k]]
    lhs <- lhs + lam * crossprod(Ak, Wk %*% Ak)
    rhs <- rhs + lam * crossprod(Ak, Wk %*% Z_list[[k]])
  }

  if (isTRUE(score_graph_lambda > 0) && !is.null(score_graph_laplacian)) {
    lhs <- lhs + score_graph_lambda * as.matrix(score_graph_laplacian)
  }

  if (!is.null(anchor_response) && !is.null(B) && as.numeric(anchor_response_alpha) > 0) {
    D <- as.numeric(anchor_response_alpha) * as.numeric(anchor_response_weights)
    if (any(D > 0)) {
      BtB <- crossprod(B)
      YB <- anchor_response %*% B
      rhs_aug <- rhs + diag(D, nrow = n_anchor) %*% YB
      sys <- kronecker(diag(K), lhs + ridge * diag(n_anchor)) +
        kronecker(BtB, diag(D, nrow = n_anchor))
      sol <- solve(sys, as.vector(rhs_aug))
      return(matrix(sol, nrow = n_anchor, ncol = K))
    }
  }

  lhs <- lhs + ridge * diag(n_anchor)
  solve(lhs, rhs)
}

.ramfa_objective <- function(Y_list,
                             Z_list,
                             B,
                             X_list,
                             V_list,
                             alpha_blocks,
                             response_alpha,
                             response_weights,
                             anchor_response = NULL,
                             anchor_response_alpha = 0,
                             anchor_response_weights = NULL,
                             S = NULL,
                             anchor_map = NULL,
                             anchor_weight = NULL,
                             coupling_lambda = NULL,
                             graph_laplacian = NULL,
                             graph_lambda = 0,
                             score_graph_laplacian = NULL,
                             score_graph_lambda = 0,
                             fg,
                             centers,
                             feature_lambda,
                             ridge = 0) {
  err_x <- 0
  err_y <- 0
  err_anchor <- 0
  err_anchor_response <- 0
  for (k in seq_along(X_list)) {
    resid_x <- X_list[[k]] - Z_list[[k]] %*% t(V_list[[k]])
    err_x <- err_x + alpha_blocks[[k]] * sum(resid_x^2)

    resid_y <- Y_list[[k]] - Z_list[[k]] %*% t(B)
    err_y <- err_y + response_alpha[[k]] * sum(response_weights[[k]] * rowSums(resid_y^2))
    if (!is.null(S) && !is.null(anchor_map) && !is.null(coupling_lambda) && coupling_lambda[[k]] > 0) {
      resid_a <- Z_list[[k]] - anchor_map[[k]] %*% S
      err_anchor <- err_anchor + coupling_lambda[[k]] * sum(anchor_weight[[k]] * rowSums(resid_a^2))
    }
  }
  if (!is.null(anchor_response) && !is.null(S) && as.numeric(anchor_response_alpha) > 0) {
    resid_anchor_response <- anchor_response - S %*% t(B)
    err_anchor_response <- as.numeric(anchor_response_alpha) *
      sum(anchor_response_weights * rowSums(resid_anchor_response^2))
  }

  pen <- 0
  if (isTRUE(graph_lambda > 0) && !is.null(graph_laplacian)) {
    V_big <- do.call(rbind, V_list)
    LV <- graph_laplacian %*% V_big
    pen <- pen + graph_lambda * sum(LV * V_big)
  }
  if (isTRUE(score_graph_lambda > 0) && !is.null(score_graph_laplacian) && !is.null(S)) {
    LS <- score_graph_laplacian %*% S
    pen <- pen + score_graph_lambda * sum(LS * S)
  }
  if (isTRUE(feature_lambda > 0) && isTRUE(fg$enabled) && !is.null(centers) && length(fg$members) > 0) {
    groups <- names(fg$members)
    if (is.null(groups)) groups <- seq_along(fg$members)
    for (gi in seq_along(fg$members)) {
      gname <- groups[gi]
      cvec <- centers$centers[match(gname, rownames(centers$centers)), , drop = TRUE]
      for (m in fg$members[[gi]]) {
        v <- V_list[[m$block]][m$feature, , drop = TRUE]
        pen <- pen + m$weight * sum((v - cvec)^2)
      }
    }
    pen <- feature_lambda * pen
  }

  ridge_pen <- ridge * (
    sum(B^2) +
      sum(vapply(Z_list, function(Z) sum(Z^2), numeric(1))) +
      sum(vapply(V_list, function(V) sum(V^2), numeric(1))) +
      if (!is.null(S)) sum(S^2) else 0
  )

  err_x + err_y + err_anchor + err_anchor_response + pen + ridge_pen
}

.ramfa_resolve_prediction_weight <- function(weight, n) {
  if (is.null(weight)) {
    return(rep(1, n))
  }
  weight <- as.numeric(weight)
  if (length(weight) == 1L) {
    weight <- rep(weight, n)
  }
  chk::chk_equal(length(weight), n)
  if (any(!is.finite(weight)) || any(weight < 0)) {
    stop("`response_weight` must be finite and non-negative.", call. = FALSE)
  }
  weight
}

.ramfa_project_scores <- function(object,
                                  new_data,
                                  block,
                                  new_response = NULL,
                                  response_weight = NULL,
                                  new_anchor_map = NULL,
                                  anchor_weight = NULL,
                                  conditional = FALSE,
                                  preprocess = TRUE) {
  resolved <- .muscal_resolve_block_id(block, names(object$V_list))
  block_name <- resolved$name
  Xnew <- as.matrix(new_data)

  if (isTRUE(preprocess)) {
    Xnew <- multivarious::transform(object$block_preproc[[block_name]], Xnew)
  }

  Vk <- object$V_list[[block_name]]
  if (ncol(Xnew) != nrow(Vk)) {
    stop(
      sprintf("Block '%s' has %d columns, expected %d.", block_name, ncol(Xnew), nrow(Vk)),
      call. = FALSE
    )
  }

  alpha_block <- as.numeric(object$alpha_blocks[[block_name]])
  VtV <- crossprod(Vk)
  XV <- Xnew %*% Vk
  K <- ncol(Vk)
  I_K <- diag(K)

  use_response <- !is.null(new_response)
  if (use_response && !isTRUE(conditional)) {
    stop(
      "`new_response` requires `conditional = TRUE`; default prediction uses only `new_data` and optional `new_anchor_map`.",
      call. = FALSE
    )
  }
  if (use_response) {
    Ynew <- as.matrix(new_response)
    chk::chk_equal(nrow(Ynew), nrow(Xnew))
    if (ncol(Ynew) != nrow(object$B)) {
      stop(
        sprintf("`new_response` has %d columns, expected %d.", ncol(Ynew), nrow(object$B)),
        call. = FALSE
      )
    }
    if (isTRUE(preprocess)) {
      Ynew <- multivarious::transform(object$response_preproc, Ynew)
    }
    response_weight <- .ramfa_resolve_prediction_weight(response_weight, nrow(Xnew))
    response_alpha <- as.numeric(object$response_alpha[[block_name]])
    BtB <- crossprod(object$B)
    YB <- Ynew %*% object$B
  } else if (!is.null(response_weight)) {
    stop("`response_weight` cannot be supplied without `new_response`.", call. = FALSE)
  }

  use_anchor <- !is.null(new_anchor_map)
  if (use_anchor) {
    if (is.null(object$S) || is.null(object$anchor_map) || object$coupling_lambda[[block_name]] <= 0) {
      stop("This fit does not have active anchor coupling for the requested block.", call. = FALSE)
    }
    anchor_info <- .ramfa_parse_prediction_anchor_map(
      new_anchor_map = new_anchor_map,
      n = nrow(Xnew),
      n_anchor = nrow(object$S),
      anchor_weight = anchor_weight
    )
    anchor_target <- anchor_info$A_new %*% object$S
  } else if (!is.null(anchor_weight)) {
    stop("`anchor_weight` cannot be supplied without `new_anchor_map`.", call. = FALSE)
  }

  scores_new <- matrix(0, nrow = nrow(Xnew), ncol = K)
  for (i in seq_len(nrow(Xnew))) {
    A <- alpha_block * VtV + object$ridge * I_K
    b <- alpha_block * XV[i, ]
    if (use_response) {
      A <- A + response_alpha * response_weight[[i]] * BtB
      b <- b + response_alpha * response_weight[[i]] * YB[i, ]
    }
    if (use_anchor) {
      cw <- as.numeric(object$coupling_lambda[[block_name]]) * anchor_info$weight[[i]]
      A <- A + cw * I_K
      b <- b + cw * anchor_target[i, ]
    }
    scores_new[i, ] <- solve(A, b)
  }
  scores_new
}

.ramfa_parse_prediction_anchor_map <- function(new_anchor_map,
                                               n,
                                               n_anchor,
                                               anchor_weight = NULL) {
  if (is.matrix(new_anchor_map) || is.data.frame(new_anchor_map)) {
    A_new <- as.matrix(new_anchor_map)
    chk::chk_equal(nrow(A_new), n)
    chk::chk_equal(ncol(A_new), n_anchor)
    if (any(!is.finite(A_new)) || any(A_new < 0)) {
      stop("Matrix `new_anchor_map` must be finite and non-negative.", call. = FALSE)
    }
    row_totals <- rowSums(A_new)
    if (any(row_totals > 0 & abs(row_totals - 1) > 1e-8)) {
      stop("Rows of matrix `new_anchor_map` with positive mass must sum to 1.", call. = FALSE)
    }
    default_weight <- as.numeric(row_totals > 0)
  } else if (is.numeric(new_anchor_map) || is.integer(new_anchor_map)) {
    idx <- as.integer(new_anchor_map)
    chk::chk_equal(length(idx), n)
    bad <- !is.na(idx) & (idx < 1L | idx > n_anchor)
    if (any(bad)) {
      stop("Integer `new_anchor_map` entries must be in 1..n_anchor or NA.", call. = FALSE)
    }
    A_new <- matrix(0, nrow = n, ncol = n_anchor)
    keep <- which(!is.na(idx))
    if (length(keep) > 0L) {
      A_new[cbind(keep, idx[keep])] <- 1
    }
    default_weight <- as.numeric(!is.na(idx))
  } else {
    stop("`new_anchor_map` must be an integer vector or numeric matrix.", call. = FALSE)
  }

  weight <- .ramfa_resolve_prediction_weight(anchor_weight, n)
  if (is.null(anchor_weight)) {
    weight <- default_weight
  }
  row_totals <- rowSums(A_new)
  if (any(weight > 0 & row_totals <= 0)) {
    stop("Rows with positive `anchor_weight` must have a non-zero `new_anchor_map` row.", call. = FALSE)
  }

  list(A_new = A_new, weight = weight)
}

#' Project New Rows into a Response-Aligned MFA Space
#'
#' @param x A fitted `response_aligned_mfa` object.
#' @param new_data Numeric matrix/data.frame of new rows for a known predictor
#'   block.
#' @param block Character or integer identifying the block used for projection.
#' @param new_response Optional observed multivariate response rows. When
#'   supplied together with `conditional = TRUE`, they refine the latent-score
#'   solve using only test-time available information.
#' @param response_weight Optional scalar or vector of non-negative row weights
#'   paired with `new_response`.
#' @param new_anchor_map Optional anchor assignments for the new rows. May be an
#'   integer vector of anchor-state ids (with `NA` for unanchored rows) or a
#'   non-negative row-stochastic matrix.
#' @param anchor_weight Optional scalar or vector of non-negative row weights
#'   paired with `new_anchor_map`.
#' @param conditional Logical; if `TRUE`, allow conditional completion with
#'   `new_response`. The default `FALSE` keeps projection pure with respect to
#'   the response target.
#' @param preprocess Logical; if `TRUE` (default), applies the fitted predictor
#'   and response preprocessing pipelines before projection.
#' @param ... Unused.
#'
#' @return A numeric matrix of latent scores.
#' @export
project.response_aligned_mfa <- function(x,
                                         new_data,
                                         block,
                                         new_response = NULL,
                                         response_weight = NULL,
                                         new_anchor_map = NULL,
                                         anchor_weight = NULL,
                                         conditional = FALSE,
                                         preprocess = TRUE,
                                         ...) {
  chk::chk_true(inherits(x, "response_aligned_mfa"))
  .ramfa_project_scores(
    object = x,
    new_data = new_data,
    block = block,
    new_response = new_response,
    response_weight = response_weight,
    new_anchor_map = new_anchor_map,
    anchor_weight = anchor_weight,
    conditional = conditional,
    preprocess = preprocess
  )
}

#' Predict from a Response-Aligned MFA Fit
#'
#' @param object A fitted `response_aligned_mfa` object.
#' @param new_data Numeric matrix/data.frame of new rows from a known predictor
#'   block.
#' @param block Character or integer identifying the block.
#' @param type One of `"response"`, `"scores"`, or `"reconstruction"`.
#' @param new_response Optional observed multivariate response rows used to
#'   refine the latent score solve when `conditional = TRUE`.
#' @param response_weight Optional scalar or vector of non-negative row weights
#'   paired with `new_response`.
#' @param new_anchor_map Optional anchor assignments for the new rows used to
#'   refine the latent score solve.
#' @param anchor_weight Optional scalar or vector of non-negative row weights
#'   paired with `new_anchor_map`.
#' @param conditional Logical; if `TRUE`, allow conditional completion with
#'   `new_response`. The default `FALSE` keeps prediction target-pure.
#' @param preprocess Logical; if `TRUE` (default), applies the fitted predictor
#'   and response preprocessing pipelines before projection.
#' @param ... Unused.
#'
#' @return A numeric matrix of predicted responses, latent scores, or
#'   reconstructed block rows.
#' @export
predict.response_aligned_mfa <- function(object,
                                         new_data,
                                         block,
                                         type = c("response", "scores", "reconstruction"),
                                         new_response = NULL,
                                         response_weight = NULL,
                                         new_anchor_map = NULL,
                                         anchor_weight = NULL,
                                         conditional = FALSE,
                                         preprocess = TRUE,
                                         ...) {
  type <- match.arg(type)
  scores_new <- .ramfa_project_scores(
    object = object,
    new_data = new_data,
    block = block,
    new_response = new_response,
    response_weight = response_weight,
    new_anchor_map = new_anchor_map,
    anchor_weight = anchor_weight,
    conditional = conditional,
    preprocess = preprocess
  )

  if (type == "scores") {
    return(scores_new)
  }

  resolved <- .muscal_resolve_block_id(block, names(object$V_list))
  block_name <- resolved$name

  if (type == "reconstruction") {
    Xhat_p <- scores_new %*% t(object$V_list[[block_name]])
    return(multivarious::inverse_transform(object$block_preproc[[block_name]], Xhat_p))
  }

  Yhat_p <- scores_new %*% t(object$B)
  multivarious::inverse_transform(object$response_preproc, Yhat_p)
}

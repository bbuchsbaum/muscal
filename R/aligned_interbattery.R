#' Aligned Interbattery Analysis
#'
#' @description
#' `aligned_interbattery()` fits a symmetric two-sided multiblock latent model
#' with subject-specific common scores `T_s`, `X`-side scores `U_s`, and
#' `Y`-side scores `V_s`. The current implementation supports:
#'
#' * one or more subjects;
#' * one or more blocks per subject on each side;
#' * nested `subject -> domain` input or flat block lists with optional
#'   `block_info`;
#' * integer row maps into a subject-level reference row set;
#' * optional subject-wise row-graph penalties on the common score space; and
#' * directional prediction from `X` to `Y`, from `Y` to `X`, or conditional
#'   completion from partially observed bundles on both sides.
#'
#' Loadings are shared across subjects when blocks on the same side share the
#' same domain/schema. This is the minimum coupling required for multi-subject
#' fits to learn a common geometry rather than decomposing into independent
#' subject-specific models. Cross-side decoding uses fitted linear coupling
#' operators from the common score space back to the side-specific score
#' systems.
#'
#' @param X A matrix/data.frame, a flat list of blocks, or a nested
#'   `X[[subject]][[domain]]` list for the `X` side.
#' @param Y A matrix/data.frame, a flat list of blocks, or a nested
#'   `Y[[subject]][[domain]]` list for the `Y` side.
#' @param x_row_map Optional row-map structure for `X`, matching the shape of
#'   `X`. Each observed block row must map to `1..N_s` for its subject. Missing
#'   row maps imply identity mapping and therefore require full subject rows.
#' @param y_row_map Optional row-map structure for `Y`, matching the shape of
#'   `Y`.
#' @param x_block_info Optional `data.frame` for flat `X` input with one row per
#'   block and columns `block`, `subject`, and optional `domain`.
#' @param y_block_info Optional `data.frame` for flat `Y` input with one row per
#'   block and columns `block`, `subject`, and optional `domain`.
#' @param ncomp Integer number of latent components.
#' @param x_preproc Single shared preprocessing pipeline for `X`-side loading
#'   types. Defaults to [multivarious::center()].
#' @param y_preproc Single shared preprocessing pipeline for `Y`-side loading
#'   types. Defaults to [multivarious::center()].
#' @param x_weight Optional non-negative scalar or vector of per-block weights
#'   for `X`. If `NULL`, MFA-style inverse-first-singular-value-squared weighting is
#'   used after preprocessing.
#' @param y_weight Optional non-negative scalar or vector of per-block weights
#'   for `Y`. If `NULL`, MFA-style inverse-first-singular-value-squared weighting is
#'   used after preprocessing.
#' @param row_graph Optional subject-wise row graph or Laplacian. For one
#'   subject, a single square matrix is allowed. For multiple subjects, supply a
#'   named or ordered list. When a graph is supplied, `row_graph_form` controls
#'   whether it is interpreted as a Laplacian or adjacency matrix.
#' @param row_graph_lambda Non-negative scalar controlling the row-graph penalty
#'   on the common scores.
#' @param row_graph_form Interpretation of `row_graph`: `"laplacian"`,
#'   `"adjacency"`, or `"normalized_laplacian"`. The normalized setting accepts
#'   either an adjacency matrix that will be normalized internally or a
#'   precomputed normalized Laplacian.
#' @param couple_x Non-negative scalar controlling the `U_s` to `T_s` coupling.
#' @param couple_y Non-negative scalar controlling the `V_s` to `T_s` coupling.
#' @param decorrelate For `decorrelate_type = "penalty"`, a non-negative scalar
#'   controlling the within-side score decorrelation penalty. For
#'   `decorrelate_type = "whiten"`, the partial-whitening exponent
#'   `gamma in [0, 1]`, where `0` is covariance-oriented and `1` is fully
#'   whitened.
#' @param decorrelate_type One of `"penalty"` or `"whiten"`. Penalty mode adds
#'   a within-side off-diagonal covariance penalty during the side-score
#'   updates. Whitening mode replaces that penalty with a fitted partial
#'   whitening operator inside the common-to-side coupling.
#' @param max_iter Maximum ALS iterations.
#' @param tol Relative convergence tolerance on the objective.
#' @param ridge Non-negative ridge stabilization for score/loadings updates.
#' @param verbose Logical; if `TRUE`, print iteration diagnostics.
#' @param use_future Logical; if `TRUE`, block-wise computations that do not
#'   depend on one another may be performed via `furrr::future_map()` when
#'   available. The main alternating-least-squares loop is intrinsically
#'   sequential and is unaffected. Accepted here primarily for interface parity
#'   with the rest of the aligned/anchored family.
#' @param ... Unused; reserved for future extensions.
#'
#' @return An object inheriting from `multivarious::multiblock_biprojector` with
#'   additional class `"aligned_interbattery"`. The common scores are stored in
#'   `s` and `common_scores`, subject-wise common scores in `T_list`,
#'   subject-wise side scores in `U_list` and `V_list`, and shared side loadings
#'   in `x_loadings` and `y_loadings`.
#' @export
aligned_interbattery <- function(X,
                                 Y,
                                 x_row_map = NULL,
                                 y_row_map = NULL,
                                 x_block_info = NULL,
                                 y_block_info = NULL,
                                 ncomp = 2,
                                 x_preproc = multivarious::center(),
                                 y_preproc = multivarious::center(),
                                 x_weight = NULL,
                                 y_weight = NULL,
                                 row_graph = NULL,
                                 row_graph_lambda = 0,
                                 row_graph_form = c("laplacian", "adjacency", "normalized_laplacian"),
                                 couple_x = 1,
                                 couple_y = 1,
                                 decorrelate = 0,
                                 decorrelate_type = c("penalty", "whiten"),
                                 max_iter = 50,
                                 tol = 1e-6,
                                 ridge = 1e-8,
                                 verbose = FALSE,
                                 use_future = FALSE,
                                 ...) {
  fit_call <- match.call(expand.dots = FALSE)
  row_graph_form <- match.arg(row_graph_form)
  decorrelate_type <- match.arg(decorrelate_type)
  chk::chk_flag(use_future)
  if (isTRUE(use_future) && !requireNamespace("furrr", quietly = TRUE)) {
    stop("use_future = TRUE requires the 'furrr' package.", call. = FALSE)
  }

  ncomp <- as.integer(ncomp)
  chk::chk_integer(ncomp)
  chk::chk_gte(ncomp, 1)
  chk::chk_numeric(row_graph_lambda)
  chk::chk_gte(row_graph_lambda, 0)
  chk::chk_numeric(couple_x)
  chk::chk_gte(couple_x, 0)
  chk::chk_numeric(couple_y)
  chk::chk_gte(couple_y, 0)
  if (couple_x <= 0 && couple_y <= 0) {
    stop("At least one of `couple_x` or `couple_y` must be strictly positive.", call. = FALSE)
  }
  chk::chk_numeric(decorrelate)
  chk::chk_gte(decorrelate, 0)
  if (identical(decorrelate_type, "whiten") && decorrelate > 1) {
    stop("`decorrelate` must lie in [0, 1] when `decorrelate_type = 'whiten'`.", call. = FALSE)
  }
  decor_penalty <- if (identical(decorrelate_type, "penalty")) as.numeric(decorrelate) else 0
  whiten_gamma <- if (identical(decorrelate_type, "whiten")) as.numeric(decorrelate) else 0
  max_iter <- as.integer(max_iter)
  chk::chk_integer(max_iter)
  chk::chk_gte(max_iter, 1)
  chk::chk_numeric(tol)
  chk::chk_gte(tol, 0)
  chk::chk_numeric(ridge)
  chk::chk_gte(ridge, 0)
  chk::chk_flag(verbose)

  x_side <- .aib_flatten_side(
    data = X,
    row_map = x_row_map,
    block_info = x_block_info,
    side = "x"
  )
  y_side <- .aib_flatten_side(
    data = Y,
    row_map = y_row_map,
    block_info = y_block_info,
    side = "y"
  )

  row_spec <- .aib_finalize_row_maps(
    x_blocks = x_side$blocks,
    y_blocks = y_side$blocks,
    x_row_map = x_side$row_map,
    y_row_map = y_side$row_map,
    x_info = x_side$block_info,
    y_info = y_side$block_info
  )

  x_info <- row_spec$x_info
  y_info <- row_spec$y_info
  x_row_map_flat <- row_spec$x_row_map
  y_row_map_flat <- row_spec$y_row_map
  subject_names <- row_spec$subject_names
  subject_sizes <- row_spec$subject_sizes
  subject_index <- .aib_make_subject_index(subject_sizes)
  total_rows <- sum(subject_sizes)
  K <- min(as.integer(ncomp), as.integer(total_rows))
  if (K < 1L) {
    stop("`ncomp` must be <= total number of subject reference rows.", call. = FALSE)
  }

  x_loading_id <- .aib_make_loading_ids(
    blocks = x_side$blocks,
    block_info = x_info,
    side = "x"
  )
  y_loading_id <- .aib_make_loading_ids(
    blocks = y_side$blocks,
    block_info = y_info,
    side = "y"
  )
  x_info$loading <- unname(x_loading_id[x_info$block])
  y_info$loading <- unname(y_loading_id[y_info$block])

  x_prep <- .aib_prepare_side_blocks(
    blocks = x_side$blocks,
    loading_id = x_loading_id,
    preproc = x_preproc,
    side = "x"
  )
  y_prep <- .aib_prepare_side_blocks(
    blocks = y_side$blocks,
    loading_id = y_loading_id,
    preproc = y_preproc,
    side = "y"
  )

  Xp <- x_prep$blocks
  Yp <- y_prep$blocks
  x_weight_parsed <- .aib_parse_block_weights(Xp, x_weight, side = "x", ridge = ridge)
  y_weight_parsed <- .aib_parse_block_weights(Yp, y_weight, side = "y", ridge = ridge)
  x_loading_weight <- .aib_make_loading_weights(x_weight_parsed, x_loading_id)
  y_loading_weight <- .aib_make_loading_weights(y_weight_parsed, y_loading_id)

  graph_spec <- .aib_parse_row_graph(
    row_graph = row_graph,
    subject_sizes = subject_sizes,
    subject_names = subject_names,
    row_graph_lambda = row_graph_lambda,
    row_graph_form = row_graph_form
  )

  init <- .aib_initialize_common_scores(
    X_list = Xp,
    Y_list = Yp,
    x_row_map = x_row_map_flat,
    y_row_map = y_row_map_flat,
    x_info = x_info,
    y_info = y_info,
    x_weight = x_weight_parsed,
    y_weight = y_weight_parsed,
    subject_sizes = subject_sizes,
    subject_names = subject_names,
    K = K
  )
  T_list <- init$T_list
  U_list <- lapply(T_list, function(x) x)
  V_list <- lapply(T_list, function(x) x)
  B_x <- diag(K)
  B_y <- diag(K)
  W_x <- .aib_identity_operator(K)
  W_y <- .aib_identity_operator(K)
  W_x_inv <- .aib_identity_operator(K)
  W_y_inv <- .aib_identity_operator(K)

  x_loadings <- .aib_update_shared_loadings(
    blocks = Xp,
    row_map = x_row_map_flat,
    block_info = x_info,
    loading_id = x_loading_id,
    block_weight = x_weight_parsed,
    subject_scores = U_list,
    ridge = ridge
  )
  y_loadings <- .aib_update_shared_loadings(
    blocks = Yp,
    row_map = y_row_map_flat,
    block_info = y_info,
    loading_id = y_loading_id,
    block_weight = y_weight_parsed,
    subject_scores = V_list,
    ridge = ridge
  )

  objective_trace <- numeric(0)
  prev_obj <- Inf

  for (iter in seq_len(max_iter)) {
    x_loadings <- .aib_update_shared_loadings(
      blocks = Xp,
      row_map = x_row_map_flat,
      block_info = x_info,
      loading_id = x_loading_id,
      block_weight = x_weight_parsed,
      subject_scores = U_list,
      ridge = ridge
    )
    y_loadings <- .aib_update_shared_loadings(
      blocks = Yp,
      row_map = y_row_map_flat,
      block_info = y_info,
      loading_id = y_loading_id,
      block_weight = y_weight_parsed,
      subject_scores = V_list,
      ridge = ridge
    )

    U_list <- .aib_update_side_scores(
      blocks = Xp,
      row_map = x_row_map_flat,
      block_info = x_info,
      loading_id = x_loading_id,
      loadings = x_loadings,
      block_weight = x_weight_parsed,
      target_scores = .aib_map_subject_scores(T_list, B_x),
      couple_transform = W_x,
      couple = as.numeric(couple_x),
      subject_sizes = subject_sizes,
      subject_names = subject_names,
      decorrelate = decor_penalty,
      ridge = ridge
    )
    V_list <- .aib_update_side_scores(
      blocks = Yp,
      row_map = y_row_map_flat,
      block_info = y_info,
      loading_id = y_loading_id,
      loadings = y_loadings,
      block_weight = y_weight_parsed,
      target_scores = .aib_map_subject_scores(T_list, B_y),
      couple_transform = W_y,
      couple = as.numeric(couple_y),
      subject_sizes = subject_sizes,
      subject_names = subject_names,
      decorrelate = decor_penalty,
      ridge = ridge
    )

    T_prev <- .aib_stack_subject_scores(T_list, subject_names)
    T_stack_curr <- .aib_stack_subject_scores(T_list, subject_names)
    U_stack_curr <- .aib_stack_subject_scores(U_list, subject_names)
    V_stack_curr <- .aib_stack_subject_scores(V_list, subject_names)
    if (identical(decorrelate_type, "whiten")) {
      x_white <- .aib_compute_whitening_pair(
        score_list = U_list,
        gamma = whiten_gamma,
        ridge = ridge
      )
      y_white <- .aib_compute_whitening_pair(
        score_list = V_list,
        gamma = whiten_gamma,
        ridge = ridge
      )
      W_x <- x_white$W
      W_y <- y_white$W
      W_x_inv <- x_white$W_inv
      W_y_inv <- y_white$W_inv
    } else {
      W_x <- .aib_identity_operator(K)
      W_y <- .aib_identity_operator(K)
      W_x_inv <- .aib_identity_operator(K)
      W_y_inv <- .aib_identity_operator(K)
    }
    B_x <- .aib_update_coupling_operator(
      T_stack = T_stack_curr,
      side_stack = U_stack_curr %*% W_x,
      ridge = ridge
    )
    B_y <- .aib_update_coupling_operator(
      T_stack = T_stack_curr,
      side_stack = V_stack_curr %*% W_y,
      ridge = ridge
    )

    T_new <- .aib_update_common_scores(
      T_prev = T_prev,
      U_list = U_list,
      V_list = V_list,
      subject_names = subject_names,
      couple_x = as.numeric(couple_x),
      couple_y = as.numeric(couple_y),
      B_x = B_x,
      B_y = B_y,
      W_x = W_x,
      W_y = W_y,
      graph_laplacian = graph_spec$L_big,
      row_graph_lambda = as.numeric(row_graph_lambda),
      ridge = ridge
    )
    T_list <- .aib_split_subject_scores(T_new, subject_index)

    obj <- .aib_objective(
      X_list = Xp,
      Y_list = Yp,
      x_row_map = x_row_map_flat,
      y_row_map = y_row_map_flat,
      x_info = x_info,
      y_info = y_info,
      x_loading_id = x_loading_id,
      y_loading_id = y_loading_id,
      x_loadings = x_loadings,
      y_loadings = y_loadings,
      x_weight = x_weight_parsed,
      y_weight = y_weight_parsed,
      T_list = T_list,
      U_list = U_list,
      V_list = V_list,
      row_graph = graph_spec$L_list,
      row_graph_lambda = as.numeric(row_graph_lambda),
      couple_x = as.numeric(couple_x),
      couple_y = as.numeric(couple_y),
      B_x = B_x,
      B_y = B_y,
      W_x = W_x,
      W_y = W_y,
      decorrelate = decor_penalty,
      ridge = ridge
    )
    objective_trace <- c(objective_trace, obj)

    rel_change <- abs(obj - prev_obj) / (abs(prev_obj) + ridge)
    if (isTRUE(verbose)) {
      message(
        sprintf(
          "aligned_interbattery iter %d: obj=%.6g, rel_change=%.3g",
          iter, obj, rel_change
        )
      )
    }
    if (is.finite(prev_obj) && rel_change < tol) break
    prev_obj <- obj
  }

  T_stack <- .aib_stack_subject_scores(T_list, subject_names)
  U_stack <- .aib_stack_subject_scores(U_list, subject_names)
  V_stack <- .aib_stack_subject_scores(V_list, subject_names)
  if (identical(decorrelate_type, "whiten")) {
    x_white <- .aib_compute_whitening_pair(
      score_list = U_list,
      gamma = whiten_gamma,
      ridge = ridge
    )
    y_white <- .aib_compute_whitening_pair(
      score_list = V_list,
      gamma = whiten_gamma,
      ridge = ridge
    )
    W_x <- x_white$W
    W_y <- y_white$W
    W_x_inv <- x_white$W_inv
    W_y_inv <- y_white$W_inv
  } else {
    W_x <- .aib_identity_operator(K)
    W_y <- .aib_identity_operator(K)
    W_x_inv <- .aib_identity_operator(K)
    W_y_inv <- .aib_identity_operator(K)
  }
  B_x <- .aib_update_coupling_operator(
    T_stack = T_stack,
    side_stack = U_stack %*% W_x,
    ridge = ridge
  )
  B_y <- .aib_update_coupling_operator(
    T_stack = T_stack,
    side_stack = V_stack %*% W_y,
    ridge = ridge
  )
  x_block_scores <- .aib_block_scores_from_subject_scores(x_info, x_row_map_flat, U_list)
  y_block_scores <- .aib_block_scores_from_subject_scores(y_info, y_row_map_flat, V_list)
  partial_scores <- c(x_block_scores, y_block_scores)
  loadings_concat <- .aib_concat_loadings(x_loadings, y_loadings)
  block_fit <- .aib_block_fit(
    X_list = Xp,
    Y_list = Yp,
    x_row_map = x_row_map_flat,
    y_row_map = y_row_map_flat,
    x_info = x_info,
    y_info = y_info,
    x_loading_id = x_loading_id,
    y_loading_id = y_loading_id,
    x_loadings = x_loadings,
    y_loadings = y_loadings,
    U_list = U_list,
    V_list = V_list
  )
  preproc_dummy <- multivarious::fit(
    multivarious::fresh(multivarious::pass()),
    matrix(0, nrow = 1, ncol = nrow(loadings_concat$v))
  )

  fit <- multivarious::multiblock_biprojector(
    v = loadings_concat$v,
    s = T_stack,
    sdev = apply(T_stack, 2, stats::sd),
    preproc = preproc_dummy,
    block_indices = loadings_concat$block_indices,
    common_scores = T_stack,
    T_list = T_list,
    U_list = U_list,
    V_list = V_list,
    x_scores = U_list,
    y_scores = V_list,
    x_loadings = x_loadings,
    y_loadings = y_loadings,
    x_loading_id = x_loading_id,
    y_loading_id = y_loading_id,
    x_block_info = x_info,
    y_block_info = y_info,
    x_block_scores = x_block_scores,
    y_block_scores = y_block_scores,
    x_loading_preproc = x_prep$preproc,
    y_loading_preproc = y_prep$preproc,
    x_weight = x_weight_parsed,
    y_weight = y_weight_parsed,
    x_loading_weight = x_loading_weight,
    y_loading_weight = y_loading_weight,
    subject_index = subject_index,
    subject_sizes = subject_sizes,
    loading_info = loadings_concat$loading_info,
    row_graph_laplacian = graph_spec$L_list,
    row_graph_lambda = as.numeric(row_graph_lambda),
    row_graph_form = row_graph_form,
    couple_x = as.numeric(couple_x),
    couple_y = as.numeric(couple_y),
    x_coupling = B_x %*% W_x_inv,
    y_coupling = B_y %*% W_y_inv,
    x_coupling_operator = B_x,
    y_coupling_operator = B_y,
    x_whitener = W_x,
    y_whitener = W_y,
    x_unwhitener = W_x_inv,
    y_unwhitener = W_y_inv,
    decorrelate = as.numeric(decorrelate),
    decorrelate_type = decorrelate_type,
    objective_trace = objective_trace,
    block_fit = block_fit,
    partial_scores = partial_scores,
    ridge = ridge,
    score_representation = "common_scores",
    classes = "aligned_interbattery"
  )

  fit <- orient_components(fit)
  fit$coupling_maps <- .aib_make_prediction_coupling_maps(
    B_x = fit$x_coupling_operator,
    B_y = fit$y_coupling_operator,
    W_x = fit$x_whitener,
    W_y = fit$y_whitener,
    W_x_inv = fit$x_unwhitener,
    W_y_inv = fit$y_unwhitener,
    ridge = ridge
  )
  fit$prediction_layer <- "native_coupling"

  .muscal_attach_fit_contract(
    fit,
    method = "aligned_interbattery",
    task = "bidirectional_prediction",
    oos_types = c("prediction", "reconstruction", "common_scores", "side_scores"),
    fit_call = fit_call,
    refit_supported = FALSE,
    prediction_target = "bidirectional"
  )
}

.aib_is_matrix_like <- function(x) {
  is.matrix(x) || is.data.frame(x)
}

.aib_align_named_list <- function(x, target_names, what) {
  if (is.null(x)) return(x)
  if (is.null(names(x))) {
    names(x) <- target_names
    return(x)
  }
  if (!setequal(names(x), target_names)) {
    stop(sprintf("%s names must match the corresponding block names exactly.", what), call. = FALSE)
  }
  x[target_names]
}

.aib_validate_flat_block_info <- function(block_info, block_names, side) {
  if (is.null(block_info)) {
    return(data.frame(
      block = block_names,
      subject = rep("S1", length(block_names)),
      domain = block_names,
      stringsAsFactors = FALSE
    ))
  }

  if (!is.data.frame(block_info)) {
    stop(sprintf("`%s_block_info` must be NULL or a data.frame.", side), call. = FALSE)
  }
  if (!("block" %in% names(block_info)) || !("subject" %in% names(block_info))) {
    stop(sprintf("`%s_block_info` must contain `block` and `subject` columns.", side), call. = FALSE)
  }
  if (!setequal(block_info$block, block_names)) {
    stop(sprintf("`%s_block_info$block` must match the flat block names exactly.", side), call. = FALSE)
  }
  info <- block_info[match(block_names, block_info$block), , drop = FALSE]
  if (!("domain" %in% names(info))) {
    info$domain <- info$block
  }
  info$block <- as.character(info$block)
  info$subject <- as.character(info$subject)
  info$domain <- as.character(info$domain)
  info
}

.aib_flatten_side <- function(data,
                              row_map = NULL,
                              block_info = NULL,
                              side = c("x", "y")) {
  side <- match.arg(side)

  if (.aib_is_matrix_like(data)) {
    data <- list(D1 = data)
    if (!is.null(row_map) && !is.list(row_map)) {
      row_map <- list(D1 = row_map)
    }
  }

  chk::chk_list(data)
  chk::chk_true(length(data) >= 1L)

  is_flat <- all(vapply(data, .aib_is_matrix_like, logical(1)))
  if (is_flat) {
    blocks <- lapply(data, as.matrix)
    if (is.null(names(blocks))) {
      names(blocks) <- paste0(toupper(side), seq_along(blocks))
    }
    info <- .aib_validate_flat_block_info(block_info, names(blocks), side = side)

    if (is.null(row_map)) {
      row_map_flat <- setNames(vector("list", length(blocks)), names(blocks))
    } else {
      if (length(blocks) == 1L && !is.list(row_map)) {
        row_map <- list(setNames(list(row_map), names(blocks))[[1]])
        names(row_map) <- names(blocks)
      }
      chk::chk_list(row_map)
      chk::chk_equal(length(row_map), length(blocks))
      row_map <- .aib_align_named_list(row_map, names(blocks), what = sprintf("%s_row_map", side))
      row_map_flat <- lapply(row_map, function(idx) if (is.null(idx)) NULL else as.integer(idx))
      names(row_map_flat) <- names(blocks)
    }

    return(list(
      blocks = blocks,
      row_map = row_map_flat,
      block_info = info
    ))
  }

  subj_names <- names(data)
  if (is.null(subj_names)) subj_names <- paste0("S", seq_along(data))
  names(data) <- subj_names
  row_map <- .aib_align_named_list(row_map, subj_names, what = sprintf("%s_row_map", side))

  blocks <- list()
  row_map_flat <- list()
  info_rows <- vector("list", 0L)

  for (si in seq_along(data)) {
    subj <- subj_names[[si]]
    subj_blocks <- data[[si]]
    subj_map <- if (is.null(row_map)) NULL else row_map[[subj]]

    if (.aib_is_matrix_like(subj_blocks)) {
      subj_blocks <- list(D1 = subj_blocks)
      if (!is.null(subj_map) && !is.list(subj_map)) {
        subj_map <- list(D1 = subj_map)
      }
    }
    if (!is.list(subj_blocks)) {
      stop(sprintf("`%s` must be flat blocks or nested subject/domain lists.", side), call. = FALSE)
    }

    dom_names <- names(subj_blocks)
    if (is.null(dom_names)) dom_names <- paste0("D", seq_along(subj_blocks))
    names(subj_blocks) <- dom_names

    if (is.null(subj_map)) {
      subj_map <- setNames(vector("list", length(subj_blocks)), dom_names)
    } else {
      chk::chk_list(subj_map)
      subj_map <- .aib_align_named_list(
        subj_map,
        dom_names,
        what = sprintf("%s_row_map for subject '%s'", side, subj)
      )
    }

    for (di in seq_along(subj_blocks)) {
      dom <- dom_names[[di]]
      if (!.aib_is_matrix_like(subj_blocks[[di]])) {
        stop(sprintf("Every observed `%s[[subject]][[domain]]` block must be matrix-like.", side), call. = FALSE)
      }
      block_name <- paste(subj, dom, sep = "__")
      blocks[[block_name]] <- as.matrix(subj_blocks[[di]])
      row_map_flat[block_name] <- list(if (is.null(subj_map[[dom]])) NULL else as.integer(subj_map[[dom]]))
      info_rows[[length(info_rows) + 1L]] <- data.frame(
        block = block_name,
        subject = subj,
        domain = dom,
        stringsAsFactors = FALSE
      )
    }
  }

  list(
    blocks = blocks,
    row_map = row_map_flat,
    block_info = do.call(rbind, info_rows)
  )
}

.aib_finalize_row_maps <- function(x_blocks,
                                   y_blocks,
                                   x_row_map,
                                   y_row_map,
                                   x_info,
                                   y_info) {
  x_subjects <- unique(as.character(x_info$subject))
  y_subjects <- unique(as.character(y_info$subject))
  if (!setequal(x_subjects, y_subjects)) {
    stop("In v1, `X` and `Y` must describe the same subject set.", call. = FALSE)
  }

  subject_names <- x_subjects
  x_subject_by_block <- setNames(as.character(x_info$subject), x_info$block)
  y_subject_by_block <- setNames(as.character(y_info$subject), y_info$block)
  subject_sizes <- setNames(integer(length(subject_names)), subject_names)

  for (subj in subject_names) {
    x_blocks_subj <- names(x_subject_by_block)[x_subject_by_block == subj]
    y_blocks_subj <- names(y_subject_by_block)[y_subject_by_block == subj]
    if (length(x_blocks_subj) == 0L || length(y_blocks_subj) == 0L) {
      stop("In v1, each subject must have at least one observed block on both sides.", call. = FALSE)
    }

    provided <- c(
      unlist(lapply(x_row_map[x_blocks_subj], function(idx) if (is.null(idx)) integer(0) else as.integer(idx)), use.names = FALSE),
      unlist(lapply(y_row_map[y_blocks_subj], function(idx) if (is.null(idx)) integer(0) else as.integer(idx)), use.names = FALSE)
    )
    if (length(provided) > 0L) {
      if (any(!is.finite(provided)) || any(is.na(provided)) || any(provided < 1L)) {
        stop("Row maps must contain positive integers with no missing values.", call. = FALSE)
      }
      N_subj <- max(as.integer(provided))
    } else {
      row_counts <- c(
        vapply(x_blocks[x_blocks_subj], nrow, integer(1)),
        vapply(y_blocks[y_blocks_subj], nrow, integer(1))
      )
      if (length(unique(row_counts)) != 1L) {
        stop(
          sprintf(
            "Subject '%s' has omitted row maps but inconsistent row counts across blocks.",
            subj
          ),
          call. = FALSE
        )
      }
      N_subj <- row_counts[[1]]
    }

    subject_sizes[[subj]] <- as.integer(N_subj)

    for (block_name in c(x_blocks_subj, y_blocks_subj)) {
      block_rows <- if (block_name %in% names(x_blocks)) nrow(x_blocks[[block_name]]) else nrow(y_blocks[[block_name]])
      row_map_ref <- if (block_name %in% names(x_row_map)) x_row_map[[block_name]] else y_row_map[[block_name]]

      if (is.null(row_map_ref)) {
        if (block_rows != N_subj) {
          stop(
            sprintf(
              "Missing row map for block '%s' implies identity mapping, but the block has %d rows and subject '%s' has %d reference rows.",
              block_name, block_rows, subj, N_subj
            ),
            call. = FALSE
          )
        }
        row_map_ref <- seq_len(N_subj)
      } else {
        chk::chk_true(is.numeric(row_map_ref) || is.integer(row_map_ref))
        row_map_ref <- as.integer(row_map_ref)
        chk::chk_equal(length(row_map_ref), block_rows)
        if (anyNA(row_map_ref) || any(row_map_ref < 1L) || any(row_map_ref > N_subj)) {
          stop(sprintf("Row map for block '%s' must map into 1..N_subject.", block_name), call. = FALSE)
        }
      }

      if (block_name %in% names(x_row_map)) {
        x_row_map[[block_name]] <- row_map_ref
      } else {
        y_row_map[[block_name]] <- row_map_ref
      }
    }
  }

  list(
    x_info = x_info,
    y_info = y_info,
    x_row_map = x_row_map,
    y_row_map = y_row_map,
    subject_names = subject_names,
    subject_sizes = subject_sizes
  )
}

.aib_block_schema_signature <- function(blocks) {
  dims <- vapply(blocks, ncol, integer(1))
  cn <- lapply(blocks, colnames)
  colnames_same <- {
    first <- cn[[1]]
    all(vapply(cn, function(x) {
      if (is.null(first) && is.null(x)) return(TRUE)
      identical(first, x)
    }, logical(1)))
  }
  list(ncol = dims[[1]], consistent_ncol = all(dims == dims[[1]]), consistent_colnames = colnames_same)
}

.aib_make_loading_ids <- function(blocks, block_info, side) {
  candidate <- as.character(block_info$domain)
  candidate[is.na(candidate) | candidate == ""] <- block_info$block[is.na(candidate) | candidate == ""]
  out <- setNames(character(nrow(block_info)), block_info$block)

  for (dom in unique(candidate)) {
    member_blocks <- block_info$block[candidate == dom]
    sig <- .aib_block_schema_signature(blocks[member_blocks])
    if (isTRUE(sig$consistent_ncol) && isTRUE(sig$consistent_colnames)) {
      out[member_blocks] <- dom
    } else {
      out[member_blocks] <- paste0(side, "__", member_blocks)
    }
  }

  out
}

.aib_prepare_side_blocks <- function(blocks, loading_id, preproc, side) {
  if (is.list(preproc) && !inherits(preproc, "pre_processor") && !inherits(preproc, "prepper")) {
    stop(sprintf("`%s_preproc` must be NULL or a single preprocessing pipeline in v1.", side), call. = FALSE)
  }

  out_blocks <- vector("list", length(blocks))
  names(out_blocks) <- names(blocks)
  preproc_fit <- vector("list", length(unique(loading_id)))
  names(preproc_fit) <- unique(loading_id)

  for (id in names(preproc_fit)) {
    members <- names(loading_id)[loading_id == id]
    pooled <- do.call(rbind, blocks[members])
    pp <- if (is.null(preproc)) {
      NULL
    } else {
      multivarious::fit(multivarious::fresh(preproc), pooled)
    }
    preproc_fit[[id]] <- pp
    for (block_name in members) {
      out_blocks[[block_name]] <- if (is.null(pp)) {
        blocks[[block_name]]
      } else {
        multivarious::transform(pp, blocks[[block_name]])
      }
    }
  }

  list(blocks = out_blocks, preproc = preproc_fit)
}

.aib_parse_block_weights <- function(blocks, weights, side, ridge = 1e-8) {
  block_names <- names(blocks)
  if (is.null(weights)) {
    first_sv <- function(mat) {
      method <- .muscal_svd_method(mat, ncomp = 1)
      tryCatch(
        multivarious::svd_wrapper(mat, ncomp = 1, method = method)$sdev[1],
        error = function(e) multivarious::svd_wrapper(mat, ncomp = 1, method = "base")$sdev[1]
      )
    }
    out <- vapply(blocks, function(mat) 1 / (first_sv(mat)^2 + ridge), numeric(1))
    return(setNames(out, block_names))
  }

  out <- as.numeric(weights)
  if (length(out) == 1L) {
    out <- rep(out, length(blocks))
  } else if (length(out) != length(blocks)) {
    stop(sprintf("`%s_weight` must be NULL, a scalar, or a vector of length equal to the number of blocks.", side), call. = FALSE)
  }
  if (any(!is.finite(out)) || any(out < 0)) {
    stop(sprintf("`%s_weight` must be finite and non-negative.", side), call. = FALSE)
  }
  setNames(out, block_names)
}

.aib_make_loading_weights <- function(block_weight, loading_id) {
  ids <- unique(unname(loading_id))
  out <- setNames(numeric(length(ids)), ids)
  for (id in ids) {
    member_blocks <- names(loading_id)[loading_id == id]
    out[[id]] <- mean(as.numeric(block_weight[member_blocks]))
  }
  out
}

.aib_make_subject_index <- function(subject_sizes) {
  out <- vector("list", length(subject_sizes))
  current <- 1L
  for (i in seq_along(subject_sizes)) {
    n_i <- as.integer(subject_sizes[[i]])
    out[[i]] <- current:(current + n_i - 1L)
    current <- current + n_i
  }
  names(out) <- names(subject_sizes)
  out
}

.aib_stack_subject_scores <- function(score_list, subject_names = names(score_list)) {
  if (is.null(subject_names)) {
    return(do.call(rbind, score_list))
  }
  do.call(rbind, score_list[subject_names])
}

.aib_split_subject_scores <- function(score_stack, subject_index) {
  out <- lapply(subject_index, function(idx) score_stack[idx, , drop = FALSE])
  names(out) <- names(subject_index)
  out
}

.aib_initialize_common_scores <- function(X_list,
                                          Y_list,
                                          x_row_map,
                                          y_row_map,
                                          x_info,
                                          y_info,
                                          x_weight,
                                          y_weight,
                                          subject_sizes,
                                          subject_names,
                                          K) {
  T_list <- vector("list", length(subject_names))
  names(T_list) <- subject_names

  x_subject_by_block <- setNames(as.character(x_info$subject), x_info$block)
  y_subject_by_block <- setNames(as.character(y_info$subject), y_info$block)

  for (subj in subject_names) {
    x_blocks_subj <- names(x_subject_by_block)[x_subject_by_block == subj]
    y_blocks_subj <- names(y_subject_by_block)[y_subject_by_block == subj]

    blocks_subj <- c(X_list[x_blocks_subj], Y_list[y_blocks_subj])
    maps_subj <- c(x_row_map[x_blocks_subj], y_row_map[y_blocks_subj])
    weight_subj <- c(x_weight[x_blocks_subj], y_weight[y_blocks_subj])
    N_subj <- as.integer(subject_sizes[[subj]])
    K_subj <- min(as.integer(K), as.integer(N_subj))

    T0 <- .amfa_initialize_scores(
      X_list = blocks_subj,
      row_index = maps_subj,
      alpha_blocks = weight_subj,
      N = N_subj,
      K = K_subj
    )

    if (K_subj < K) {
      extra <- matrix(stats::rnorm(N_subj * (K - K_subj), sd = 0.05), nrow = N_subj)
      T0 <- cbind(T0, extra)
    }
    T_list[[subj]] <- T0
  }

  T_stack <- .aib_stack_subject_scores(T_list, subject_names)
  T_stack <- .muscal_stiefel_retract(T_stack)
  list(T_list = .aib_split_subject_scores(T_stack, .aib_make_subject_index(subject_sizes)))
}

.aib_update_shared_loadings <- function(blocks,
                                        row_map,
                                        block_info,
                                        loading_id,
                                        block_weight,
                                        subject_scores,
                                        ridge = 1e-8) {
  block_subject <- setNames(as.character(block_info$subject), block_info$block)
  out <- vector("list", length(unique(loading_id)))
  names(out) <- unique(loading_id)

  for (id in names(out)) {
    member_blocks <- names(loading_id)[loading_id == id]
    S_stack <- do.call(rbind, lapply(member_blocks, function(block_name) {
      subj <- block_subject[[block_name]]
      sqrt(block_weight[[block_name]]) * subject_scores[[subj]][row_map[[block_name]], , drop = FALSE]
    }))
    X_stack <- do.call(rbind, lapply(member_blocks, function(block_name) {
      sqrt(block_weight[[block_name]]) * blocks[[block_name]]
    }))
    out[[id]] <- .lmfa_update_loadings(S = S_stack, X = X_stack, alpha = 1, ridge = ridge)$V
  }

  out
}

.aib_update_side_scores <- function(blocks,
                                    row_map,
                                    block_info,
                                    loading_id,
                                    loadings,
                                    block_weight,
                                    target_scores,
                                    couple_transform = .aib_identity_operator(ncol(target_scores[[1L]])),
                                    couple,
                                    subject_sizes,
                                    subject_names,
                                    decorrelate = 0,
                                    ridge = 1e-8) {
  out <- vector("list", length(subject_names))
  names(out) <- subject_names
  block_subject <- setNames(as.character(block_info$subject), block_info$block)
  K <- ncol(target_scores[[subject_names[[1]]]])
  I_K <- diag(K)
  couple_gram <- couple_transform %*% t(couple_transform)

  for (subj in subject_names) {
    N_subj <- as.integer(subject_sizes[[subj]])
    subj_blocks <- names(block_subject)[block_subject == subj]
    target_subj <- target_scores[[subj]]

    if (length(subj_blocks) == 0L) {
      out[[subj]] <- if (couple > 0) target_subj else matrix(0, nrow = N_subj, ncol = K)
      next
    }

    block_contrib <- lapply(subj_blocks, function(block_name) {
      L <- loadings[[loading_id[[block_name]]]]
      agg <- .muscal_rowsum_counts(
        blocks[[block_name]] %*% L,
        row_map[[block_name]],
        N_subj
      )
      list(
        PtP = crossprod(L),
        sums = agg$sums,
        counts = agg$counts,
        weight = as.numeric(block_weight[[block_name]])
      )
    })

    S_subj <- matrix(0, nrow = N_subj, ncol = K)
    for (i in seq_len(N_subj)) {
      A_i <- ridge * I_K
      b_i <- rep(0, K)
      if (couple > 0) {
        A_i <- A_i + couple * couple_gram
        b_i <- b_i + couple * drop(target_subj[i, , drop = FALSE] %*% t(couple_transform))
      }

      for (bc in block_contrib) {
        cnt <- bc$counts[[i]]
        if (cnt == 0L) next
        A_i <- A_i + bc$weight * cnt * bc$PtP
        b_i <- b_i + bc$weight * bc$sums[i, ]
      }

      if (all(abs(b_i) < 1e-14)) {
        S_subj[i, ] <- 0
      } else {
        S_subj[i, ] <- solve(A_i, b_i)
      }
    }
    if (decorrelate > 0) {
      S_subj <- .aib_refine_side_scores_decorrelation(
        S = S_subj,
        blocks = blocks[subj_blocks],
        row_map = row_map[subj_blocks],
        loading_id = loading_id[subj_blocks],
        loadings = loadings,
        block_weight = block_weight[subj_blocks],
        target_scores = target_subj,
        couple_transform = couple_transform,
        couple = couple,
        decorrelate = decorrelate,
        ridge = ridge
      )
    }
    out[[subj]] <- S_subj
  }

  out
}

.aib_parse_row_graph <- function(row_graph,
                                 subject_sizes,
                                 subject_names,
                                 row_graph_lambda = 0,
                                 row_graph_form = c("laplacian", "adjacency", "normalized_laplacian")) {
  row_graph_form <- match.arg(row_graph_form)
  if (is.null(row_graph) || row_graph_lambda <= 0) {
    return(list(
      L_list = setNames(lapply(subject_sizes, function(n) NULL), subject_names),
      L_big = NULL
    ))
  }

  if (!(is.list(row_graph) && !.aib_is_matrix_like(row_graph))) {
    if (length(subject_names) == 1L) {
      row_graph <- list(row_graph)
      names(row_graph) <- subject_names
    } else {
      same_n <- length(unique(subject_sizes)) == 1L
      if (!same_n) {
        stop("A single `row_graph` matrix can only be reused across subjects when every subject has the same reference-row count.", call. = FALSE)
      }
      row_graph <- setNames(rep(list(row_graph), length(subject_names)), subject_names)
    }
  }

  row_graph <- .aib_align_named_list(row_graph, subject_names, what = "row_graph")
  L_list <- vector("list", length(subject_names))
  names(L_list) <- subject_names

  for (subj in subject_names) {
    G <- row_graph[[subj]]
    if (is.null(G)) {
      L_list[[subj]] <- NULL
      next
    }
    n_subj <- as.integer(subject_sizes[[subj]])
    if (nrow(G) != ncol(G) || nrow(G) != n_subj) {
      stop(sprintf("`row_graph` for subject '%s' must be square with dimension %d.", subj, n_subj), call. = FALSE)
    }
    L_list[[subj]] <- .aib_as_laplacian_matrix(
      G = G,
      form = row_graph_form,
      what = sprintf("`row_graph` for subject '%s'", subj)
    )
  }

  L_big <- NULL
  if (length(L_list) > 0L) {
    blocks <- lapply(seq_along(subject_names), function(i) {
      subj <- subject_names[[i]]
      if (is.null(L_list[[subj]])) {
        Matrix::Diagonal(n = as.integer(subject_sizes[[subj]]), x = 0)
      } else {
        L_list[[subj]]
      }
    })
    L_big <- Matrix::bdiag(blocks)
  }

  list(L_list = L_list, L_big = L_big)
}

.aib_as_laplacian_matrix <- function(G, form, what) {
  form <- match.arg(form, c("laplacian", "adjacency", "normalized_laplacian"))

  if (identical(form, "laplacian")) {
    return(.gamfa_validate_graph_matrix(G, "laplacian", what))
  }
  if (identical(form, "adjacency")) {
    A <- .gamfa_validate_graph_matrix(G, "adjacency", what)
    return(.gamfa_laplacian_from_adjacency(A, normalized = FALSE))
  }

  A_try <- tryCatch(
    .gamfa_validate_graph_matrix(G, "adjacency", what),
    error = function(e) NULL
  )
  if (!is.null(A_try)) {
    return(.gamfa_laplacian_from_adjacency(A_try, normalized = TRUE))
  }

  .gamfa_validate_graph_matrix(G, "normalized_laplacian", what)
}

.aib_map_subject_scores <- function(score_list, operator) {
  lapply(score_list, function(S) S %*% operator)
}

.aib_identity_operator <- function(K) {
  diag(K)
}

.aib_matrix_power_sym <- function(A, power, tol = 1e-10) {
  A <- (A + t(A)) / 2
  eig <- eigen(A, symmetric = TRUE)
  vals <- pmax(eig$values, tol)
  eig$vectors %*% diag(vals^power, nrow = length(vals)) %*% t(eig$vectors)
}

.aib_compute_whitening_pair <- function(score_list,
                                        gamma = 0,
                                        ridge = 1e-8) {
  K <- ncol(score_list[[1L]])
  if (gamma <= 0) {
    I_K <- .aib_identity_operator(K)
    return(list(W = I_K, W_inv = I_K))
  }

  S_stack <- do.call(rbind, score_list)
  Sigma <- crossprod(S_stack) / max(1, nrow(S_stack))
  Sigma_r <- Sigma + ridge * diag(K)

  list(
    W = .aib_matrix_power_sym(Sigma_r, -gamma / 2),
    W_inv = .aib_matrix_power_sym(Sigma_r, gamma / 2)
  )
}

.aib_update_coupling_operator <- function(T_stack,
                                          side_stack,
                                          ridge = 1e-8) {
  K <- ncol(T_stack)
  solve(crossprod(T_stack) + ridge * diag(K), crossprod(T_stack, side_stack))
}

.aib_update_common_scores <- function(T_prev,
                                      U_list,
                                      V_list,
                                      subject_names,
                                      couple_x,
                                      couple_y,
                                      B_x,
                                      B_y,
                                      W_x = .aib_identity_operator(ncol(T_prev)),
                                      W_y = .aib_identity_operator(ncol(T_prev)),
                                      graph_laplacian = NULL,
                                      row_graph_lambda = 0,
                                      ridge = 1e-8) {
  U_stack <- .aib_stack_subject_scores(U_list, subject_names)
  V_stack <- .aib_stack_subject_scores(V_list, subject_names)
  rhs <- couple_x * (U_stack %*% W_x) %*% t(B_x) + couple_y * (V_stack %*% W_y) %*% t(B_y)
  A_i <- couple_x * (B_x %*% t(B_x)) + couple_y * (B_y %*% t(B_y))
  A_list <- rep(list(A_i), nrow(rhs))

  .muscal_stiefel_mm(
    S = T_prev,
    A_list = A_list,
    rhs = rhs,
    graph_laplacian = graph_laplacian,
    graph_lambda = row_graph_lambda,
    ridge = ridge
  )$S
}

.aib_score_offdiag_penalty <- function(S) {
  n <- nrow(S)
  if (n <= 0L || ncol(S) <= 1L) return(0)
  C <- crossprod(S) / n
  diag(C) <- 0
  sum(C^2)
}

.aib_score_offdiag_gradient <- function(S) {
  n <- nrow(S)
  if (n <= 0L || ncol(S) <= 1L) {
    return(matrix(0, nrow = nrow(S), ncol = ncol(S)))
  }
  C <- crossprod(S) / n
  diag(C) <- 0
  (4 / n) * S %*% C
}

.aib_side_subject_objective <- function(S,
                                        blocks,
                                        row_map,
                                        loading_id,
                                        loadings,
                                        block_weight,
                                        target_scores,
                                        couple_transform = .aib_identity_operator(ncol(S)),
                                        couple,
                                        decorrelate = 0,
                                        ridge = 0) {
  obj <- 0
  for (block_name in names(blocks)) {
    L <- loadings[[loading_id[[block_name]]]]
    resid <- blocks[[block_name]] - S[row_map[[block_name]], , drop = FALSE] %*% t(L)
    obj <- obj + as.numeric(block_weight[[block_name]]) * sum(resid^2)
  }
  if (couple > 0) {
    obj <- obj + couple * sum((S %*% couple_transform - target_scores)^2)
  }
  if (decorrelate > 0) {
    obj <- obj + decorrelate * .aib_score_offdiag_penalty(S)
  }
  if (ridge > 0) {
    obj <- obj + ridge * sum(S^2)
  }
  obj
}

.aib_side_subject_gradient <- function(S,
                                       blocks,
                                       row_map,
                                       loading_id,
                                       loadings,
                                       block_weight,
                                       target_scores,
                                       couple_transform = .aib_identity_operator(ncol(S)),
                                       couple,
                                       decorrelate = 0,
                                       ridge = 0) {
  grad <- matrix(0, nrow = nrow(S), ncol = ncol(S))
  N <- nrow(S)

  for (block_name in names(blocks)) {
    L <- loadings[[loading_id[[block_name]]]]
    resid <- S[row_map[[block_name]], , drop = FALSE] %*% t(L) - blocks[[block_name]]
    grad_obs <- 2 * as.numeric(block_weight[[block_name]]) * resid %*% L
    grad <- grad + .muscal_rowsum_counts(grad_obs, row_map[[block_name]], N)$sums
  }
  if (couple > 0) {
    grad <- grad + 2 * couple * (S %*% couple_transform - target_scores) %*% t(couple_transform)
  }
  if (decorrelate > 0) {
    grad <- grad + decorrelate * .aib_score_offdiag_gradient(S)
  }
  if (ridge > 0) {
    grad <- grad + 2 * ridge * S
  }

  grad
}

.aib_refine_side_scores_decorrelation <- function(S,
                                                  blocks,
                                                  row_map,
                                                  loading_id,
                                                  loadings,
                                                  block_weight,
                                                  target_scores,
                                                  couple_transform = .aib_identity_operator(ncol(S)),
                                                  couple,
                                                  decorrelate,
                                                  ridge = 0,
                                                  max_iter = 8,
                                                  tol = 1e-8) {
  if (decorrelate <= 0 || ncol(S) <= 1L) return(S)

  objective_fn <- function(Sx) {
    .aib_side_subject_objective(
      S = Sx,
      blocks = blocks,
      row_map = row_map,
      loading_id = loading_id,
      loadings = loadings,
      block_weight = block_weight,
      target_scores = target_scores,
      couple_transform = couple_transform,
      couple = couple,
      decorrelate = decorrelate,
      ridge = ridge
    )
  }
  gradient_fn <- function(Sx) {
    .aib_side_subject_gradient(
      S = Sx,
      blocks = blocks,
      row_map = row_map,
      loading_id = loading_id,
      loadings = loadings,
      block_weight = block_weight,
      target_scores = target_scores,
      couple_transform = couple_transform,
      couple = couple,
      decorrelate = decorrelate,
      ridge = ridge
    )
  }

  obj <- objective_fn(S)
  if (!is.finite(obj)) return(S)

  for (iter in seq_len(max_iter)) {
    G <- gradient_fn(S)
    grad_norm2 <- sum(G^2)
    if (!is.finite(grad_norm2) || grad_norm2 <= tol^2) break

    step <- 1
    accepted <- FALSE
    while (step > 1e-8) {
      S_try <- S - step * G
      obj_try <- objective_fn(S_try)
      if (is.finite(obj_try) && obj_try <= obj - 1e-4 * step * grad_norm2) {
        accepted <- TRUE
        break
      }
      step <- step / 2
    }
    if (!accepted) break
    if (abs(obj - obj_try) / (abs(obj) + 1e-12) < tol) {
      S <- S_try
      break
    }
    S <- S_try
    obj <- obj_try
  }

  S
}

.aib_graph_penalty <- function(T_list, row_graph, row_graph_lambda) {
  if (isTRUE(row_graph_lambda <= 0) || is.null(row_graph)) return(0)
  pen <- 0
  for (subj in names(T_list)) {
    L <- row_graph[[subj]]
    if (is.null(L)) next
    T <- T_list[[subj]]
    LT <- L %*% T
    pen <- pen + row_graph_lambda * sum(LT * T)
  }
  pen
}

.aib_score_decorrelation_penalty <- function(score_list) {
  sum(vapply(score_list, .aib_score_offdiag_penalty, numeric(1)))
}

.aib_objective <- function(X_list,
                           Y_list,
                           x_row_map,
                           y_row_map,
                           x_info,
                           y_info,
                           x_loading_id,
                           y_loading_id,
                           x_loadings,
                           y_loadings,
                           x_weight,
                           y_weight,
                           T_list,
                           U_list,
                           V_list,
                           row_graph = NULL,
                           row_graph_lambda = 0,
                           couple_x = 1,
                           couple_y = 1,
                           B_x,
                           B_y,
                           W_x = .aib_identity_operator(ncol(T_list[[1L]])),
                           W_y = .aib_identity_operator(ncol(T_list[[1L]])),
                           decorrelate = 0,
                           ridge = 0) {
  err_x <- 0
  err_y <- 0

  for (i in seq_len(nrow(x_info))) {
    block_name <- x_info$block[[i]]
    subj <- x_info$subject[[i]]
    resid <- X_list[[block_name]] - U_list[[subj]][x_row_map[[block_name]], , drop = FALSE] %*% t(x_loadings[[x_loading_id[[block_name]]]])
    err_x <- err_x + x_weight[[block_name]] * sum(resid^2)
  }
  for (i in seq_len(nrow(y_info))) {
    block_name <- y_info$block[[i]]
    subj <- y_info$subject[[i]]
    resid <- Y_list[[block_name]] - V_list[[subj]][y_row_map[[block_name]], , drop = FALSE] %*% t(y_loadings[[y_loading_id[[block_name]]]])
    err_y <- err_y + y_weight[[block_name]] * sum(resid^2)
  }

  couple_pen <- 0
  for (subj in names(T_list)) {
    couple_pen <- couple_pen +
      couple_x * sum((U_list[[subj]] %*% W_x - T_list[[subj]] %*% B_x)^2) +
      couple_y * sum((V_list[[subj]] %*% W_y - T_list[[subj]] %*% B_y)^2)
  }

  ridge_pen <- ridge * (
    sum(vapply(T_list, function(x) sum(x^2), numeric(1))) +
      sum(vapply(U_list, function(x) sum(x^2), numeric(1))) +
      sum(vapply(V_list, function(x) sum(x^2), numeric(1))) +
      sum(vapply(x_loadings, function(x) sum(x^2), numeric(1))) +
      sum(vapply(y_loadings, function(x) sum(x^2), numeric(1)))
  )

  decor_pen <- decorrelate * (
    .aib_score_decorrelation_penalty(U_list) +
      .aib_score_decorrelation_penalty(V_list)
  )

  err_x + err_y + couple_pen + decor_pen + .aib_graph_penalty(T_list, row_graph, row_graph_lambda) + ridge_pen
}

.aib_make_prediction_coupling_maps <- function(B_x,
                                               B_y,
                                               W_x = .aib_identity_operator(ncol(B_x)),
                                               W_y = .aib_identity_operator(ncol(B_y)),
                                               W_x_inv = .aib_identity_operator(ncol(B_x)),
                                               W_y_inv = .aib_identity_operator(ncol(B_y)),
                                               ridge = 1e-8) {
  K <- ncol(B_x)
  I_K <- diag(K)
  list(
    x_to_common = W_x %*% t(B_x) %*% solve(B_x %*% t(B_x) + ridge * I_K),
    y_to_common = W_y %*% t(B_y) %*% solve(B_y %*% t(B_y) + ridge * I_K),
    common_to_x = B_x %*% W_x_inv,
    common_to_y = B_y %*% W_y_inv
  )
}

.aib_block_scores_from_subject_scores <- function(block_info, row_map, subject_scores) {
  out <- vector("list", nrow(block_info))
  names(out) <- block_info$block
  for (i in seq_len(nrow(block_info))) {
    block_name <- block_info$block[[i]]
    subj <- block_info$subject[[i]]
    out[[block_name]] <- subject_scores[[subj]][row_map[[block_name]], , drop = FALSE]
  }
  out
}

.aib_concat_loadings <- function(x_loadings, y_loadings) {
  block_indices <- list()
  loading_info <- vector("list", 0L)
  current <- 1L
  mats <- list()

  for (side in c("x", "y")) {
    side_loadings <- if (side == "x") x_loadings else y_loadings
    for (id in names(side_loadings)) {
      mat <- side_loadings[[id]]
      key <- paste(side, id, sep = "::")
      block_indices[[key]] <- current:(current + nrow(mat) - 1L)
      mats[[key]] <- mat
      loading_info[[length(loading_info) + 1L]] <- data.frame(
        side = side,
        loading = id,
        n_features = nrow(mat),
        stringsAsFactors = FALSE
      )
      current <- current + nrow(mat)
    }
  }

  list(
    v = do.call(rbind, mats),
    block_indices = block_indices,
    loading_info = do.call(rbind, loading_info)
  )
}

.aib_block_fit <- function(X_list,
                           Y_list,
                           x_row_map,
                           y_row_map,
                           x_info,
                           y_info,
                           x_loading_id,
                           y_loading_id,
                           x_loadings,
                           y_loadings,
                           U_list,
                           V_list) {
  rows <- vector("list", 0L)

  for (i in seq_len(nrow(x_info))) {
    block_name <- x_info$block[[i]]
    subj <- x_info$subject[[i]]
    Xhat <- U_list[[subj]][x_row_map[[block_name]], , drop = FALSE] %*% t(x_loadings[[x_loading_id[[block_name]]]])
    sse <- sum((X_list[[block_name]] - Xhat)^2)
    tss <- sum(X_list[[block_name]]^2)
    rows[[length(rows) + 1L]] <- data.frame(
      side = "x",
      block = block_name,
      subject = subj,
      domain = x_info$domain[[i]],
      loading = x_loading_id[[block_name]],
      n = nrow(X_list[[block_name]]),
      p = ncol(X_list[[block_name]]),
      sse = sse,
      tss = tss,
      r2 = if (tss > 0) 1 - sse / tss else NA_real_,
      stringsAsFactors = FALSE
    )
  }
  for (i in seq_len(nrow(y_info))) {
    block_name <- y_info$block[[i]]
    subj <- y_info$subject[[i]]
    Yhat <- V_list[[subj]][y_row_map[[block_name]], , drop = FALSE] %*% t(y_loadings[[y_loading_id[[block_name]]]])
    sse <- sum((Y_list[[block_name]] - Yhat)^2)
    tss <- sum(Y_list[[block_name]]^2)
    rows[[length(rows) + 1L]] <- data.frame(
      side = "y",
      block = block_name,
      subject = subj,
      domain = y_info$domain[[i]],
      loading = y_loading_id[[block_name]],
      n = nrow(Y_list[[block_name]]),
      p = ncol(Y_list[[block_name]]),
      sse = sse,
      tss = tss,
      r2 = if (tss > 0) 1 - sse / tss else NA_real_,
      stringsAsFactors = FALSE
    )
  }

  do.call(rbind, rows)
}

.aib_inverse_transform <- function(preproc, x) {
  if (is.null(preproc)) return(x)
  multivarious::inverse_transform(preproc, x)
}

.aib_resolve_loading_query <- function(query,
                                       loading_id,
                                       loadings,
                                       side) {
  loading_names <- names(loadings)
  if (missing(query) || is.null(query)) {
    stop("`block` must identify a source loading type or training block.", call. = FALSE)
  }

  if (is.numeric(query)) {
    if (length(query) != 1L || is.na(query) || query != as.integer(query)) {
      stop("Numeric block identifiers must be single integer indices.", call. = FALSE)
    }
    idx <- as.integer(query)
    if (idx < 1L || idx > length(loading_names)) {
      stop(sprintf("`block` index %d is outside the valid range 1:%d.", idx, length(loading_names)), call. = FALSE)
    }
    return(loading_names[[idx]])
  }

  query <- as.character(query)
  if (length(query) != 1L || is.na(query) || query == "") {
    stop("`block` must be a single non-empty block identifier.", call. = FALSE)
  }
  if (query %in% loading_names) {
    return(query)
  }
  if (query %in% names(loading_id)) {
    return(unname(loading_id[[query]]))
  }
  stop(
    sprintf(
      "Unknown %s-side block/loading '%s'. Available loading ids: %s.",
      side,
      query,
      paste(loading_names, collapse = ", ")
    ),
    call. = FALSE
  )
}

.aib_empty_new_side_bundle <- function() {
  list(
    blocks = list(),
    loading_id = setNames(character(0), character(0)),
    row_map = list(),
    nrow = NA_integer_
  )
}

.aib_align_new_side_row_map <- function(row_map_new,
                                        block_names,
                                        side,
                                        allow_empty = FALSE) {
  if (is.null(row_map_new)) {
    return(setNames(vector("list", length(block_names)), block_names))
  }

  if ((is.numeric(row_map_new) || is.integer(row_map_new)) && length(block_names) == 1L) {
    out <- setNames(list(as.integer(row_map_new)), block_names)
    return(out)
  }

  if (!is.list(row_map_new)) {
    stop(
      sprintf(
        "`%s_row_map_new` must be NULL, an integer vector for a single observed block, or a list matching the observed bundle.",
        side
      ),
      call. = FALSE
    )
  }
  if (length(row_map_new) == 0L && allow_empty) {
    return(setNames(vector("list", length(block_names)), block_names))
  }

  .aib_align_named_list(
    row_map_new,
    block_names,
    what = sprintf("%s_row_map_new", side)
  )
}

.aib_finalize_new_side_row_maps <- function(blocks,
                                            row_map_new,
                                            side) {
  block_names <- names(blocks)
  row_map_list <- .aib_align_new_side_row_map(
    row_map_new = row_map_new,
    block_names = block_names,
    side = side
  )

  provided <- unlist(lapply(row_map_list, function(idx) {
    if (is.null(idx)) integer(0) else as.integer(idx)
  }), use.names = FALSE)

  if (length(provided) > 0L) {
    if (any(!is.finite(provided)) || any(is.na(provided)) || any(provided < 1L)) {
      stop(sprintf("`%s_row_map_new` must contain positive integers with no missing values.", side), call. = FALSE)
    }
    N_ref <- max(as.integer(provided))
  } else {
    row_counts <- vapply(blocks, nrow, integer(1))
    if (length(unique(row_counts)) != 1L) {
      stop(
        sprintf(
          paste(
            "Observed %s-side blocks have inconsistent row counts with no",
            "row map supplied; provide `%s_row_map_new` to map them into a",
            "shared reference row set."
          ),
          side,
          side
        ),
        call. = FALSE
      )
    }
    N_ref <- row_counts[[1]]
  }

  for (block_name in block_names) {
    block_rows <- nrow(blocks[[block_name]])
    idx <- row_map_list[[block_name]]

    if (is.null(idx)) {
      if (block_rows != N_ref) {
        stop(
          sprintf(
            paste(
              "Missing row map for new block '%s' implies identity mapping,",
              "but the block has %d rows and the inferred reference row set",
              "has %d rows."
            ),
            block_name,
            block_rows,
            N_ref
          ),
          call. = FALSE
        )
      }
      idx <- seq_len(N_ref)
    } else {
      chk::chk_true(is.numeric(idx) || is.integer(idx))
      idx <- as.integer(idx)
      chk::chk_equal(length(idx), block_rows)
      if (anyNA(idx) || any(idx < 1L) || any(idx > N_ref)) {
        stop(sprintf("New row map for block '%s' must map into 1..N_reference.", block_name), call. = FALSE)
      }
    }

    row_map_list[[block_name]] <- idx
  }

  list(
    row_map = row_map_list,
    nrow = as.integer(N_ref)
  )
}

.aib_prepare_new_side_bundle <- function(object,
                                         new_data,
                                         block = NULL,
                                         row_map_new = NULL,
                                         side = c("x", "y"),
                                         preprocess = TRUE,
                                         allow_empty = FALSE) {
  side <- match.arg(side)
  if (is.null(new_data)) {
    if (allow_empty) return(.aib_empty_new_side_bundle())
    stop("`new_data` must not be NULL for the requested source side.", call. = FALSE)
  }

  loadings <- if (side == "x") object$x_loadings else object$y_loadings
  loading_id_map <- if (side == "x") object$x_loading_id else object$y_loading_id
  preproc_map <- if (side == "x") object$x_loading_preproc else object$y_loading_preproc

  if (.aib_is_matrix_like(new_data)) {
    if (is.null(block)) {
      stop(
        paste(
          "A single matrix/data.frame `new_data` requires `block` to identify",
          "the source loading."
        ),
        call. = FALSE
      )
    }
    block_list <- list(new_data)
    block_queries <- as.character(block)
    if (length(block_queries) != 1L) {
      stop("Single-block projection requires a single `block` identifier.", call. = FALSE)
    }
    block_names <- "block1"
  } else if (is.list(new_data)) {
    if (length(new_data) == 0L) {
      if (allow_empty) return(.aib_empty_new_side_bundle())
      stop("`new_data` must contain at least one observed block.", call. = FALSE)
    }
    if (!all(vapply(new_data, .aib_is_matrix_like, logical(1)))) {
      stop("Observed block bundles must contain only matrices or data.frames.", call. = FALSE)
    }
    block_list <- new_data
    if (!is.null(names(block_list)) && all(!(is.na(names(block_list)) | names(block_list) == ""))) {
      block_queries <- names(block_list)
      block_names <- make.unique(block_queries)
    } else {
      if (is.null(block)) {
        stop(
          paste(
            "Unnamed block bundles require `block` to identify each observed",
            "loading."
          ),
          call. = FALSE
        )
      }
      block_queries <- as.character(block)
      if (length(block_queries) != length(block_list)) {
        stop("`block` must have one entry per observed block in `new_data`.", call. = FALSE)
      }
      block_names <- make.unique(block_queries)
    }
  } else {
    stop(
      paste(
        "`new_data` must be a matrix/data.frame, a block bundle, or for",
        "`from = 'both'` a list with `x`/`y` bundles."
      ),
      call. = FALSE
    )
  }

  out_blocks <- vector("list", length(block_list))
  out_loading_id <- setNames(character(length(block_list)), block_names)
  n_rows <- integer(length(block_list))

  for (i in seq_along(block_list)) {
    id <- .aib_resolve_loading_query(block_queries[[i]], loading_id_map, loadings, side = side)
    Xi <- as.matrix(block_list[[i]])
    if (isTRUE(preprocess) && !is.null(preproc_map[[id]])) {
      Xi <- multivarious::transform(preproc_map[[id]], Xi)
    }
    if (ncol(Xi) != nrow(loadings[[id]])) {
      stop(
        sprintf(
          "Source block '%s' has %d columns after preprocessing, expected %d.",
          id, ncol(Xi), nrow(loadings[[id]])
        ),
        call. = FALSE
      )
    }
    out_blocks[[i]] <- Xi
    out_loading_id[[i]] <- id
    n_rows[[i]] <- nrow(Xi)
  }

  names(out_blocks) <- block_names
  names(out_loading_id) <- block_names
  row_map_spec <- .aib_finalize_new_side_row_maps(
    blocks = out_blocks,
    row_map_new = row_map_new,
    side = side
  )

  list(
    blocks = out_blocks,
    loading_id = out_loading_id,
    row_map = row_map_spec$row_map,
    nrow = row_map_spec$nrow
  )
}

.aib_project_side_bundle <- function(object, bundle, side = c("x", "y")) {
  side <- match.arg(side)
  if (length(bundle$blocks) == 0L) return(NULL)

  loadings <- if (side == "x") object$x_loadings else object$y_loadings
  loading_weight <- if (side == "x") object$x_loading_weight else object$y_loading_weight
  K <- ncol(loadings[[1]])
  I_K <- diag(K)
  N_ref <- as.integer(bundle$nrow)

  block_contrib <- lapply(names(bundle$blocks), function(block_name) {
    id <- bundle$loading_id[[block_name]]
    L <- loadings[[id]]
    agg <- .muscal_rowsum_counts(
      bundle$blocks[[block_name]] %*% L,
      bundle$row_map[[block_name]],
      N_ref
    )
    list(
      PtP = crossprod(L),
      sums = agg$sums,
      counts = agg$counts,
      weight = as.numeric(loading_weight[[id]])
    )
  })

  S <- matrix(0, nrow = N_ref, ncol = K)
  for (i in seq_len(N_ref)) {
    A_i <- object$ridge * I_K
    b_i <- rep(0, K)
    for (bc in block_contrib) {
      cnt <- bc$counts[[i]]
      if (cnt == 0L) next
      A_i <- A_i + bc$weight * cnt * bc$PtP
      b_i <- b_i + bc$weight * bc$sums[i, ]
    }
    if (all(abs(b_i) < 1e-14)) {
      S[i, ] <- 0
    } else {
      S[i, ] <- solve(A_i, b_i)
    }
  }

  S
}

.aib_parse_new_row_graph <- function(object,
                                     row_graph_new,
                                     n,
                                     row_graph_form = NULL) {
  if (is.null(row_graph_new) || object$row_graph_lambda <= 0) {
    return(NULL)
  }

  if (is.null(row_graph_form)) {
    row_graph_form <- object$row_graph_form
  }
  row_graph_form <- match.arg(row_graph_form, c("laplacian", "adjacency", "normalized_laplacian"))
  if (!all(dim(row_graph_new) == c(n, n))) {
    stop(sprintf("`row_graph_new` must be %d x %d.", n, n), call. = FALSE)
  }
  .aib_as_laplacian_matrix(
    G = row_graph_new,
    form = row_graph_form,
    what = "`row_graph_new`"
  )
}

.aib_solve_common_scores_oos <- function(object,
                                         x_scores = NULL,
                                         y_scores = NULL,
                                         row_graph_new = NULL,
                                         row_graph_form = NULL) {
  if (is.null(x_scores) && is.null(y_scores)) {
    stop("At least one of `x_scores` or `y_scores` must be supplied.", call. = FALSE)
  }

  n <- if (!is.null(x_scores)) nrow(x_scores) else nrow(y_scores)
  K <- if (!is.null(x_scores)) ncol(x_scores) else ncol(y_scores)
  rhs <- matrix(0, nrow = n, ncol = K)
  A <- object$ridge * diag(K)

  if (!is.null(x_scores)) {
    rhs <- rhs + object$couple_x * (x_scores %*% object$x_whitener) %*% t(object$x_coupling_operator)
    A <- A + object$couple_x * (object$x_coupling_operator %*% t(object$x_coupling_operator))
  }
  if (!is.null(y_scores)) {
    rhs <- rhs + object$couple_y * (y_scores %*% object$y_whitener) %*% t(object$y_coupling_operator)
    A <- A + object$couple_y * (object$y_coupling_operator %*% t(object$y_coupling_operator))
  }
  if (sum(abs(A)) <= 0) {
    stop("At least one observed side must have positive coupling weight.", call. = FALSE)
  }

  L_new <- .aib_parse_new_row_graph(
    object = object,
    row_graph_new = row_graph_new,
    n = n,
    row_graph_form = row_graph_form
  )

  if (is.null(L_new)) {
    return(rhs %*% solve(A))
  }

  X0 <- rhs %*% solve(A)
  .muscal_matrix_cg(
    rhs = rhs,
    apply_A = function(S) S %*% A + object$row_graph_lambda * (L_new %*% S),
    X0 = X0,
    tol = 1e-8,
    max_iter = max(50L, 10L * K)
  )
}

.aib_decode_side_scores <- function(object,
                                    side_scores,
                                    side = c("x", "y"),
                                    target_ids) {
  side <- match.arg(side)
  loadings <- if (side == "x") object$x_loadings else object$y_loadings
  preproc_map <- if (side == "x") object$x_loading_preproc else object$y_loading_preproc

  preds <- lapply(target_ids, function(id) {
    pred_p <- side_scores %*% t(loadings[[id]])
    .aib_inverse_transform(preproc_map[[id]], pred_p)
  })
  names(preds) <- target_ids

  if (length(preds) == 1L) {
    return(preds[[1]])
  }
  preds
}

.aib_decode_side_scores_list <- function(object,
                                         side_scores,
                                         side = c("x", "y"),
                                         target_ids) {
  side <- match.arg(side)
  loadings <- if (side == "x") object$x_loadings else object$y_loadings
  preproc_map <- if (side == "x") object$x_loading_preproc else object$y_loading_preproc

  preds <- lapply(target_ids, function(id) {
    pred_p <- side_scores %*% t(loadings[[id]])
    .aib_inverse_transform(preproc_map[[id]], pred_p)
  })
  names(preds) <- target_ids
  preds
}

.aib_reconstruct_bundle <- function(object,
                                    side_scores,
                                    bundle,
                                    side = c("x", "y"),
                                    target_block = NULL) {
  side <- match.arg(side)
  if (length(bundle$blocks) == 0L || is.null(side_scores)) {
    stop(
      sprintf(
        "Reconstruction on the %s side requires observed %s-side blocks in `new_data`.",
        side,
        side
      ),
      call. = FALSE
    )
  }

  loadings <- if (side == "x") object$x_loadings else object$y_loadings
  preproc_map <- if (side == "x") object$x_loading_preproc else object$y_loading_preproc
  loading_id_map <- if (side == "x") object$x_loading_id else object$y_loading_id

  target_blocks <- if (is.null(target_block)) {
    names(bundle$blocks)
  } else if (target_block %in% names(bundle$blocks)) {
    target_block
  } else {
    target_id <- .aib_resolve_loading_query(
      target_block,
      loading_id_map,
      loadings,
      side = side
    )
    names(bundle$loading_id)[bundle$loading_id == target_id]
  }

  if (length(target_blocks) == 0L) {
    stop(
      sprintf(
        "Reconstruction on the %s side is only defined for observed blocks in the new bundle.",
        side
      ),
      call. = FALSE
    )
  }

  preds <- lapply(target_blocks, function(block_name) {
    id <- bundle$loading_id[[block_name]]
    pred_p <- side_scores[bundle$row_map[[block_name]], , drop = FALSE] %*% t(loadings[[id]])
    .aib_inverse_transform(preproc_map[[id]], pred_p)
  })
  names(preds) <- target_blocks

  if (length(preds) == 1L) {
    return(preds[[1]])
  }
  preds
}

.aib_resolve_joint_target <- function(object, target_block) {
  if (is.null(target_block)) return(NULL)
  if (is.numeric(target_block)) {
    stop(
      paste(
        "Numeric `target_block` is not supported when `from = 'both'`; use a",
        "block/loading name, optionally prefixed by `x:` or `y:`."
      ),
      call. = FALSE
    )
  }

  query <- as.character(target_block)
  if (length(query) != 1L || is.na(query) || query == "") {
    stop("`target_block` must be a single non-empty identifier.", call. = FALSE)
  }

  if (startsWith(query, "x:")) {
    return(list(
      side = "x",
      id = .aib_resolve_loading_query(
        substring(query, 3L),
        object$x_loading_id,
        object$x_loadings,
        side = "x"
      )
    ))
  }
  if (startsWith(query, "y:")) {
    return(list(
      side = "y",
      id = .aib_resolve_loading_query(
        substring(query, 3L),
        object$y_loading_id,
        object$y_loadings,
        side = "y"
      )
    ))
  }

  hit_x <- query %in% names(object$x_loadings) || query %in% names(object$x_loading_id)
  hit_y <- query %in% names(object$y_loadings) || query %in% names(object$y_loading_id)
  if (hit_x && hit_y) {
    stop(
      sprintf(
        "`target_block = '%s'` is ambiguous across sides; use `x:%s` or `y:%s`.",
        query, query, query
      ),
      call. = FALSE
    )
  }
  if (hit_x) {
    return(list(
      side = "x",
      id = .aib_resolve_loading_query(query, object$x_loading_id, object$x_loadings, side = "x")
    ))
  }
  if (hit_y) {
    return(list(
      side = "y",
      id = .aib_resolve_loading_query(query, object$y_loading_id, object$y_loadings, side = "y")
    ))
  }

  stop(sprintf("Unknown target block/loading '%s'.", query), call. = FALSE)
}

.aib_project_new_bundle <- function(object,
                                    new_data,
                                    block = NULL,
                                    from = c("x", "y", "both"),
                                    preprocess = TRUE,
                                    x_row_map_new = NULL,
                                    y_row_map_new = NULL,
                                    row_graph_new = NULL,
                                    row_graph_form = NULL) {
  from <- match.arg(from)

  if (from == "both") {
    if (!is.list(new_data) || .aib_is_matrix_like(new_data)) {
      stop(
        paste(
          "For `from = 'both'`, `new_data` must be a list with optional",
          "elements `x` and `y`."
        ),
        call. = FALSE
      )
    }

    block_x <- NULL
    block_y <- NULL
    if (!is.null(block)) {
      if (!is.list(block)) {
        stop(
          paste(
            "When `from = 'both'`, `block` must be NULL or a list with",
            "optional elements `x` and `y`."
          ),
          call. = FALSE
        )
      }
      block_x <- block$x
      block_y <- block$y
    }

    x_bundle <- .aib_prepare_new_side_bundle(
      object = object,
      new_data = new_data$x,
      block = block_x,
      row_map_new = x_row_map_new,
      side = "x",
      preprocess = preprocess,
      allow_empty = TRUE
    )
    y_bundle <- .aib_prepare_new_side_bundle(
      object = object,
      new_data = new_data$y,
      block = block_y,
      row_map_new = y_row_map_new,
      side = "y",
      preprocess = preprocess,
      allow_empty = TRUE
    )

    n_candidates <- c(x_bundle$nrow, y_bundle$nrow)
    n_candidates <- n_candidates[is.finite(n_candidates)]
    if (length(n_candidates) == 0L) {
      stop("`new_data` must contain at least one observed block bundle.", call. = FALSE)
    }
    if (length(unique(n_candidates)) != 1L) {
      stop("Observed `x` and `y` bundles must have the same number of rows.", call. = FALSE)
    }

    x_scores <- .aib_project_side_bundle(object, x_bundle, side = "x")
    y_scores <- .aib_project_side_bundle(object, y_bundle, side = "y")
    common_scores <- .aib_solve_common_scores_oos(
      object = object,
      x_scores = x_scores,
      y_scores = y_scores,
      row_graph_new = row_graph_new,
      row_graph_form = row_graph_form
    )

    return(list(
      common_scores = common_scores,
      x_scores = x_scores,
      y_scores = y_scores,
      x_bundle = x_bundle,
      y_bundle = y_bundle
    ))
  }

  bundle <- .aib_prepare_new_side_bundle(
    object = object,
    new_data = new_data,
    block = block,
    row_map_new = if (from == "x") x_row_map_new else y_row_map_new,
    side = from,
    preprocess = preprocess,
    allow_empty = FALSE
  )
  side_scores <- .aib_project_side_bundle(object, bundle, side = from)
  common_scores <- .aib_solve_common_scores_oos(
    object = object,
    x_scores = if (from == "x") side_scores else NULL,
    y_scores = if (from == "y") side_scores else NULL,
    row_graph_new = row_graph_new,
    row_graph_form = row_graph_form
  )

  list(
    side_scores = side_scores,
    common_scores = common_scores,
    loading_id = unique(unname(bundle$loading_id)),
    bundle = bundle
  )
}

#' Project New Rows into an Aligned Interbattery Space
#'
#' @param x A fitted `aligned_interbattery` object.
#' @param new_data Observed out-of-sample data. For `from = "x"` or `"y"`, this
#'   may be a single matrix/data.frame or a list of observed blocks from that
#'   side. For `from = "both"`, supply a list with optional elements `x` and
#'   `y`, each containing a single block or block bundle for that side.
#' @param block Optional source block identifier(s). For single-matrix input,
#'   this identifies the source loading type or training block. For unnamed
#'   block bundles, provide one identifier per observed block. For
#'   `from = "both"`, `block` may be `NULL` or a list with optional `x`/`y`
#'   entries matching the supplied side bundles.
#' @param from One of `"x"`, `"y"`, or `"both"`.
#' @param side One of `"common"`, `"side"`, `"x"`, `"y"`, or `"both"`.
#' @param preprocess Logical; if `TRUE` (default), apply the fitted shared
#'   preprocessing pipeline(s) before projection.
#' @param x_row_map_new Optional integer row-map structure for new `x`-side
#'   blocks. Missing maps imply identity mapping and therefore require aligned
#'   row counts within the supplied `x` bundle.
#' @param y_row_map_new Optional integer row-map structure for new `y`-side
#'   blocks.
#' @param row_graph_new Optional graph/Laplacian on the new rows used to refine
#'   the common-score solve. This is never required.
#' @param row_graph_form Interpretation of `row_graph_new`. When `NULL`
#'   (default), the fitted object's `row_graph_form` is used.
#' @param ... Unused.
#'
#' @return A numeric score matrix or, for `side = "both"` / `from = "both"`, a
#'   named list of score matrices.
#' @export
project.aligned_interbattery <- function(x,
                                         new_data,
                                         block = NULL,
                                         from = c("x", "y", "both"),
                                         side = c("common", "side", "x", "y", "both"),
                                         preprocess = TRUE,
                                         x_row_map_new = NULL,
                                         y_row_map_new = NULL,
                                         row_graph_new = NULL,
                                         row_graph_form = NULL,
                                         ...) {
  chk::chk_true(inherits(x, "aligned_interbattery"))
  from <- match.arg(from)
  side <- match.arg(side)
  proj <- .aib_project_new_bundle(
    object = x,
    new_data = new_data,
    block = block,
    from = from,
    preprocess = preprocess,
    x_row_map_new = x_row_map_new,
    y_row_map_new = y_row_map_new,
    row_graph_new = row_graph_new,
    row_graph_form = row_graph_form
  )

  if (from == "both") {
    if (side == "side") {
      stop("`side = 'side'` is only defined for one-sided projection.", call. = FALSE)
    }
    if (side == "common") return(proj$common_scores)
    if (side == "x") return(proj$x_scores)
    if (side == "y") return(proj$y_scores)
    return(list(common = proj$common_scores, x = proj$x_scores, y = proj$y_scores))
  }

  if (side == "common") {
    return(proj$common_scores)
  }
  if (side == "side") {
    return(proj$side_scores)
  }
  if (side == from) {
    return(proj$side_scores)
  }
  if (side == "both") {
    return(setNames(list(proj$common_scores, proj$side_scores), c("common", from)))
  }
  stop(sprintf("`side = '%s'` is incompatible with `from = '%s'`.", side, from), call. = FALSE)
}

#' Predict from an Aligned Interbattery Fit
#'
#' @param object A fitted `aligned_interbattery` object.
#' @inheritParams project.aligned_interbattery
#' @param type One of `"prediction"`, `"reconstruction"`, `"common_scores"`, or
#'   `"side_scores"`.
#' @param target_block Optional target loading type or training block name. For
#'   one-sided prediction this targets the opposite side. For `from = "both"`,
#'   it may identify either side; ambiguous names can be disambiguated via
#'   `x:<name>` or `y:<name>`. When omitted, multi-target calls return named
#'   lists.
#' @return A numeric matrix or a named list of matrices for multi-target
#'   prediction / completion.
#' @export
predict.aligned_interbattery <- function(object,
                                         new_data,
                                         block = NULL,
                                         from = c("x", "y", "both"),
                                         type = c("prediction", "reconstruction", "common_scores", "side_scores"),
                                         target_block = NULL,
                                         preprocess = TRUE,
                                         x_row_map_new = NULL,
                                         y_row_map_new = NULL,
                                         row_graph_new = NULL,
                                         row_graph_form = NULL,
                                         ...) {
  chk::chk_true(inherits(object, "aligned_interbattery"))
  from <- match.arg(from)
  type <- match.arg(type)
  proj <- .aib_project_new_bundle(
    object = object,
    new_data = new_data,
    block = block,
    from = from,
    preprocess = preprocess,
    x_row_map_new = x_row_map_new,
    y_row_map_new = y_row_map_new,
    row_graph_new = row_graph_new,
    row_graph_form = row_graph_form
  )

  if (type == "common_scores") {
    return(proj$common_scores)
  }
  if (type == "side_scores") {
    if (from == "both") {
      return(list(x = proj$x_scores, y = proj$y_scores))
    }
    return(proj$side_scores)
  }

  if (type == "reconstruction") {
    if (from == "both") {
      if (is.null(target_block)) {
        out <- list()
        if (!is.null(proj$x_scores)) {
          out$x <- .aib_reconstruct_bundle(object, proj$x_scores, proj$x_bundle, side = "x")
        }
        if (!is.null(proj$y_scores)) {
          out$y <- .aib_reconstruct_bundle(object, proj$y_scores, proj$y_bundle, side = "y")
        }
        return(out)
      }
      target <- .aib_resolve_joint_target(object, target_block)
      if (target$side == "x") {
        return(.aib_reconstruct_bundle(object, proj$x_scores, proj$x_bundle, side = "x", target_block = target$id))
      }
      return(.aib_reconstruct_bundle(object, proj$y_scores, proj$y_bundle, side = "y", target_block = target$id))
    }

    return(.aib_reconstruct_bundle(
      object = object,
      side_scores = proj$side_scores,
      bundle = proj$bundle,
      side = from,
      target_block = target_block
    ))
  }

  if (from == "both") {
    x_target_scores <- proj$common_scores %*% object$x_coupling
    y_target_scores <- proj$common_scores %*% object$y_coupling

    if (is.null(target_block)) {
      return(list(
        x = .aib_decode_side_scores_list(object, x_target_scores, side = "x", target_ids = names(object$x_loadings)),
        y = .aib_decode_side_scores_list(object, y_target_scores, side = "y", target_ids = names(object$y_loadings))
      ))
    }

    target <- .aib_resolve_joint_target(object, target_block)
    side_scores <- if (target$side == "x") x_target_scores else y_target_scores
    return(.aib_decode_side_scores(object, side_scores, side = target$side, target_ids = target$id))
  }

  to_side <- if (from == "x") "y" else "x"
  target_scores <- if (to_side == "x") {
    proj$common_scores %*% object$x_coupling
  } else {
    proj$common_scores %*% object$y_coupling
  }
  target_ids <- if (!is.null(target_block)) {
    .aib_resolve_loading_query(
      target_block,
      if (to_side == "x") object$x_loading_id else object$y_loading_id,
      if (to_side == "x") object$x_loadings else object$y_loadings,
      side = to_side
    )
  } else {
    if (to_side == "x") names(object$x_loadings) else names(object$y_loadings)
  }

  .aib_decode_side_scores(object, target_scores, side = to_side, target_ids = target_ids)
}

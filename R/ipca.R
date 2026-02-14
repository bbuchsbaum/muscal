#' @md
#' @rdname ipca
#' @details
#' `ipca.list()` converts a list of blocks to a `multiblock` object and then
#' dispatches to [ipca.multiblock()].
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' X <- list(
#'   X1 = matrix(rnorm(60 * 40), 60, 40),
#'   X2 = matrix(rnorm(60 * 80), 60, 80),
#'   X3 = matrix(rnorm(60 * 50), 60, 50)
#' )
#' fit <- ipca(X, ncomp = 3, lambda = 1)
#' stopifnot(ncol(multivarious::scores(fit)) == 3)
#' }
#' @export
ipca.list <- function(data,
                      preproc = multivarious::center(),
                      ncomp = 2,
                      lambda = 1,
                      method = c("auto", "gram", "dense"),
                      max_iter = 100,
                      tol = 1e-6,
                      normalize_trace = TRUE,
                      use_future = FALSE,
                      eig_solver = c("auto", "full", "truncated"),
                      eig_rank = NULL,
                      eig_trunc_min_n = 400,
                      ...) {
  data <- multidesign::multiblock(data)
  ipca.multiblock(
    data = data,
    preproc = preproc,
    ncomp = ncomp,
    lambda = lambda,
    method = method,
    max_iter = max_iter,
    tol = tol,
    normalize_trace = normalize_trace,
    use_future = use_future,
    eig_solver = eig_solver,
    eig_rank = eig_rank,
    eig_trunc_min_n = eig_trunc_min_n,
    ...
  )
}


#' @md
#' @rdname ipca
#' @details
#' `ipca.multiblock()` implements multiplicative Frobenius iPCA from Tang & Allen
#' (2021) using a Flip-Flop algorithm. In `method = "gram"` mode, fitting is
#' performed in sample space (`n x n`) and avoids `p_k x p_k` eigendecompositions.
#' This is exact for the multiplicative Frobenius updates.
#'
#' @export
ipca.multiblock <- function(data,
                            preproc = multivarious::center(),
                            ncomp = 2,
                            lambda = 1,
                            method = c("auto", "gram", "dense"),
                            max_iter = 100,
                            tol = 1e-6,
                            normalize_trace = TRUE,
                            use_future = FALSE,
                            eig_solver = c("auto", "full", "truncated"),
                            eig_rank = NULL,
                            eig_trunc_min_n = 400,
                            .init_state = NULL,
                            .return_state = FALSE,
                            ...) {
  chk::chk_true(length(data) > 1)
  for (i in seq_along(data)) {
    chk::chkor_vld(chk::vld_matrix(data[[i]]), chk::vld_s4_class(data[[i]], "Matrix"))
  }

  nrs <- vapply(data, nrow, integer(1))
  chk::chk_true(all(nrs == nrs[1]))
  n <- nrs[1]
  K <- length(data)

  method <- match.arg(method)
  eig_solver <- match.arg(eig_solver)
  chk::chk_flag(normalize_trace)
  chk::chk_flag(use_future)
  chk::chk_flag(.return_state)
  chk::chk_number(max_iter)
  chk::chk_gte(max_iter, 1)
  max_iter <- as.integer(max_iter)
  chk::chk_number(tol)
  chk::chk_gt(tol, 0)

  chk::chk_number(ncomp)
  chk::chk_gte(ncomp, 1)
  if (abs(ncomp - round(ncomp)) > 1e-8) {
    stop("ncomp must be an integer-like value.", call. = FALSE)
  }
  ncomp_eff <- min(as.integer(round(ncomp)), n)
  if (!is.null(eig_rank)) {
    chk::chk_number(eig_rank)
    chk::chk_gte(eig_rank, ncomp_eff)
    eig_rank <- as.integer(round(eig_rank))
  }
  chk::chk_number(eig_trunc_min_n)
  chk::chk_gte(eig_trunc_min_n, 2)
  eig_trunc_min_n <- as.integer(round(eig_trunc_min_n))

  if (is.null(names(data))) {
    names(data) <- paste0("B", seq_len(K))
  }

  prep <- prepare_block_preprocessors(data, preproc, check_consistent_ncol = FALSE)
  Xp <- lapply(prep$Xp, function(x) {
    if (inherits(x, "Matrix")) as.matrix(x) else x
  })
  proclist <- prep$proclist

  p_vec <- vapply(Xp, ncol, integer(1))
  p_total <- sum(p_vec)

  if (length(lambda) == 1) {
    lambda <- rep(as.numeric(lambda), K)
  } else {
    lambda <- as.numeric(lambda)
  }
  chk::chk_true(length(lambda) == K)
  chk::chk_true(all(is.finite(lambda)))
  chk::chk_true(all(lambda > 0))

  method_used <- method
  if (method == "auto") {
    method_used <- if (any(p_vec > n)) "gram" else "dense"
  }

  block_indices <- .ipca_block_indices(Xp)
  proc <- multivarious::concat_pre_processors(proclist, block_indices)

  fit <- if (identical(method_used, "gram")) {
    .ipca_fit_gram(
      X = Xp,
      p_vec = p_vec,
      lambda = lambda,
      ncomp = ncomp_eff,
      max_iter = max_iter,
      tol = tol,
      normalize_trace = normalize_trace,
      use_future = use_future,
      eig_solver = eig_solver,
      eig_rank = eig_rank,
      eig_trunc_min_n = eig_trunc_min_n,
      init_state = .init_state,
      return_state = .return_state
    )
  } else {
    .ipca_fit_dense(
      X = Xp,
      p_vec = p_vec,
      lambda = lambda,
      ncomp = ncomp_eff,
      max_iter = max_iter,
      tol = tol,
      normalize_trace = normalize_trace,
      eig_solver = eig_solver,
      eig_rank = eig_rank,
      eig_trunc_min_n = eig_trunc_min_n,
      init_state = .init_state,
      return_state = .return_state
    )
  }

  U <- fit$U
  phi <- fit$phi

  sdev <- sqrt(pmax(phi[seq_len(ncomp_eff)], 0))
  S_scores <- U[, seq_len(ncomp_eff), drop = FALSE] %*%
    diag(sdev, nrow = ncomp_eff, ncol = ncomp_eff)

  v_concat <- matrix(0, nrow = p_total, ncol = ncomp_eff)
  for (k in seq_along(Xp)) {
    v_concat[block_indices[[k]], ] <- fit$loadings[[k]]
  }

  partial_scores <- lapply(seq_along(Xp), function(k) Xp[[k]] %*% fit$loadings[[k]])
  names(partial_scores) <- names(Xp)
  projection_map <- .ipca_fit_projection_map(Xp, S_scores)

  out <- multivarious::multiblock_biprojector(
    v = v_concat,
    s = S_scores,
    sdev = sdev,
    preproc = proc,
    block_indices = block_indices,
    Sigma_eigenvalues = phi,
    Sigma_eigenvectors = U,
    Delta_leading_eigenvalues = fit$delta_leading,
    lambda = lambda,
    method_used = method_used,
    converged = fit$converged,
    iter = fit$iter,
    rel_change_phi = fit$rel_change_phi,
    rel_change_delta = fit$rel_change_delta,
    eig_solver = eig_solver,
    eig_solver_used = fit$eig_solver_used,
    normalize_trace = normalize_trace,
    partial_scores = partial_scores,
    projection_map = projection_map,
    names = names(Xp),
    classes = "ipca"
  )

  if (isTRUE(.return_state) && !is.null(fit$warm_state)) {
    out$warm_state <- fit$warm_state
  }
  out
}


.ipca_block_indices <- function(blocks) {
  out <- vector("list", length(blocks))
  idx <- 1L
  for (k in seq_along(blocks)) {
    p_k <- ncol(blocks[[k]])
    out[[k]] <- idx:(idx + p_k - 1L)
    idx <- idx + p_k
  }
  names(out) <- names(blocks)
  out
}


.ipca_sym <- function(M) {
  (M + t(M)) / 2
}


.ipca_rel_change <- function(new, old, eps = 1e-12) {
  if (length(new) == 0L || length(old) == 0L) return(Inf)
  max(abs(new - old) / (abs(old) + eps))
}


.ipca_hypot_const <- function(d, cst) {
  if (!is.finite(cst) || cst < 0) {
    stop("Invalid shrinkage constant encountered in iPCA updates.", call. = FALSE)
  }
  if (cst == 0) return(abs(d))
  b <- sqrt(cst)
  a <- abs(d)
  m <- pmax(a, b)
  out <- m * sqrt((a / m)^2 + (b / m)^2)
  out[m == 0] <- 0
  out
}


.ipca_shrink_eigs <- function(d, cst, denom) {
  s <- .ipca_hypot_const(d, cst)
  eig <- (d + s) / denom
  inv_eig <- denom / (d + s)
  list(eig = eig, inv_eig = inv_eig)
}


.ipca_choose_eig_solver <- function(eig_solver, n, eig_trunc_min_n) {
  if (identical(eig_solver, "full")) return("full")
  if (identical(eig_solver, "truncated")) return("truncated")
  if (n >= eig_trunc_min_n && requireNamespace("RSpectra", quietly = TRUE)) {
    return("truncated")
  }
  "full"
}


.ipca_effective_rank <- function(eig_rank, n, ncomp) {
  if (n <= 2) return(NA_integer_)
  if (is.null(eig_rank)) {
    r <- max(ncomp + 5L, as.integer(ceiling(0.20 * n)))
  } else {
    r <- as.integer(round(eig_rank))
  }
  r <- max(r, ncomp)
  r <- min(r, n - 1L)
  if (r < ncomp) stop("eig_rank must be >= ncomp for truncated eigensolver.", call. = FALSE)
  r
}


.ipca_pad_cols <- function(M, ncol_target) {
  if (ncol(M) == ncol_target) return(M)
  if (ncol(M) > ncol_target) return(M[, seq_len(ncol_target), drop = FALSE])
  cbind(M, matrix(0, nrow = nrow(M), ncol = ncol_target - ncol(M)))
}


.ipca_pad_vec <- function(x, n_target) {
  if (length(x) == n_target) return(x)
  if (length(x) > n_target) return(x[seq_len(n_target)])
  c(x, rep(0, n_target - length(x)))
}


.ipca_fit_projection_map <- function(X_blocks, S_scores, ridge = 1e-8) {
  if (length(X_blocks) == 0L) {
    return(list(type = "identity", coef_blocks = list()))
  }

  XXt <- Reduce(`+`, lapply(X_blocks, tcrossprod))
  tr <- mean(diag(XXt))
  ridge_eff <- ridge * if (is.finite(tr) && tr > 0) tr else 1
  XXt_r <- XXt + ridge_eff * diag(nrow(XXt))

  A <- tryCatch(
    solve(XXt_r, S_scores),
    error = function(e) MASS::ginv(XXt_r) %*% S_scores
  )

  coef_blocks <- lapply(X_blocks, function(Xk) crossprod(Xk, A))
  list(
    type = "dual_ridge",
    n_blocks = length(X_blocks),
    ridge = ridge_eff,
    coef_blocks = coef_blocks
  )
}


.ipca_bind_blocks <- function(new_data, block_indices, block_names = NULL) {
  if (is.matrix(new_data) || is.data.frame(new_data)) {
    return(as.matrix(new_data))
  }

  if (!is.list(new_data)) {
    stop("new_data must be a matrix/data.frame or list of blocks.", call. = FALSE)
  }

  K <- length(block_indices)
  if (length(new_data) != K) {
    stop(sprintf("Expected %d blocks in new_data, got %d.", K, length(new_data)), call. = FALSE)
  }

  if (!is.null(names(new_data)) && !is.null(block_names) && all(block_names %in% names(new_data))) {
    new_data <- new_data[block_names]
  }

  blocks <- lapply(seq_along(new_data), function(k) {
    Xk <- new_data[[k]]
    chk::chk_true(is.matrix(Xk) || is.data.frame(Xk))
    Xk <- as.matrix(Xk)
    expected_p <- length(block_indices[[k]])
    if (ncol(Xk) != expected_p) {
      stop(
        sprintf("Block %d has %d columns, expected %d.", k, ncol(Xk), expected_p),
        call. = FALSE
      )
    }
    Xk
  })

  nrs <- vapply(blocks, nrow, integer(1))
  if (!all(nrs == nrs[1])) {
    stop("All new_data blocks must have the same number of rows.", call. = FALSE)
  }

  do.call(cbind, blocks)
}


#' @rdname project
#' @export
project.ipca <- function(x, new_data, ...) {
  X_raw <- .ipca_bind_blocks(new_data, x$block_indices, x$names)
  Xp <- reprocess(x, X_raw)
  Xp <- if (inherits(Xp, "Matrix")) as.matrix(Xp) else Xp

  n_new <- nrow(Xp)
  ncomp <- ncol(x$s)
  map <- x$projection_map
  if (!is.null(map$coef_blocks) && length(map$coef_blocks) == length(x$block_indices)) {
    out <- matrix(0, nrow = n_new, ncol = ncomp)
    for (k in seq_along(x$block_indices)) {
      idx <- x$block_indices[[k]]
      B_k <- map$coef_blocks[[k]]
      out <- out + Xp[, idx, drop = FALSE] %*% B_k
    }
    return(out)
  }

  Xp %*% x$v
}


.ipca_fit_dense <- function(X,
                            p_vec,
                            lambda,
                            ncomp,
                            max_iter,
                            tol,
                            normalize_trace,
                            eig_solver,
                            eig_rank,
                            eig_trunc_min_n,
                            init_state = NULL,
                            return_state = FALSE) {
  K <- length(X)
  n <- nrow(X[[1]])
  p_total <- sum(p_vec)

  Sigma_inv <- diag(n)
  Delta_inv <- lapply(p_vec, diag)
  delta_inv_norm2 <- as.numeric(p_vec)

  phi <- rep(1, n)
  U <- diag(n)
  loadings <- lapply(p_vec, function(p) matrix(0, p, ncomp))
  delta_leading <- lapply(seq_along(p_vec), function(...) rep(0, ncomp))

  converged <- FALSE
  rel_phi <- Inf
  rel_delta <- Inf
  iter <- 0L
  eig_solver_mode <- .ipca_choose_eig_solver(eig_solver, n, eig_trunc_min_n)
  eig_rank_eff <- NULL
  if (identical(eig_solver_mode, "truncated")) {
    eig_rank_eff <- .ipca_effective_rank(eig_rank, n, ncomp)
    if (is.na(eig_rank_eff)) eig_solver_mode <- "full"
  }

  if (is.list(init_state)) {
    if (!is.null(init_state$Sigma_inv) &&
      is.matrix(init_state$Sigma_inv) &&
      all(dim(init_state$Sigma_inv) == c(n, n))) {
      Sigma_inv <- .ipca_sym(init_state$Sigma_inv)
    }
    if (!is.null(init_state$Delta_inv) &&
      is.list(init_state$Delta_inv) &&
      length(init_state$Delta_inv) == K) {
      ok_delta <- all(vapply(seq_len(K), function(k) {
        Dk <- init_state$Delta_inv[[k]]
        is.matrix(Dk) && all(dim(Dk) == c(p_vec[k], p_vec[k]))
      }, logical(1)))
      if (ok_delta) {
        Delta_inv <- lapply(init_state$Delta_inv, .ipca_sym)
      }
    }
    if (!is.null(init_state$delta_inv_norm2) && length(init_state$delta_inv_norm2) == K) {
      delta_inv_norm2 <- as.numeric(init_state$delta_inv_norm2)
    } else if (!is.null(init_state$Delta_inv)) {
      delta_inv_norm2 <- vapply(Delta_inv, function(D) sum(D^2), numeric(1))
    }
    if (!is.null(init_state$phi) && length(init_state$phi) == n) {
      phi <- as.numeric(init_state$phi)
    }
    if (!is.null(init_state$U) && is.matrix(init_state$U) &&
      nrow(init_state$U) == n && ncol(init_state$U) >= ncomp) {
      U <- init_state$U
    }
  }

  eig_solver_used <- NULL

  for (it in seq_len(max_iter)) {
    iter <- it
    phi_prev <- phi
    delta_prev <- delta_inv_norm2

    S <- matrix(0, n, n)
    for (k in seq_len(K)) {
      S <- S + X[[k]] %*% Delta_inv[[k]] %*% t(X[[k]])
    }
    S <- .ipca_sym(S)

    c_sigma <- 8 * p_total * sum(lambda * delta_inv_norm2)
    use_full <- !identical(eig_solver_mode, "truncated")
    if (!use_full) {
      eigS_tr <- tryCatch(
        RSpectra::eigs_sym(S, k = eig_rank_eff, which = "LM"),
        error = function(e) NULL
      )
      if (is.null(eigS_tr) || any(!is.finite(eigS_tr$values))) {
        use_full <- TRUE
      } else {
        ord <- order(eigS_tr$values, decreasing = TRUE)
        U_r <- eigS_tr$vectors[, ord, drop = FALSE]
        gamma_r <- pmax(eigS_tr$values[ord], 0)
        tail_n <- n - length(gamma_r)
        gamma0 <- if (tail_n > 0) {
          pmax((sum(diag(S)) - sum(gamma_r)) / tail_n, 0)
        } else {
          0
        }

        shr_r <- .ipca_shrink_eigs(gamma_r, c_sigma, 2 * p_total)
        if (tail_n > 0) {
          shr0 <- .ipca_shrink_eigs(gamma0, c_sigma, 2 * p_total)
          phi0 <- shr0$eig
          inv_phi0 <- shr0$inv_eig
          phi <- c(shr_r$eig, rep(phi0, tail_n))
          Sigma_inv <- inv_phi0 * diag(n) + U_r %*% ((shr_r$inv_eig - inv_phi0) * t(U_r))
          sigma_inv_norm2 <- sum(shr_r$inv_eig^2) + tail_n * (inv_phi0^2)
        } else {
          phi <- shr_r$eig
          Sigma_inv <- U_r %*% (shr_r$inv_eig * t(U_r))
          sigma_inv_norm2 <- sum(shr_r$inv_eig^2)
        }
        Sigma_inv <- .ipca_sym(Sigma_inv)
        U <- U_r
      }
    }

    if (use_full) {
      eigS <- eigen(S, symmetric = TRUE)
      U <- eigS$vectors
      gamma <- pmax(eigS$values, 0)
      shr_sigma <- .ipca_shrink_eigs(gamma, c_sigma, 2 * p_total)
      phi <- shr_sigma$eig
      inv_phi <- shr_sigma$inv_eig
      Sigma_inv <- U %*% (inv_phi * t(U))
      Sigma_inv <- .ipca_sym(Sigma_inv)
      sigma_inv_norm2 <- sum(inv_phi^2)
    }

    eig_iter <- if (use_full) "full" else "truncated"
    if (is.null(eig_solver_used)) {
      eig_solver_used <- eig_iter
    } else if (!identical(eig_solver_used, eig_iter)) {
      eig_solver_used <- "mixed"
    }

    for (k in seq_len(K)) {
      T_k <- crossprod(X[[k]], Sigma_inv %*% X[[k]])
      T_k <- .ipca_sym(T_k)
      eigT <- eigen(T_k, symmetric = TRUE)
      V_k <- eigT$vectors
      psi <- pmax(eigT$values, 0)

      c_k <- 8 * n * lambda[k] * sigma_inv_norm2
      shr_delta <- .ipca_shrink_eigs(psi, c_k, 2 * n)
      delta_eig <- shr_delta$eig
      inv_delta <- shr_delta$inv_eig

      Delta_inv[[k]] <- V_k %*% (inv_delta * t(V_k))
      Delta_inv[[k]] <- .ipca_sym(Delta_inv[[k]])
      delta_inv_norm2[k] <- sum(inv_delta^2)

      m <- min(ncomp, ncol(V_k))
      loadings[[k]] <- .ipca_pad_cols(V_k[, seq_len(m), drop = FALSE], ncomp)
      delta_leading[[k]] <- .ipca_pad_vec(delta_eig[seq_len(m)], ncomp)
    }

    if (normalize_trace) {
      c_scale <- mean(phi)
      if (is.finite(c_scale) && c_scale > 0) {
        phi <- phi / c_scale
        Sigma_inv <- Sigma_inv * c_scale
        for (k in seq_len(K)) {
          Delta_inv[[k]] <- Delta_inv[[k]] / c_scale
          delta_inv_norm2[k] <- delta_inv_norm2[k] / (c_scale^2)
          delta_leading[[k]] <- delta_leading[[k]] * c_scale
        }
      }
    }

    rel_phi <- .ipca_rel_change(phi, phi_prev)
    rel_delta <- .ipca_rel_change(delta_inv_norm2, delta_prev)
    if (rel_phi < tol && rel_delta < tol) {
      converged <- TRUE
      break
    }
  }

  list(
    U = U,
    phi = phi,
    loadings = loadings,
    delta_leading = delta_leading,
    converged = converged,
    iter = iter,
    rel_change_phi = rel_phi,
    rel_change_delta = rel_delta,
    eig_solver_used = eig_solver_used %||% "full",
    warm_state = if (isTRUE(return_state)) {
      list(
        Sigma_inv = Sigma_inv,
        Delta_inv = Delta_inv,
        delta_inv_norm2 = delta_inv_norm2,
        phi = phi,
        U = U
      )
    } else {
      NULL
    }
  )
}


.ipca_update_gram_block <- function(G_k,
                                    p_k,
                                    U,
                                    phi,
                                    lambda_k,
                                    sigma_inv_norm2,
                                    n,
                                    ncomp) {
  sqrt_phi <- sqrt(phi)
  inv_sqrt_phi <- 1 / sqrt_phi

  B_k <- crossprod(U, G_k %*% U)
  M_k <- (inv_sqrt_phi %o% inv_sqrt_phi) * B_k
  M_k <- .ipca_sym(M_k)

  eigM <- eigen(M_k, symmetric = TRUE)
  q_tilde <- eigM$vectors
  s <- pmax(eigM$values, 0)

  c_k <- 8 * n * lambda_k * sigma_inv_norm2
  shr_all <- .ipca_shrink_eigs(s, c_k, 2 * n)
  inv_delta_all <- shr_all$inv_eig

  # S_k = Sigma^{1/2} Q diag(s * inv_delta) Q^T Sigma^{1/2}
  weights <- s * inv_delta_all
  Qw <- sweep(q_tilde, 2, sqrt(pmax(weights, 0)), `*`)
  C_k <- tcrossprod(Qw)
  S_k <- U %*% (((sqrt_phi %o% sqrt_phi) * C_k) %*% t(U))
  S_k <- .ipca_sym(S_k)

  m_k <- min(n, p_k)
  shr_m <- .ipca_shrink_eigs(s[seq_len(m_k)], c_k, 2 * n)
  inv_m <- shr_m$inv_eig

  extra_zeros <- max(p_k - n, 0L)
  inv_delta0 <- 2 * n / sqrt(c_k)
  delta_inv_norm2 <- sum(inv_m^2) + extra_zeros * (inv_delta0^2)

  nlead <- min(ncomp, m_k)
  delta_leading <- .ipca_pad_vec(shr_m$eig[seq_len(nlead)], ncomp)

  list(
    S_k = S_k,
    delta_inv_norm2 = delta_inv_norm2,
    q_tilde = q_tilde,
    s = s,
    delta_leading = delta_leading
  )
}


.ipca_extract_loadings_gram <- function(X, U, phi, q_tilde, s, ncomp) {
  inv_sqrt_phi <- 1 / sqrt(phi)
  m <- min(ncomp, ncol(X), length(s), nrow(q_tilde))
  if (m < 1) return(matrix(0, ncol(X), ncomp))

  s_top <- pmax(s[seq_len(m)], 0)
  q_top <- q_tilde[, seq_len(m), drop = FALSE]

  B <- U %*% sweep(q_top, 1, inv_sqrt_phi, `*`)
  V <- crossprod(X, B)

  for (j in seq_len(m)) {
    if (s_top[j] > 1e-12) {
      V[, j] <- V[, j] / sqrt(s_top[j])
    } else {
      V[, j] <- 0
    }
  }

  .ipca_pad_cols(V, ncomp)
}


.ipca_fit_gram <- function(X,
                           p_vec,
                           lambda,
                           ncomp,
                           max_iter,
                           tol,
                           normalize_trace,
                           use_future,
                           eig_solver,
                           eig_rank,
                           eig_trunc_min_n,
                           init_state = NULL,
                           return_state = FALSE) {
  K <- length(X)
  n <- nrow(X[[1]])
  p_total <- sum(p_vec)

  G <- lapply(X, function(x) {
    g <- tcrossprod(x)
    .ipca_sym(g)
  })

  S_k <- G
  delta_inv_norm2 <- as.numeric(p_vec)
  phi <- rep(1, n)
  U <- diag(n)

  if (is.list(init_state)) {
    if (!is.null(init_state$S_k) &&
      is.list(init_state$S_k) &&
      length(init_state$S_k) == K) {
      ok_sk <- all(vapply(init_state$S_k, function(Si) {
        is.matrix(Si) && all(dim(Si) == c(n, n))
      }, logical(1)))
      if (ok_sk) S_k <- lapply(init_state$S_k, .ipca_sym)
    }
    if (!is.null(init_state$delta_inv_norm2) && length(init_state$delta_inv_norm2) == K) {
      delta_inv_norm2 <- as.numeric(init_state$delta_inv_norm2)
    }
    if (!is.null(init_state$phi) && length(init_state$phi) == n) {
      phi <- as.numeric(init_state$phi)
    }
    if (!is.null(init_state$U) &&
      is.matrix(init_state$U) &&
      nrow(init_state$U) == n &&
      ncol(init_state$U) == n) {
      U <- init_state$U
    }
  }

  q_tilde_last <- lapply(seq_len(K), function(...) diag(n))
  s_last <- lapply(seq_len(K), function(...) rep(0, n))
  delta_leading <- lapply(seq_len(K), function(...) rep(0, ncomp))
  loadings <- lapply(p_vec, function(p) matrix(0, p, ncomp))

  map_fun <- if (isTRUE(use_future)) {
    function(.x, .f) {
      furrr::future_map(.x, .f, .options = furrr::furrr_options(seed = TRUE))
    }
  } else {
    function(.x, .f) lapply(.x, .f)
  }

  converged <- FALSE
  rel_phi <- Inf
  rel_delta <- Inf
  iter <- 0L

  for (it in seq_len(max_iter)) {
    iter <- it
    phi_prev <- phi
    delta_prev <- delta_inv_norm2

    S <- Reduce(`+`, S_k)
    S <- .ipca_sym(S)

    eigS <- eigen(S, symmetric = TRUE)
    U <- eigS$vectors
    gamma <- pmax(eigS$values, 0)

    c_sigma <- 8 * p_total * sum(lambda * delta_inv_norm2)
    shr_sigma <- .ipca_shrink_eigs(gamma, c_sigma, 2 * p_total)
    phi <- shr_sigma$eig
    inv_phi <- shr_sigma$inv_eig
    sigma_inv_norm2 <- sum(inv_phi^2)

    updates <- map_fun(seq_len(K), function(k) {
      .ipca_update_gram_block(
        G_k = G[[k]],
        p_k = p_vec[k],
        U = U,
        phi = phi,
        lambda_k = lambda[k],
        sigma_inv_norm2 = sigma_inv_norm2,
        n = n,
        ncomp = ncomp
      )
    })

    for (k in seq_len(K)) {
      S_k[[k]] <- updates[[k]]$S_k
      delta_inv_norm2[k] <- updates[[k]]$delta_inv_norm2
      q_tilde_last[[k]] <- updates[[k]]$q_tilde
      s_last[[k]] <- updates[[k]]$s
      delta_leading[[k]] <- updates[[k]]$delta_leading
    }

    if (normalize_trace) {
      c_scale <- mean(phi)
      if (is.finite(c_scale) && c_scale > 0) {
        phi <- phi / c_scale
        for (k in seq_len(K)) {
          S_k[[k]] <- S_k[[k]] / c_scale
          delta_inv_norm2[k] <- delta_inv_norm2[k] / (c_scale^2)
          delta_leading[[k]] <- delta_leading[[k]] * c_scale
        }
      }
    }

    rel_phi <- .ipca_rel_change(phi, phi_prev)
    rel_delta <- .ipca_rel_change(delta_inv_norm2, delta_prev)
    if (rel_phi < tol && rel_delta < tol) {
      converged <- TRUE
      break
    }
  }

  for (k in seq_len(K)) {
    loadings[[k]] <- .ipca_extract_loadings_gram(
      X = X[[k]],
      U = U,
      phi = phi,
      q_tilde = q_tilde_last[[k]],
      s = s_last[[k]],
      ncomp = ncomp
    )
  }

  list(
    U = U,
    phi = phi,
    loadings = loadings,
    delta_leading = delta_leading,
    converged = converged,
    iter = iter,
    rel_change_phi = rel_phi,
    rel_change_delta = rel_delta,
    eig_solver_used = "full",
    warm_state = if (isTRUE(return_state)) {
      list(
        S_k = S_k,
        delta_inv_norm2 = delta_inv_norm2,
        phi = phi,
        U = U
      )
    } else {
      NULL
    }
  )
}


.ipca_lambda_from_alpha <- function(alpha, p_vec, tie = c("inv_p", "inv_pbar", "equal")) {
  tie <- match.arg(tie)
  chk::chk_number(alpha)
  chk::chk_true(is.finite(alpha) && alpha > 0)
  p_vec <- as.numeric(p_vec)
  chk::chk_true(all(is.finite(p_vec)))
  chk::chk_true(all(p_vec > 0))

  if (tie == "inv_p") {
    alpha / p_vec
  } else if (tie == "inv_pbar") {
    alpha * mean(p_vec) / p_vec
  } else {
    rep(alpha, length(p_vec))
  }
}


.ipca_make_masked_training <- function(X, mask) {
  Xt <- X
  for (j in seq_len(ncol(Xt))) {
    mj <- mask[, j]
    if (any(mj)) {
      obs <- !mj
      mu <- if (any(obs)) mean(Xt[obs, j]) else 0
      Xt[mj, j] <- mu
    }
  }
  Xt
}


#' Tune Tied Penalties for iPCA via 1D Alpha Search
#'
#' @description
#' Performs a lightweight 1D search over `alpha` and ties block penalties as
#' `lambda_k = alpha / p_k` (or related schemes), then selects the value with
#' the lowest held-out entry reconstruction MSE.
#'
#' @param data A list of matrices/data.frames or a `multiblock` object. Blocks
#'   must share rows.
#' @param preproc A `multivarious` preprocessing pipeline.
#' @param ncomp Integer; number of joint components to evaluate.
#' @param alpha_grid Positive numeric vector of candidate alpha values.
#' @param tie Penalty tying rule: `"inv_p"` (`alpha / p_k`), `"inv_pbar"`
#'   (`alpha * mean(p) / p_k`), or `"equal"` (`alpha` for all blocks).
#' @param holdout_frac Fraction of entries to hold out in each block.
#' @param warm_start Logical; if `TRUE`, alpha candidates are fit in ascending
#'   order and each fit is initialized from the previous alpha within a mask.
#' @param method One of `"auto"`, `"gram"`, or `"dense"` for iPCA fits.
#' @param max_iter Maximum number of iPCA iterations per alpha.
#' @param tol Convergence tolerance for iPCA fits.
#' @param normalize_trace Logical; passed to [ipca()].
#' @param use_future Logical; passed to [ipca()].
#' @param eig_solver Eigensolver policy passed to [ipca()].
#' @param eig_rank Optional truncated eigensolver rank passed to [ipca()].
#' @param eig_trunc_min_n Minimum `n` at which `eig_solver = "auto"` can switch
#'   to truncated eigendecomposition in dense mode.
#' @param seed Integer random seed for holdout mask generation.
#' @param verbose Logical; print per-alpha progress.
#' @param ... Additional arguments passed to [ipca()].
#'
#' @return A list with:
#' \describe{
#'   \item{best_alpha}{Selected alpha value.}
#'   \item{best_lambda}{Selected lambda vector of length K.}
#'   \item{results}{Data frame of alpha candidates and held-out MSE.}
#'   \item{fit}{iPCA fit refit on full data with `best_alpha`.}
#' }
#'
#' @export
ipca_tune_alpha <- function(data,
                            preproc = multivarious::center(),
                            ncomp = 2,
                            alpha_grid = 10^seq(-3, 3, by = 1),
                            tie = c("inv_p", "inv_pbar", "equal"),
                            holdout_frac = 0.05,
                            n_masks = 1,
                            warm_start = TRUE,
                            method = c("auto", "gram", "dense"),
                            max_iter = 100,
                            tol = 1e-6,
                            normalize_trace = TRUE,
                            use_future = FALSE,
                            eig_solver = c("auto", "full", "truncated"),
                            eig_rank = NULL,
                            eig_trunc_min_n = 400,
                            seed = 1,
                            verbose = FALSE,
                            ...) {
  tie <- match.arg(tie)
  method <- match.arg(method)
  eig_solver <- match.arg(eig_solver)
  chk::chk_flag(verbose)
  chk::chk_flag(warm_start)
  chk::chk_number(holdout_frac)
  chk::chk_true(is.finite(holdout_frac) && holdout_frac > 0 && holdout_frac < 1)
  chk::chk_number(n_masks)
  chk::chk_gte(n_masks, 1)
  if (abs(n_masks - round(n_masks)) > 1e-8) {
    stop("n_masks must be an integer-like value.", call. = FALSE)
  }
  n_masks <- as.integer(round(n_masks))
  if (!is.null(eig_rank)) {
    chk::chk_number(eig_rank)
    chk::chk_gte(eig_rank, ncomp)
    eig_rank <- as.integer(round(eig_rank))
  }
  chk::chk_number(eig_trunc_min_n)
  chk::chk_gte(eig_trunc_min_n, 2)
  eig_trunc_min_n <- as.integer(round(eig_trunc_min_n))

  alpha_grid <- as.numeric(alpha_grid)
  chk::chk_true(length(alpha_grid) >= 1)
  chk::chk_true(all(is.finite(alpha_grid)))
  chk::chk_true(all(alpha_grid > 0))

  X <- lapply(data, function(x) {
    chk::chk_true(is.matrix(x) || is.data.frame(x))
    as.matrix(x)
  })
  if (length(X) <= 1) stop("ipca_tune_alpha requires at least two blocks.", call. = FALSE)

  nrs <- vapply(X, nrow, integer(1))
  if (!all(nrs == nrs[1])) stop("All blocks must have the same number of rows.", call. = FALSE)
  p_vec <- vapply(X, ncol, integer(1))

  lambda_grid <- t(vapply(alpha_grid, function(a) .ipca_lambda_from_alpha(a, p_vec, tie), numeric(length(p_vec))))
  colnames(lambda_grid) <- names(X) %||% paste0("B", seq_along(X))

  results <- data.frame(
    alpha = alpha_grid,
    holdout_mse = Inf,
    n_success = 0L,
    converged_rate = NA_real_,
    converged = FALSE,
    iter = NA_integer_,
    method_used = NA_character_,
    stringsAsFactors = FALSE
  )
  holdout_sum <- rep(0, length(alpha_grid))
  iter_sum <- rep(0, length(alpha_grid))
  converged_sum <- rep(0, length(alpha_grid))
  method_first <- rep(NA_character_, length(alpha_grid))

  for (m in seq_len(n_masks)) {
    set.seed(seed + m - 1L)
    masks <- lapply(X, function(x) {
      mm <- matrix(stats::runif(length(x)) < holdout_frac, nrow = nrow(x), ncol = ncol(x))
      if (!any(mm)) {
        mm[sample.int(length(x), 1)] <- TRUE
      }
      mm
    })
    X_train <- Map(.ipca_make_masked_training, X, masks)

    prep_cv <- prepare_block_preprocessors(X_train, preproc, check_consistent_ncol = FALSE)
    X_train_fit <- lapply(prep_cv$Xp, function(x) if (inherits(x, "Matrix")) as.matrix(x) else x)
    if (is.null(prep_cv$proclist)) {
      X_eval <- X
    } else {
      X_eval <- lapply(seq_along(X), function(k) {
        multivarious::transform(prep_cv$proclist[[k]], X[[k]])
      })
    }

    alpha_order <- if (isTRUE(warm_start)) order(alpha_grid) else seq_along(alpha_grid)
    warm_state <- NULL
    for (ii in seq_along(alpha_order)) {
      i <- alpha_order[ii]
      lam <- lambda_grid[i, ]
      fit_i <- tryCatch(
        ipca(
          data = X_train_fit,
          preproc = multivarious::pass(),
          ncomp = ncomp,
          lambda = lam,
          method = method,
          max_iter = max_iter,
          tol = tol,
          normalize_trace = normalize_trace,
          use_future = use_future,
          eig_solver = eig_solver,
          eig_rank = eig_rank,
          eig_trunc_min_n = eig_trunc_min_n,
          .init_state = if (isTRUE(warm_start)) warm_state else NULL,
          .return_state = isTRUE(warm_start),
          ...
        ),
        error = function(e) e
      )

      if (inherits(fit_i, "error")) {
        if (isTRUE(verbose)) {
          message(sprintf("mask=%d alpha=%.4g failed: %s", m, alpha_grid[i], conditionMessage(fit_i)))
        }
        next
      }
      if (isTRUE(warm_start) && !is.null(fit_i$warm_state)) {
        warm_state <- fit_i$warm_state
      }

      S <- multivarious::scores(fit_i)
      mse_num <- 0
      mse_den <- 0
      for (k in seq_along(X)) {
        idx <- fit_i$block_indices[[k]]
        V_k <- fit_i$v[idx, , drop = FALSE]
        Xhat_k <- S %*% t(V_k)
        mk <- masks[[k]]
        diff <- X_eval[[k]][mk] - Xhat_k[mk]
        mse_num <- mse_num + sum(diff^2)
        mse_den <- mse_den + length(diff)
      }
      mse <- mse_num / max(mse_den, 1)

      holdout_sum[i] <- holdout_sum[i] + mse
      results$n_success[i] <- results$n_success[i] + 1L
      converged_sum[i] <- converged_sum[i] + as.integer(isTRUE(fit_i$converged))
      iter_sum[i] <- iter_sum[i] + fit_i$iter
      if (is.na(method_first[i])) {
        method_first[i] <- fit_i$method_used
      } else if (!identical(method_first[i], fit_i$method_used)) {
        method_first[i] <- "mixed"
      }

      if (isTRUE(verbose)) {
        message(sprintf("mask=%d alpha=%.4g mse=%.6f iter=%d", m, alpha_grid[i], mse, fit_i$iter))
      }
    }
  }

  success <- results$n_success > 0L
  results$holdout_mse[success] <- holdout_sum[success] / results$n_success[success]
  results$converged_rate[success] <- converged_sum[success] / results$n_success[success]
  results$converged[success] <- converged_sum[success] == results$n_success[success]
  results$iter[success] <- as.integer(round(iter_sum[success] / results$n_success[success]))
  results$method_used <- method_first

  finite_idx <- which(success & is.finite(results$holdout_mse))
  if (length(finite_idx) == 0L) {
    stop("All alpha candidates failed during tuning.", call. = FALSE)
  }

  best_idx <- which.min(results$holdout_mse)
  best_alpha <- results$alpha[best_idx]
  best_lambda <- as.numeric(lambda_grid[best_idx, ])
  names(best_lambda) <- colnames(lambda_grid)

  best_fit <- ipca(
    data = X,
    preproc = preproc,
    ncomp = ncomp,
    lambda = best_lambda,
    method = method,
    max_iter = max_iter,
    tol = tol,
    normalize_trace = normalize_trace,
    use_future = use_future,
    eig_solver = eig_solver,
    eig_rank = eig_rank,
    eig_trunc_min_n = eig_trunc_min_n,
    ...
  )

  out <- list(
    best_alpha = best_alpha,
    best_lambda = best_lambda,
    tie = tie,
    holdout_frac = holdout_frac,
    warm_start = warm_start,
    results = results,
    lambda_grid = lambda_grid,
    fit = best_fit
  )
  class(out) <- "ipca_alpha_tuning"
  out
}

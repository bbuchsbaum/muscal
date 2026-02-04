## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 6,
  fig.height = 4
)

## ----basic-example, eval = FALSE----------------------------------------------
# library(muscal)
# library(multivarious)
# 
# # Simulate three blocks of data on the same 50 observations
# set.seed(42)
# n <- 50
# block1 <- matrix(rnorm(n * 20), n, 20)  # 20 variables
# block2 <- matrix(rnorm(n * 15), n, 15)  # 15 variables
# block3 <- matrix(rnorm(n * 25), n, 25)  # 25 variables
# 
# # Fit MFA with 3 components
# fit <- mfa(list(block1, block2, block3), ncomp = 3)
# 
# # Extract global scores
# S <- scores(fit)
# dim(S)  # 50 x 3

## ----normalization, eval = FALSE----------------------------------------------
# # RV-based normalization
# fit_rv <- mfa(list(block1, block2, block3), ncomp = 3, normalization = "RV")
# 
# # Custom weights
# ncols <- c(20, 15, 25)
# custom_weights <- rep(1/ncols, ncols)  # inverse-column-count weighting
# fit_custom <- mfa(list(block1, block2, block3),
#                   ncomp = 3,
#                   normalization = "custom",
#                   A = custom_weights)

## ----results, eval = FALSE----------------------------------------------------
# # Global scores (observations in the compromise space)
# S <- scores(fit)
# 
# # Loadings in the concatenated variable space (sum(p_k) x ncomp)
# V <- fit$v
# 
# # Block-specific loadings are slices of V using block_indices
# V_block1 <- V[fit$block_indices[[1]], , drop = FALSE]
# 
# # Project a single block onto the fitted model
# proj <- project_block(fit, block1, block = 1)
# 
# # Partial factor scores: how each block "sees" the observations
# # (useful for assessing block agreement)
# partial_scores <- project_block(fit, block2, block = 2)

## ----alpha, eval = FALSE------------------------------------------------------
# fit$alpha


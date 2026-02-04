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
# set.seed(123)
# 
# # Reference block: 100 subjects, 20 variables
# N <- 100
# Y <- matrix(rnorm(N * 20), N, 20)
# 
# # Auxiliary block 1: 60 subjects (subset), 30 variables
# X1 <- matrix(rnorm(60 * 30), 60, 30)
# idx1 <- sample.int(N, 60, replace = FALSE)  # which reference subjects
# 
# # Auxiliary block 2: 45 subjects (different subset), 25 variables
# X2 <- matrix(rnorm(45 * 25), 45, 25)
# idx2 <- sample.int(N, 45, replace = FALSE)
# 
# # Fit linked MFA
# fit <- linked_mfa(
#   Y = Y,
#   X = list(X1 = X1, X2 = X2),
#   row_index = list(X1 = idx1, X2 = idx2),
#   ncomp = 3
# )
# 
# # Scores are defined for all N reference subjects
# S <- scores(fit)
# dim(S)  # 100 x 3

## ----mapping, eval = FALSE----------------------------------------------------
# # row_index[[k]][i] says: "row i of X[[k]] corresponds to row row_index[[k]][i] of Y"
# 
# # Example: X1 has 60 rows, mapping to specific Y rows
# idx1 <- sample.int(N, 60, replace = FALSE)  # length 60, values in 1..N
# 
# # Subjects can appear multiple times (repeated measures)
# idx_repeated <- c(rep(1L, 3), rep(2L, 2), 3L)  # subject 1 measured 3 times

## ----feature-groups, eval = FALSE---------------------------------------------
# # Automatic grouping by column name
# # Features with the same name across blocks are grouped together
# fit <- linked_mfa(
#   Y = Y,
#   X = list(X1 = X1, X2 = X2),
#   row_index = list(X1 = idx1, X2 = idx2),
#   ncomp = 3,
#   feature_groups = "colnames",
#   feature_lambda = 0.1
# )
# 
# # Manual grouping via data frame
# groups_df <- data.frame(
#   block = c("X1", "X1", "X2", "X2"),
#   feature = c(1, 5, 1, 3),       # feature indices within each block
#   group = c("ROI_A", "ROI_B", "ROI_A", "ROI_B"),
#   weight = c(1, 1, 1, 1)
# )
# 
# fit <- linked_mfa(
#   Y = Y,
#   X = list(X1 = X1, X2 = X2),
#   row_index = list(X1 = idx1, X2 = idx2),
#   ncomp = 3,
#   feature_groups = groups_df,
#   feature_lambda = 0.5
# )

## ----normalization, eval = FALSE----------------------------------------------
# # Custom weights: emphasize the reference block
# fit <- linked_mfa(
#   Y = Y,
#   X = list(X1 = X1, X2 = X2),
#   row_index = list(X1 = idx1, X2 = idx2),
#   ncomp = 3,
#   normalization = "custom",
#   alpha = c(2, 1, 1)  # Y weight, X1 weight, X2 weight
# )

## ----verbose, eval = FALSE----------------------------------------------------
# fit <- linked_mfa(
#   Y = Y,
#   X = list(X1 = X1, X2 = X2),
#   row_index = list(X1 = idx1, X2 = idx2),
#   ncomp = 3,
#   max_iter = 100,
#   tol = 1e-6,
#   verbose = TRUE
# )
# 
# # Inspect convergence
# plot(fit$objective_trace, type = "l",
#      xlab = "Iteration", ylab = "Objective")

## ----results, eval = FALSE----------------------------------------------------
# # Global scores (N x ncomp)
# S <- scores(fit)
# 
# # Concatenated loadings (rows correspond to Y then each X block)
# V <- fit$v
# 
# # Block-specific loadings
# fit$B        # Y loadings (q x ncomp)
# fit$V_list   # list of X_k loadings
# 
# # Block indices in concatenated space
# fit$block_indices
# 
# # Row mappings
# fit$row_index
# 
# # Convergence trace
# fit$objective_trace

## ----neuroimaging, eval = FALSE-----------------------------------------------
# # Y: structural MRI on 200 subjects (all have this)
# # X1: task fMRI session 1 (150 subjects completed)
# # X2: task fMRI session 2 (120 subjects completed)
# # X3: resting state (180 subjects completed)
# 
# fit <- linked_mfa(
#   Y = structural_data,
#   X = list(task1 = task1_data, task2 = task2_data, rest = rest_data),
#   row_index = list(task1 = task1_subjects,
#                    task2 = task2_subjects,
#                    rest = rest_subjects),
#   ncomp = 5,
#   feature_groups = "colnames",  # same ROIs across sessions
# 
#   feature_lambda = 0.2
# )
# 
# # Scores are defined for all 200 subjects,
# # integrating information from whatever sessions each completed


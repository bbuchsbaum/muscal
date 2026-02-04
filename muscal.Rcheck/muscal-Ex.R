pkgname <- "muscal"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('muscal')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("bamfa")
### * bamfa

flush(stderr()); flush(stdout())

### Name: bamfa
### Title: Barycentric Multiple Factor Analysis (BaMFA)
### Aliases: bamfa bamfa.default bamfa.list bamfa.multiblock
###   bamfa.multidesign

### ** Examples

# Generate example multi-block data (e.g., 3 subjects, 10 conditions, 50 features)
set.seed(123)
n_obs <- 10
n_features <- 50
n_subjects <- 3
data_list <- lapply(1:n_subjects, function(i) {
  matrix(rnorm(n_obs * n_features), n_obs, n_features) +
  matrix(rnorm(n_obs * 1, mean=i), n_obs, n_features) # Add subject offset
})
names(data_list) <- paste0("Subject_", 1:n_subjects)

# Run BaMFA with k_g=3 global, k_l=2 local components
result <- bamfa(data_list, k_g = 3, k_l = 2, niter = 10)
print(result)




cleanEx()
nameEx("check_duality")
### * check_duality

flush(stderr()); flush(stdout())

### Name: check_duality
### Title: Verify duality between ROI- and RV-space coordinates
### Aliases: check_duality

### ** Examples

# Create a list of correlation matrices
Xs <- lapply(1:5, function(i) matrix(rnorm(10*10), 10, 10))
Xs <- lapply(Xs, cor)
# Apply COVSTATIS
res <- covstatis(Xs, ncomp=3)
# Check duality
check_duality(res)



cleanEx()
nameEx("covstatis")
### * covstatis

flush(stderr()); flush(stdout())

### Name: covstatis
### Title: STATIS for Covariance Matrices (generic)
### Aliases: covstatis covstatis.list

### ** Examples

# Create a list of correlation matrices
Xs <- lapply(1:5, function(i) matrix(rnorm(10*10), 10, 10))
Xs <- lapply(Xs, cor)

# Apply COVSTATIS
res <- covstatis(Xs, ncomp=3)

# Project a new correlation matrix
new_mat <- cor(matrix(rnorm(10*10), 10, 10))
proj <- project_cov(res, new_mat)




cleanEx()
nameEx("linked_mfa")
### * linked_mfa

flush(stderr()); flush(stdout())

### Name: linked_mfa
### Title: Linked Multiple Factor Analysis (Linked MFA)
### Aliases: linked_mfa

### ** Examples




cleanEx()
nameEx("mfa")
### * mfa

flush(stderr()); flush(stdout())

### Name: mfa
### Title: Multiple Factor Analysis (generic)
### Aliases: mfa mfa.list mfa.multiblock

### ** Examples




cleanEx()
nameEx("penalized_mfa")
### * penalized_mfa

flush(stderr()); flush(stdout())

### Name: penalized_mfa
### Title: Penalized Multiple Factor Analysis (MFA)
### Aliases: penalized_mfa penalized_mfa.list penalized_mfa.multiblock
###   penalized_mfa.multidesign

### ** Examples

## Not run: 
##D # Example 1: Basic usage with simulated data
##D set.seed(123)
##D data_list <- lapply(1:3, function(i) {
##D   matrix(rnorm(100), 10, 10)
##D })
##D res <- penalized_mfa(data_list, ncomp=2, lambda=1, penalty_method="projection",
##D                      optimizer="adam", max_iter=50, verbose=TRUE)
##D print(res)
##D 
##D # Access block-specific loadings
##D V_list <- attr(res, "V_list")
##D print(V_list[[1]])  # Loadings for first block
##D 
##D # Plot convergence
##D obj_vals <- attr(res, "obj_values")
##D plot(obj_vals, type='b', xlab='Iteration', ylab='Objective',
##D      main='Convergence of Penalized MFA')
##D 
##D # Example 2: With consensus and custom preprocessing
##D library(multivarious)
##D res2 <- penalized_mfa(data_list, ncomp=3, lambda=2,
##D                       preproc=standardize(),  # Center and scale
##D                       compute_consensus=TRUE,
##D                       verbose=FALSE)
##D consensus_loadings <- attr(res2, "consensus")
##D 
##D # Example 3: Comparing penalty methods
##D res_proj <- penalized_mfa(data_list, ncomp=2, lambda=1, penalty_method="projection")
##D res_mean <- penalized_mfa(data_list, ncomp=2, lambda=1, penalty_method="global_mean")
##D res_pair <- penalized_mfa(data_list, ncomp=2, lambda=1, penalty_method="pairwise")
##D 
##D # Example 4: Lambda selection via objective values
##D lambdas <- c(0, 0.1, 0.5, 1, 2, 5, 10)
##D results <- lapply(lambdas, function(lam) {
##D   fit <- penalized_mfa(data_list, ncomp=2, lambda=lam, verbose=FALSE)
##D   list(lambda=lam, final_obj=tail(attr(fit, "obj_values"), 1))
##D })
##D 
##D # Example 5: Using with multiblock object
##D mb <- multiblock(data_list)
##D res_mb <- penalized_mfa(mb, ncomp=2, lambda=1)
##D 
##D # Example 6: Different preprocessors per block
##D preproc_list <- list(
##D   center(),
##D   standardize(),
##D   pass()  # No preprocessing for third block
##D )
##D res_custom <- penalized_mfa(data_list, ncomp=2, lambda=1,
##D                             preproc=preproc_list)
##D 
##D # Example 7: Gradient descent instead of Adam
##D res_gd <- penalized_mfa(data_list, ncomp=2, lambda=1,
##D                         optimizer="gradient",
##D                         learning_rate=0.001,
##D                         max_iter=100)
## End(Not run)




cleanEx()
nameEx("penalized_mfa_clusterwise")
### * penalized_mfa_clusterwise

flush(stderr()); flush(stdout())

### Name: penalized_mfa_clusterwise
### Title: Penalized MFA with Clusterwise Spatial Smoothness Constraints
### Aliases: penalized_mfa_clusterwise pmfa_cluster

### ** Examples

## Not run: 
##D # Example 1: Basic usage with simulated spatial cluster data
##D set.seed(123)
##D S <- 3  # 3 subjects
##D 
##D # Generate data with spatial structure
##D data_list <- lapply(1:S, function(s) {
##D   n <- 50  # 50 observations
##D   k <- 20  # 20 clusters
##D   matrix(rnorm(n * k), n, k)
##D })
##D 
##D # Generate 3D spatial coordinates for clusters
##D coords_list <- lapply(1:S, function(s) {
##D   matrix(runif(20 * 3, 0, 10), 20, 3)  # Random 3D positions
##D })
##D 
##D # Fit model with spatial smoothness
##D res <- penalized_mfa_clusterwise(
##D   data_list, coords_list,
##D   ncomp = 3,
##D   lambda = 1,
##D   max_iter = 20,
##D   verbose = TRUE
##D )
##D print(res)
##D 
##D # Plot convergence
##D plot(res$obj_values, type = 'b', xlab = 'Iteration', ylab = 'Objective',
##D      main = 'Convergence of Spatially-Regularized MFA')
##D 
##D # Example 2: Using the shorter alias
##D res2 <- pmfa_cluster(data_list, coords_list, ncomp = 2, lambda = 0.5)
##D 
##D # Example 3: Compare different lambda values
##D lambdas <- c(0, 0.1, 0.5, 1, 2, 5)
##D results <- lapply(lambdas, function(lam) {
##D   fit <- pmfa_cluster(data_list, coords_list, ncomp = 2, lambda = lam,
##D                       verbose = FALSE)
##D   list(
##D     lambda = lam,
##D     final_obj = tail(fit$obj_values, 1),
##D     iterations = fit$iterations_run
##D   )
##D })
##D 
##D # Example 4: Using Adam optimizer for faster convergence
##D res_adam <- pmfa_cluster(
##D   data_list, coords_list,
##D   ncomp = 3,
##D   lambda = 1,
##D   optimizer = "adam",
##D   learning_rate = 0.05,
##D   max_iter = 30,
##D   verbose = TRUE
##D )
##D 
##D # Example 5: Controlling k-NN graph construction
##D res_dense <- pmfa_cluster(
##D   data_list, coords_list,
##D   ncomp = 2,
##D   lambda = 1,
##D   adjacency_opts = list(k_nn = 10),  # More neighbors = denser graph
##D   verbose = TRUE
##D )
##D 
##D # Example 6: Using unnormalized Laplacian
##D res_unnorm <- pmfa_cluster(
##D   data_list, coords_list,
##D   ncomp = 2,
##D   lambda = 1,
##D   normalized_laplacian = FALSE
##D )
##D 
##D # Example 7: Memory-constrained settings
##D # For large datasets, reduce memory budget
##D res_mem <- pmfa_cluster(
##D   data_list, coords_list,
##D   ncomp = 2,
##D   lambda = 1,
##D   memory_budget_mb = 256,  # Limit memory per block
##D   verbose = TRUE
##D )
##D # Check which blocks used precomputed gradients
##D print(res_mem$precompute_info)
##D 
##D # Example 8: Extracting block-specific loadings
##D V_list <- res$V_list
##D 
##D # Loadings for first subject
##D V1 <- V_list[[1]]
##D dim(V1)  # k_s x ncomp
##D 
##D # Visualize spatial smoothness of first loading vector
##D library(scatterplot3d)
##D scatterplot3d(
##D   coords_list[[1]],
##D   color = rank(V1[, 1]),
##D   main = "Spatial Pattern of First Loading (Subject 1)",
##D   xlab = "X", ylab = "Y", zlab = "Z"
##D )
##D 
##D # Example 9: Examining the spatial penalty contribution
##D # Extract Laplacian and loadings
##D L <- res$Sadj
##D LV <- res$LV
##D v <- res$v
##D 
##D # Compute spatial penalty: tr(V' L V)
##D spatial_penalty <- sum(LV * v)
##D cat("Spatial penalty term:", spatial_penalty, "\n")
##D 
##D # Example 10: Variable-rank loadings (blocks with different sizes)
##D # Simulate data with varying cluster numbers
##D data_var <- list(
##D   matrix(rnorm(50 * 15), 50, 15),  # 15 clusters
##D   matrix(rnorm(50 * 20), 50, 20),  # 20 clusters
##D   matrix(rnorm(50 * 10), 50, 10)   # 10 clusters
##D )
##D 
##D coords_var <- list(
##D   matrix(runif(15 * 3), 15, 3),
##D   matrix(runif(20 * 3), 20, 3),
##D   matrix(runif(10 * 3), 10, 3)
##D )
##D 
##D res_var <- pmfa_cluster(
##D   data_var, coords_var,
##D   ncomp = 8,  # Will be capped at 10 (smallest block)
##D   lambda = 1,
##D   verbose = TRUE
##D )
##D # Check effective ranks per block
##D print(res_var$ncomp_block)
##D 
##D # Example 11: Custom preprocessing per block
##D library(multivarious)
##D 
##D preproc_list <- list(
##D   center(),        # Just center first block
##D   standardize(),   # Center and scale second block
##D   center()         # Just center third block
##D )
##D 
##D res_preproc <- pmfa_cluster(
##D   data_list, coords_list,
##D   ncomp = 2,
##D   lambda = 1,
##D   preproc = preproc_list
##D )
##D 
##D # Example 12: Lambda selection via cross-validation style approach
##D # (Simplified - would need proper CV in practice)
##D lambda_grid <- 10^seq(-1, 1, length.out = 10)
##D cv_results <- data.frame(
##D   lambda = lambda_grid,
##D   recon_error = NA,
##D   spatial_penalty = NA,
##D   total_obj = NA
##D )
##D 
##D for (i in seq_along(lambda_grid)) {
##D   fit <- pmfa_cluster(data_list, coords_list, ncomp = 2,
##D                       lambda = lambda_grid[i], verbose = FALSE)
##D   cv_results$total_obj[i] <- tail(fit$obj_values, 1)
##D   # Could decompose into reconstruction and spatial components
##D }
##D 
##D # Plot objective vs lambda
##D plot(cv_results$lambda, cv_results$total_obj, log = "x",
##D      type = "b", xlab = "Lambda (log scale)", ylab = "Final Objective",
##D      main = "Objective Function vs. Smoothness Penalty")
## End(Not run)




cleanEx()
nameEx("project_covariate")
### * project_covariate

flush(stderr()); flush(stdout())

### Name: project_covariate
### Title: Project a Subject-Level Covariate (generic)
### Aliases: project_covariate project_covariate.covstatis

### ** Examples

# Create a list of correlation matrices
Xs <- lapply(1:5, function(i) matrix(rnorm(10*10), 10, 10))
Xs <- lapply(Xs, cor)

# Apply COVSTATIS
res <- covstatis(Xs, ncomp=3)

# Create a random covariate vector (e.g., episodic memory scores)
y <- rnorm(length(Xs))

# Project the covariate to get dimension-wise coordinates
dim_cos <- project_covariate(res, y, what = "dimension", scale = "cosine")

# Project the covariate to get an ROI-wise pattern
roi_beta <- project_covariate(res, y, what = "observation", scale = "beta")




### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')

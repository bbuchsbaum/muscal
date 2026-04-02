# Penalized MFA with Clusterwise Spatial Smoothness Constraints

This function implements a spatially-regularized Multiple Factor
Analysis for multi-subject neuroimaging or spatial cluster data. Unlike
standard MFA which encourages similarity among block loadings globally,
this method uses spatial coordinates of clusters (e.g., brain regions,
spatial locations) to construct a graph and enforces smoothness of
loadings among spatially adjacent clusters via a graph Laplacian
penalty.

## Usage

``` r
penalized_mfa_clusterwise(
  data_list,
  coords_list,
  ncomp = 2L,
  lambda = 1,
  adjacency_opts = list(),
  max_iter = 10L,
  nsteps_inner = 1L,
  learning_rate = 0.01,
  optimizer = c("gradient", "adam"),
  beta1 = 0.9,
  beta2 = 0.999,
  adam_epsilon = 1e-08,
  tol_obj = 1e-07,
  tol_inner = NULL,
  verbose = FALSE,
  preproc = NULL,
  memory_budget_mb = 1024,
  normalized_laplacian = TRUE
)

pmfa_cluster(
  data_list,
  coords_list,
  ncomp = 2L,
  lambda = 1,
  adjacency_opts = list(),
  max_iter = 10L,
  nsteps_inner = 1L,
  learning_rate = 0.01,
  optimizer = c("gradient", "adam"),
  beta1 = 0.9,
  beta2 = 0.999,
  adam_epsilon = 1e-08,
  tol_obj = 1e-07,
  tol_inner = NULL,
  verbose = FALSE,
  preproc = NULL,
  memory_budget_mb = 1024,
  normalized_laplacian = TRUE
)
```

## Arguments

- data_list:

  A list of length \\S\\ containing numeric matrices. Each element
  \\\mathbf{X}\_s\\ is an \\n_s \times k_s\\ matrix where \\n_s\\ is the
  number of observations and \\k_s\\ is the number of clusters/features
  for subject \\s\\. Data should be column-centered (mean zero per
  cluster) for proper reconstruction.

- coords_list:

  A list of length \\S\\ containing coordinate matrices. Each element is
  a \\k_s \times 3\\ matrix of (x, y, z) spatial coordinates for the
  clusters in the corresponding data block. Must have exactly 3 columns.

- ncomp:

  Integer number of latent components to extract per block (default: 2).
  When a block has fewer than \`ncomp\` clusters, the effective rank is
  automatically reduced to match the number of clusters.

- lambda:

  Non-negative scalar controlling the spatial smoothness penalty
  strength (default: 1). Larger values enforce stronger smoothness among
  spatially adjacent clusters. When \\\lambda = 0\\, the method reduces
  to independent PCA per block. Typical range: 0.1 to 10.

- adjacency_opts:

  A named list of options passed to the internal \`spatial_constraints\`
  function for adjacency matrix construction (default: empty list).
  Supported options:

  - `k_nn`: Number of nearest neighbors (default: 6)

- max_iter:

  Maximum number of outer block-coordinate descent iterations (default:
  10). More iterations allow better convergence but increase computation
  time. Typical range: 10 to 100.

- nsteps_inner:

  Number of gradient update steps per block in each outer iteration
  (default: 1). Higher values perform more thorough optimization of each
  block before moving to the next. For stability, starting with 1 is
  recommended.

- learning_rate:

  Step size for the optimizer (default: 0.01). Controls the magnitude of
  updates in the gradient direction. Too large may cause divergence; too
  small slows convergence. Typical range: 0.001 to 0.1.

- optimizer:

  Character string specifying the optimization algorithm (default:
  "gradient"). Options are:

  - `"gradient"`: Standard gradient descent with fixed learning rate

  - `"adam"`: Adaptive Moment Estimation with momentum and adaptive
    learning rates

  Adam is generally more robust but gradient descent offers more
  control.

- beta1:

  Adam hyperparameter for first moment decay (default: 0.9). Controls
  the exponential decay rate for gradient moving average. Typical range:
  0.8 to 0.95. Only used when \`optimizer = "adam"\`.

- beta2:

  Adam hyperparameter for second moment decay (default: 0.999). Controls
  the exponential decay rate for squared gradient moving average.
  Typical range: 0.99 to 0.9999. Only used when \`optimizer = "adam"\`.

- adam_epsilon:

  Small constant added to denominator in Adam for numerical stability
  (default: 1e-8). Prevents division by zero. Only used when \`optimizer
  = "adam"\`.

- tol_obj:

  Numeric tolerance for outer loop convergence based on relative change
  in objective function (default: 1e-7). The algorithm stops when
  \\\|f^{(t+1)} - f^{(t)}\| / (\|f^{(t)}\| + \epsilon) \<
  \text{tol_obj}\\. Smaller values require more precise convergence.
  Typical range: 1e-8 to 1e-5.

- tol_inner:

  Optional numeric tolerance for inner loop convergence based on
  Frobenius norm of change in loading matrix (default: NULL, no early
  stopping). When specified, inner loop stops if
  \\\\\mathbf{V}\_s^{\text{new}} - \mathbf{V}\_s\\\_F \<
  \text{tol_inner}\\.

- verbose:

  Logical indicating whether to print iteration progress (default:
  FALSE). When TRUE, displays detailed information about objective
  values, convergence, memory usage, and potential issues using the
  \`cli\` package.

- preproc:

  Preprocessing specification (default: NULL, will force centering). Can
  be:

  - `NULL`: Automatically applies centering to each block

  - A single \`pre_processor\` object: Applied independently to each
    block

  - A list of \`pre_processor\` objects: One per block

  \*\*Critical:\*\* Preprocessing must preserve the number of columns to
  maintain spatial correspondence.

- memory_budget_mb:

  Numeric specifying the maximum memory in megabytes allocated per block
  for precomputing \\\mathbf{X}\_s^\top \mathbf{X}\_s\\ (default: 1024).
  If a block's \\k_s^2 \times 8 / (1024^2)\\ exceeds this budget,
  gradients are computed on-the-fly instead. Increase for better speed
  with large memory; decrease if memory is constrained.

- normalized_laplacian:

  Logical indicating whether to use the normalized graph Laplacian
  (default: TRUE). When TRUE, uses \\\mathbf{L}\_{\text{sym}} =
  \mathbf{D}^{-1/2} \mathbf{L} \mathbf{D}^{-1/2}\\ for scale-free
  smoothness that is independent of node degree. When FALSE, uses the
  unnormalized Laplacian \\\mathbf{L} = \mathbf{D} - \mathbf{A}\\.
  Normalized is recommended for most applications.

## Value

A \`multiblock_projector\` object of class
\`"penalized_mfa_clusterwise"\` with the following components:

- `v`:

  Concatenated loading matrix (total_clusters × ncomp) formed by
  vertically stacking all block-specific loading matrices. Padded with
  zeros for variable-rank blocks.

- `preproc`:

  A pass-through preprocessor (identity transformation), since this
  function expects pre-processed input.

- `block_indices`:

  Named list indicating which rows of \`v\` correspond to each
  block/subject. Each element is an integer vector of row indices.

Additional information stored as named elements (accessible via \`\$\`):

- `V_list`:

  List of length S containing the final orthonormal loading matrices for
  each block. Each element is a \\k_s \times r_s\\ matrix where \\r_s =
  \min(\text{ncomp}, k_s)\\.

- `Sadj`:

  The (normalized or unnormalized) graph Laplacian matrix used for the
  spatial smoothness penalty. Sparse matrix of dimension \\\sum_s k_s
  \times \sum_s k_s\\.

- `LV`:

  The matrix product \\\mathbf{L} \mathbf{V}\\ from the final iteration,
  useful for computing spatial penalty contributions.

- `obj_values`:

  Numeric vector of objective function values at each outer iteration,
  including the initial value. Length is (iterations_run + 1).

- `lambda`:

  The spatial smoothness penalty weight used.

- `precompute_info`:

  Logical vector of length S indicating which blocks used precomputed
  \\\mathbf{X}\_s^\top \mathbf{X}\_s\\ (TRUE) vs. on-the-fly computation
  (FALSE).

- `iterations_run`:

  Integer indicating how many outer iterations were completed before
  convergence or reaching max_iter.

- `ncomp_block`:

  Integer vector of length S containing the effective number of
  components extracted for each block (may differ from \`ncomp\` for
  blocks with fewer clusters).

## Details

\## Optimization Problem

For \\S\\ subjects/blocks with data matrices \\\mathbf{X}\_s \in
\mathbb{R}^{n_s \times k_s}\\ (where \\k_s\\ is the number of clusters
for subject \\s\\), we estimate orthonormal loading matrices
\\\mathbf{V}\_s \in \mathbb{R}^{k_s \times r_s}\\ by minimizing:

\$\$ \min\_{\\\mathbf{V}\_s\\} \sum\_{s=1}^S \\\mathbf{X}\_s -
\mathbf{X}\_s \mathbf{V}\_s \mathbf{V}\_s^\top\\\_F^2 + \lambda
\\\text{tr}(\mathbf{V}^\top \mathbf{L} \mathbf{V}) \$\$

where:

- The first term is the reconstruction error (sum of squared residuals)

- \\\mathbf{V}\\ is the vertical concatenation of all \\\mathbf{V}\_s\\

- \\\mathbf{L} = \mathbf{D} - \mathbf{A}\\ is the graph Laplacian

- \\\mathbf{A}\\ is the adjacency matrix constructed from spatial
  coordinates

- \\\mathbf{D}\\ is the degree matrix (\\D\_{ii} = \sum_j A\_{ij}\\)

- \\\lambda \geq 0\\ controls the spatial smoothness strength

\## Spatial Smoothness Penalty

The Laplacian penalty \\\text{tr}(\mathbf{V}^\top \mathbf{L}
\mathbf{V})\\ can be rewritten as: \$\$ \frac{1}{2} \sum\_{i,j} A\_{ij}
\\\mathbf{v}\_i - \mathbf{v}\_j\\\_2^2 \$\$

This penalizes the squared Euclidean distance between loadings of
adjacent clusters. When two clusters \\i\\ and \\j\\ are connected
(\\A\_{ij} = 1\\), their loading vectors \\\mathbf{v}\_i\\ and
\\\mathbf{v}\_j\\ are encouraged to be similar.

\*\*Normalized vs. Unnormalized Laplacian:\*\*

- **Unnormalized** (`normalized_laplacian = FALSE`):

  \\\mathbf{L} = \mathbf{D} - \mathbf{A}\\. Smoothness is scaled by node
  degree; high-degree nodes have stronger smoothness constraints.

- **Normalized** (`normalized_laplacian = TRUE`, default):

  \\\mathbf{L}\_{\text{sym}} = \mathbf{D}^{-1/2} \mathbf{L}
  \mathbf{D}^{-1/2}\\. Scale-free smoothness that is independent of node
  degree. Recommended for most applications.

\## Graph Construction

The spatial adjacency graph is built using k-nearest neighbors (k-NN) in
3D space:

1.  Cluster coordinates from all subjects are pooled

2.  For each cluster, find its k nearest neighbors (default: k=6)

3.  Create edges between each cluster and its neighbors

4.  Symmetrize the adjacency matrix

The number of neighbors can be controlled via \`adjacency_opts =
list(k_nn = k)\`.

\## Optimization Algorithm

The method uses block-coordinate descent (BCD) with Riemannian
optimization:

1\. \*\*Initialization\*\*: Initialize each \\\mathbf{V}\_s\\ via SVD of
\\\mathbf{X}\_s\\ 2. \*\*Outer loop\*\* (max_iter iterations): Cycle
through all subjects/blocks 3. \*\*Inner loop\*\* (nsteps_inner steps
per block): - Compute reconstruction gradient: \\\nabla\_{\mathbf{V}\_s}
= -2 \mathbf{X}\_s^\top \mathbf{X}\_s \mathbf{V}\_s\\ - Compute spatial
gradient: \\2\lambda \mathbf{L}\_{ss} \mathbf{V}\_s + 2\lambda
\mathbf{L}\_{s,-s} \mathbf{V}\_{-s}\\ - Project combined gradient onto
Stiefel manifold tangent space - Update via gradient descent or Adam -
Retract to manifold via QR decomposition 4. \*\*Convergence check\*\*:
Stop when relative change in objective \< \`tol_obj\`

\## Memory Management

For large datasets, computing and storing \\\mathbf{X}\_s^\top
\mathbf{X}\_s\\ for all blocks may exceed available memory. The
\`memory_budget_mb\` parameter controls this tradeoff:

- \*\*Precomputed mode\*\* (if \\k_s^2 \times 8 / (1024^2)\\ \<
  \`memory_budget_mb\`): Store \\\mathbf{X}\_s^\top \mathbf{X}\_s\\ for
  faster gradient computation

- \*\*On-the-fly mode\*\* (otherwise): Compute gradients as needed,
  saving memory at the cost of computation time

\## Variable-Rank Loadings

When subjects have different numbers of clusters, each block can have a
different effective rank. The algorithm automatically sets \\r_s =
\min(\text{ncomp}, k_s)\\ for each block and pads the concatenated
loading matrix with zeros for consistency.

\## Preprocessing

Data preprocessing is crucial for reconstruction-based methods. The
function defaults to centering each block if no preprocessing is
specified. \*\*Important:\*\* Preprocessing must preserve the number of
columns (clusters) to maintain correspondence with spatial coordinates.

\## When Lambda = 0

When \\\lambda = 0\\, the spatial penalty is disabled and the method
reduces to independent PCA on each block. This provides a useful
baseline for comparison.

\## Practical Considerations

- **Choosing lambda**:

  Start with \\\lambda = 0\\ (no smoothness) and gradually increase.
  Typical range: 0.1 to 10. Monitor the objective function components
  (reconstruction vs. smoothness) to balance the tradeoff.

- **Number of neighbors**:

  More neighbors (larger k-NN) create a denser graph with stronger
  smoothness constraints. Default k=6 works well for 3D spatial data.

- **Convergence**:

  The method may converge slowly for large \\\lambda\\. Increase
  \`max_iter\` or adjust \`learning_rate\` if needed.

- **Optimizer choice**:

  Adam is generally more robust and converges faster. Use gradient
  descent for better control or when Adam's adaptive behavior is
  undesired.

## Engineering Improvements

This implementation includes several optimizations:

- \*\*Mathematical\*\*: Normalized Laplacian default, \\\lambda=0\\
  optimization skips

- \*\*Numerical\*\*: Sparse matrix operations, in-block QR retraction,
  fast orthogonalization

- \*\*Stability\*\*: Variable-rank loadings, robust parameter
  validation, Adam state management

- \*\*Architecture\*\*: Modular helper functions, comprehensive logging,
  clear error messages

- \*\*Memory\*\*: Adaptive memory budget, efficient gradient computation
  strategies

- \*\*API\*\*: Shorter \`pmfa_cluster()\` alias, consistent naming
  conventions

## Examples

``` r
if (FALSE) { # \dontrun{
# Example 1: Basic usage with simulated spatial cluster data
set.seed(123)
S <- 3  # 3 subjects

# Generate data with spatial structure
data_list <- lapply(1:S, function(s) {
  n <- 50  # 50 observations
  k <- 20  # 20 clusters
  matrix(rnorm(n * k), n, k)
})

# Generate 3D spatial coordinates for clusters
coords_list <- lapply(1:S, function(s) {
  matrix(runif(20 * 3, 0, 10), 20, 3)  # Random 3D positions
})

# Fit model with spatial smoothness
res <- penalized_mfa_clusterwise(
  data_list, coords_list,
  ncomp = 3,
  lambda = 1,
  max_iter = 20,
  verbose = TRUE
)
print(res)

# Plot convergence
plot(res$obj_values, type = 'b', xlab = 'Iteration', ylab = 'Objective',
     main = 'Convergence of Spatially-Regularized MFA')

# Example 2: Using the shorter alias
res2 <- pmfa_cluster(data_list, coords_list, ncomp = 2, lambda = 0.5)

# Example 3: Compare different lambda values
lambdas <- c(0, 0.1, 0.5, 1, 2, 5)
results <- lapply(lambdas, function(lam) {
  fit <- pmfa_cluster(data_list, coords_list, ncomp = 2, lambda = lam,
                      verbose = FALSE)
  list(
    lambda = lam,
    final_obj = tail(fit$obj_values, 1),
    iterations = fit$iterations_run
  )
})

# Example 4: Using Adam optimizer for faster convergence
res_adam <- pmfa_cluster(
  data_list, coords_list,
  ncomp = 3,
  lambda = 1,
  optimizer = "adam",
  learning_rate = 0.05,
  max_iter = 30,
  verbose = TRUE
)

# Example 5: Controlling k-NN graph construction
res_dense <- pmfa_cluster(
  data_list, coords_list,
  ncomp = 2,
  lambda = 1,
  adjacency_opts = list(k_nn = 10),  # More neighbors = denser graph
  verbose = TRUE
)

# Example 6: Using unnormalized Laplacian
res_unnorm <- pmfa_cluster(
  data_list, coords_list,
  ncomp = 2,
  lambda = 1,
  normalized_laplacian = FALSE
)

# Example 7: Memory-constrained settings
# For large datasets, reduce memory budget
res_mem <- pmfa_cluster(
  data_list, coords_list,
  ncomp = 2,
  lambda = 1,
  memory_budget_mb = 256,  # Limit memory per block
  verbose = TRUE
)
# Check which blocks used precomputed gradients
print(res_mem$precompute_info)

# Example 8: Extracting block-specific loadings
V_list <- res$V_list

# Loadings for first subject
V1 <- V_list[[1]]
dim(V1)  # k_s x ncomp

# Visualize spatial smoothness of first loading vector
library(scatterplot3d)
scatterplot3d(
  coords_list[[1]],
  color = rank(V1[, 1]),
  main = "Spatial Pattern of First Loading (Subject 1)",
  xlab = "X", ylab = "Y", zlab = "Z"
)

# Example 9: Examining the spatial penalty contribution
# Extract Laplacian and loadings
L <- res$Sadj
LV <- res$LV
v <- res$v

# Compute spatial penalty: tr(V' L V)
spatial_penalty <- sum(LV * v)
cat("Spatial penalty term:", spatial_penalty, "\n")

# Example 10: Variable-rank loadings (blocks with different sizes)
# Simulate data with varying cluster numbers
data_var <- list(
  matrix(rnorm(50 * 15), 50, 15),  # 15 clusters
  matrix(rnorm(50 * 20), 50, 20),  # 20 clusters
  matrix(rnorm(50 * 10), 50, 10)   # 10 clusters
)

coords_var <- list(
  matrix(runif(15 * 3), 15, 3),
  matrix(runif(20 * 3), 20, 3),
  matrix(runif(10 * 3), 10, 3)
)

res_var <- pmfa_cluster(
  data_var, coords_var,
  ncomp = 8,  # Will be capped at 10 (smallest block)
  lambda = 1,
  verbose = TRUE
)
# Check effective ranks per block
print(res_var$ncomp_block)

# Example 11: Custom preprocessing per block
library(multivarious)

preproc_list <- list(
  center(),        # Just center first block
  standardize(),   # Center and scale second block
  center()         # Just center third block
)

res_preproc <- pmfa_cluster(
  data_list, coords_list,
  ncomp = 2,
  lambda = 1,
  preproc = preproc_list
)

# Example 12: Lambda selection via cross-validation style approach
# (Simplified - would need proper CV in practice)
lambda_grid <- 10^seq(-1, 1, length.out = 10)
cv_results <- data.frame(
  lambda = lambda_grid,
  recon_error = NA,
  spatial_penalty = NA,
  total_obj = NA
)

for (i in seq_along(lambda_grid)) {
  fit <- pmfa_cluster(data_list, coords_list, ncomp = 2,
                      lambda = lambda_grid[i], verbose = FALSE)
  cv_results$total_obj[i] <- tail(fit$obj_values, 1)
  # Could decompose into reconstruction and spatial components
}

# Plot objective vs lambda
plot(cv_results$lambda, cv_results$total_obj, log = "x",
     type = "b", xlab = "Lambda (log scale)", ylab = "Final Objective",
     main = "Objective Function vs. Smoothness Penalty")
} # }
```

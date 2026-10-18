# `penalized_mfa_clusterwise()` Recovery Plan

## Status and decision

This is a repair plan for the existing `penalized_mfa_clusterwise()` /
`pmfa_cluster()` API. The current `lambda = 0` path has a valid independent-PCA
subspace interpretation, but the spatial path is not ready for scientific use.

Until Milestone 5 is complete:

- document `lambda > 0` as experimental;
- do not use this method as a backend for `design_mfa()`;
- preserve the current public function and alias, but avoid adding new wrappers
  around the existing optimizer;
- treat every confirmed defect below as a required regression test.

The recovery should be a contract-first rewrite of the graph and optimization
internals, not a sequence of learning-rate adjustments.

## Scope

### In scope

- Define the precise penalized objective and its weighting.
- Replace dense, order-dependent graph construction with sparse graph
  primitives.
- Correct the combinatorial and normalized Laplacians.
- Replace the variable-rank fiction with a coherent common-rank contract.
- Replace the current fixed-step/Adam update with a convergent product-Stiefel
  optimizer.
- Repair preprocessing, projection, diagnostics, and fit metadata.
- Establish correctness, robustness, and performance gates.
- Preserve the public S3 class and the `pmfa_cluster()` alias where possible.

### Out of scope

- Factorial design decomposition or group-level design inference. That belongs
  in `design_mfa()`.
- Dense all-pairs coordinate distances, even as an "exact" large-data mode.
- C++ optimization before the sparse R implementation is correct and profiled.
- Variable numbers of identifiable components across blocks. A future masked
  or hierarchical component model would need a different objective.

## Mathematical contract

For preprocessed block `s`, let

```math
X_s \in \mathbb{R}^{n_s \times p_s}, \qquad
V_s \in \mathbb{R}^{p_s \times K}, \qquad V_s^\top V_s = I_K.
```

Let `V = rbind(V_1, ..., V_S)` and let `L` be a symmetric positive
semidefinite Laplacian over all block features. The repaired method minimizes

```math
F(V_1,\ldots,V_S) =
\sum_{s=1}^S \alpha_s
\left(\lVert X_s\rVert_F^2 - \lVert X_sV_s\rVert_F^2\right)
+ \lambda\operatorname{tr}(V^\top L V)
```

subject to the blockwise Stiefel constraints.

Consequences of this contract:

- `K <= min(p_s)` is required. All blocks share the same component columns.
- After preprocessing, `K` must not exceed the common numerically supported
  data rank unless an explicitly documented spatial-only component model is
  introduced later.
- `lambda = 0` is exactly independent rank-`K` PCA for every block.
- Component columns are coupled by the spatial penalty, but only the joint
  subspace is identifiable when roots are tied.
- `L` must be symmetric positive semidefinite and use the same concatenated
  feature order as `V`.
- The optimizer must report both objective convergence and first-order
  stationarity on the product Stiefel manifold.

### Block weighting

Add an explicit `block_weight` argument. The objective always records the
resolved numeric `alpha_s` values.

- `NULL` retains `alpha_s = 1` for one compatibility release.
- Numeric weights provide an exact user-controlled objective.
- A documented equal-subject option uses `alpha_s = 1 / n_s`.
- An MFA-style option uses the inverse squared leading singular value after
  preprocessing.

When trial counts differ and compatibility weights are used, print a concise
warning that the reconstruction term is row-count weighted. The new
`design_mfa()` method must not inherit this compatibility default.

## Confirmed defects to encode as regressions

1. Symmetrizing a directed k-nearest-neighbor matrix with
   `forceSymmetric(..., uplo = "U")` discards lower-triangle-only edges.
2. Reordering features changes the resulting graph after mapping it back to
   the original feature identities.
3. The current normalized-Laplacian path is identical to the unnormalized path.
4. `ncomp > p_s` for one block produces incompatible padded and unpadded
   gradients.
5. `k_nn >= P` can produce `NA` sparse-matrix indices.
6. Penalized objective traces can increase, and some fits end worse than their
   initialization.
7. Relative objective change can declare apparent convergence while the
   projected gradient remains large.
8. The returned projector does not contain the fitted centering preprocessors.
9. Non-finite data and coordinates are not rejected at the API boundary.
10. Dense distance construction is quadratic in the pooled feature count.

## Target public API

Keep the existing arguments where they remain meaningful, and evolve toward:

```r
penalized_mfa_clusterwise(
  data_list,
  coords_list,
  ncomp = 2L,
  rank_action = c("error", "reduce"),
  lambda = 1,
  block_weight = NULL,
  graph = NULL,
  adjacency_opts = list(
    mode = "within_plus_exact",
    within_k = 6L,
    between_k = 1L,
    symmetrize = "union",
    normalized = TRUE
  ),
  optimizer = "riemannian",
  max_iter = 100L,
  tol_obj = 1e-8,
  tol_grad = 1e-5,
  step_control = list(method = "armijo"),
  preproc = multivarious::center(),
  memory_budget_mb = 1024,
  verbose = FALSE
)
```

Compatibility policy:

- Accept `normalized_laplacian` temporarily and translate it into graph
  construction options.
- Accept `optimizer = "gradient"` as a deprecated alias for
  `"riemannian"`.
- Deprecate Adam for this constrained quadratic. It may remain as an explicitly
  non-monotone experimental option for one release, but it is not a release
  gate and must not be the default.
- Reject heterogeneous block ranks rather than silently padding them.

## Milestone 0: Quarantine and executable baselines

### Work

- Add a lifecycle warning to documentation for `lambda > 0`.
- Save deterministic minimal reproductions for every confirmed defect.
- Record the current `lambda = 0` PCA subspace result as a mathematical oracle,
  not as a golden serialization.
- Add a development benchmark script with fixed seeds and session information.
- Split the current monolithic source into graph, optimizer, and public-wrapper
  internals only when doing so does not obscure `git blame` for the regressions.

### Acceptance gate

- Every defect above has a failing test or benchmark demonstrating the old
  behavior.
- The existing API remains callable.
- The `lambda = 0` projector-distance oracle is below a tolerance based on
  matrix conditioning, normally `1e-10` for the canonical fixtures.

## Milestone 1: Input and rank contract

### Work

- Validate that every data block is numeric, two-dimensional, nonempty, and
  finite.
- Validate finite numeric coordinates with a common coordinate dimension.
- Require `nrow(coords_list[[s]]) == ncol(data_list[[s]])`.
- Require unique block names and preserve them through all returned objects.
- Require a common rank satisfying `1 <= K <= min(p_s)`.
- Determine the effective rank of every preprocessed block using the same
  scale-aware threshold used by initialization. If `K` exceeds the minimum
  effective rank, either fail or reduce the common rank according to
  `rank_action`; never fill unsupported columns randomly.
- Validate nonnegative finite block weights with at least one positive entry.
- Define behavior for constant columns, rank-deficient blocks, duplicate
  coordinates, disconnected graphs, and zero-degree vertices.
- Fail before preprocessing or SVD if the contract is violated.

### Acceptance gate

- Invalid inputs fail with stable, specific messages.
- Minimal valid blocks (`p_s = K`, small `n_s`) fit without special-case crashes.
- Rank-deficient data follow the documented error-or-reduce policy and report
  each block's effective rank; they do not create random hidden components.

## Milestone 2: Sparse graph rewrite

### Graph representation

Introduce an internal graph object containing:

- symmetric sparse adjacency `A`;
- combinatorial and optional normalized Laplacians;
- degree vector computed from `A`, never from `L`;
- block offsets and feature identities;
- graph-construction metadata and edge counts by edge type.

The normalized Laplacian is

```math
L_{sym} = D^{-1/2}(D-A)D^{-1/2},
```

with zero rows for isolated vertices.

### Construction modes

- Build within-block spatial neighbors separately using `RANN::nn2()` or an
  equivalent sparse spatial index.
- Represent exact cross-subject coordinate matches as a separate edge class.
- Support sparse user-provided correspondence edges for softly aligned grids.
- If approximate cross-block kNN is requested, query block pairs or a named
  reference grid rather than silently mixing it with within-block smoothing.
- Symmetrize directed neighbor edges with an explicit union, mutual, or average
  rule. Never use one triangle as the source of truth.
- Cap neighbor counts safely and define a deterministic policy for distance
  ties and duplicate coordinates.
- Preserve a custom sparse `graph` input after validating symmetry,
  nonnegative edges, dimensions, zero diagonal, and finite values.

### Required tests

- Hand-calculated line, square, disconnected, and duplicate-coordinate graphs.
- Feature-permutation equivariance after restoring original feature identities.
- Block-order equivariance when names and all associated inputs are permuted.
- Union and mutual symmetrization edge oracles.
- Exact-match and custom-edge equivalence.
- Combinatorial and normalized Laplacian analytic oracles.
- Symmetry and positive-semidefinite quadratic-form properties.
- `k = 0`, `k = 1`, `k >= p_s`, isolated nodes, and one-feature blocks.
- A guard proving that no dense `P x P` distance object is allocated.

### Acceptance gate

- All graph tests pass under repeated feature and block permutations.
- `max(abs(L - t(L)))` is within sparse floating-point tolerance.
- Random vectors satisfy `x' L x >= -tol`.
- Graph memory grows with stored edges, not with `P^2`.

## Milestone 3: Product-Stiefel optimizer

### Algorithm

Implement block-coordinate Riemannian gradient descent with Armijo
backtracking:

1. Compute the exact Euclidean gradient using the current values of every
   block, not a stale outer-iteration copy.
2. Project each block gradient onto the tangent space:

   ```math
   G_s^R = G_s - V_s\operatorname{sym}(V_s^\top G_s).
   ```

3. Propose a step and retract with a thin QR or polar retraction.
4. Evaluate the complete objective and backtrack until the Armijo condition
   holds.
5. Cycle blocks until both objective and projected-gradient criteria are met.

Expose the reconstruction, spatial, and total objectives independently. Store
per-iteration step sizes, backtracking counts, objective changes, and projected
gradient norms.

### Numerical implementation

- Use `crossprod(X_s, X_s %*% V_s)` when a Gram matrix exceeds the memory
  budget.
- Use a Gram matrix only when the measured cost model favors it.
- Use `irlba` or a small-Gram eigenpath for initialization when truncated SVD
  is materially cheaper; retain exact `svd()` for small matrices.
- Apply a deterministic sign convention only after convergence.
- Compare tied components through projectors or principal angles, not raw
  columns.

### Required tests

- Central finite-difference checks of the Euclidean and projected gradients.
- Differential comparison against exhaustive or trusted manifold solutions on
  tiny `K = 1` and `K = 2` problems.
- Monotone total objective under the default optimizer.
- Final orthonormality and projected-gradient residual checks.
- `lambda = 0` equivalence to independent PCA.
- Continuity as `lambda -> 0` and increasing smoothness as `lambda` grows.
- Equivalent results from precomputed-Gram and on-the-fly paths.
- Randomized fuzzing over dimensions, condition numbers, weights, graphs, and
  penalty strengths.
- Multi-start subspace agreement where the optimum is separated; honest
  non-uniqueness diagnostics otherwise.

### Acceptance gate

- The accepted objective trace is nonincreasing up to a scale-aware tolerance.
- The final projected-gradient norm is below `tol_grad`, or the fit returns
  `converged = FALSE` with a specific termination reason.
- No canonical penalized fixture ends worse than its initialization.
- The two Gram strategies agree in objective and fitted projectors.

## Milestone 4: Preprocessing and fitted-object contract

### Work

- Store the fitted block preprocessors returned by
  `prepare_block_preprocessors()`.
- Build the concatenated projector preprocessor with
  `multivarious::concat_pre_processors()`.
- Store authoritative `V_list`, block indices, graph metadata, resolved block
  weights, rank, objective components, convergence diagnostics, and the fitted
  call.
- Attach a `muscal_fit_contract` and a complete refit specification.
- Define `project()` for one named block and for a complete block list.
- Ensure projection transforms raw new data exactly once.
- Define whether prediction means scores, reconstruction, or both; do not let
  inherited behavior decide this accidentally.

### Required tests

- Raw training projections equal explicit
  `transform(preprocessor, X) %*% V_s` calculations.
- Centering, scaling, custom preprocessor lists, and `preproc = NULL` follow the
  documented semantics.
- Named block reordering does not alter projections.
- Refit metadata reproduces the original deterministic fit.
- Serialization and reload preserve projection behavior.

### Acceptance gate

- Projection discrepancies are below numerical tolerance for every supported
  preprocessing path.
- `infer_muscal()` can resample/refit the object without private reconstruction
  of arguments.

## Milestone 5: Statistical behavior and performance

### Statistical simulations

- Independent blocks with no spatial correspondence.
- Shared smooth loading fields with known cross-block correspondence.
- Misregistered coordinates and deliberately wrong graphs.
- Unequal trial counts, unequal noise levels, and explicit block weights.
- Tied eigenvalues, weak signal, rank deficiency, and disconnected graphs.

Evaluate subspace recovery, reconstruction error, graph roughness, sensitivity
to block replication, and failure signaling. Do not claim group inference from
this method; these simulations characterize the estimator only.

### Performance work

- Benchmark graph construction independently from optimization.
- Record elapsed time, peak resident memory, edge count, iterations, and sparse
  matrix sizes.
- Benchmark exact versus truncated SVD and Gram versus on-the-fly gradients.
- Add chunked/operator paths before considering RcppEigen.

Provisional reference-runner gates, to be calibrated and then frozen:

- `P = 100,000`, `k = 6` graph construction uses less than 1.5 GB peak memory
  and completes within 60 seconds.
- A representative `4 x (80 trials x 5,000 features)`, `K = 5` fit uses less
  than 2 GB and completes within 60 seconds.
- Scheduled million-node graph tests demonstrate sparse growth and never
  allocate a dense distance matrix.
- After baselines are established, PR performance fails only for a robust
  regression such as more than 2x the rolling reference runtime; nightly jobs
  retain the stricter size envelope.

### Acceptance gate

- Canonical recovery simulations meet predeclared subspace and roughness
  thresholds.
- Performance budgets are recorded with hardware and package versions.
- `R CMD check` and the full package test suite pass.
- Documentation no longer labels the repaired implementation experimental.

## Test and benchmark layout

Recommended files:

- `tests/testthat/test-cluster-graph.R`
- `tests/testthat/test-penalized-mfa-clusterwise-contract.R`
- `tests/testthat/test-penalized-mfa-clusterwise-oracles.R`
- `tests/testthat/test-penalized-mfa-clusterwise-optimizer.R`
- `tests/testthat/test-penalized-mfa-clusterwise-projection.R`
- `tests/testthat/helper-clusterwise-sim.R`
- `dev/benchmark_penalized_mfa_clusterwise.R`

CI pyramid:

- Every PR: contracts, regressions, analytic oracles, finite differences, and
  small metamorphic tests.
- Nightly: randomized fuzzing, large sparse graphs, optimizer stress tests, and
  performance envelopes.
- Pre-release: full simulation matrix, serialized projection checks, reverse
  dependency checks if applicable, and rendered-vignette review.

## Completion definition

The recovery is complete only when:

1. the objective and graph semantics are documented precisely;
2. every confirmed defect has a passing regression;
3. the default optimizer is monotone and meets a stationarity criterion;
4. raw-data projection uses the actual fitted preprocessing pipeline;
5. graph construction is sparse and permutation equivariant;
6. the common-rank contract is enforced;
7. statistical simulations and performance gates pass; and
8. no part of `design_mfa()` needs to rely on legacy clusterwise behavior.

# Codebase Concerns

**Analysis Date:** 2026-01-22

## Tech Debt

**Missing Package Dependencies and Private API Access:**
- Issue: The package relies on triple-colon operator (`:::`) to access private/internal functions from the `multivarious` package
- Files: `R/bada.R` (lines 161, 237, 239)
- Impact: These private APIs (`multivarious:::fresh`, `multivarious:::concat_pre_processors`, `multivarious:::discriminant_projector`) could break without warning if the parent package's internal structure changes. No guarantee of API stability
- Fix approach: Coordinate with `multivarious` maintainers to export these functions officially, or implement equivalent functionality directly in muscal

**Attribute-Based Metadata Storage (Anti-pattern):**
- Issue: Both `penalized_mfa` and `penalized_mfa_clusterwise` store auxiliary results as object attributes instead of list elements within the multiblock_projector
- Files: `R/penalized_mfa.R` (lines 466-472 accessing via `attr()`) and `R/penalized_mfa_clusterwise.R` (similar pattern)
- Impact: Attributes are fragile - they don't serialize reliably, may be stripped during subsetting operations, and make code harder to understand. Inconsistent with the design goal mentioned in `/dev/todo.md`
- Fix approach: Refactor return value construction to use `...` mechanism for storing auxiliary results (e.g., `obj_values`, `lambda`, `consensus`, `V_list`) as named list elements within the `multiblock_projector` instead of `attr()`. Update all accessor code to use `$` instead of `attr()`

**Inconsistent S3 Method Signatures:**
- Issue: Multiple S3 generic functions have argument mismatches between their generic definition and method implementations
- Files: `R/all_generic.R`, `R/penalized_mfa.R`, `R/mfa.R`, `R/covstatis.R`, `R/bada.R`
- Examples documented in `check_todo.md` (lines 26-33):
  - `project_cov` vs `project_cov.covstatis` argument mismatch
  - `mfa` vs `mfa.list` vs `mfa.multiblock` argument mismatch
  - `covstatis` vs `covstatis.list` argument mismatch
  - `reconstruct` vs `reconstruct.covstatis` argument mismatch
  - `bootstrap` vs `bootstrap.bada` argument mismatch
  - `reprocess` vs `reprocess.bada` argument mismatch
  - `project` vs `project.bada` argument mismatch
- Impact: R CMD check will complain and cause dispatch confusion. Users can pass arguments that silently get ignored
- Fix approach: Ensure all method signatures match their generic counterparts or use `...` to capture extra arguments explicitly

**Unused Imports in DESCRIPTION:**
- Issue: `deflist` and possibly `irlba` are listed in DESCRIPTION Imports but may not be actively used
- Files: `/DESCRIPTION` (line 40)
- Impact: Unnecessary dependencies increase installation time and surface area for bugs
- Fix approach: Audit actual usage and move to Suggests or remove entirely

**Heavy Dependency on External Packages:**
- Issue: Package imports 22 dependencies including `MatrixCorrelation`, `RANN`, `genpca`, `pracma`, `furrr` - some with conditional usage
- Files: `/DESCRIPTION` (lines 14-41)
- Impact: Large dependency footprint creates maintenance burden and potential version conflicts. Some dependencies (MatrixCorrelation, RANN, genpca) are used conditionally or in specific functions
- Fix approach: Consider moving conditionally-used packages to Suggests instead of Imports; review pracma usage for potential reimplementation using base R

---

## Known Bugs & Missing Functionality

**Examples Failing in MFA Function:**
- Issue: `mfa()` examples fail because they depend on `genpca` package which must be installed
- Files: `R/mfa.R` (line 177 checks for genpca)
- Status: ERROR per `check_todo.md` (line 4)
- Impact: R CMD check will fail; package cannot be released without this fix
- Trigger: Running `mfa()` examples without genpca installed
- Workaround: Install genpca separately; for release, either make genpca a hard dependency or wrap examples in `\dontrun{}`

**Undocumented Functions and Missing S3 Generics:**
- Issue: `significant_components` function has no S3 generic and minimal documentation
- Files: `R/components_utils.R`, needs entry in `R/all_generic.R`
- Status: Documented in `check_todo.md` (lines 15, 38)
- Impact: Function is exported but lacks proper documentation in man/ folder; S3 methods for this function cannot be created
- Fix approach: Create `significant_components.R` generic wrapper and document properly in roxygen

**Unresolved Cross-Reference Links:**
- Issue: Documentation references functions from multivarious package (e.g., `multidesign::split_by`, `multivarious::between_class_scatter`) that don't exist or aren't exported
- Files: Various R files with roxygen comments
- Status: Multiple unresolved references per `check_todo.md` (lines 15-19)
- Impact: HTML documentation will have broken links; violates R documentation standards
- Fix approach: Verify actual function names in dependency packages and fix references or adjust documentation approach

---

## Performance Bottlenecks

**Inefficient Nested Loop in spatial_constraints Function:**
- Issue: k-NN graph construction uses inefficient distance-based search with full distance matrix computation
- Files: `R/penalized_mfa_clusterwise.R` (lines 42-68)
- Problem: `dist()` computes all pairwise distances (O(n²) space), then loops to find k-NN (O(n² log n) time). For 1000+ clusters, this becomes slow
- Current solution: Comment suggests RANN should be used but falls back to dist()
- Cause: Dependency on RANN is conditional; `dist()` is a fallback
- Improvement path: Always use RANN for k-NN (faster than O(n²)), or implement kd-tree approach. Add early exit for small graphs

**Memory-Critical: Large Cross-Product Matrices:**
- Issue: Computing and storing X^T X matrices for large blocks can exceed available memory
- Files: `R/penalized_mfa_clusterwise.R` (lines 314-325, `should_precompute()` function)
- Current mitigation: `memory_budget_mb` parameter controls precomputation vs on-the-fly computation
- Impact: On-the-fly computation is 2-3x slower due to repeated matrix multiplications (lines 87 in same file)
- Cause: Large neuroimaging datasets (e.g., k_s > 1000) generate massive XtX matrices
- Improvement path: Implement block-wise gradient computation; consider sparse representations for spatial structure

**Redundant Laplacian Computation When Lambda = 0:**
- Issue: Even when lambda=0 (no spatial penalty), sparse Laplacian matrix is constructed and stored
- Files: `R/penalized_mfa_clusterwise.R` (lines 882-887 added optimization but initial computation at line 755+ still happens)
- Impact: Wastes memory and computation for users exploring without spatial constraints
- Cause: Optimization partially added but Laplacian construction happens earlier
- Improvement path: Move all Laplacian operations behind `if (lambda > 0)` guard; defer construction

**SVD in Loop During Initialization:**
- Issue: SVD computed for each block independently during initialization without any caching
- Files: `R/penalized_mfa_clusterwise.R` (lines 815-842)
- Impact: If initialization is called multiple times (e.g., in cross-validation), SVD is recomputed from scratch
- Cause: No memoization or SVD caching mechanism
- Improvement path: Cache SVD results if initialization called multiple times; accept optional pre-computed SVD as parameter

**Adam Optimizer Storage Overhead:**
- Issue: When using Adam optimizer, maintains two additional moment matrices (M and V2) per block per iteration
- Files: `R/penalized_mfa_clusterwise.R` (lines 845-850, 962-968)
- Impact: Storage scales as O(S * k_max * ncomp * num_moments) = roughly 2x memory overhead
- Cause: Adam maintains full moment matrices instead of scalar learning rates
- Improvement path: Consider gradient accumulation buffer instead; or reduce precision to float32 for moments

---

## Fragile Areas

**penalized_mfa_clusterwise Function (1191 lines):**
- Files: `R/penalized_mfa_clusterwise.R`
- Why fragile: Monolithic function with complex control flow, many branching conditions, multiple optimization strategies (gradient/adam), and special cases for empty blocks, variable-rank loadings, and lambda=0
- Safe modification: Add comprehensive test cases for edge cases (empty blocks, single-component, lambda=0, very small memory budget). Test SVD failure recovery
- Test coverage gaps:
  - Empty blocks (k_s = 0) recovery
  - Single component extraction (ncomp=1)
  - Numerical instability when Laplacian becomes singular
  - Memory budget transitions between precomputed and on-the-fly modes
  - Adam convergence vs gradient descent convergence

**Objective Function Calculation in penalized_mfa_clusterwise:**
- Files: `R/penalized_mfa_clusterwise.R` (lines 890-922)
- Why fragile: Handles both precomputed and on-the-fly gradient modes, skips empty blocks, contains checks for numerical validity
- Fragility: Line 915 silently clamps negative reconstruction costs (`max(0, recon_term)`) due to numerical precision - masks real issues
- Risk: Silent clamping could hide convergence problems or numerical instability
- Safe modification: Add assertions that reconstruction should be positive-semidefinite; log warnings when clamping occurs; investigate root cause

**Block Update Function with Riemannian Manifold Operations:**
- Files: `R/penalized_mfa_clusterwise.R` (lines 127-227, `block_update_cluster` function)
- Why fragile: QR decomposition for retraction can fail; multiple projection steps onto Stiefel manifold with potential numerical issues
- Risk: Singular or near-singular matrices in QR decomposition could cause NaN propagation
- Safe modification: Add error handling around QR decomposition; use SVD-based retraction as fallback; add assertions about orthonormality before/after retraction

**Preprocessing Pipeline Consistency:**
- Files: `R/penalized_mfa.R` (lines 43-76), `R/penalized_mfa_clusterwise.R` (preprocessing section)
- Why fragile: Preprocessing must preserve column count for spatial coordinates to match. If preprocessing drops or rearranges columns, entire spatial penalty becomes meaningless
- Risk: Silent dimension mismatch between data and coordinates
- Current safeguard: Line 73 in penalized_mfa.R has check; penalized_mfa_clusterwise relies on user to ensure preprocessor preserves columns
- Safe modification: Add explicit validation that preprocessing doesn't change column count; raise error if mismatch detected

**Variable-Rank Loadings Logic:**
- Files: `R/penalized_mfa_clusterwise.R` (lines 327-331 in docs, 852-876 in implementation)
- Why fragile: When blocks have different cluster counts, algorithm must handle variable-rank matrices with zero-padding. Combines inconsistent dimensions which complicates indexing
- Risk: Off-by-one errors when mapping between padded and unpadded representations
- Safe modification: Add explicit tests for variable-rank scenarios (blocks with 5, 10, 20 clusters with ncomp=8). Verify block_indices correctly map to padded matrix positions

---

## Scaling Limits

**Spatial Coordinates Memory:**
- Current capacity: All coordinates must fit in memory as dense matrices
- Limit: ~100,000 clusters (3D coordinates per cluster = 2.4 MB at 8 bytes/double)
- Scaling path: For brain imaging, typical datasets are 5,000-50,000 voxels per subject - acceptable. For larger spatial domains, implement out-of-core k-NN computation

**Laplacian Matrix Sparsity:**
- Current capacity: Sparse matrices stored as Matrix package dgCMatrix
- Limit: ~1,000,000 non-zero entries before manipulation becomes slow
- Scaling path: For k=6 nearest neighbors, expected non-zeros ≈ 12*num_clusters. At 100k clusters = 1.2M non-zeros - at limit. Consider hierarchical approaches or approximate Laplacians

**Cross-Product Matrix (XtX):**
- Current capacity: Dense precomputed XtX matrices limited by memory_budget_mb (default 1024 MB per block)
- Limit: k_s ≤ ~11,000 columns (11000² * 8 bytes ≈ 968 MB)
- Scaling path: Implement incremental QR for on-the-fly XtX products; use sketching techniques; move to iterative methods that don't need explicit XtX

**Optimization Iterations with Verbose Output:**
- Current capacity: CLI logging scales linearly with iterations
- Limit: >10,000 iterations produces massive output. At 10 subjects × max_iter=100, ~1000 log lines
- Scaling path: Implement summary logging instead of per-iteration; use progress bars instead of detailed output

---

## Dependencies at Risk

**multivarious Package (Critical Private API Dependency):**
- Risk: Reliance on `:::` private functions with no stability guarantee
- Impact: If multivarious refactors internals, bada.R breaks completely
- Current usage:
  - `multivarious:::fresh` - used to copy preprocessor state (line 161, 49)
  - `multivarious:::concat_pre_processors` - used to combine preprocessor lists (line 237)
  - `multivarious:::discriminant_projector` - used to construct final result object (line 239)
- Migration plan: Contact multivarious maintainers to export these functions officially, or implement equivalent functionality in muscal

**genpca Package (Conditional, but blocks MFA functionality):**
- Risk: genpca is required for `mfa()` examples but not in Imports (only implicitly)
- Impact: R CMD check fails because examples reference undefined objects
- Current status: Checked at runtime with `requireNamespace()` (line 177 in mfa.R)
- Migration plan: Either (1) add genpca as hard dependency, (2) wrap examples in `\dontrun{}`, or (3) implement MFA normalization without genpca

**RANN Package (Conditional, affects performance):**
- Risk: Package falls back to dist() when RANN unavailable, causing O(n²) performance penalty
- Impact: Users without RANN installed get poor performance on large spatial datasets
- Current status: Checked at runtime in `spatial_constraints` (line 46) but no fallback gracefully implemented
- Migration plan: Either add RANN as hard dependency or implement pure-R k-NN construction; document performance difference

**pracma Package (Orthogonalization fallback):**
- Risk: `pracma::orth()` used as fallback for QR-based orthogonalization
- Impact: If pracma becomes unavailable, initialization and retraction fail
- Current usage: Lines 200, 822 in penalized_mfa_clusterwise.R
- Migration plan: Replace with base R `qr.Q()` approach exclusively; remove pracma dependency

---

## Security Considerations

**No Input Validation on Spatial Coordinates:**
- Risk: If coords_list contains NaN, Inf, or malformed values, k-NN graph construction silently fails or produces garbage
- Files: `R/penalized_mfa_clusterwise.R` (line 705 checks row count but not value validity)
- Current mitigation: None
- Recommendations: Add explicit checks for finite values in coordinates; reject NaN/Inf; validate coordinate dimensions match data dimensions

**No Bounds Checking on Hyperparameters:**
- Risk: User can pass extreme values (e.g., learning_rate = 1e6, lambda = 1e10) causing numerical overflow/underflow
- Files: `R/penalized_mfa_clusterwise.R` validation at lines 3-33 checks basic ranges but not practical bounds
- Current mitigation: Validation uses `chk_range()` for beta1/beta2 but learning_rate only checked > 0
- Recommendations: Add practical bounds (e.g., 0.0001 ≤ learning_rate ≤ 1.0); warn if parameters outside typical ranges

**Attribute Stripping in Serialization:**
- Risk: If user saves penalized_mfa object and reloads it, attributes (V_list, obj_values, lambda) may be lost
- Files: `R/penalized_mfa.R` and `R/penalized_mfa_clusterwise.R` (all attr() calls)
- Impact: Print methods fail; auxiliary information vanishes
- Recommendations: Use proper object structure with list elements instead of attributes; ensure serialize/deserialize preserves data

---

## Test Coverage Gaps

**Empty and Edge-Case Blocks:**
- What's not tested: Blocks with zero rows, blocks with one column, blocks with more columns than rows
- Files: `R/penalized_mfa_clusterwise.R` (lines 810-813 attempt to handle k_s==0 but untested)
- Risk: Silent failures or dimension mismatch errors not caught until production
- Priority: HIGH (blocks are real in spatial clustering scenarios with missing subjects)

**Numerical Stability Under Extreme Conditions:**
- What's not tested: Nearly-collinear columns, extremely small singular values, nearly-singular Laplacian matrices
- Files: Objective calculation and SVD operations throughout
- Risk: Convergence may fail silently; NaN values may propagate without warning
- Priority: HIGH (real data often has numerical issues)

**Lambda Grid Selection (No Cross-Validation):**
- What's not tested: Systematic lambda selection, cross-validation for hyperparameter tuning
- Files: Example 12 in penalized_mfa_clusterwise.R (lines 641-659) is pseudo-code, not tested
- Risk: Users must manually choose lambda with no guidance; no automated selection implemented
- Priority: MEDIUM (affects practical usability but not correctness)

**Memory Budget Transitions:**
- What's not tested: Behavior when transitioning between precomputed and on-the-fly modes, boundary cases where memory_budget_mb equals exact computation cost
- Files: `R/penalized_mfa_clusterwise.R` (lines 73-77, 792-801)
- Risk: Could cause memory exhaustion or silent correctness issues at boundary values
- Priority: MEDIUM (affects large-scale deployments)

**Adam vs Gradient Descent Convergence:**
- What's not tested: Comparative convergence behavior; whether Adam actually outperforms gradient descent on real problems
- Files: `R/penalized_mfa_clusterwise.R` (optimizer parameter, lines 958-977)
- Risk: Adam may be less stable on manifold constraints; no evidence it's better than gradient descent
- Priority: MEDIUM (affects algorithmic soundness)

**Error Recovery in Optimization:**
- What's not tested: Behavior when SVD fails during initialization or retraction; when objective function becomes NA/Inf
- Files: `R/penalized_mfa_clusterwise.R` (lines 815-818, 927-935, 995-999 have try-catch but untested)
- Risk: Silent failures or unexpected error messages when numerical instability occurs
- Priority: HIGH (affects robustness in production)

---

## Missing Critical Features

**No Automated Lambda Selection:**
- Problem: Users must manually choose spatial penalty weight; no cross-validation or information criterion available
- Blocks: `penalized_mfa_clusterwise()` and `penalized_mfa()`
- Current state: Example 12 in penalized_mfa_clusterwise.R shows manual grid search but no implementation

**No Convergence Diagnostics:**
- Problem: No easy way to assess whether optimization converged or got stuck; convergence plots not generated
- Available: obj_values output from penalized_mfa_clusterwise can be plotted manually
- Missing: Automated convergence assessment, stopping criteria recommendations

**No Block-Specific Component Selection:**
- Problem: Some blocks might need different numbers of components; currently all blocks forced to same ncomp
- Files: `R/penalized_mfa_clusterwise.R` (ncomp_block computed as minimum at lines 758-761)
- Impact: Information lost when different blocks have different intrinsic dimensionality

**Missing Documentation for Internal Functions:**
- Problem: Many helper functions lack proper roxygen documentation
- Files: `R/penalized_mfa_clusterwise.R` has many @keywords internal functions with incomplete docs (spatial_constraints, block_update_cluster, etc.)
- Impact: Difficult to maintain or extend; unclear intent and contracts

---

## Missing Imports

**TODO Author Update:**
- Location: `/DESCRIPTION` (line 5)
- Issue: Comment says "# TODO: UPDATE AUTHORS" but is in the actual DESCRIPTION file
- Fix: Either remove comment or actually update to include additional authors if any

---

*Concerns audit: 2026-01-22*

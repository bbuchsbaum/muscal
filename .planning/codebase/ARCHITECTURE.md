# Architecture

**Analysis Date:** 2026-01-22

## Pattern Overview

**Overall:** S3 Generic-Based Multivariate Analysis Framework

This is an R package implementing "French school" multivariate statistical analysis methods. The architecture uses S3 generics as the primary dispatch mechanism, with pluggable preprocessing pipelines and standardized return types based on the `multivarious` package's projector classes.

**Key Characteristics:**
- S3 generic-based method dispatch with multiple input type support (list, multiblock, multidesign)
- Standardized preprocessing pipeline integration via `multivarious::prepper` objects
- Return values as `multivarious::multiblock_projector` or `multivarious::multiblock_biprojector` subclasses
- Auxiliary outputs stored as named list elements within returned objects
- Block-aware computations with per-block loading matrices and indices
- C++ acceleration via Rcpp for critical numerical operations (SVD, ridge solve)

## Layers

**Analysis Layer:**
- Purpose: User-facing analysis functions implementing multivariate methods
- Location: `R/mfa.R`, `R/bamfa.R`, `R/penalized_mfa.R`, `R/penalized_mfa_clusterwise.R`, `R/bada.R`, `R/covstatis.R`
- Contains: S3 generics and concrete methods for MFA, BaMFA, penalized MFA, BaDA, STATIS
- Depends on: Preprocessing layer, computation layer, utility layer
- Used by: End users and downstream analysis pipelines

**Generics Layer:**
- Purpose: Define S3 generic interfaces and re-export key generics from `multivarious`
- Location: `R/all_generic.R`
- Contains: Function definitions for `bada()`, `mfa()`, `penalized_mfa()`, `covstatis()`, `project_cov()`, `project_subjects()`, `project_covariate()`, `project()`, `reprocess()`
- Depends on: None (base S3)
- Used by: Analysis methods via `UseMethod()` dispatch

**Preprocessing Layer:**
- Purpose: Handle flexible block preprocessing before analysis
- Location: `R/utils.R`, `R/penalized_mfa.R`
- Contains: `prepare_block_preprocessors()` function, handling NULL/single/list preprocessor arguments
- Depends on: `multivarious::prep()`, `multivarious::fresh()`, `multivarious::init_transform()`
- Used by: All analysis functions for preprocessing setup

**Computation Layer:**
- Purpose: Core numerical algorithms and helper computations
- Location: `R/mfa.R`, `R/bamfa.R`, `R/penalized_mfa.R`, `R/covstatis.R`, `R/bada.R`
- Contains: Block normalization (MFA, RV, RV2), SVD computations, alternating optimization, ridge solvers
- Depends on: `RSpectra` (sparse SVD), `Matrix` (linear algebra), Rcpp/C++ kernels
- Used by: Analysis methods for numerical operations

**Utility/Helper Layer:**
- Purpose: General-purpose utilities for matrix operations, data synthesis, and component analysis
- Location: `R/utils.R`, `R/components_utils.R`, `R/synthdat.R`
- Contains: Matrix helpers, bootstrap summarization, synthetic data generation, component retention, loading reliability
- Depends on: Base R, `abind`, `dplyr`, `purrr`, `furrr`
- Used by: Analysis layer, testing, user-facing analysis workflows

**C++/Rcpp Layer:**
- Purpose: High-performance numerical operations
- Location: `src/ridge_solve.cpp`, `src/RcppExports.cpp`
- Contains: Ridge regression solver, C++ templates via RcppArmadillo and RcppEigen
- Depends on: Rcpp, RcppArmadillo, RcppEigen
- Used by: `ls_ridge()` function in bamfa.R, SVD operations

## Data Flow

**MFA Analysis Flow:**

1. User calls `mfa(data, preproc=..., ncomp=..., normalization=...)`
2. S3 dispatch routes to appropriate method (`mfa.list()`, `mfa.multiblock()`, or `mfa.multidesign()`)
3. Input conversion: `mfa.list()` converts list to multiblock object
4. Preprocessing: `prepare_block_preprocessors()` fits preprocessors independently to each block
5. Block normalization: `normalization_factors()` computes per-block weights (MFA, RV, RV2, Frob, or custom)
6. Concatenation: All preprocessed blocks concatenated column-wise into single matrix
7. SVD: `multivarious::svd_wrapper()` with optional weighting applied
8. Return: `multiblock_projector` object with `v` (loadings), `s` (scores), `sdev` (singular values), `block_indices`, and `alpha` (normalization factors)

**BaMFA Analysis Flow:**

1. User calls `bamfa(data, k_g=..., k_l=..., niter=..., preproc=...)`
2. S3 dispatch determines input type
3. Preprocessing: `prep_bamfa_blocks()` preprocesses each block independently
4. Concatenation: Blocks joined into global matrix
5. Alternating optimization loop (niter iterations):
   - Global step: SVD on concatenated matrix → global loadings `G`
   - Block-specific steps: For each block, ridge regression with optional sparsity penalty
   - Adam optimizer updates block loadings `B_i`
6. Return: `bamfa`/`multiblock_projector` with `B_list` (block-specific loadings), `G` (global loadings), `s` (scores)

**Penalized MFA Analysis Flow:**

1. User calls `penalized_mfa(data, ncomp=..., lambda=..., preproc=...)`
2. S3 dispatch to appropriate method
3. Preprocessing: Handled via `prepare_block_preprocessors()`
4. Alternating optimization loop:
   - SVD on concatenated matrix
   - Per-block ridge regression with L2/L1 penalties (opt. clusterwise penalization)
   - Adam optimizer for gradient-based updates
5. Return: `penalized_mfa`/`multiblock_projector` with penalty path and convergence details

**STATIS (Covariance) Flow:**

1. User calls `covstatis(data, ncomp=..., normalize=..., dcenter=...)`
2. Input: List of covariance matrices
3. Optional preprocessing: Double-centering (Gower transformation), Frobenius normalization
4. Similarity computation: Inner products between preprocessed matrices
5. SVD on similarity matrix → compromise structure
6. Return: `covstatis` with compromise loadings, component importance

**State Management:**
- Per-block preprocessing state: Stored in `multivarious::concat_pre_processors()` output, accessible via `$preproc` slot
- Block boundaries: Maintained via `block_indices` (list of column ranges for each block)
- Block normalization weights: Stored as `alpha` vector in analysis results
- Convergence state: For iterative methods, stored as `obj_values`, iteration counts, penalty values in result object via `...` mechanism

## Key Abstractions

**S3 Generic Methods:**
- Purpose: Enable polymorphism across input types without if-else chains
- Examples: `mfa()`, `bamfa()`, `penalized_mfa()`, `bada()`, `covstatis()`
- Pattern: Generic function calls `UseMethod("mfa")`, which dispatches to `mfa.list()`, `mfa.multiblock()`, or `mfa.multidesign()`

**Block Preprocessor Pipeline:**
- Purpose: Apply consistent preprocessing independently to each data block
- Examples: `prepare_block_preprocessors()` in `R/utils.R`, `prepare_block_preprocessors()` in `R/penalized_mfa.R`
- Pattern: Accepts NULL (no preprocessing), single `prepper` object (replicate per block), or list of `prepper` objects (one per block). Returns fitted preprocessors and preprocessed data.

**Multiblock Projector:**
- Purpose: Standardized return type representing a fitted multivariate analysis
- Examples: Results from `mfa()`, `bamfa()`, `penalized_mfa()` inherit from `multiblock_projector` or `multiblock_biprojector`
- Pattern: Contains `v` (loadings), `s` (scores), `sdev` (singular values), `preproc` (fitted preprocessors), `block_indices` (per-block column ranges). Named list elements store method-specific outputs.

**Normalization Factors:**
- Purpose: Compute per-block weighting schemes for MFA-style analyses
- Examples: `normalization_factors()` in `R/mfa.R`
- Pattern: Takes list of preprocessed blocks and normalization type string, returns numeric vector of weights (one per block)

**Block Index Mapping:**
- Purpose: Track which columns belong to which original data block after concatenation
- Examples: `block_indices` field in returned objects
- Pattern: List where element i contains integer vector of column indices for block i in concatenated matrix

## Entry Points

**Analysis Functions:**
- Location: `R/mfa.R`, `R/bamfa.R`, `R/penalized_mfa.R`, `R/bada.R`, `R/covstatis.R`, `R/penalized_mfa_clusterwise.R`
- Triggers: Direct user calls to `mfa()`, `bamfa()`, `penalized_mfa()`, etc.
- Responsibilities: Input validation, S3 dispatch, result assembly

**Core Method Implementations:**
- Locations:
  - `mfa.list()` and `mfa.multiblock()` in `R/mfa.R`
  - `bamfa.list()`, `bamfa.multiblock()`, `bamfa.multidesign()` in `R/bamfa.R`
  - `penalized_mfa.list()`, `penalized_mfa.multiblock()` in `R/penalized_mfa.R`
  - `covstatis.list()` in `R/covstatis.R`
- Triggers: S3 dispatch from generic functions
- Responsibilities: Argument parsing, preprocessing, core algorithm execution, result object construction

**Helper Functions:**
- Locations: `R/utils.R`, `R/components_utils.R`, `R/synthdat.R`
- Triggers: Called by analysis functions and user code
- Responsibilities: Preprocessing, matrix operations, data synthesis, component analysis

## Error Handling

**Strategy:** Validation-first with informative messages

**Patterns:**

- Input validation via `chk::chk_*()` functions for type, dimension, and range checks
- S3 method dispatch with `.default` methods catching unsupported input types
- Dimensionality checks: Ensure all blocks have same number of rows, consistent column counts after preprocessing
- Ridge solver fallback: If Rcpp implementation fails, fall back to R-based `solve()` with warning
- SVD method selection: Automatically choose between fast (`svds`) and base SVD based on matrix dimensions (fallback if svds fails)

Example from `R/bamfa.R`:
```r
coefs <- tryCatch(
  ridge_solve(Z, X, effective_lambda),
  error = function(e) {
    warning("Rcpp ridge_solve failed: ", e$message, ". Falling back to R.", call. = FALSE)
    # R fallback implementation
  }
)
```

## Cross-Cutting Concerns

**Logging:**
- Framework: `cli` package (colored console output)
- Pattern: Use `cli::cli_h1()`, `cli::cli_h2()`, `cli::cli_alert_*()` for structured messages
- Examples: Status messages during preprocessing, normalization type notifications, convergence details

**Validation:**
- Framework: `chk` package for assertions, `assertthat` for custom validation
- Pattern: Check input dimensions, types, and constraints at function entry
- Examples: `chk::chk_matrix()`, `chk::chk_gte()`, `chk::chk_list()`

**Authentication/Authorization:**
- Not applicable (analytical package with no auth layer)

**Block Awareness:**
- Pattern: All analyses maintain `block_indices` mapping to support per-block extraction and manipulation
- Examples: `project()` methods can work with per-block scores, loadings can be extracted per block via `block_indices`

**Preprocessing Integration:**
- Pattern: Every analysis function accepts flexible `preproc` argument integrated via `multivarious::prepper` interface
- Examples: Centering, scaling, normalization applied independently to each block before concatenation

**Auxiliary Output Storage:**
- Pattern: Method-specific outputs (convergence metrics, penalty values, intermediate matrices) stored as named list elements in projector object
- Examples: `obj_values` for objective function, `B_list` for block-specific loadings, `lambda` for penalty value


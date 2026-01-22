# Coding Conventions

**Analysis Date:** 2026-01-22

## Naming Patterns

**Files:**
- Lowercase with underscores: `mfa.R`, `bada.R`, `penalized_mfa.R`, `components_utils.R`
- Test files: `test-<topic>.R` format (e.g., `test-mfa.R`, `test-bada-multidesign.R`)
- Dashed separator in test names for clarity

**Functions:**
- Exported functions: snake_case (e.g., `normalization_factors()`, `prepare_block_preprocessors()`, `estimate_components()`)
- Internal functions prefixed with dot: `.extract_V_list()`, `.infer_n()`, `.pre_process_new_cov()`, `.get_G_scores()`
- S3 method dispatch: function.class pattern (e.g., `mfa.list()`, `mfa.multiblock()`, `bada.multidesign()`, `bootstrap.bada()`)
- Helper functions: descriptive snake_case names

**Variables:**
- Local variables: snake_case (e.g., `block_names`, `p_post`, `proclist`, `Xp`)
- Matrix/data objects: capital letters for mathematical objects (e.g., `M`, `V`, `G`, `X1`, `X2`, `Z`, `B`)
- Loop counters: single lowercase letters (`i`, `j`, `s`)
- Abbreviated: `ncomp` (number of components), `niter` (number of iterations), `nboot` (number of bootstrap samples)

**Types:**
- S3 classes define object types: `"mfa"`, `"bada"`, `"penalized_mfa"`, `"covstatis"`, `"bamfa"`
- Class hierarchy via inheritance: inherits from `"multiblock_projector"`

## Code Style

**Formatting:**
- Roxygen2 documentation with `#'` comments
- Roxygen2 configuration in DESCRIPTION: `RoxygenNote: 7.3.3`
- Indentation: 2 spaces (standard R convention)
- Line length: practical limits (some functions are long for mathematical clarity)

**Linting:**
- Package uses `chk` package for input validation (see `assertthat` and `chk` in imports)
- No ESLint or Prettier (R package, not JavaScript)
- Standards follow R package conventions

## Import Organization

**Order:**
1. Package documentation comments: `@useDynLib`, `@importFrom Rcpp`, `@import` for bulk imports
2. External package imports listed separately: listed in alphabetical order in roxygen blocks
3. Selective imports for specific functions: `@importFrom chk chk_list chk_integer...`

**Path Aliases:**
- No formal path aliases in R (not applicable)
- File organization: Related functions grouped in single files (e.g., `mfa.R` contains `mfa()` generic, `mfa.list()`, `mfa.multiblock()`)

**Example import block from `bamfa.R` (lines 1-22):**
```r
#' @useDynLib muscal, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @import RcppArmadillo
#' @import RcppEigen
#' @importFrom stats rnorm sd
#' @importFrom chk chk_list chk_integer chk_numeric chk_gte chk_flag chk_matrix chk_not_empty
#' @importFrom chk chk_true
#' @importFrom Matrix tcrossprod Diagonal sparseMatrix rankMatrix
#' @importFrom MASS ginv
#' @importFrom crayon bold magenta green blue yellow
#' @importFrom purrr map map2
#' @importFrom multivarious init_transform prep fresh center concat_pre_processors pass
#' @importFrom multidesign multiblock
#' @importFrom RSpectra svds
#' @importFrom rlang enquo quo_is_missing quo_get_expr
#' @importFrom dplyr select pull
#' @importFrom multivarious multiblock_projector
```

## Error Handling

**Patterns:**
- Input validation using `chk` package: `chk::chk_integer()`, `chk::chk_gte()`, `chk::chk_matrix()`
  - Example from `bamfa.R` (lines 191-202):
  ```r
  chk::chk_not_empty(data)
  chk::chk_list(data)
  chk::chk_integer(k_g)
  chk::chk_gte(lambda_l, 0)
  chk::chk_gte(tol, 0)
  ```

- Error messages with `stop()`: Include context and formatting with `sprintf()`
  - Example from `utils.R` (lines 35, 52-53, 59, 66):
  ```r
  stop("Package 'multivarious' needed for preprocessing. Please install it.", call. = FALSE)
  stop(sprintf("If 'preproc' is a list, its length (%d) must match the number of data blocks (%d).",
               length(preproc_arg), S), call. = FALSE)
  ```

- Warnings with `warning()`: Include context about what condition was detected and its impact
  - Example from `utils.R` (line 86):
  ```r
  warning("Preprocessing resulted in blocks with different numbers of columns. Subsequent steps might require consistent dimensions.", call.=FALSE)
  ```

- Messages with `message()`: Informational status updates
  - Example from `mfa.R` (line 39):
  ```r
  message("normalization type:", type)
  ```

- Optional logging with `cli` package: For verbose operations
  - Used in `penalized_mfa_clusterwise.R` with `cli::cli_alert_info()`, `cli::cli_alert_danger()`, `cli::cli_alert_success()`
  - Controlled by `verbose` parameter in function signatures

- Logical operators for fallback: `%||%` operator for NULL coalescing
  - Example from `components_utils.R` (line 39):
  ```r
  sdev <- sdev %||% fit$sdev %||% stop("No singular values available.")
  ```

## Logging

**Framework:** `cli` package (via `@importFrom cli cli_alert_info cli_alert_danger...`)

**Patterns:**
- Status messages: `cli::cli_alert_info()`
- Error/warning alerts: `cli::cli_alert_danger()`
- Success messages: `cli::cli_alert_success()`
- Section headers: `cli::cli_h1()`, `cli::cli_h2()`
- Bulleted lists: `cli::cli_ul()`

Example from `penalized_mfa_clusterwise.R`:
```r
if (verbose) cli::cli_alert_danger("Iter {iter}, Block {s}: Skipping update due to inconsistent V dimensions.")
if (verbose) cli::cli_alert_info("Iter {iter}: obj={format(new_obj, digits=6)}, rel_change={format(rel_change, scientific=TRUE, digits=2)}")
cli::cli_alert_success("Converged after {iter} iterations.")
```

## Comments

**When to Comment:**
- Complex algorithms: Roxygen `@details` sections for extensive explanation
- Non-obvious mathematical operations: Inline comments explaining math logic
- Helper function purpose: Brief comment above function body
- Example from `penalized_mfa.R` (lines 20-22):
```r
# Update biased first and second moment estimates
M <- beta1 * M + (1 - beta1) * G
V2 <- beta2 * V2 + (1 - beta2) * (G * G)
```

**JSDoc/TSDoc:**
- Not applicable (R, not TypeScript/JavaScript)
- Roxygen2 documentation instead: `#'` comments with `@param`, `@return`, `@details`, `@examples`

**Example roxygen block from `mfa.R` (lines 1-11):**
```r
#' Compute a similarity matrix from blocks of data
#'
#' Creates a symmetric matrix where each element [i,j] is the similarity between
#' blocks i and j, calculated using the supplied function.
#'
#' @param blocks A list of numeric matrices or data frames
#' @param FUN Function to compute similarity between two blocks
#' @param ... Additional arguments passed to FUN
#' @return A symmetric similarity matrix with dimensions length(blocks) × length(blocks)
#' @importFrom utils combn
#' @noRd
#' @keywords internal
```

## Function Design

**Size:**
- Functions typically 20-100 lines for core logic
- Helper functions: 5-40 lines
- Large complex functions: `penalized_mfa_clusterwise.R` functions can be 200+ lines due to optimization loops

**Parameters:**
- Required parameters first, defaults last
- Type specification through `chk::` validation
- Example from `mfa.list()` (line 82):
```r
mfa.list <- function(data, preproc=center(), ncomp=2,
                     normalization=c("MFA", "RV", "None", "Frob", "custom"),
                     A=NULL, M=NULL, ...)
```

**Return Values:**
- Named lists for complex returns: `list(Xp = Xp, proclist = proclist, p_post = p_post)`
- S3 objects with class attribute: `structure(result, class = c("mfa", "multiblock_projector"))`
- Vectors/matrices for simple returns

## Module Design

**Exports:**
- Declared with `@export` roxygen tag
- S3 methods for dispatching: generic + method implementations
- Example from `all_generic.R`: Generic functions like `mfa()`, `bada()`, `covstatis()` that dispatch to methods

**Barrel Files:**
- `all_generic.R` acts as export hub for all generic functions and methods
- `RcppExports.R` auto-generated from C++ sources (do not edit)
- `tibble-imports.R` minimal, re-exports tibble utilities

---

*Convention analysis: 2026-01-22*

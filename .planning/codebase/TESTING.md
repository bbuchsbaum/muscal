# Testing Patterns

**Analysis Date:** 2026-01-22

## Test Framework

**Runner:**
- `testthat` (version >= 3.0.0, listed in Suggests in DESCRIPTION)
- Config: Tests located in `/Users/bbuchsbaum/code/muscal/tests/testthat/`
- No dedicated config file; uses standard testthat convention

**Assertion Library:**
- `testthat::expect_*` family of functions (built into testthat)

**Run Commands:**
```bash
# Run all tests
devtools::test()
# or in R:
testthat::test_file("tests/testthat/test-mfa.R")

# Watch mode
devtools::test(filter = "mfa")

# Coverage
devtools::test_coverage()
# or
covr::package_coverage()
```

## Test File Organization

**Location:**
- Co-located pattern: Tests in separate `tests/testthat/` directory (standard R package structure)
- Source code in `/Users/bbuchsbaum/code/muscal/R/`
- Tests in `/Users/bbuchsbaum/code/muscal/tests/testthat/`

**Naming:**
- Pattern: `test-<topic>.R`
- Examples: `test-mfa.R`, `test-bada-multidesign.R`, `test-penalized-mfa.R`, `test-bamfa.R`, `test-utils-misc.R`
- Multiple test files for related functionality (16 test files total)

**Structure:**
```
tests/
└── testthat/
    ├── test-bada-multidesign.R
    ├── test-bamfa-recon.R
    ├── test-bamfa.R
    ├── test-covstatis.R
    ├── test-duality.R
    ├── test-internal-utils.R
    ├── test-mfa.R
    ├── test-penalized-mfa-clusterwise.R
    ├── test-penalized-mfa-opts.R
    ├── test-penalized-mfa.R
    ├── test-preproc-utils.R
    ├── test-project-covariate.R
    ├── test-summary-utils.R
    ├── test-synthdat.R
    ├── test-utils-extra.R
    └── test-utils-misc.R
```

## Test Structure

**Suite Organization:**
From `test-mfa.R`:
```r
library(testthat)
library(muscal)

test_that("mfa.list produces expected structure", {
  set.seed(1)
  blocks <- replicate(3, matrix(rnorm(20), nrow = 5), simplify = FALSE)
  res <- mfa(blocks, ncomp = 2)
  expect_s3_class(res, "mfa")
  expect_equal(length(res$block_indices), 3)
  # ... more assertions
})

test_that("mfa.multiblock validates input", {
  x1 <- matrix(rnorm(10), nrow = 5)
  expect_error(mfa(list(x1)))
  expect_error(mfa(list(x1, x2)))
  expect_error(mfa(list(x1, x1), normalization = "custom"))
})
```

**Patterns:**

1. **Setup Pattern:**
   - Load libraries: `library(testthat)` and `library(muscal)` at file top
   - Set seed for reproducibility: `set.seed(123)` before data generation
   - Create test data with `matrix()`, `data.frame()`, `replicate()`
   - Optional: Create helper functions (e.g., `make_cov()` in `test-internal-utils.R`)

2. **Test Organization Pattern:**
   - One test per logical unit: `test_that("description", { ... })`
   - Group related tests with comments: `# --- .function_name -------`
   - Multiple assertions per test allowed

3. **Assertion Pattern:**
   - Class checking: `expect_s3_class(obj, "class_name")`
   - Type checking: `expect_equal()`, `expect_true()`, `expect_false()`, `expect_null()`
   - Length/dimension checking: `expect_length()`, `expect_equal(dim())`
   - Error checking: `expect_error(expr, pattern)`
   - Warning checking: `expect_warning(expr, pattern)`

From `test-utils-misc.R` (lines 4-7):
```r
# Tests for %||% operator
expect_equal(NULL %||% 5, 5)
expect_equal(3 %||% 4, 3)
```

## Mocking

**Framework:**
- No formal mocking library used (unittest.mock, mockery, etc. not in imports)
- Manual test data generation using `matrix()`, `data.frame()`, seed-controlled randomness

**Patterns:**
From `test-bada-multidesign.R`:
```r
# Test data setup
set.seed(123)
x <- matrix(rnorm(20), nrow = 10)
design <- data.frame(
  subj_id = rep(c("S1", "S2"), each = 5),
  y = rep(c("A", "B"), 5)
)
md <- multidesign::multidesign(x, design)

res <- bada(md, y = y, subject = subj_id, ncomp = 1)

expect_equal(levels(res$subjects), c("S1", "S2"))
```

**What to Mock:**
- External package results: Create minimal reproducible data structures
- RNG-dependent operations: Use `set.seed()` for deterministic tests
- Package dependencies: Use `skip_if_not_installed()` for optional packages

From `test-bamfa.R`:
```r
skip_if_not_installed("multidesign")
mb <- multidesign::multiblock(list(A = matrix(0, 2, 2), B = matrix(0, 3, 2)))
expect_error(bamfa(mb, k_g = 1L, k_l = 1L, niter = 1L))
```

**What NOT to Mock:**
- Core functions under test: Test actual behavior, not mocked behavior
- Validation functions: Let `chk::` validation run in tests to verify error messages
- Mathematical operations: Use exact equality checks with small tolerances when needed

## Fixtures and Factories

**Test Data:**
From `test-internal-utils.R` (lines 4-12):
```r
set.seed(1)
make_cov <- function(seed) {
  set.seed(seed)
  A <- matrix(rnorm(9), 3, 3)
  tcrossprod(A)
}
subject_data <- lapply(1:3, make_cov)
res_none <- covstatis(subject_data, ncomp = 2, norm_method = "none", dcenter = FALSE)
res_mfa  <- covstatis(subject_data, ncomp = 2, norm_method = "mfa", dcenter = FALSE)
```

**Location:**
- Fixtures created inline in test files using base R functions
- No separate fixture files or factory libraries
- Use `set.seed()` and `replicate()` for generating test matrices

Example from `test-penalized-mfa.R` (lines 4-9):
```r
set.seed(123)
X1 <- scale(matrix(rnorm(100 * 4), 100, 4), scale = FALSE)
X2 <- scale(matrix(rnorm(100 * 16), 100, 16), scale = FALSE)
data_list_uneven <- list(B1 = X1, B2 = X2)
data_list_even <- list(B1 = X1, B2 = X1)
```

## Coverage

**Requirements:**
- No explicit coverage target enforced
- Coverage testing available through `devtools::test_coverage()` or `covr::package_coverage()`

**View Coverage:**
```bash
# In R console:
devtools::test_coverage()

# Or using covr package:
coverage <- covr::package_coverage()
print(coverage)
```

## Test Types

**Unit Tests:**
- Scope: Individual functions and methods
- Approach: Test function inputs, outputs, and error conditions
- Example from `test-mfa.R`:
```r
test_that("normalization_factors MFA matches svd", {
  X1 <- matrix(1:4, nrow = 2)
  X2 <- matrix(2:5, nrow = 2)
  nf <- muscal:::normalization_factors(list(X1, X2), type = "MFA")
  expect_equal(nf[1], 1/(svd(X1)$d[1]^2))
  expect_equal(nf[2], 1/(svd(X2)$d[1]^2))
})
```

**Integration Tests:**
- Scope: Multiple functions working together
- Approach: Test full pipeline from data input to result output
- Example from `test-bamfa.R`:
```r
test_that("bamfa.list returns projector with expected structure", {
  set.seed(123)
  blocks <- list(A = matrix(rnorm(20), 5, 4), B = matrix(rnorm(20), 5, 4))
  res <- bamfa(blocks, k_g = 2L, k_l = 1L, niter = 2L, preproc = multivarious::pass())
  expect_s3_class(res, "bamfa")
  expect_s3_class(res, "multiblock_projector")
  expect_equal(res$k_g, 2)
  expect_equal(res$k_l, 1)
  expect_length(res$B_list, 2)
  expect_length(res$block_indices, 2)
})
```

**E2E Tests:**
- Framework: Not explicitly used
- Some integration tests approach E2E scope (testing full algorithm pipelines)

## Common Patterns

**Async Testing:**
- Not applicable to R (no async/await pattern)
- Parallel processing using `furrr::future_map()` tested through deterministic output verification

**Error Testing:**
From `test-mfa.R` (lines 23-31):
```r
test_that("mfa.multiblock validates input", {
  x1 <- matrix(rnorm(10), nrow = 5)
  expect_error(mfa(list(x1)))

  x2 <- matrix(rnorm(12), nrow = 6)
  expect_error(mfa(list(x1, x2)))

  expect_error(mfa(list(x1, x1), normalization = "custom"))
})
```

**Warning Testing:**
From `test-penalized-mfa.R` (lines 12-18):
```r
test_that("penalized_mfa.list handles uneven and even blocks correctly", {
  # Test with uneven blocks: should warn and produce no consensus matrix
  expect_warning(
    res_uneven <- penalized_mfa.list(data_list_uneven, ncomp = 2, lambda = 1,
                                      compute_consensus = TRUE, verbose = FALSE,
                                      penalty_method = "projection"),
    "Cannot compute consensus"
  )
  expect_s3_class(res_uneven, "penalized_mfa")
  expect_null(attr(res_uneven, "consensus"))
})
```

**Private Function Testing:**
- Use `:::` operator to access internal functions (not exported)
- Example from `test-mfa.R` (line 38):
```r
nf <- muscal:::normalization_factors(list(X1, X2), type = "MFA")
```

**Skipping Tests:**
From `test-bamfa.R` (line 16):
```r
skip_if_not_installed("multidesign")
mb <- multidesign::multiblock(list(A = matrix(0, 2, 2), B = matrix(0, 3, 2)))
```

---

*Testing analysis: 2026-01-22*

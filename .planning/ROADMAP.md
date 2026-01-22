# Roadmap: muscal

**Created:** 2026-01-22
**Depth:** Quick (3-5 phases)
**Core Value:** Confidence in correctness for CRAN submission

## Overview

This roadmap delivers test coverage and code quality improvements for the muscal R package. The approach is tests-first: fix existing issues, add comprehensive tests to core modules, then tackle the complex penalized MFA modules, and finally refactor with the safety net of tests in place. Focus is on the 1191-line `penalized_mfa_clusterwise.R` which needs most attention.

## Phases

### Phase 1: Foundation Fix

**Goal:** R CMD check passes and all existing tests pass

**Dependencies:** None (starting phase)

**Requirements:**
- TEST-01: All existing tests pass (currently 1 failing)
- TEST-03: R CMD check passes with no errors
- QUAL-01: Fix deprecated chkor() calls (replace with chkor_vld())
- QUAL-02: Fix S3 method signature mismatches
- QUAL-03: Fix or wrap failing examples (genpca dependency)
- QUAL-04: Remove/fix unused imports in DESCRIPTION

**Success Criteria:**
1. Running `devtools::test()` reports 0 failures (all 150 tests pass)
2. Running `R CMD check` reports 0 errors
3. MFA examples run without requiring genpca (wrapped in dontrun or conditional)
4. No S3 method signature warnings from R CMD check

---

### Phase 2: Core Module Tests

**Goal:** Core analysis modules have comprehensive test coverage

**Dependencies:** Phase 1 (clean baseline needed)

**Requirements:**
- CORE-01: mfa.R has comprehensive tests (input validation, normalization types, edge cases)
- CORE-02: bamfa.R has comprehensive tests (k_g/k_l combinations, convergence, edge cases)
- CORE-03: covstatis.R has comprehensive tests (normalization, projection methods)
- CORE-04: bada.R has comprehensive tests (multidesign input, bootstrap)
- CORE-05: utils.R has comprehensive tests (preprocessing, helpers)

**Success Criteria:**
1. mfa.R tests cover all normalization types (MFA, RV, RV2, Frob, none)
2. bamfa.R tests verify convergence for various k_g/k_l combinations
3. bada.R tests verify bootstrap resampling produces stable results
4. Each core module has at least one edge case test (empty input, single block, single component)
5. Running test suite for core modules shows no failures

---

### Phase 3: Penalized MFA Tests

**Goal:** Penalized MFA modules have tests covering core functionality and edge cases

**Dependencies:** Phase 2 (test patterns established)

**Requirements:**
- PMFA-01: penalized_mfa.R has comprehensive tests (penalty methods, convergence, edge cases)
- PMFA-02: penalized_mfa_clusterwise.R has tests for core functionality
- PMFA-03: penalized_mfa_clusterwise.R has tests for edge cases (empty blocks, lambda=0, variable-rank)
- PMFA-04: penalized_mfa_clusterwise.R has tests for optimization paths (gradient vs Adam)
- TEST-02: Test coverage reaches 75% overall

**Success Criteria:**
1. penalized_mfa.R tests cover L1 and L2 penalty methods
2. penalized_mfa_clusterwise.R tests verify convergence with spatial constraints
3. Tests verify lambda=0 produces same results as penalized_mfa without spatial penalty
4. Tests verify gradient and Adam optimizers both converge
5. `covr::package_coverage()` reports >= 75% overall coverage

---

### Phase 4: Refactoring

**Goal:** penalized_mfa_clusterwise.R is decomposed into testable units

**Dependencies:** Phase 3 (test safety net in place)

**Requirements:**
- REFAC-01: Extract helper functions from penalized_mfa_clusterwise
- REFAC-02: Separate optimization logic from data preparation
- REFAC-03: Create testable units for Riemannian operations

**Success Criteria:**
1. penalized_mfa_clusterwise.R main function is under 500 lines
2. Extracted helper functions have their own unit tests
3. Riemannian operations (retraction, gradient) are in separate, documented functions
4. All existing tests still pass after refactoring (no regression)

---

## Progress

| Phase | Status | Requirements | Completed |
|-------|--------|--------------|-----------|
| 1 - Foundation Fix | Not Started | 6 | 0 |
| 2 - Core Module Tests | Not Started | 5 | 0 |
| 3 - Penalized MFA Tests | Not Started | 5 | 0 |
| 4 - Refactoring | Not Started | 3 | 0 |

**Total:** 0/19 requirements completed

---

## Coverage Map

| Requirement | Phase | Description |
|-------------|-------|-------------|
| TEST-01 | 1 | All existing tests pass |
| TEST-03 | 1 | R CMD check passes |
| QUAL-01 | 1 | Fix deprecated chkor() calls |
| QUAL-02 | 1 | Fix S3 method signature mismatches |
| QUAL-03 | 1 | Fix failing examples |
| QUAL-04 | 1 | Fix unused imports |
| CORE-01 | 2 | mfa.R comprehensive tests |
| CORE-02 | 2 | bamfa.R comprehensive tests |
| CORE-03 | 2 | covstatis.R comprehensive tests |
| CORE-04 | 2 | bada.R comprehensive tests |
| CORE-05 | 2 | utils.R comprehensive tests |
| PMFA-01 | 3 | penalized_mfa.R tests |
| PMFA-02 | 3 | penalized_mfa_clusterwise core tests |
| PMFA-03 | 3 | penalized_mfa_clusterwise edge cases |
| PMFA-04 | 3 | penalized_mfa_clusterwise optimization tests |
| TEST-02 | 3 | 75% coverage target |
| REFAC-01 | 4 | Extract helper functions |
| REFAC-02 | 4 | Separate optimization from data prep |
| REFAC-03 | 4 | Testable Riemannian operations |

**Coverage:** 19/19 requirements mapped

---

*Roadmap created: 2026-01-22*

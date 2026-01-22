# muscal

## What This Is

An R package implementing "French school" multivariate statistical analysis methods (MFA, BaMFA, penalized MFA, BADA, STATIS) for multi-block data analysis. Built on the `multivarious` framework with S3 generics, preprocessing pipelines, and C++ acceleration via Rcpp.

## Core Value

**Confidence in correctness.** Tests verify the algorithms work as intended, enabling safe CRAN submission and future refactoring.

## Requirements

### Validated

<!-- Existing functionality that works -->

- ✓ MFA (Multiple Factor Analysis) with list/multiblock/multidesign inputs — existing
- ✓ BaMFA (Barycentric MFA) with global/block-specific loadings — existing
- ✓ Penalized MFA with L2/L1 penalties and Adam optimizer — existing
- ✓ Penalized MFA Clusterwise with spatial constraints — existing
- ✓ BADA (Barycentric Discriminant Analysis) — existing
- ✓ STATIS for covariance matrices — existing
- ✓ Block normalization (MFA, RV, RV2, Frob) — existing
- ✓ S3 generic dispatch for multiple input types — existing
- ✓ Preprocessing pipeline integration via multivarious — existing
- ✓ C++ ridge solver via Rcpp — existing
- ✓ Bootstrap resampling for BADA — existing
- ✓ Synthetic data generation utilities — existing

### Active

<!-- Current goals -->

- [ ] Achieve 75% test coverage across package
- [ ] Fix R CMD check errors (example failures, documentation issues)
- [ ] Fix bugs discovered during testing
- [ ] Refactor penalized_mfa_clusterwise (1191 lines → smaller functions)
- [ ] Fix S3 method signature mismatches flagged by R CMD check
- [ ] Address private API usage (multivarious:::) for CRAN compliance

### Out of Scope

- Performance optimization of penalized_mfa_clusterwise — defer until after refactor
- New analysis methods — focus on quality of existing methods
- Vignettes/tutorials — CRAN submission doesn't require them

## Context

**Codebase state:**
- 16 test files in `tests/testthat/`
- Current coverage unknown (needs measurement)
- `penalized_mfa_clusterwise.R` is 1191 lines — monolithic, hard to test
- CONCERNS.md documents tech debt, fragile areas, test coverage gaps

**CRAN submission blockers (from CONCERNS.md):**
- MFA examples fail without genpca installed
- S3 method signature mismatches
- Private API access via `:::`
- Undocumented exported functions

**Key files to test:**
- `R/mfa.R` (220 lines)
- `R/bamfa.R` (525 lines)
- `R/penalized_mfa.R` (693 lines)
- `R/penalized_mfa_clusterwise.R` (1191 lines) — highest priority
- `R/bada.R` (345 lines)
- `R/covstatis.R` (517 lines)
- `R/utils.R` (381 lines)

## Constraints

- **Timeline**: Weeks — need to move quickly toward CRAN submission
- **Approach**: Tests first, refactor after — build confidence before changing code
- **Dependencies**: Must work with multivarious, genpca, multidesign ecosystem

## Key Decisions

| Decision | Rationale | Outcome |
|----------|-----------|---------|
| 75% coverage target | Personal threshold for confidence, not CRAN requirement | — Pending |
| Tests before refactor | Safety net needed before touching penalized_mfa_clusterwise | — Pending |
| Focus on penalized_mfa_clusterwise | Largest, most complex, least tested module | — Pending |

---
*Last updated: 2026-01-22 after initialization*

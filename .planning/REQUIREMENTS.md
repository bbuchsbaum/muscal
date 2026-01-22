# Requirements: muscal

**Defined:** 2026-01-22
**Core Value:** Confidence in correctness for CRAN submission

## v1 Requirements

Requirements for this milestone. Each maps to roadmap phases.

### Testing Infrastructure

- [ ] **TEST-01**: All existing tests pass (currently 1 failing)
- [ ] **TEST-02**: Test coverage reaches 75% overall
- [ ] **TEST-03**: R CMD check passes with no errors

### Core Module Tests

- [ ] **CORE-01**: mfa.R has comprehensive tests (input validation, normalization types, edge cases)
- [ ] **CORE-02**: bamfa.R has comprehensive tests (k_g/k_l combinations, convergence, edge cases)
- [ ] **CORE-03**: covstatis.R has comprehensive tests (normalization, projection methods)
- [ ] **CORE-04**: bada.R has comprehensive tests (multidesign input, bootstrap)
- [ ] **CORE-05**: utils.R has comprehensive tests (preprocessing, helpers)

### Penalized MFA Tests

- [ ] **PMFA-01**: penalized_mfa.R has comprehensive tests (penalty methods, convergence, edge cases)
- [ ] **PMFA-02**: penalized_mfa_clusterwise.R has tests for core functionality
- [ ] **PMFA-03**: penalized_mfa_clusterwise.R has tests for edge cases (empty blocks, lambda=0, variable-rank)
- [ ] **PMFA-04**: penalized_mfa_clusterwise.R has tests for optimization paths (gradient vs Adam)

### Code Quality

- [ ] **QUAL-01**: Fix deprecated chkor() calls (replace with chkor_vld())
- [ ] **QUAL-02**: Fix S3 method signature mismatches
- [ ] **QUAL-03**: Fix or wrap failing examples (genpca dependency)
- [ ] **QUAL-04**: Remove/fix unused imports in DESCRIPTION

### Refactoring

- [ ] **REFAC-01**: Extract helper functions from penalized_mfa_clusterwise
- [ ] **REFAC-02**: Separate optimization logic from data preparation
- [ ] **REFAC-03**: Create testable units for Riemannian operations

## v2 Requirements

Deferred to future work.

- **PERF-01**: Optimize k-NN graph construction (use RANN consistently)
- **PERF-02**: Improve memory efficiency for large datasets
- **DOC-01**: Add vignettes with usage examples
- **DOC-02**: Improve roxygen documentation completeness

## Out of Scope

| Feature | Reason |
|---------|--------|
| New analysis methods | Focus on quality of existing methods |
| Performance optimization | Defer until after refactor confirms correctness |
| 100% test coverage | 75% is sufficient for confidence |
| Comprehensive vignettes | Not required for CRAN submission |

## Traceability

| Requirement | Phase | Status |
|-------------|-------|--------|
| TEST-01 | Phase 1 | Pending |
| TEST-03 | Phase 1 | Pending |
| QUAL-01 | Phase 1 | Pending |
| QUAL-02 | Phase 1 | Pending |
| QUAL-03 | Phase 1 | Pending |
| QUAL-04 | Phase 1 | Pending |
| CORE-01 | Phase 2 | Pending |
| CORE-02 | Phase 2 | Pending |
| CORE-03 | Phase 2 | Pending |
| CORE-04 | Phase 2 | Pending |
| CORE-05 | Phase 2 | Pending |
| PMFA-01 | Phase 3 | Pending |
| PMFA-02 | Phase 3 | Pending |
| PMFA-03 | Phase 3 | Pending |
| PMFA-04 | Phase 3 | Pending |
| TEST-02 | Phase 3 | Pending |
| REFAC-01 | Phase 4 | Pending |
| REFAC-02 | Phase 4 | Pending |
| REFAC-03 | Phase 4 | Pending |

**Coverage:**
- v1 requirements: 19 total
- Mapped to phases: 19
- Unmapped: 0

---

*Requirements defined: 2026-01-22*
*Last updated: 2026-01-22 after roadmap creation*

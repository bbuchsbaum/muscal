---
phase: 01-foundation-fix
plan: 02
subsystem: core
tags: [s3-methods, api-compliance, r-package]

# Dependency graph
requires:
  - phase: none
    provides: n/a
provides:
  - Clean S3 method signatures matching multivarious generics
  - Public API usage only (no ::: calls)
  - CRAN-compliant bada.R
affects: [01-03, 01-04, testing, cran-submission]

# Tech tracking
tech-stack:
  added: []
  patterns:
    - S3 method signatures must match generic exactly
    - Use imported functions directly or with :: (never :::)
    - Extra method args passed via ... dots

key-files:
  created: []
  modified:
    - R/bada.R
    - NAMESPACE
    - man/bada.Rd

key-decisions:
  - "bootstrap.bada: data, alpha, verbose moved to ... to match generic bootstrap(x, nboot, ...)"
  - "reprocess.bada: block moved to ... to match generic reprocess(x, new_data, colind, ...)"
  - "project.bada: block moved to ... to match generic project(x, new_data, ...)"

patterns-established:
  - "S3 method args beyond generic signature go via ... dots"
  - "Use rlang %||% operator for default values from dots"

# Metrics
duration: 3min
completed: 2026-01-22
---

# Phase 01 Plan 02: S3 Signature Fix Summary

**Removed private API (:::) calls and fixed S3 method signature mismatches to match multivarious generics**

## Performance

- **Duration:** 3 min
- **Started:** 2026-01-22T21:22:21Z
- **Completed:** 2026-01-22T21:25:05Z
- **Tasks:** 2
- **Files modified:** 3 (R/bada.R, NAMESPACE, man/bada.Rd)

## Accomplishments

- Eliminated all triple-colon (:::) private API calls from bada.R
- Fixed bootstrap.bada signature to match multivarious::bootstrap generic
- Fixed reprocess.bada signature to match multivarious::reprocess generic
- Fixed project.bada signature to match multivarious::project generic
- R CMD check now passes without S3 signature warnings

## Task Commits

Each task was committed atomically:

1. **Task 1: Fix private API usage in bada.R** - `f4768b9` (fix)
2. **Task 2: Fix S3 method signature mismatches** - `849956a` (fix)

## Files Created/Modified

- `R/bada.R` - Core bada methods with corrected signatures and public API calls
- `NAMESPACE` - Added importFrom(rlang, %||%) for null-coalescing operator
- `man/bada.Rd` - Updated documentation reflecting signature changes

## Decisions Made

- **bootstrap.bada signature change:** The multivarious::bootstrap generic is `bootstrap(x, nboot, ...)`. The method-specific `data`, `alpha`, and `verbose` args now come via `...` and are extracted with `dots$arg %||% default`.
- **reprocess.bada signature change:** The generic is `reprocess(x, new_data, colind, ...)`. The `block` argument specific to bada now comes via `...`.
- **project.bada signature change:** The generic is `project(x, new_data, ...)`. The `block` argument specific to bada now comes via `...`.
- **Added rlang %||% import:** Used for providing defaults when extracting args from dots.

## Deviations from Plan

None - plan executed exactly as written.

## Issues Encountered

None - all changes were straightforward.

## User Setup Required

None - no external service configuration required.

## Next Phase Readiness

- S3 method signatures are now compliant with generics
- No private API usage remaining
- Ready for 01-03 (chkor deprecation fix) and subsequent phases
- R CMD check passes signature-related checks

---
*Phase: 01-foundation-fix*
*Completed: 2026-01-22*

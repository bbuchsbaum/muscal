---
phase: 01-foundation-fix
plan: 03
subsystem: testing
tags: [chk, mfa, genpca, testthat, R-CMD-check]

# Dependency graph
requires:
  - phase: 01-foundation-fix/01
    provides: "deprecated chkor() replaced with chkor_vld()"
  - phase: 01-foundation-fix/02
    provides: "S3 signature fixes"
provides:
  - "All 150 tests passing"
  - "MFA examples wrapped for R CMD check"
  - "chkor_vld validation fixed"
affects: [01-foundation-fix/04, phase-02]

# Tech tracking
tech-stack:
  added: []
  patterns:
    - "Use direct validation instead of chkor_vld for simple null checks"
    - "Wrap examples depending on Suggested packages in donttest{}"

key-files:
  created: []
  modified:
    - "R/mfa.R"
    - "tests/testthat/test-penalized-mfa.R"
    - "man/mfa.Rd"

key-decisions:
  - "Replaced chkor_vld() with direct is.null() check due to quosure evaluation issues in namespace context"
  - "Used donttest{} instead of dontrun{} so examples remain manually runnable"

patterns-established:
  - "Simple validation: Use direct is.null() checks instead of chkor_vld() when validating NULL arguments"
  - "Examples with Suggested deps: Wrap in donttest{} for R CMD check safety"

# Metrics
duration: 7min
completed: 2026-01-22
---

# Phase 01-03: Fix Failing Tests Summary

**Fixed 2 failing tests and wrapped MFA examples - all 150 tests now pass with clean R CMD check**

## Performance

- **Duration:** 7 min
- **Started:** 2026-01-22T21:27:21Z
- **Completed:** 2026-01-22T21:34:10Z
- **Tasks:** 2
- **Files modified:** 3

## Accomplishments
- Fixed mfa.R custom normalization validation that was silently passing when both A and M were NULL
- Updated test-penalized-mfa.R to match actual warning message for projection penalty with uneven blocks
- Wrapped MFA examples in \donttest{} blocks so they don't fail R CMD check when genpca is unavailable
- All 150 tests now pass (0 failures)

## Task Commits

Each task was committed atomically:

1. **Task 1: Identify and fix the failing test** - `ceb0d0a` (fix)
2. **Task 2: Fix or wrap MFA examples for genpca dependency** - `dfad867` (docs)

## Files Created/Modified
- `R/mfa.R` - Fixed custom normalization validation, wrapped examples in donttest{}
- `tests/testthat/test-penalized-mfa.R` - Updated expected warning message
- `man/mfa.Rd` - Regenerated with donttest{} example wrappers

## Decisions Made

1. **Replaced chkor_vld() with direct is.null() check** - The chkor_vld() function was not erroring when called from within the muscal namespace due to quosure evaluation context issues. Using a simple `if (is.null(A) && is.null(M)) stop(...)` pattern is more reliable.

2. **Used donttest{} instead of dontrun{}** - The \donttest{} wrapper skips examples during R CMD check but still allows users to run them manually when genpca is installed. This is preferable to \dontrun{} which marks examples as "never run."

3. **Updated test expectation instead of changing code** - The penalized_mfa warning message was intentionally changed to be more informative ("Cannot use 'projection' penalty with blocks of different dimensions; disabling penalty"). The test expectation was outdated.

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 1 - Bug] Fixed chkor_vld namespace evaluation issue**
- **Found during:** Task 1 (investigating test failure)
- **Issue:** Plan 01-01 changed `chkor(chk_not_null(...))` to `chkor_vld(vld_not_null(...))` but the latter doesn't error correctly when called from within the muscal namespace due to rlang quosure evaluation context
- **Fix:** Replaced with simple `if (is.null(A) && is.null(M)) stop(...)` pattern
- **Files modified:** R/mfa.R
- **Verification:** Test now passes, mfa(list(x1, x1), normalization="custom") correctly errors
- **Committed in:** ceb0d0a

---

**Total deviations:** 1 auto-fixed (1 bug)
**Impact on plan:** Bug fix was necessary for correct validation behavior. No scope creep.

## Issues Encountered

- Initially found 2 failing tests instead of the 1 mentioned in STATE.md
- The chkor_vld() function behaves differently when called from namespace context vs global environment - discovered through systematic debugging
- Both test failures were related to code/test expectation mismatches from plan 01-01 changes

## User Setup Required

None - no external service configuration required.

## Next Phase Readiness

- All tests pass (150/150)
- MFA examples properly wrapped for R CMD check
- Ready for plan 01-04 (final foundation cleanup)

---
*Phase: 01-foundation-fix*
*Completed: 2026-01-22*

---
type: quick-fix
target: R CMD check warnings
created: 2026-01-22
files_modified:
  - R/bada.R
  - R/penalized_mfa.R
  - R/utils.R
  - R/components_utils.R
  - R/synthdat.R
  - R/all_generic.R
  - man/*.Rd (via roxygen2)
autonomous: true
---

<objective>
Fix R CMD check warnings to achieve clean CRAN submission status.

Current state: 5 WARNINGs, 4 NOTEs
Target state: 0 WARNINGs (NOTEs about .claude/.planning are acceptable)

Purpose: Clean R CMD check is required for CRAN submission
Output: Updated R source files with regenerated documentation
</objective>

<context>
@.planning/STATE.md
@DESCRIPTION
</context>

<tasks>

<task type="auto">
  <name>Task 1: Fix missing/unexported dependency functions</name>
  <files>R/bada.R, R/penalized_mfa.R</files>
  <action>
Two functions are used but not exported from dependency packages:
1. `multidesign::split_by` - used in penalized_mfa.R line 438
2. `multivarious::coef.projector` - used in bada.R line 366

For `split_by`: Replace `multidesign::split_by(data, !!subject_quo)` with the base R equivalent using `split()` method that IS exported: `split(data, ...)` on the multidesign object.

For `coef.projector`: This is an S3 method, not directly exported. Replace `multivarious::coef.projector(x)` with `coef(x)` which will dispatch correctly. The coef generic is exported and the method is registered.
  </action>
  <verify>
Run `R -e "devtools::check(args = c('--no-manual', '--no-tests'), error_on = 'never')" 2>&1 | grep -E "split_by|coef.projector"` - should return empty (no warnings)
  </verify>
  <done>No warnings about missing/unexported objects from multidesign or multivarious</done>
</task>

<task type="auto">
  <name>Task 2: Fix Rd cross-reference links</name>
  <files>R/all_generic.R, R/bada.R, R/mfa.R</files>
  <action>
The link `[multivarious]{preprocessing pipeline}` is invalid Rd syntax.

In R/all_generic.R (bada generic) and R/mfa.R, replace:
- `\link[multivarious]{preprocessing pipeline}` with `\code{\link[multivarious]{center}}` or simply remove the link and use plain text like "a preprocessing pipeline from multivarious"

The preprocessing pipeline concept doesn't have a single help page - it's a collection of functions. Use `center()` as the representative link since that's the most commonly used preprocessor.
  </action>
  <verify>
Run `R -e "devtools::check(args = c('--no-manual', '--no-tests'), error_on = 'never')" 2>&1 | grep "preprocessing pipeline"` - should return empty
  </verify>
  <done>No warnings about missing cross-reference links</done>
</task>

<task type="auto">
  <name>Task 3: Document significant_components function</name>
  <files>R/utils.R</files>
  <action>
The `significant_components` function is exported (`@export`) but lacks documentation.

Add roxygen2 documentation block above the function (line 108):
```r
#' Determine Significant Components via RMT and Inter-block Coherence
#'
#' Identifies statistically significant components using Random Matrix Theory
#' (Marchenko-Pastur edge) and inter-block coherence tests.
#'
#' @param fit A fitted multiblock model containing `V_list` (list of loading matrices)
#'   and optionally `sdev` (singular values).
#' @param n Integer; number of observations (rows) in the original data.
#' @param k_vec Optional integer vector; number of columns in each block. If NULL,
#'   derived from V_list dimensions.
#' @param alpha Numeric; significance level for the ICC test (default 0.05).
#' @param check_rmt Logical; whether to apply the Marchenko-Pastur test (default TRUE).
#' @param tail_frac Numeric; fraction of eigenvalues to use for noise estimation (default 0.3).
#'
#' @return A list with elements:
#'   \item{keep}{Integer vector of component indices passing both tests}
#'   \item{rmt_pass}{Logical vector; which components pass the RMT test}
#'   \item{icc_pass}{Logical vector; which components pass the ICC test}
#'   \item{icc}{Numeric vector; inter-block coherence values}
#'   \item{icc_pvalue}{Numeric vector; p-values for ICC test}
#'   \item{mp_edge}{Marchenko-Pastur edge threshold}
#'   \item{sigma2_est}{Estimated noise variance}
#'   \item{lambda}{Eigenvalues (squared singular values)}
#'   \item{n, k_vec, alpha}{Input parameters}
#'
#' @export
```
  </action>
  <verify>
Run `R -e "devtools::document(); devtools::check(args = c('--no-manual', '--no-tests'), error_on = 'never')" 2>&1 | grep "significant_components"` - should return empty (no undocumented warning)
  </verify>
  <done>significant_components is documented and exported</done>
</task>

<task type="auto">
  <name>Task 4: Sync documentation with function signatures</name>
  <files>R/components_utils.R, R/synthdat.R, R/all_generic.R, R/covstatis.R</files>
  <action>
Multiple documentation/signature mismatches need fixing:

1. **estimate_components.Rd** - Missing `tail_q` and `V_ref`:
   - `tail_q` is in the function signature (line 36) but not documented
   - `V_ref` is a parameter for `loading_reliability`, add to @param

2. **prepare_block_preprocessors.Rd** - Wrong param names:
   - Change `@param data` to `@param data_list`
   - Change `@param preproc` to `@param preproc_arg`

3. **project.Rd** - Missing `x`, `new_data`, `...`, `colind`:
   - Add @param entries for all generic arguments

4. **reconstruct.covstatis.Rd** - Missing `...`:
   - Add `@param ... Additional arguments (currently unused)`

5. **synthetic_multiblock.Rd** - OLD vs NEW signature:
   - The roxygen comments document OLD params (p_list, ncomp, g_list, u, orthogonal_space, A, spatial)
   - The actual function has NEW params (S, n, p, r, sigma, sphere, k_nn, seed)
   - Replace the old @param entries with documentation for the actual parameters

After fixing roxygen comments, regenerate documentation with `devtools::document()`.
  </action>
  <verify>
Run `R -e "devtools::document(); devtools::check(args = c('--no-manual', '--no-tests'), error_on = 'never')" 2>&1 | grep -E "Undocumented arguments|Documented arguments not in"` - should return empty
  </verify>
  <done>All function arguments documented and matching signatures</done>
</task>

<task type="auto">
  <name>Task 5: Fix estimate_components call to significant_components</name>
  <files>R/components_utils.R</files>
  <action>
R CMD check detected: "possible error in significant_components(..., tail_quantile = tail_q): unused argument"

The issue is in estimate_components() at line 51:
```r
res <- significant_components(
  ...
  tail_quantile = tail_q  # WRONG - param is named tail_frac
)
```

The `significant_components` function uses `tail_frac`, not `tail_quantile`.

Change `tail_quantile = tail_q` to `tail_frac = tail_q`.
  </action>
  <verify>
Run `R -e "devtools::check(args = c('--no-manual', '--no-tests'), error_on = 'never')" 2>&1 | grep "tail_quantile"` - should return empty
  </verify>
  <done>Function call uses correct parameter name</done>
</task>

</tasks>

<verification>
After all tasks:
1. Run `devtools::document()` to regenerate all Rd files
2. Run `devtools::check()`
3. Verify: 0 WARNINGS related to documentation/exports
4. Acceptable NOTEs: hidden files (.claude, .planning), non-standard top-level files
</verification>

<success_criteria>
- R CMD check shows 0 WARNINGS (down from 5)
- All functions have matching documentation
- No missing cross-references
- No unexported dependency usage
</success_criteria>

<output>
After completion, commit changes with message:
"fix: resolve R CMD check warnings for CRAN submission"
</output>

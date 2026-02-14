# iPCA Implementation Plan (Tang & Allen 2021)

Date: 2026-02-13
Owner: muscal package
Status: In Progress

## Objective

Implement a production-ready multiplicative Frobenius iPCA estimator that is:
- fast in the common high-dimensional setting (`p_k >> n`),
- numerically stable,
- compatible with existing muscal API conventions,
- test-covered for correctness and regressions.

## Scope

In scope:
- New S3 generic `ipca()` plus `ipca.list()` and `ipca.multiblock()`.
- Multiplicative Frobenius Flip-Flop updates (Algorithm 1).
- Efficient sample-space/Gram update path using only `n x n` matrices during fitting.
- Optional dense update path for moderate dimensions and regression testing.
- Basic tuning entrypoint via explicit `lambda` vector.
- Output as `multivarious::multiblock_biprojector` subclass `"ipca"`.
- Unit tests for structure, numerical behavior, and dense-vs-gram consistency.

Out of scope (v1):
- Missing-data imputation tuning (Algorithms 7/9).
- L1/additive Frobenius penalties.
- Robust t-likelihood variant.
- Full plotting suite.

## Design Decisions

1. Default estimator is multiplicative Frobenius iPCA with strictly positive `lambda_k`.
2. Default solver mode is `method = "auto"`:
   - use Gram/sample-space updates when any `p_k > n`,
   - otherwise allow dense updates.
3. Fit using eigenvalue/eigenvector representations; avoid explicit matrix inverses.
4. Stabilize scale non-identifiability with optional trace normalization:
   - enforce `mean(diag(Sigma)) = 1` each iteration.
5. Convergence criterion:
   - max relative change in `Sigma` eigenvalues and in per-block `||Delta_k^{-1}||_F^2`.

## Work Breakdown

1. API and docs
- Add `ipca` generic in `R/all_generic.R`.
- Add roxygen docs and exports.

2. Core algorithm
- Add `R/ipca.R` with:
  - input checks/preprocessing,
  - shrinkage/eigen helpers,
  - `MF_iPCA_gram` core loop,
  - `MF_iPCA_dense` loop,
  - post-fit top-`ncomp` loading extraction.

3. Tests
- Add `tests/testthat/test-ipca.R`:
  - object structure + dimensions,
  - high-dimensional (`p>n`) fit success,
  - dense/gram subspace agreement on small synthetic data,
  - validation for invalid lambdas and mismatched rows.

4. Validation
- Run targeted test file and fix issues.

## Acceptance Criteria

1. `ipca()` is exported and documented.
2. Fitting succeeds for mixed block sizes including `p_k >> n` without constructing `p_k x p_k` matrices in Gram mode.
3. Returned object includes:
   - `s`, `sdev`, `v`, `block_indices`, block-aware `preproc`,
   - `Sigma_eigenvalues`, `Sigma_eigenvectors`,
   - per-block loading/eigen diagnostics.
4. Unit tests pass for new behavior.

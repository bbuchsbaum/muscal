# Orthonormal Score Refactor Plan

## Status

This note describes the intended long-term repair for the MFA-style models in
`muscal`:

- `anchored_mfa()` / `linked_mfa()`
- `aligned_mfa()`
- `graph_anchored_mfa()`
- `coupled_graph_anchored_mfa()`

The goal is to make orthonormal scores part of the fitted model rather than a
post-hoc QR reparameterization. The mathematically clean formulation is:

\[
\min \; f(S, \theta) \quad \text{s.t. } S^\top S = I_K,
\]

where `theta` collects the loading parameters (`B`, `V_k`, and in the coupled
model also `Z_k`).

At the moment, the worktree is in a transitional state:

- the older code paths historically used an unconstrained Euclidean score update
  followed by in-loop QR reparameterization,
- a first-pass Stiefel projected-gradient replacement has been started,
- but that replacement is not yet considered stable enough to merge as the only
  score-update path.

The immediate recovery plan is:

1. keep the stable unconstrained path available,
2. add the orthonormal path behind an explicit switch,
3. replace the generic manifold solver with a problem-specific orthonormal-score
   update,
4. only then decide whether low-level optimization is needed.


## Why The Old QR Step Is Wrong

The historical code updated `S` in Euclidean space and then orthonormalized it
inside the ALS loop:

- [R/linked_mfa.R](../R/linked_mfa.R)
- [R/aligned_mfa.R](../R/aligned_mfa.R)
- [R/graph_anchored_mfa.R](../R/graph_anchored_mfa.R)

For a pure reconstruction objective, that is just a change of basis. But once
the objective includes substantive penalties on the factors,

- ridge on `B` and `V_k`,
- graph penalties such as `tr(V^T L_V V)`,
- grouped-loading penalties such as `sum ||v_j - c_g||^2`,
- score-graph penalties such as `tr(S^T L_S S)`,

the non-orthogonal QR change of basis is no longer penalty-invariant.

So the old story is mathematically inconsistent:

- unconstrained update in one coordinate system,
- then non-orthogonal reparameterization to a different penalty geometry.

The correct repair is not "never orthonormalize". It is:

- choose `S^T S = I` as the identification scheme,
- optimize that constrained model directly.


## API Plan

Introduce a score constraint switch on the MFA-style models:

```r
score_constraint = c("none", "orthonormal")
```

Temporary default:

- `"none"` for the stable branch while the orthonormal solver is being finished.

Long-term intended default:

- `"orthonormal"` once the constrained path is correct, fast enough, and fully
  covered by tests.

This switch should be added consistently to:

- `anchored_mfa()`
- `aligned_mfa()`
- `graph_anchored_mfa()`
- `coupled_graph_anchored_mfa()`

Interpretation:

- `"none"` means the historical Euclidean score path,
- `"orthonormal"` means `S` is updated under the explicit constraint
  `S^T S = I`.


## Unified Score Subproblem

For all four MFA-style models, the score subproblem can be written as a
constrained quadratic on the Stiefel manifold:

\[
\min_{S^\top S = I} \;
\operatorname{tr}(S^\top \mathcal{A}(S)) - 2 \operatorname{tr}(C^\top S).
\]

The operator `A` is self-adjoint positive semidefinite, and `C` is an `N x K`
matrix.

This is the core observation that makes a generic manifold optimizer
unnecessary.


## Operator Form By Model

### 1. Anchored MFA / Linked MFA

Model:

\[
Y \approx S B^\top,\qquad X_k \approx S[idx_k, ] V_k^\top.
\]

Ignoring constants, the score objective is:

\[
f(S) =
\alpha_Y \|Y - S B^\top\|_F^2
 + \sum_k \alpha_k \|X_k - M_k S V_k^\top\|_F^2
 + \rho \|S\|_F^2.
\]

This gives:

\[
\mathcal{A}(S)
=
\alpha_Y S(B^\top B)
 + \sum_k \alpha_k D_k S (V_k^\top V_k)
 + \rho S
\]

and

\[
C
=
\alpha_Y Y B
 + \sum_k \alpha_k M_k^\top X_k V_k.
\]

Definitions:

- `M_k` is the row-selector / gather matrix implied by `row_index[[k]]`,
- `D_k = M_k^\top M_k` is diagonal with anchor-row repetition counts.


### 2. Aligned MFA

Model:

\[
X_k \approx S[idx_k, ] V_k^\top.
\]

The score objective is:

\[
f(S)
=
\sum_k \alpha_k \|X_k - M_k S V_k^\top\|_F^2.
\]

So:

\[
\mathcal{A}(S)
=
\sum_k \alpha_k D_k S (V_k^\top V_k),
\qquad
C
=
\sum_k \alpha_k M_k^\top X_k V_k.
\]


### 3. Graph-Anchored MFA

Add a score-graph penalty:

\[
\lambda_S \operatorname{tr}(S^\top L_S S).
\]

Then:

\[
\mathcal{A}(S)
=
\alpha_Y S(B^\top B)
 + \sum_k \alpha_k D_k S (V_k^\top V_k)
 + \lambda_S L_S S
 + \rho S,
\]

\[
C
=
\alpha_Y Y B
 + \sum_k \alpha_k M_k^\top X_k V_k.
\]


### 4. Coupled Graph-Anchored MFA

Model:

\[
Y \approx S B^\top,\qquad
X_k \approx Z_k V_k^\top,
\]

with coupling:

\[
\mu \sum_k \|Z_k - M_k S\|_F^2.
\]

Then:

\[
\mathcal{A}(S)
=
\alpha_Y S(B^\top B)
 + \mu \sum_k D_k S
 + \lambda_S L_S S
 + \rho S,
\]

\[
C
=
\alpha_Y Y B
 + \mu \sum_k M_k^\top Z_k.
\]

This is exactly the same type of operator as the previous models, which is why
one orthonormal-score update strategy should cover all of them.


## Recommended Constrained Score Update

Do **not** use a generic projected-gradient solver as the final implementation.

Use a problem-specific MM / generalized power iteration update for

\[
\min_{S^\top S = I} \operatorname{tr}(S^\top \mathcal{A}(S)) -
2 \operatorname{tr}(C^\top S).
\]

Given a current iterate `S_t`, let `beta` be an upper bound on the operator norm
of `A`. Form:

\[
G_t = C + \beta S_t - \mathcal{A}(S_t).
\]

Then update:

\[
S_{t+1} = \operatorname{polar}(G_t).
\]

Implementation notes:

- `polar(G)` can be computed as the orthonormal factor of `G`,
- practically: thin QR followed by a tiny `K x K` SVD if needed,
- this avoids line search and generic closure-driven optimization,
- this should give monotone descent for the constrained quadratic objective.


## Warm Start

Keep the old Euclidean score solve as a warm start only:

1. compute the old unconstrained `S_euc`,
2. set `S_0 = polar(S_euc)`,
3. run a small number of MM/GPI iterations,
4. accept the constrained update only if the constrained objective decreases.

This is clean because the Euclidean solve is merely an initializer, not a hidden
coordinate change inside the constrained model.


## Replace The Kronecker Solves Too

The graph code still builds large Kronecker systems. Even aside from the
orthonormal-score work, those are the wrong asymptotic objects.

Relevant existing code:

- `.gamfa_update_V_graph()` in
  [R/graph_anchored_mfa.R](../R/graph_anchored_mfa.R)
- `.gamfa_update_scores_graph()` in
  [R/graph_anchored_mfa.R](../R/graph_anchored_mfa.R)

These should move to matrix-operator CG.

### Graph V update

Solve for `V` in matrix form:

\[
\lambda_G L_G V + \text{blockwise}(V_k G_k) + \rho V = R.
\]

### Unconstrained graph S update

\[
\lambda_S L_S S + \alpha_Y S(B^\top B) + \sum_k \alpha_k D_k S G_k + \rho S = C.
\]

### Constrained graph S update

The same operator becomes the core of the MM/GPI score update.

So the implementation should revolve around:

- `apply_A_S(S, state)` for each model family,
- `apply_A_V(V, state)` for graph-smoothed loadings.


## Aligned Initializer

The aligned model should use a deterministic, permutation-equivariant
initializer tied to the actual objective.

Recommended initializer:

1. initialize `V_k` from each block's top `K` right singular vectors,
2. compute block scores `T_k = X_k V_k`,
3. aggregate to latent rows using `row_index`,
4. sum across blocks using block weights,
5. set `S_0 = polar(C_0)`.

This should replace any ad hoc random/jittered initialization strategy.


## Testing Contract

The test suite must distinguish the unconstrained and orthonormal paths.

### Tests that should apply to both paths

- reconstruction quality
- graph-off equivalence
- name-based `row_index` alignment
- edge-list vs adjacency equivalence
- zero-weight semantics
- prediction contract
- synthetic regime checks (shared vs coupled)

### Tests only for `score_constraint = "orthonormal"`

- `crossprod(S) ≈ I`
- constrained objective trace is nonincreasing
- projected-gradient / KKT residual is small:

\[
\|\nabla f(S) - S \, \mathrm{sym}(S^\top \nabla f(S))\|_F
\]

- permutation invariance of the constrained fit up to numerical tolerance

### Tests only for `score_constraint = "none"`

- do not assert orthonormality,
- compare subspaces/reconstructions rather than raw scores.


## Immediate Implementation Order

1. Restore one stable branch with `score_constraint = "none"` as default.
2. Keep the current noncontroversial fixes:
   - name alignment,
   - zero-weight handling,
   - adjacency validation,
   - removed over-hard `min p_k` caps where appropriate.
3. Remove the generic projected-gradient Stiefel solver as the final plan.
4. Implement operator-based MM/GPI for orthonormal scores.
5. Replace graph Kronecker solves with matrix-CG.
6. Re-profile.
7. Only then decide whether any RcppEigen kernels are warranted.


## Rcpp Guidance

Rcpp is **not** the first move.

Only after the algorithm is correct and operator-based should we consider
RcppEigen for:

- gather/scatter kernels for `M_k^T (X_k V_k)` and repetition counts,
- operator application `S -> A(S)`,
- matrix-CG kernels if pure R overhead remains measurable.

Do **not** optimize the current generic manifold backtracking path in C++; that
would speed up the wrong architecture.

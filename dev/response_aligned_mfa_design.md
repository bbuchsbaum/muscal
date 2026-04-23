# Response-Aligned MFA Design Note

## Summary

The `anchor_mfa` family should be centered on a shared, response-informed
latent space for multiblock `X`, with optional anchor structure when genuine
shared response states exist. Exact repeated rows of `Y` are one instantiation
of anchor structure, not the definition of the family.

This note records the target architecture and the staged implementation plan.
The current implementation slice in this branch adds `response_aligned_mfa()`
as the missing no-anchor supervised member of the family.

## Unifying Objective

The intended internal objective is

```math
\eta_0 \|R_0^{1/2}(Y^{(a)} - S B^\top)\|_F^2
+ \sum_{k=1}^K \Big[
  \alpha_k \|X_k - Z_k V_k^\top\|_F^2
  + \eta_k \|R_k^{1/2}(Y_k - Z_k B^\top)\|_F^2
  + \mu_k \|M_k^{1/2}(Z_k - A_k S)\|_F^2
\Big]
+ \lambda_V \mathcal P_V(V)
+ \lambda_S \mathcal P_S(S).
```

When all `\mu_k = 0` and there is no anchor-level response `Y^(a)`, the
engine should drop `S`, `A_k`, `M_k`, and `\lambda_S \mathcal P_S(S)`
entirely rather than carry degenerate zero-weight scaffolding.

Objects:

- `X_k`: block `k`
- `Y_k`: optional paired multivariate response for block `k`
- `Y^(a)`: optional anchor-level response table
- `Z_k`: block-specific latent scores
- `S`: anchor-state latent scores
- `B`: shared response loading matrix
- `V_k`: block loading matrix for block `k`
- `A_k`: anchor map from block rows to anchor states
- `M_k`: rowwise anchor-confidence weights
- `R_k`: rowwise response weights

## Why This Is the Right Center

This separates three things that should not be conflated:

1. Prediction of `Y` from `X`
   via `Y_k \approx Z_k B^\top`
2. Alignment of the `X` blocks in a common latent space
   via shared `B` and optional feature coupling on `V_k`
3. Extra borrowing from genuine shared anchor states
   via `Z_k \approx A_k S`

That is more principled than stretching "repeated rows of `Y`" into cases
where `Y` is simply a multivariate regression target.

## Special Cases

- `aligned_mfa()`
  No `Y`, hard row linkage, `Z_k = A_k S`
- `anchored_mfa()`
  Anchor response table only, hard `Z_k = A_k S`
- `graph_anchored_mfa()`
  Anchored MFA plus graph penalties on `V` and/or `S`
- `coupled_graph_anchored_mfa()`
  Anchored MFA with soft block-specific scores `Z_k`
- `response_aligned_mfa()`
  Blockwise paired `(X_k, Y_k)`, no anchor structure, `\mu_k = 0`
- Mixed repeated-plus-novel `Y`
  Supported by `A_k`, `M_k`, and `R_k` without inventing a separate method

## Primitive Inputs

- `anchor_map` should be primitive
  `row_index` is the one-hot special case
- `anchor_weight` should be primitive
  unanchored rows must have zero coupling weight rather than a zero anchor
  target
- `response_weight` should be primitive
  partial, missing, or low-confidence `Y_k` should be represented directly

In the soft-anchor case, row-stochastic `A_k` should be the default: rows with
positive anchor weight should be non-negative and sum to one. The geometry of
the association belongs in `A_k`; confidence belongs in `M_k`.

The `R_k` objects are rowwise in v1. They support missing rows and row-level
confidence, but they do not yet support missing components within a multivariate
response row. Entrywise response masks are a later extension.

## Identifiability

Identifiability should be stated as part of the model.

- If anchor structure is active, constrain the reference score basis directly,
  e.g. `S^T S = I_r` or a weighted analogue.
- If there is no anchor structure, impose a pooled score normalization across
  the block-specific scores, e.g. `sum_k w_k Z_k^T Z_k = I_r`.
- Then apply deterministic component order and sign conventions.

This is more robust than relying on `B^T B = I`, which is unavailable whenever
the response dimension or effective response rank is smaller than the latent
rank.

Order/sign conventions only resolve separated components. Tied or nearly tied
roots remain identifiable only at the subspace level.

## Prediction Contract

Prediction should be defined at test time, not inferred from the training
objective. Pure prediction and conditional completion should be kept separate.

For a new row in block `k`, the default `predict()` path is

```math
x_\text{new} \mapsto \hat z_k(x_\text{new}) \mapsto \hat y_\text{new},
```

with

```math
\hat z_k
=
\arg\min_z\;
\alpha_k \|x_\text{new} - z V_k^\top\|^2
+ \mu_k \|m^{1/2}(z - a S)\|^2
+ \rho \|z\|^2,
```

using only terms for information actually available at prediction time.

If `anchor_map_new` is unavailable, that term is omitted rather than
approximated with training-time coupling structure.

Conditional completion is a separate path: when partial `Y_new` is explicitly
provided, the latent solve may add a response-refinement term

```math
\eta_k \|r^{1/2}(y_\text{obs} - z B^\top)\|^2,
```

but that should not be the default semantics of `predict()`. In v1 this should
be exposed only behind an explicit flag such as `conditional = TRUE`.

## Object Contract

In the no-anchor case there is no shared observation-level `S`, only the
block-specific score matrices `Z_k`.

For fit-object compatibility:

- `Z_list` is the authoritative score representation;
- `s` may store a concatenated stack of the `Z_k` scores;
- `score_index` maps rows of `s` back to blocks.
- the fit should mark this explicitly, e.g. with `score_representation =
  "stacked_block_scores"`.

Any plotting, projection, or wrapper logic for this family should treat
`Z_list` rather than `s` as the primary score object.

## Weighting and Scaling Defaults

The defaults need to be part of the method specification, not left implicit.

- `\alpha_k` defaults to the existing MFA-style predictor scaling on `X_k`.
- `\eta_k` defaults to unit block weight after response preprocessing, so the
  supervision term operates in the shared preprocessed response space.
- `\mu_k` is interpreted relative to row-stochastic `A_k` and rowwise
  `M_k`; it should be documented as a per-block anchor-coupling strength, not
  as a raw scale on `A_k`.
- `Y_k` should be preprocessed in a shared response space before fitting. In
  v1, centering is the default and stronger standardization is explicit.

## Alignment Caveat

Shared `B` aligns blocks only to the extent that `Y` is informative.

If the response has low effective rank, response coupling can align at most
that many latent directions. Remaining directions must be aligned through
anchor coupling, `P_V`, or be treated explicitly as residual dimensions.

## Implementation Plan

1. Generalize the coupled engine around `Z_k`, `B`, `V_k`, optional `S`,
   optional `A_k`, optional `M_k`, optional `R_k`, optional `Y^(a)`, and the
   penalty objects `P_V`, `P_S`.
2. Expose `response_aligned_mfa()` first as the no-anchor supervised case.
3. Keep `anchored_mfa()` strict and honest about exact anchor linkage.
4. Re-express the current anchored methods as wrappers or special cases of the
   generalized engine.
5. Add `aligned_rrr()` later as a pure reduced-rank regression baseline.

## Current Slice

This implementation slice adds:

- this design note;
- `response_aligned_mfa()` as the supervised member of the family, with optional
  soft anchor coupling via `anchor_map`, `anchor_weight`,
  `coupling_lambda`, and optional anchor-level responses `Y^(a)`;
- a hard-anchor core used by `anchored_mfa()`, so strict anchored fits now
  expose the same family-level objects (`S`, `Z_list`, `score_index`) rather
  than only a bare shared score matrix;
- an explicit pooled score normalization rule;
- an explicit test-time latent projection contract for prediction.
- `aligned_rrr()` as the pure shared reduced-rank regression baseline, with the
  same pooled response preprocessing and stacked-score object contract but no
  `X` reconstruction term or anchor machinery.

The private engine now unifies the anchored and response-aligned wrappers,
including optional score/row-graph support on the anchored side. The baseline
comparison model `aligned_rrr()` remains a separate, narrower method by
design.

## Naming Note

If this architecture becomes the true family center, then `anchor_mfa` is no
longer a good umbrella name. A broader family label such as linked MFA or
response-aligned MFA would be more faithful, while `anchored_mfa()` remains the
exact-anchor wrapper.

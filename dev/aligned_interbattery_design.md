# Aligned Interbattery / Bidirectional MFA Design Note

## Summary

This note sketches a symmetric sibling to `response_aligned_mfa()`: a
bidirectional encoding-decoding model for two multiblock sides, `X` and `Y`,
with a shared latent space in the middle.

The target use case is:

- one or more subjects;
- each subject may contain one or more blocks / domains on each side;
- one side may have block-specific feature sets, missing blocks, or block-sparse
  structure;
- the other side may have a common feature schema across subjects;
- we want to predict `Y` from `X`, predict `X` from `Y`, and recover a common
  latent geometry that aligns subjects / rows / blocks.

The defining principle should be:

> not “repeated rows of `Y`”, and not “one side is response by definition”,
> but a shared latent space linking two multiblock sides, with optional row-
> alignment structure and an optional CCA-like decorrelation knob.

Working public name: `aligned_interbattery()`.
If the package naming should stay closer to existing MFA wrappers,
`bidirectional_aligned_mfa()` is a reasonable alternative. Conceptually,
this is the symmetric sibling of `response_aligned_mfa()`.

## Setting and Notation

Let subjects be indexed by `s = 1, ..., S`.

For each subject `s`:

- `X_sd` is block / domain `d` on the `X` side;
- `Y_sg` is block / domain `g` on the `Y` side.

Allow blocks to be missing. Allow different `X` blocks to have different column
schemas. If some block types share the same columns across subjects, their
loading matrices may be shared by block type rather than re-estimated per
subject.

To handle non-identical row sets across blocks, let

- `A^X_sd` map rows of `X_sd` into a subject-level reference row set;
- `A^Y_sg` map rows of `Y_sg` into a subject-level reference row set.

Exact row sharing is the identity-map special case. Soft or sparse row maps are
allowed when blocks only partially align at the row level.

Latent objects:

- `U_s`: `X`-side subject scores;
- `V_s`: `Y`-side subject scores;
- `T_s`: shared subject compromise scores;
- `P_b`: `X`-side loading matrix for block type `b`;
- `Q_c`: `Y`-side loading matrix for block type `c`.

Here `b = type_X(s, d)` and `c = type_Y(s, g)` optionally let blocks of the same
schema share loadings across subjects.

## Unifying Objective

A canonical internal objective is

```math
\begin{aligned}
\mathcal L
=&\; \sum_{s,d} \alpha_{sd}
\left\| (W^X_{sd})^{1/2}\big(X_{sd} - A^X_{sd} U_s P_{b(sd)}^\top\big) \right\|_F^2 \\
&\; + \sum_{s,g} \beta_{sg}
\left\| (W^Y_{sg})^{1/2}\big(Y_{sg} - A^Y_{sg} V_s Q_{c(sg)}^\top\big) \right\|_F^2 \\
&\; + \tau_X \sum_s \|U_s - T_s\|_F^2
    + \tau_Y \sum_s \|V_s - T_s\|_F^2 \\
&\; + \lambda_{\mathrm{row}} \sum_s \operatorname{tr}(T_s^\top L_s T_s) \\
&\; + \omega\, \Omega_{\mathrm{decorr}}(U, V)
    + \lambda_P \, \mathcal P_X(P)
    + \lambda_Q \, \mathcal P_Y(Q).
\end{aligned}
```

Objects:

- `W^X_sd`, `W^Y_sg`: rowwise weights / masks for missingness or confidence;
- `L_s`: optional row-graph Laplacian on the subject reference rows;
- `Omega_decorr(U, V)`: within-side decorrelation / whitening penalty;
- `P_X(P)`, `P_Y(Q)`: graph, group, block-sparse, or other loading penalties.

Interpretation:

1. the first line reconstructs `X` from `X`-side scores;
2. the second line reconstructs `Y` from `Y`-side scores;
3. the coupling terms tie both sides to a common middle space `T_s`;
4. the row-graph term says nearby rows should stay nearby in the shared space;
5. the decorrelation term gives a CCA-like knob;
6. the loading penalties encode block structure or prior feature similarity.

When rows that are similar in `Y` should be close in latent space, the natural
place to encode that is the graph term on `T_s`, with `L_s` derived from a
similarity graph on the relevant `Y`-side rows. More generally, `L_s` can come
from `X`, `Y`, or external metadata; the motivating case is `Y`-derived row
geometry.

## Why This Is the Right Center

This separates four things that should not be conflated:

1. `X <-> Y` prediction;
2. bidirectional reconstruction of both sides;
3. alignment of subjects / blocks in a common latent space;
4. within-side decorrelation in the CCA sense.

A single shared score matrix can express a symmetric joint factor model, but it
cannot cleanly represent the distinction between

- side-specific latent coordinates,
- a shared compromise space,
- and a user-controlled CCA-like whitening pressure.

That is why the minimal rigorous architecture uses `U_s`, `V_s`, and `T_s`
rather than only one score matrix.

## Decorrelate / CCA Continuum

The user-facing parameter should be `decorrelate`, with default `0`.

Two implementation options are compatible with the same public API.

### v1 default: score decorrelation penalty

A practical first implementation is

```math
\Omega_{\mathrm{decorr}}(U, V)
=
\sum_s \left\|\operatorname{offdiag}(n_s^{-1} U_s^\top U_s)\right\|_F^2
+
\sum_s \left\|\operatorname{offdiag}(n_s^{-1} V_s^\top V_s)\right\|_F^2.
```

Then:

- `decorrelate = 0` means purely reconstruction / covariance oriented;
- larger `decorrelate` pushes within-side latent coordinates toward
  canonical-style decorrelation.

This is easy to add without changing the rest of the model.

### v2 refinement: partial whitening continuum

If a more literal PLS-to-CCA continuum is desired, keep the same outer model but
replace the coupling metric with partially whitened scores:

```math
\tilde U_s = U_s (\Sigma_U + \rho I)^{-\gamma/2},
\qquad
\tilde V_s = V_s (\Sigma_V + \rho I)^{-\gamma/2},
```

where `gamma in [0, 1]` and `Sigma_U`, `Sigma_V` are pooled score covariances.
Then couple `tilde U_s` and `tilde V_s` to `T_s` instead of `U_s` and `V_s`.

- `gamma = 0`: covariance / PLS-like;
- `gamma = 1`: fully whitened / CCA-like.

The rest of the architecture stays unchanged.

## Prediction Contract

Prediction must be defined using only information available at test time.

### `Y` from `X`

For a new subject or row bundle with observed `X` blocks only:

1. estimate `hat U_s` from the observed `X`-side blocks;
2. form `hat T_s` from `hat U_s` via the fitted coupling rule;
3. decode the requested `Y` blocks via the `Q` loadings.

### `X` from `Y`

For observed `Y` blocks only:

1. estimate `hat V_s` from the observed `Y`-side blocks;
2. form `hat T_s` from `hat V_s`;
3. decode the requested `X` blocks via the `P` loadings.

### Joint completion

If both sides are partially observed, estimate `hat U_s` and `hat V_s` from the
available blocks and solve for a compromise `hat T_s` using only those observed
terms. This is a completion / fusion mode, not the default meaning of pure
out-of-side prediction.

Optional `row_graph_new`, `x_row_map_new`, or `y_row_map_new` may refine the
projection when available, but they must not be hidden requirements of
`predict()`.

## Identifiability

Identifiability must be stated as part of the model.

A robust default is to normalize the compromise basis directly:

```math
\sum_s w_s T_s^\top T_s = I_r,
```

with non-negative subject weights `w_s`, followed by deterministic component
ordering and sign conventions.

Then align `U_s` and `V_s` to `T_s` through the coupling terms. This is cleaner
than trying to identify the model through the loading matrices alone.

As usual, order/sign conventions only identify separated components. If roots are
tied or nearly tied, only the corresponding latent subspace is identifiable.

## Primitive Inputs

The primitive inputs should be:

- nested `X` and `Y` block structures;
- optional `x_row_map` and `y_row_map`;
- optional block-type metadata for sharing loadings across repeated schemas;
- rowwise weights on both sides;
- an optional `row_graph` (or directly a Laplacian);
- an optional `decorrelate` parameter;
- side-specific loading penalties.

Notes:

- identity row maps are the exact shared-row special case;
- sparse row maps are the soft-alignment case;
- if a side has a common feature set across subjects, loadings for that block
  type should be shareable by default;
- if a block type has no shared schema, it simply gets its own loading matrix.

## Proposed Public API

A reasonable first wrapper is:

```r
aligned_interbattery(
  X,
  Y,
  x_row_map = NULL,
  y_row_map = NULL,
  x_block_type = NULL,
  y_block_type = NULL,
  ncomp = 2,
  x_preproc = multivarious::center(),
  y_preproc = multivarious::center(),
  x_weight = NULL,
  y_weight = NULL,
  row_graph = NULL,
  row_graph_lambda = 0,
  couple_x = 1,
  couple_y = 1,
  decorrelate = 0,
  decorrelate_type = c("penalty", "whiten"),
  x_loading_graph = NULL,
  y_loading_graph = NULL,
  x_loading_lambda = 0,
  y_loading_lambda = 0,
  share_loadings = c("by_schema", "none"),
  use_future = FALSE,
  ...
)
```

Suggested semantics:

- `X`, `Y`: nested lists `subject -> block`, or flat block lists with explicit
  `subject` metadata;
- `x_row_map`, `y_row_map`: nested lists matching `X` and `Y`;
- `x_block_type`, `y_block_type`: metadata for sharing loadings across repeated
  schemas;
- `row_graph`: optional graph / kernel / Laplacian on subject reference rows;
- `couple_x`, `couple_y`: strength of the side-to-compromise coupling;
- `decorrelate`: CCA-like decorrelation strength;
- `decorrelate_type`: v1 penalty mode vs. v2 whitening mode.

The fit object should expose at least:

- `common_scores`: `T_s` by subject;
- `x_scores`: `U_s` by subject;
- `y_scores`: `V_s` by subject;
- `x_loadings`, `y_loadings` by block type / block;
- `partial_fits` or per-block reconstructions;
- `predict()` for `Y_from_X`, `X_from_Y`, and joint completion;
- `project()` / `scores()` methods with `side = c("common", "x", "y")`.

## Special Cases

- **Single common score system**
  If `couple_x` and `couple_y` are effectively infinite, then `U_s = V_s = T_s`
  and the model reduces to a symmetric shared-score joint factor model.

- **One block per side, shared rows, no penalties**
  This becomes a classical interbattery / joint latent model, with prediction in
  both directions.

- **`decorrelate = 0`**
  Reconstruction / covariance-oriented behavior, closer to multiblock PLS or
  joint factor modeling.

- **Strong decorrelation**
  More canonical / CCA-like latent coordinates.

- **One-sided supervision**
  If only one side is reconstructed / predicted, the model collapses toward
  `response_aligned_mfa()`.

- **Identity row maps**
  Standard shared-row multiblock setting.

- **Sparse row maps**
  Soft row alignment across blocks within subject.

- **Add private factors (v2)**
  Adding `X`-private and `Y`-private score subspaces yields a JIVE / IBFA-style
  extension when only part of the variation is cross-side.

## What This Method Is and Is Not

This method should be the symmetric, bidirectional sibling of the
response-informed MFA family. It is not just `response_aligned_mfa()` with the
letters swapped.

The defining feature is the separation of

- side-specific score systems,
- a shared compromise space,
- optional row-geometry priors,
- and optional CCA-like decorrelation.

That separation is what makes the family broad enough to cover joint factor
models, PLS-like behavior, and CCA-like behavior under one roof.

## Implementation Plan

1. Build the general engine around `U_s`, `V_s`, `T_s`, row maps, optional row
   graphs, and block-type loading sharing.
2. Start with `decorrelate_type = "penalty"` and add partial whitening later if
   needed.
3. Expose `aligned_interbattery()` as the first public wrapper.
4. Make `predict()` explicit about direction: `from = "x"`, `from = "y"`, or
   `from = "both"`.
5. Keep private-factor structure out of v1 unless reconstruction quality clearly
   demands it.

## Relationship to Existing Family

- `response_aligned_mfa()` remains the asymmetric, `X -> Y`-centered sibling.
- `aligned_interbattery()` is the symmetric `X <-> Y` sibling.
- `aligned_mcca()` remains the correlation-first counterpart when the focus is
  on canonical compromise rather than bidirectional reconstruction.

So the architectural picture becomes:

- `aligned_mfa()` / anchored family: shared-row or anchor-centered factor models;
- `response_aligned_mfa()`: asymmetric response-informed common-space model;
- `aligned_interbattery()`: symmetric bidirectional common-space model;
- `aligned_mcca()`: correlation-first compromise model.

# Q06 discovery consolidation

The former nonstandard specification file is preserved as Q06-original-spec.R.
Its twelve cases map to discoverable coverage as follows:

| Original property | Discoverable coverage |
|---|---|
| Anchored low-rank recovery | test-linked-mfa.R: strong shared low-rank fit |
| Linked/anchored alias | test-aligned-models-synthetic.R: exact alias |
| Matched loading shrinkage | test-linked-mfa.R: feature grouping |
| Aligned latent recovery | test-aligned-models-synthetic.R: synthetic shared scores; meaningful canonical overlap threshold |
| Named row indices | test-aligned-models-synthetic.R: named indices; same RNG, tighter stopping, 1e-8 threshold |
| Rank above narrow block width | test-aligned-models-synthetic.R: no forced narrowest-width cap |
| Zero-weight feature-prior isolation | test-aligned-models-synthetic.R: zero-weight blocks |
| MCCA latent recovery | test-aligned-models-synthetic.R: synthetic aligned data |
| Anchored/aligned MCCA equivalence | test-aligned-models-synthetic.R: identical None normalization |
| All-zero weights rejected | test-aligned-models-synthetic.R: rejects all-zero weights |
| Anchored objective descent | test-linked-mfa.R: objective trace |
| Aligned objective descent | test-aligned-mfa.R: objective trace/relabeling |

Helper names are prefixed to avoid collisions with graph-model fixtures. These
fixtures test declared software behavior and synthetic recovery only.

# Musca TODO List

This file tracks tasks related to improving consistency and functionality based on the principles in `design.md`.

## Consistency Refactoring

- **`penalized_mfa`:**
    - [ ] Refactor return value construction (`penalized_mfa.list`) to use the `...` mechanism for storing auxiliary results (e.g., `obj_values`, `lambda`, `consensus`, `V_list`) as named list elements within the `multiblock_projector`, instead of using `attr()`.
    - [ ] Update `print.penalized_mfa` to access these auxiliary results using `$` instead of `attr()`.

- **`penalized_mfa_clusterwise`:**
    - [ ] Create an S3 generic `penalized_mfa_clusterwise`.
    - [ ] Move primary documentation from the current function definition to the new generic.
    - [ ] Rename the current function to `penalized_mfa_clusterwise.list` (or similar appropriate input type) and use `@rdname penalized_mfa_clusterwise`.
    - [ ] Create a `penalized_mfa_clusterwise.default` method with an informative error message for unsupported types.
    - [ ] Refactor return value construction to use the `...` mechanism for storing auxiliary results (e.g., `Sadj`, `LV`, `obj_values`, `lambda`, `V_list`) as named list elements within the `multiblock_projector`, instead of using `attr()`.
    - [ ] Update `print.penalized_mfa_clusterwise` to access these auxiliary results using `$` instead of `attr()`.

## Potential Enhancements

- [ ] Consider adding calculation of scores (`s = Xp %*% v`) to `penalized_mfa` and `penalized_mfa_clusterwise` to potentially enable returning a `multiblock_biprojector` instead of just a `projector`, if the computational cost is acceptable and conceptually aligned.
- [ ] Add usage examples for the `significant_components` utility function in the documentation of relevant analysis methods (e.g., `mfa`, `penalized_mfa`).

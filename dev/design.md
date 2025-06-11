# Musca Package Design Principles

This document outlines key design principles for ensuring consistency and maintainability across the multi-block analysis methods within the `musca` package.

## 1. S3 Generics and Methods

- Each core analysis method (e.g., `mfa`, `bamfa`, `penalized_mfa`) should be implemented as an S3 generic function.
- Provide specific S3 methods for common input data structures (`list`, `multiblock`, `multidesign`, `hyperdesign`) *where the method's logic and input requirements are applicable* to that structure. For example, `gpca_align` specifically requires `hyperdesign` because it needs per-block design tables containing the response variable `y`.
- Implement a `.default` method for each generic. This method should catch unsupported input types and issue an informative error message, guiding the user towards supported input classes (e.g., directing them to use `hyperdesign` or `multiblock`).

## 2. Documentation

- The primary user-facing documentation, including `@description`, `@details`, `@param` descriptions for common parameters, `@return` value overview, `@examples`, and `@seealso` should reside with the S3 **generic** function definition.
- Documentation for specific S3 **methods** should be minimal. Use `@rdname generic_function_name` to link back to the main generic documentation.
- Method-specific documentation should only cover deviations, specific requirements (e.g., why `gpca_align.hyperdesign` is the main method), or implementation details unique to that method.

## 3. Preprocessing Argument (`preproc`)

- Analysis functions that perform internal preprocessing should accept a `preproc` argument.
- This argument should consistently accept:
    - `NULL`: No preprocessing is performed by the function.
    - A single `multivarious::prepper` object: This definition is applied independently to each input data block.
    - A `list` of `multivarious::prepper` objects: One `prepper` definition per input block. The list length must match the number of blocks.
- Utilize the `prepare_block_preprocessors` utility function (defined in `R/utils.R`) internally to handle this logic consistently.

## 4. Return Values

- **Prefer Standard Projector Classes:** Where the analysis conceptually results in a projection (i.e., computes loading vectors `v` defining a subspace), the function should return a standard `multivarious` projector class (`multivarious::multiblock_projector` or `multivarious::multiblock_biprojector` if scores `s` and singular values `sdev` are also naturally computed).
- **Correct Preprocessor Slot:** The `preproc` slot within the returned projector object **must** contain the block-aware preprocessor created using `multivarious::concat_pre_processors`. This ensures that subsequent projection methods applied to the result object use the correct, potentially block-specific, preprocessing pipeline.
- **Auxiliary Results via `...`:** Method-specific outputs beyond the core projector components (e.g., convergence details, penalty values, intermediate results like `B_list` in `bamfa`, the original `V_list`, etc.) should be stored as **named list elements** within the returned object. This is achieved by passing them as named arguments via the `...` mechanism to the `multiblock_projector` or `multiblock_biprojector` constructor functions. **Do not use `attr()`** for storing these results.
- **Method-Specific Class:** Add a method-specific class name (e.g., "mfa", "bamfa", "penalized_mfa") to the class list of the returned object. This allows for specific method dispatch (e.g., for `print` methods) while retaining compatibility with the `multivarious` projector ecosystem.

## 5. Print Methods

- Provide a dedicated `print.method_name` S3 method for each core analysis function.
- This method should be tailored to the specific return structure (typically a `multiblock_projector` or `biprojector` subclass).
- Access auxiliary results stored within the object using the standard list accessor `$` (e.g., `x$lambda`, `x$obj_values`).
- Provide a clear, informative summary suitable for the specific analysis performed.

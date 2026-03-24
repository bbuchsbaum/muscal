# Perform a single Adam update step on a block's loadings.

This function updates the loadings matrix `V` using one step of the Adam
optimizer. It maintains and returns the updated first (`M`) and second
(`V2`) moment estimates. The update is performed in the ambient space.

## Usage

``` r
adam_update_block(
  V,
  G,
  M,
  V2,
  step_count,
  beta1,
  beta2,
  adam_epsilon,
  learning_rate
)
```

## Arguments

- V:

  Current loading matrix (\\p \times k\\).

- G:

  The gradient of the loss with respect to V (\\p \times k\\).

- M, V2:

  Current first and second moment estimate matrices.

- step_count:

  The global step number, used for bias correction.

- beta1, beta2:

  Adam hyperparameters (e.g., 0.9, 0.999).

- adam_epsilon:

  A small constant to prevent division by zero (e.g., 1e-8).

- learning_rate:

  The base step size for the update.

## Value

A list containing the updated \`V\`, \`M\`, and \`V2\` matrices.

# Construct a Small Deterministic Design-MFA Example

Returns a balanced six-subject trialwise field dataset used by the
runnable examples for the Design-MFA API. The fixture is deterministic
and small enough for documentation and teaching.

## Usage

``` r
design_mfa_example_data()
```

## Value

A list with named \`Y\` and \`design\` subject blocks and
\`feature_names\`.

## Examples

``` r
example_data <- design_mfa_example_data()
vapply(example_data$Y, nrow, integer(1))
#> S1 S2 S3 S4 S5 S6 
#> 30 30 30 30 30 30 
```

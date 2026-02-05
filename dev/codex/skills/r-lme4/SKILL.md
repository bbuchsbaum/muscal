---
name: r-lme4
description: Use the R package lme4 for its core workflows (Linear Mixed-Effects Models using 'Eigen' and S4). Trigger on requests to write, review, or debug R code that uses lme4, including selecting functions, interpreting outputs, and troubleshooting common errors.
---

# lme4 (R)

## Quick start

1. Confirm install/version: `packageVersion("lme4")`
2. Load: `library(lme4)`
3. Find entrypoints: see `references/help-index.txt` and `references/exports.txt`

## Common workflows

### Fit a linear mixed model (`lmer`)

```r
library(lme4)

m <- lmer(reaction ~ days + (1 + days | subject), data = sleepstudy)
summary(m)
```

Use `REML = FALSE` when comparing fixed effects via likelihood-ratio tests:

```r
m0 <- lmer(reaction ~ 1 + (1 + days | subject), data = sleepstudy, REML = FALSE)
m1 <- lmer(reaction ~ days + (1 + days | subject), data = sleepstudy, REML = FALSE)
anova(m0, m1)
```

### Fit a generalized linear mixed model (`glmer`)

```r
library(lme4)

m <- glmer(
  cbind(incidence, size - incidence) ~ period + (1 | herd),
  data = cbpp,
  family = binomial
)
summary(m)
```

### Specify common random-effects structures

- Random intercept: `(1 | group)`
- Random intercept + slope: `(1 + x | group)`
- Remove intercept–slope correlation: `(1 + x || group)`
- Nested: `(1 | g1 / g2)` (equivalent to `(1 | g1) + (1 | g1:g2)`)
- Crossed: `(1 | subject) + (1 | item)`

### Extract effects, variance components, and model matrices

```r
fixef(m)          # fixed effects
ranef(m)          # random effects (BLUPs)
coef(m)           # per-group conditional coefficients
VarCorr(m)        # variance/covariance of random effects
sigma(m)          # residual SD (LMMs)

getME(m, "X")     # fixed-effects model matrix
getME(m, "Z")     # random-effects model matrix
getME(m, "theta") # variance parameters (internal scale)
```

### Predict with/without random effects

```r
predict(m)                         # include fitted random effects
predict(m, re.form = NA)           # population-level (fixed effects only)
predict(m, newdata = nd, re.form = NA)
predict(m, newdata = nd, allow.new.levels = TRUE) # new groups in nd
```

### Uncertainty and model comparison

- Prefer confidence intervals via `confint(m, method = "profile")` or bootstrapping (`bootMer`) over “Wald SEs only” when it matters.
- `lme4` does not provide default p-values for `lmer()`; for tests consider likelihood-ratio tests (`anova()` with `REML=FALSE`), parametric bootstrap, or companion packages like `lmerTest` (if appropriate for your workflow).

## Troubleshooting

### "boundary (singular) fit"

1. Confirm with `isSingular(m, tol = 1e-4)`.
2. Simplify random effects: drop slopes/correlations (use `||`), remove weak terms, or use a more parsimonious grouping structure.
3. Center/scale predictors (especially for random slopes) and check for near-zero variance components.

### "Model failed to converge" / gradient warnings

1. Center/scale predictors; check collinearity and separation (GLMMs).
2. Increase iterations / change optimizer:

```r
ctrl <- lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
m2 <- update(m, control = ctrl)
```

3. Inspect `summary(m2)$optinfo$conv` and consider simplifying the random-effects structure.

### GLMM fit issues (binomial/Poisson)

- Check for (quasi-)separation, zero-inflation, and overly complex random-effects terms relative to data size.
- Consider alternative families/links or specialized packages (e.g., zero-inflation) when diagnostics indicate misspecification.

## References

- `references/package-description.txt`: package metadata (installed)
- `references/exports.txt`: exported symbols
- `references/help-index.txt`: help index dump
- `references/vignettes.tsv`: vignette names + titles (when present)

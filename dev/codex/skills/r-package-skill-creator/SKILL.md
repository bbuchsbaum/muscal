---
name: r-package-skill-creator
description: Create Codex skills for specific R packages/libraries (e.g., lme4, dplyr, ggplot2) by scaffolding a new skill folder and optionally extracting lightweight reference dumps (exports, help index, vignettes) from the installed package. Use when asked to “make a skill for R package X”, “scaffold a Codex skill for an R library”, or to standardize a repeatable template for R-package skills.
---

# R Package Skill Creator

## Quick start

1. Pick a package name (e.g., `lme4`) and an output folder (e.g., `dev/codex/skills/`).
2. Run `python scripts/scaffold_r_pkg_skill.py lme4 --out dev/codex/skills`.
3. Edit the generated skill’s `SKILL.md` to add:
   - 3–6 “common workflows” users actually ask for
   - 3–6 troubleshooting bullets for common errors
4. Validate the generated skill folder (see “Validation” below).

## Workflow

### 1) Scaffold the skill

Use the bundled script to generate a consistent skeleton (`SKILL.md`, `agents/openai.yaml`, and `references/` dumps when the package is installed).

Recommended naming convention for generated skills:
- `r-<pkg>` (example: `r-lme4`)

### 2) Draft the trigger description

In the generated skill’s YAML frontmatter `description`, include:
- The package name, and the *actual tasks* users ask for (not package history)
- A few key function names (when useful)
- Typical contexts: debugging `.R` / `.Rmd` / Quarto code, interpreting output, choosing options

Keep it short but trigger-friendly; the description is the main trigger.

### 3) Curate the generated `SKILL.md`

Use the skeleton as a starting point, then add:
- **Common workflows**: 3–6 canonical tasks, each with a minimal code snippet (≤ ~10 lines)
- **Troubleshooting**: common errors/warnings and what to check first
- **Decision points**: “If you need X, use function Y; otherwise use Z”

Do *not* paste whole help pages into `SKILL.md`; put bulk material under `references/` and link to it.

### 4) Add references (progressive disclosure)

The scaffold script can write:
- `references/package-description.txt`
- `references/exports.txt`
- `references/help-index.txt`
- `references/vignettes.tsv` (when present)

Keep these as “raw dumps” to consult when writing the curated sections.

### 5) Validation

Run the system validator if available in your environment:

`python $CODEX_HOME/skills/.system/skill-creator/scripts/quick_validate.py <path/to/generated-skill>`

If `$CODEX_HOME` isn’t set, run the validator you have access to (or just ensure the YAML frontmatter has only `name` and `description`).

## Script

### `scripts/scaffold_r_pkg_skill.py`

Create a new `r-<pkg>` skill folder.

Example:

`python scripts/scaffold_r_pkg_skill.py lme4 --out dev/codex/skills`

Options:
- `--skill-name r-lme4` to override the default `r-<pkg>` naming
- `--force` to overwrite an existing skill folder
- `--no-references` to skip extracting package docs (still scaffolds the skill)

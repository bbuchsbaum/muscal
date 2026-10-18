---
name: r-package-skill-creator
description: "Create a Codex skill for a specific R package, with optional lightweight references from an installed package."
---

# R Package Skill Creator

Create a focused `r-<package>` skill that helps Codex choose and use the
package's real public API. Use the bundled scaffold when it saves work; do not
overwrite an existing skill unless the user asked to update it.

## Scaffold

```sh
python scripts/scaffold_r_pkg_skill.py <package> --out <skills-directory>
```

The script creates `SKILL.md`, `agents/openai.yaml`, and, when the package is
installed, lightweight reference dumps. Use `--skill-name` only when the normal
`r-<package>` name is unsuitable, `--no-references` when installed-package
metadata is not useful, and `--force` only for an authorized replacement.

## Curate the result

Keep the trigger description short: name the package and the concrete tasks the
skill owns. In the root `SKILL.md`, retain only the package-specific decisions,
common workflows, pitfalls, and completion checks that materially improve the
model's work.

Put bulk exports, help indexes, package metadata, and vignette inventories under
`references/`. Link them from the root with a condition that tells the model
when to read each file. Do not copy whole help pages into the entrypoint or turn
every exported function into a workflow.

Use examples when they clarify a package-specific choice, not to fill a quota.
Prefer short public-API examples that can be checked against the installed
package. Preserve uncertainty when the local package version or optional
dependencies are unavailable.

## Validate

Run the system skill validator when available:

```sh
python "${CODEX_HOME:-$HOME/.codex}/skills/.system/skill-creator/scripts/quick_validate.py" <skill-directory>
```

Also inspect the generated description, reference links, and examples. For a
substantial package skill, test a representative request or runnable example so
validation covers behavior rather than only YAML structure.

Finish with a usable, discoverable skill and report what was generated,
validated, and left unverified.

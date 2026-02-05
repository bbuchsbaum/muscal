#!/usr/bin/env python3
from __future__ import annotations

import argparse
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path


PKG_RE = re.compile(r"^[A-Za-z][A-Za-z0-9.]*$")
SKILL_RE = re.compile(r"^[a-z0-9]+(?:-[a-z0-9]+)*$")


def _die(message: str, exit_code: int = 2) -> "None":
    print(f"error: {message}", file=sys.stderr)
    raise SystemExit(exit_code)


def _write_text(path: Path, content: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content, encoding="utf-8")


def _run(cmd: list[str]) -> subprocess.CompletedProcess[str]:
    return subprocess.run(cmd, check=False, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)


def _rscript_exists() -> bool:
    return shutil.which("Rscript") is not None


def _extract_references(pkg: str, refs_dir: Path) -> None:
    if not _rscript_exists():
        _die("Rscript not found; rerun with --no-references or install R")

    refs_dir.mkdir(parents=True, exist_ok=True)

    r_code = f"""
pkg <- "{pkg}"
out_dir <- "{refs_dir.as_posix()}"

if (!requireNamespace(pkg, quietly=TRUE)) {{
  cat("Package not installed:", pkg, "\\n", file=stderr())
  quit(status=3)
}}

desc <- tryCatch(utils::packageDescription(pkg), error=function(e) NULL)
if (is.null(desc)) {{
  cat("Unable to read packageDescription() for:", pkg, "\\n", file=stderr())
  quit(status=4)
}}

writeLines(capture.output(print(desc)), con=file.path(out_dir, "package-description.txt"), useBytes=TRUE)

exports <- tryCatch(getNamespaceExports(pkg), error=function(e) character())
if (length(exports) > 0) {{
  writeLines(sort(unique(exports)), con=file.path(out_dir, "exports.txt"), useBytes=TRUE)
}}

h <- tryCatch(utils::help(package=I(pkg)), error=function(e) NULL)
if (!is.null(h) && !is.null(h$info) && length(h$info) >= 2) {{
  writeLines(h$info[[2]], con=file.path(out_dir, "help-index.txt"), useBytes=TRUE)
}}

if (!is.null(h) && !is.null(h$info) && length(h$info) >= 3 && is.matrix(h$info[[3]])) {{
  v <- h$info[[3]]
  tsv <- paste(v[,1], v[,2], sep="\\t")
  writeLines(tsv, con=file.path(out_dir, "vignettes.tsv"), useBytes=TRUE)
}}
"""

    proc = _run(["Rscript", "--vanilla", "-e", r_code])
    if proc.returncode != 0:
        details = proc.stderr.strip() or proc.stdout.strip()
        raise SystemExit(f"Rscript failed ({proc.returncode}): {details}")


def _infer_title_from_pkg_desc(refs_dir: Path) -> str | None:
    path = refs_dir / "package-description.txt"
    if not path.exists():
        return None

    # The printed packageDescription() includes a line like: "Title: ..."
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        m = re.match(r"^Title:\s*(.+?)\s*$", line)
        if m:
            return m.group(1).strip().strip('"')
    return None


def _default_skill_md(skill_name: str, pkg: str, title: str | None) -> str:
    title_bit = f" ({title})" if title else ""
    return f"""---
name: {skill_name}
description: Use the R package {pkg} for its core workflows{title_bit}. Trigger on requests to write, review, or debug R code that uses {pkg}, including selecting functions, interpreting outputs, and troubleshooting common errors.
---

# {pkg} (R)

## Quick start

1. Confirm install/version: `packageVersion("{pkg}")`
2. Load: `library({pkg})`
3. Find entrypoints: see `references/exports.txt` and `references/help-index.txt`

## Common workflows (fill in)

- [TODO] Add 3–6 workflows users ask for (each with a short snippet).

## Troubleshooting (fill in)

- [TODO] Add common errors/warnings and what to check first.

## References

- `references/package-description.txt`: package metadata (installed)
- `references/exports.txt`: exported symbols
- `references/help-index.txt`: help index dump
- `references/vignettes.tsv`: vignette names + titles (when present)
"""


def _default_openai_yaml(skill_name: str, pkg: str) -> str:
    display = f"R: {pkg}"
    short = f"Use the R package {pkg}"
    prompt = f"Use ${skill_name} to help me use {pkg} in R."
    return f"""interface:
  display_name: "{display}"
  short_description: "{short}"
  default_prompt: "{prompt}"
"""


def _parse_args(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Scaffold a Codex skill for an R package (r-<pkg>) with optional reference dumps.",
    )
    parser.add_argument("package", help="R package name (e.g., lme4)")
    parser.add_argument(
        "--out",
        required=True,
        help="Output directory to create the skill folder in (the skill name becomes a subfolder).",
    )
    parser.add_argument(
        "--skill-name",
        help="Override skill folder/name (default: r-<pkg>). Must be lowercase hyphen-case.",
    )
    parser.add_argument("--force", action="store_true", help="Overwrite existing skill folder.")
    parser.add_argument(
        "--no-references",
        action="store_true",
        help="Skip extracting references via Rscript (still scaffolds the skill).",
    )
    return parser.parse_args(argv)


def main(argv: list[str]) -> int:
    args = _parse_args(argv)

    pkg = args.package.strip()
    if not PKG_RE.match(pkg):
        _die(f"Invalid package name: {pkg!r}")

    skill_name = args.skill_name.strip() if args.skill_name else f"r-{pkg.lower()}"
    if not SKILL_RE.match(skill_name):
        _die(f"Invalid skill name (must be lowercase hyphen-case): {skill_name!r}")

    out_root = Path(args.out).expanduser().resolve()
    skill_dir = out_root / skill_name
    if skill_dir.exists():
        if not args.force:
            _die(f"Skill folder already exists: {skill_dir} (use --force to overwrite)")
        shutil.rmtree(skill_dir)

    (skill_dir / "agents").mkdir(parents=True, exist_ok=True)
    (skill_dir / "references").mkdir(parents=True, exist_ok=True)

    if not args.no_references:
        _extract_references(pkg=pkg, refs_dir=skill_dir / "references")

    title = _infer_title_from_pkg_desc(skill_dir / "references")
    _write_text(skill_dir / "SKILL.md", _default_skill_md(skill_name=skill_name, pkg=pkg, title=title))
    _write_text(skill_dir / "agents" / "openai.yaml", _default_openai_yaml(skill_name=skill_name, pkg=pkg))

    print(f"[OK] Created skill: {skill_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))

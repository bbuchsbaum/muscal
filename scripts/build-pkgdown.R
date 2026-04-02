#!/usr/bin/env Rscript

# pkgdown's default callr/processx build path is fragile in this repository's
# current local setup. Build the site in-process, including article renders.

cache_root <- file.path(tempdir(), "muscal-pkgdown-cache")
dir.create(file.path(cache_root, "xdg"), recursive = TRUE, showWarnings = FALSE)

Sys.setenv(
  R_USER_CACHE_DIR = cache_root,
  XDG_CACHE_HOME = file.path(cache_root, "xdg")
)

override <- list()
pkg <- pkgdown:::section_init(".", override = override)

pkgdown:::pkgdown_sitrep(pkg)
pkgdown:::init_site(pkg, override)
pkgdown:::build_home(pkg, override = override, quiet = TRUE, preview = FALSE)
pkgdown:::build_reference(
  pkg,
  lazy = FALSE,
  examples = TRUE,
  run_dont_run = FALSE,
  seed = 1014L,
  override = override,
  preview = FALSE,
  devel = FALSE
)

pkgdown:::build_articles_index(pkg)
for (name in pkg$vignettes$name[pkg$vignettes$type == "rmd"]) {
  pkgdown:::build_article(
    name,
    pkg = pkg,
    lazy = FALSE,
    seed = 1014L,
    new_process = FALSE,
    quiet = TRUE
  )
}
pkgdown:::build_quarto_articles(pkg, quiet = TRUE)
pkgdown:::build_tutorials(pkg, override = override, preview = FALSE)
pkgdown:::build_news(pkg, override = override, preview = FALSE)
pkgdown:::build_sitemap(pkg)

if (pkg$bs_version > 3) {
  pkgdown:::build_llm_docs(pkg)
}

pkgdown:::build_redirects(pkg, override = override)

if (pkg$bs_version == 3) {
  pkgdown:::build_docsearch_json(pkg)
} else {
  pkgdown:::build_search(pkg, override = override)
}

pkgdown:::check_built_site(pkg)

# Contributing

Thanks for looking. `tidymarkers` is an experimental package, but it is built and
gated like a normal R package: the merge gate below is what decides whether a
change lands.

## Quick development loop

```r
# environment: Imports + dev/test tooling, without the heavy Suggests
# (faster than renv::restore(): no historical 262-package dev lock)
# Rscript dev/install-deps.R                  # dependencies + the package
# Rscript dev/install-deps.R --with-suggests  # full set, for R CMD check

# load the package from source, without installing
pkgload::load_all()

# fast tests; --full adds the external cross-checks (BLAST/Rsubread/annotations/SOCK)
# Rscript dev/test.R --full
testthat::test_local()

# documentation from roxygen -> NAMESPACE + man/
roxygen2::roxygenise()

# full check
#   Rscript dev/check.R            # the repository gate (see merge gate below)
#   R CMD build . && R CMD check --no-manual tidymarkers_*.tar.gz
```

## Merge gate

`Rscript dev/check.R` must pass: docs and artifacts, `styler`, `lintr`, roxygen
sync, figures and the install smoke test. `--full` adds
`R CMD check --as-cran`, which runs the complete test suite. `--quick` is the
routine loop (roxygen, style, lint, criterion fixtures, docs/artifacts).

CI (`.github/workflows/R-CMD-check.yaml`) runs `R CMD check --as-cran` on push
and pull request, with OligoArrayAux installed from conda. `dev/check.R` is the
local equivalent and the one that also checks figures and artifact hashes, which
CI does not render.

## Style

- **tidyverse style**: `|>` (base pipe), backticks, `snake_case`, `snake_case()`
  for verbs; no `*_dfr()` — use `purrr::list_rbind()`.
- Inside `dplyr`/`ggplot2` use the `.data$col` pronoun (not bare names).
- Comments explain **why**, not what.
- Formatting: `styler::style_pkg()`; lint: `lintr::lint_package()` (config in
  `.lintr`).

## Documentation

- Every exported function: `@param` for all arguments, `@return`, `@family` (the
  group drives "See also" and pkgdown), `@examples` (heavy ones in `\dontrun{}`).
- New behaviour → a new test in `tests/testthat/`.
- A meaningful change → an entry in `NEWS.md` and in `ROADMAP.md`.

All public-facing prose is English; `ROADMAP.md`, `dev/README.md` and
`docs/archive/inventory-ru.md` are internal working documents and stay Russian.

## What lives where

Only the files R needs actually ship. Everything below is development material
and is excluded from the tarball by `.Rbuildignore`:

| Path | What it is |
|---|---|
| `R/`, `man/`, `vignettes/`, `tests/` | the package |
| `data-raw/` | generators for `inst/extdata/` and the source reads behind the rules (e.g. `polymarker-source.R`) |
| `dev/` | the gate (`check.R`, `test.R`), renderers, benchmarks, `dev/lecture/` slides, and the vendored ELK/libavoid layout engine in `.tools/elk` |
| `docs/` | `INDEX.md`, `artifacts.yml` (figure manifest), generated artifact graph, `docs/archive/` |
| `.agents/`, `.cursor/` | agent instructions and repo skills; not part of the product |
| `VALIDATION.md` | the evidence accounting behind the design rules |

## Code review

No independent human code review is recorded in this repository. What actually
runs on every change is the gate above plus the test suite, including the
falsifiability test (`test-falsifiable.R`), which fails if a published acceptance
criterion has no fixture that can drive it to `fail`. Treat that as a floor, not
as a substitute for review.

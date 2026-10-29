# Contributing

Thanks for looking. `tidymarkers` designs, checks and audits plant SNP markers.
It is experimental, but it is built and gated like a normal R package: the merge
gate below is what decides whether a change lands.

## What a contribution is worth here

This repository does not value volume. It values a claim that can be checked.
Before you write code, know which of four kinds of change you are making,
because each has a different gate:

| Kind | Example | What proves it |
|---|---|---|
| **Rule** — a design or acceptance rule | a new mismatch severity, a changed Tm window | a measurement read from a source you opened, and a `data-raw/<paper>.R` script that re-derives it and fails on disagreement |
| **Fixture** — an input the suite runs against | a panel, a haplotypes table, a reference slice | the file plus the id in `dev/check-criteria.R` that names a real test driving it |
| **Wiring** — existing code reaches a caller | an exported function that nothing calls | a call site, or the function is unexported and its `@keywords internal` says why |
| **Prose** — README, vignette, roadmap | a corrected sentence | it matches what the code and the tests now do |

A rule with no measurement is a default. Say so in the roxygen comment and in
`NEWS.md`; a default written down honestly is welcome, a default dressed as a
finding is not.

## Scope

The package answers one question: given a target variant and a genome, which
marker assay should be designed for it, and how far can you trust the answer?
Everything in `R/` serves that path — marker design (KASP, TaqMan, LAMP, RPA,
HRM, CRISPR), specificity against a genome, decision ladders, and the audit
trail behind a verdict. A change that does not serve that path belongs in a
downstream package or in the dossier, not here.

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

Two gates are specific to this project and will reject a change that looks fine:

- **`dev/check-criteria.R`** — every criterion must be able to fail. A new
  `gate` or `check` id needs either a fixture naming a real test that drives it
  to `fail`, or an exemption with a non-empty reason (a genome, a BLAST
  database or a population the in-memory suite does not carry).
- **`tests/testthat/test-falsifiable.R`** — a published acceptance criterion
  with no fixture that can drive it to `fail` fails the suite.

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

All public-facing prose is English, except the changelog: `NEWS.md` is Russian
by the repository convention (see `AGENTS.md`). `ROADMAP.md`, `dev/README.md`
and `docs/archive/inventory-ru.md` are internal working documents and stay
Russian.

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

## Reviewer protocol

A change under review is asked exactly two questions, in this order. A proposal
that anyone can check against the repository is treated as a contribution
regardless of who raised it; the credits below record who raised what.

1. **Is the claim measured?** Point at the file, the script in `data-raw/` and
   the test that drives the criterion to `fail`. "The paper says so" starts the
   answer, it does not finish it.
2. **Does it move the marker path?** If it does, which assay, which decision
   rung, and what does the caller do differently?

Review is recorded in the pull request or in `ROADMAP.md`, with the section id
the change closes.

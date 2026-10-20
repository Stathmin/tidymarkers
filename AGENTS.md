# Agent instructions — tidymarkers

Project-specific rules for agents working in this repository. Read these before
changing anything, then read `README.md`, `CONTRIBUTING.md`, `ROADMAP.md` (the open
queue) and `INVENTORY.md` (the archived inventory of what already exists).

## Running anything (read first — the environment is NOT the default library)

The lockfile targets **R 4.6.1 / Bioconductor 3.23**. On this machine that R is
`/opt/R-4.6.1`; the default `/usr/bin/R` (4.5) does **not** match the lockfile and
cannot load the renv library. No `.Rprofile`/`renv/` is committed on purpose
(activation would make `R CMD build` run against an empty library; see README
§Воспроизводимость). Run `Rscript` from the project root; activate renv manually
only when you need the pinned snapshot. The external tools (OligoArrayAux) live in
`.tools/` and are found by the package at run time.

Always prefix commands with:

```sh
export PATH=/opt/R-4.6.1/bin:/root/development/tidymarkers/.tools/bin:$PATH
export UNAFOLDDAT=/root/development/tidymarkers/.tools/share/oligoarrayaux
```

Then the usual commands:

```sh
# load the package from source
Rscript -e 'pkgload::load_all(".")'
# fast tests (core contracts; external cross-checks skip)
Rscript dev/test.R
# the same, one file/scheme, and the full external suite
Rscript dev/test.R crispr
Rscript dev/test.R --full          # Opt-in; R CMD check always runs it
# the internal gate (styler, lintr, roxygen sync, figures, install smoke)
Rscript dev/check.R              # add --full for R CMD check, --no-render to skip figures
```

If `pkgload::load_all()` cannot find a dependency, run `Rscript dev/install-deps.R`.
It installs the hard `Imports`, the small dev/test tools and the package itself; add
`--with-suggests` only when you need the heavy annotation/BLAST/Rsubread cross-checks.
Tests never install anything: `dev/test.R` stops if the package is not installed.
`renv::restore()` also works but restores the historical 262-package snapshot and is
much slower. Do **not** conclude the code or deps are missing.

## Merge gate

`Rscript dev/check.R` must pass (with the env above). It is the local stand-in for
CI, which this repo does not ship. `--full` adds `R CMD check --as-cran`.

## Conventions

- tidyverse style (`|>`, `.data$`, `snake_case`); `styler` + `lintr` clean.
- New behaviour → a test in `tests/testthat/`; significant change → `NEWS.md` +
  `ROADMAP.md`.
- Keep the tree minimal; prefer editing existing files.

## Current work

The decision DAGs are being rebuilt by a 10-agent effort; the plan, the invariants
(C1–C14) and each agent's scope/brief are in `dev/dag-rebuild-plan.md`.

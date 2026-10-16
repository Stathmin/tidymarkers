# Agent instructions — tidykasp

Project-specific rules for agents working in this repository. Read these before
changing anything, then read `README.md`, `CONTRIBUTING.md`, `ROADMAP.md` (the open
queue) and `INVENTORY.md` (the archived inventory of what already exists).

## Running anything (read first — the environment is NOT the default library)

The package dependencies live in a **renv cache library**, not the default R
library. Without `R_LIBS` set, `pkgload::load_all()` fails at import
(`Biostrings`, `DECIPHER`, `TmCalculator`, `rBLAST`, `tidygraph`, … are "required"
but not found) and **nothing can be run or tested**. The project also keeps the
external tools (OligoArrayAux) in `.tools/`.

Always prefix commands with:

```sh
export R_LIBS=/home/stathmin/.cache/R/renv/library/kasp-example-cee55197/linux-debian-forky/R-4.6/x86_64-pc-linux-gnu
export PATH=/home/stathmin/Documents/devel/tidykasp/.tools/bin:$PATH
export UNAFOLDDAT=/home/stathmin/Documents/devel/tidykasp/.tools/share/oligoarrayaux
```

Then the usual commands:

```sh
# load the package from source
Rscript -e 'pkgload::load_all(".")'
# tests (one file per scheme, so filter by scheme name)
Rscript -e 'testthat::test_local(filter = "crispr")'
# the internal gate (styler, lintr, roxygen sync, figures, install smoke)
Rscript dev/check.R              # add --full for R CMD check, --no-render to skip figures
```

If these are not set, do **not** conclude the code or deps are missing — set them.

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

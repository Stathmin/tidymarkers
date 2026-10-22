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
export R_LIBS_SITE=/opt/R-4.6.1/lib/R/library
```

`R_LIBS_USER` is not the project library, and the `.tmp-rt/lib` path in old handoffs
is dead. `R_LIBS_SITE` puts the R 4.6.1 site library after the two system ones; check
it with `vapply(c("GenomicRanges", "Rsamtools", "TmCalculator", "tidygraph",
"DECIPHER"), requireNamespace, logical(1), quietly = TRUE)` → all `TRUE`.
Do **not** set `R_LIBS` to that path: R then puts it first and loads system-R's
`rlang.so`, and every package needing rlang fails with
`undefined symbol: R_MakeMissingBinding`. `ggforce` is the tell: it resolves and
`GenomicRanges` does not (or the reverse) means the library order is wrong.
`/opt/R-4.6.1/lib/R/library/ggforce` is a broken 0.5.0 build
(`undefined symbol: R_getRegisteredNamespace`); the working copy is 0.4.2 from
`/usr/lib/R/site-library`, and the package only uses `facet_zoom` and `geom_arc`.

### R behaviour on this build that has bitten real work

These are measured on R 4.6.1 here, not general R folklore. Each one silently
returned a wrong answer before it was found, so treat the symptoms as a signal
to check the idiom, not to lengthen the code.

- **`x %in% code` with a multi-base code is element-wise, not a set test.**
  `"N" %in% c("A", "C", "G", "T")` is a length-4 logical; `isTRUE()` or a
  one-element `vapply` drops it and the answer becomes `FALSE`. Use
  `!is.na(match(x, code))`.
- **`which()` treats `NA` as unset**, and `x[NA]` reads as out of range. A
  logical with an `NA` therefore gives a short (or empty) index, not a
  skip-the-`NA`s filter. Resolve `NA` before the lookup.
- **`substr(x, from, stop)` is not vectorised over `x` here**: it truncates to the
  first element. Split instead of slicing.
- **`c()` on a matrix** unlists it, so `c(a_matrix, TRUE, FALSE)` counts the
  matrix as its number of columns, not its length. `as.vector()` first.
- **`x[cbind(i, j)]` with repeated `i` returns the cross product**, not pairwise
  lookups. Build a linear index.
- `pkgload::load_all()` caches a compiled/lazyloaded copy: if a function you
  edited still behaves like the old body, `recompile = TRUE`.

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

## When something is slow, profile before editing

A straight pipeline with nested helpers hides its cost. Profile it, do not guess:

```sh
Rscript -e 'Rprof("x.out", interval = 0.01); <the call>; Rprof(NULL)'
Rscript -e 'summaryRprof("x.out")'
```

- `summaryRprof()` aggregates by function name and **hides callers**. When the
  question is "who repeats this?", parse the raw `x.out` lines and read the
  innermost frame together with its parent; `rle()` the parent chain to find the
  repeated run.
- A profile dominated by `validObject` / `initialize` / `getClassDef` /
  `standardGeneric` / `getClassesFromCache` is **S4 object construction**, not
  the algorithm. The fix is to build fewer objects, not to rewrite the maths.
- A grid search multiplies every axis. Count the objects a call builds (strands x
  lengths x positions x bases) before optimising any one of them.
- Fix the cheap doomed path first: when the whole grid is built and then
discarded, answer the discard condition before the grid exists.
- RDF export is not an R debugging standard (it is an interchange format, and R's
  call graph knows nothing about it). `Rprof` is the standard here.

## Reading the literature (never from memory)

The scientific sources behind `inst/extdata/scientific_sources.csv` live in the
operator's Zotero library, collection `tidymarkers — marker design sources`
(`4RSE6KC2`, tag `merovingen:tidymarkers`). 30 of the 72 records carry a full-text
attachment; 4 rows are `DOI-UNRESOLVED`. Read the full text through Merovingen's
read-only MCP lane before you rely on a claim:

```sh
cd /root/development/Merovingen
ZOTERO_MCP_ENABLED=true ZOTERO_MCP_MODE=live .venv/bin/python \
  data/scratch/reports/scripts/zotero_read_items.py <ITEM_KEY>   # -> /tmp/read/<KEY>.md
ZOTERO_MCP_ENABLED=true ZOTERO_MCP_MODE=live .venv/bin/python \
  data/scratch/reports/scripts/zotero_quote.py --grep "phrase"  # phrase across the library
```

From Python, the tools are `zotero_get_collection_items`, `zotero_get_item_metadata`,
`zotero_get_item_fulltext`, `zotero_get_pdf_outline` and `zotero_read_pdf_pages` in
`src.tools.readonly.zotero_bridge` (allowlist `zotero/allowlist.v1.yaml`).
Writes to the library are an operator action (ADR-0037) — an agent cites, it does
not add records. Check the `claim` column of `scientific_sources.csv` against the
sentence in the fulltext, not against the title: a title can neither confirm nor
refute a number. A source without an attachment is cited only for what is visible
in its abstract, and the gap is noted in `access_detail`.

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

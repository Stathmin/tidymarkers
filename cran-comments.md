## R CMD check results

0 errors | 0 warnings | 1 note

- This is a new submission.

The note is the standard `New submission` note. No other note is expected: the
scoped check of the built tarball (`--no-manual --as-cran`, `_R_CHECK_FORCE_SUGGESTS_=false`) returned `Status: OK`.

## Test environments

- local: Ubuntu 24.04, R 4.6.1 (2026-06-24), 64-bit
- GitHub Actions (`R-CMD-check`), each leg `--no-manual --as-cran` with
  `error-on: warning`:
    - Ubuntu latest, R release / oldrel-1 / devel
    - macOS latest, R release
    - Windows latest, R release

## Optional external engines

`SystemRequirements` names two optional command-line engines and one optional
runtime:

- **OligoArrayAux** (`hybrid-min`, `hybrid-ss-min`) — nearest-neighbour
  thermodynamics for hairpins, dimers and amplicon melting.
- **NCBI BLAST+** — one of the four measurement engines reported in the
  provenance stamp.
- **Node.js** — layout of the vendored decision-graph engine.

None is a build or check dependency. Every call is guarded, and the test suite
skips the assertion that needs an absent engine rather than failing, so a check
machine without them still runs the whole suite. The CI workflow installs
OligoArrayAux and BLAST+ on Linux and BLAST+ on macOS; the Windows leg runs
without them by design.

`rprimer`, `Rsubread`, `rBLAST`, `VariantAnnotation`, `BSgenome`,
`TxDb.Athaliana.BioMart.plantsmart28` and `crisprVerse` are in `Suggests` and
are used only for cross-checks against a second implementation. Their absence
removes one cross-check, not a capability.

## Downstream dependencies

There are none: this is a new submission and the package is not on CRAN. See
`revdep/README.md`.

## Package size

The tarball is 2.3 MB; the installed size is 10.3 MB, of which `inst/extdata`
is 4.8 MB. That data is the evidence base: curated tables transcribed from the
published sources named in `inst/extdata/scientific_sources.csv`, each one
readable and checkable against its source. It is data, not a build artifact, so
it cannot be generated at install time. No file is over 1.5 MB.

## Spelling

The package passes the package-wide spell check (`spelling::spell_check_files()`
over `man/`, `README.md` and `vignettes/`, gated against `inst/WORDLIST`).
`NEWS.md` and `ROADMAP.md` are written in Russian for the maintainer's work log;
they are not part of the English documentation and are not spell-checked.

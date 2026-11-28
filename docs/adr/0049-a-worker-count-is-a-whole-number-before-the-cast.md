# ADR 0049 — A worker count is a positive whole number before the cast

- Status: accepted
- Date: 2026-11-28
- Builds on: ADR 0044 (one whole-number predicate), ADR 0047 (a search radius
  is checked before the cast), ADR 0048 (a base-pair threshold is checked
  before the cast)

## Context

Four public entry points take a worker count that ends up in an external tool:

| path | before the fix |
|---|---|
| `genome_index(threads=)` | `as.integer(threads)` then stored in `index$threads` and passed to every index build |
| `marker_bowtie_index(threads=)` | passed raw to `Rbowtie::bowtie_build()` |
| `specificity_bowtie(nthreads=)` | passed raw to `Rbowtie::bowtie(p=)` |
| `specificity_fm(nthreads=)` | passed raw to `Rsubread::align()` |

Measured on R 4.6.1: `as.integer(0.5)` is `0` and `as.integer(Inf)` is `NA`
with a warning. A zero-worker build can stall or fall back to an engine
default the caller did not ask for; `NA` reaches `bowtie-build --threads` as
the string `NA` and fails inside the tool with a message that never names the
argument. The two specificity paths never cast at all, so the fraction reaches
the engine and its error, if any, is the engine's.

This is the same defect class as ADR 0047 and ADR 0048: a caller-supplied
number crosses the package boundary before the package has said what a valid
number is.

## Decision

A worker count is one positive whole number. `check_worker_count(x, arg)` in
`R/marker_validate.R` runs `is_whole_number(x, min = 1L)` on the value the
caller gave, aborts with
`` `{arg}` must be one positive whole number of workers. `` otherwise, and
returns `as.integer(x)`. All four entry points call it before any genome,
index or engine work, and the returned integer is what flows onwards.

`is_whole_number()` runs first, so `NA`, `NaN`, `Inf`, a character, a vector
and `0` all reach the package message instead of base R's condition error or
an engine error.

## Alternatives rejected

- **Clamp the value with `max(1L, ...)`.** A caller who asked for half a
  worker has a bug; silently running one hides it.
- **Leave `specificity_bowtie()` / `specificity_fm()` to the engine.** Their
  error names the engine, not the argument, and it only appears when the
  engine is installed, so a test cannot pin the contract offline.
- **One gate only in `genome_index()`.** `marker_bowtie_index()`,
  `specificity_bowtie()` and `specificity_fm()` are exported and are called
  directly, so each one is a trust boundary of its own.

## Consequences

- A fractional, zero, negative, infinite or vector worker count is refused
  with one package message before any file or engine is touched.
- The accepted values behave exactly as before: a whole number `>= 1` is cast
  once and used everywhere.
- `check_worker_count()` is internal; no exported signature changes.
- The rule is narrower than "any count": it says workers, so the message can
  name the unit the caller chose.

## Verification

`tests/testthat/test-genome-index.R` covers `genome_index(threads=)` (bad
values refused, `4L` stored as `4L`); `test-specificity.R` covers
`specificity_bowtie()` and `specificity_fm()`; `test-bowtie.R` covers
`marker_bowtie_index()`. The focused run
`Rscript dev/test.R "genome-index|^specificity$|^bowtie$|^fm$"` finished with
`FAIL 0 | WARN 0 | SKIP 10 | PASS 427`.

## Amendment (2026-11-28) — three more `cores` paths reach an engine

The first sweep named the four entry points above. A second read of `cores=`
found three more exported paths where the count crossed into a backend before
any check:

| path | before the amendment |
|---|---|
| `pcr_efficiency(cores=)` | `processors = max(1L, as.integer(cores))` |
| `assemble_kasp_assay(cores=)` and `assemble_primer_assay(cores=)` | `marker_bpparam()` did `n <- max(1L, as.integer(cores))` |
| `render_design_search(cores=)` | `ncores <- min(as.integer(cores), length(uniq))` |

Measured on R 4.6.1: `cores = 1.9` ran one worker, `Inf` became `NA` and
reached DECIPHER / the plot pool as a missing value, and `c(1, 2)` was
silently resolved by `max()` or `min()` instead of being refused.

`marker_bpparam()` is the single place `cores` becomes a `BPPARAM`, so the
check sits there rather than in each constructor; `pcr_efficiency()` and
`render_design_search()` each call `check_worker_count(cores, "cores")` on
their own value before any engine, plot or frame work. The accepted values
behave exactly as before, and no exported signature changes.

Verification: `test-marker-assemble.R` refuses `1.5`, `0`, `Inf`, `c(1, 2)`
and `"2"` at both constructors; `test-thermo.R` at `pcr_efficiency()` and
keeps `cores = 1L` reaching the engine check; `test-search-anim.R` at
`render_design_search()` before the first frame. `dev/test.R
"marker-assemble|^thermo$|search-anim"` — 83 pass, 0 fail.

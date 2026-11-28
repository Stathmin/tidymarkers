# 0006 — The panel pipeline is an orchestrator with falsifiable seams

- Status: accepted
- Date: 2026-11-26
- Supersedes: none

## Context

`ZY1` landed on 2026-11-19: `R/marker_panel.R` (280 lines) composes the bricks
that already exist — `marker_plan()`, `marker_solve()`, `marker_one_pot()` and
`marker_artefact()` — behind one call, `marker_design_panel()`. The committed
fixture is `tests/testthat/fixtures/panel/{genome.fa,haplotypes.csv}`, and
`tests/testthat/test-panel.R` (132 lines) already checks that every target ends
in one row: a design, or the cause chain of every chemistry tried.

The composition is correct but its seams are implicit:

- **GWAS loci.** `targets` and `region` both feed `marker_plan()`, but duplicate
  targets are not rejected; two rows with the same `(chrom, pos, ref, alt)` can
  produce two report rows, and `targets.csv` does not say which one is the
  duplicate.
- **LD.** There is no LD step. The criterion `"design relies on linkage"` in
  `R/marker_schemes.R` is a per-assay boolean about ALT-primer coverage, not a
  population LD prune. A GWAS input with two perfectly linked signals designs
  both.
- **Phased haplotypes.** `normalise_haplotypes()` checks the six columns and
  that `freq` is one value in `[0, 1]` per haplotype. It does not check that a
  haplotype has at most one allele at a position, that a row's `ref`/`alt`
  matches the target, or that the frequencies sum to at most one per chromosome.
  Absence is treated as reference by the design code, but the panel neither
  states that convention nor counts how often it was used.
- **Marker selection.** `marker_solve()` writes one `report` row per target and
  one `attempts` row per (target, scheme); the panel copies the report and folds
  the attempts into `cause`. The contract is real but it is asserted in one
  test, not stated as a seam.
- **Chemistry choice.** The solver walks the allowed schemes cheapest-first and
  stops at the first success. The panel has no `costs=` argument, so a caller
  cannot exercise the ordering; the test cannot show that the cheapest success
  won, only that a success happened.
- **Assembly and checks.** `assemble_kasp_assay()` and `marker_design_checks()`
  are called through `marker_solve()`; a half-built assay or an unknown
  criterion id would surface as a downstream `NA`, not at the seam.
- **Artefact.** `marker_artefact()` writes `artefact.json`, `decisions.csv`,
  `checks.csv`, `why.csv`, `targets.csv` and `solved.rds`. `artefact_id` already
  excludes `created` and is stable, but `created` itself is `Sys.time()`, so a
  rerun cannot be byte-identical. `solved.rds` carries an igraph pointer and a
  random id (the existing test says so); it cannot be compared byte-for-byte.

The blocker named in `ZY1` was "no committed phased-haplotype fixture". That
blocker is removed: the fixture is committed and hand-written, and the test
fails if a target is dropped. It was not removed by a stub. The remaining work
is to make each seam falsifiable and to make the artefact reproducible.

## Decision

### 1. One orchestrator, eight named seams

`marker_design_panel()` stays the only public entry point. It does not grow a
second pipeline object: `marker_plan` and `marker_solved` are already the
objects. Each seam is either an existing exported brick or one new internal
function with a typed input and output.

| Seam              | Input                                  | Output                         | Falsifiable criterion                                                                 |
| ----------------- | -------------------------------------- | ------------------------------ | ------------------------------------------------------------------------------------- |
| GWAS loci         | `targets` table or `region`            | `plan$targets`                 | every distinct target appears once; a duplicate aborts with `tidymarkers_duplicate_target` |
| LD                | `ld = NULL`, a table, or a threshold   | filtered target table          | with `ld = 0.8` a perfect-LD pair leaves one target; with `ld = NULL` both stay and the manifest records `ld = "not_applied"` |
| Phased haplotypes | `haplotypes` table or CSV path         | normalised haplotypes          | no duplicate `(haplotype, chrom, pos)`; a row at a target position matches the target's `ref`/`alt` (other rows are checked against the variant table when supplied); absence is reference and the count is recorded; `sum(freq) <= 1` per chromosome; a conflict aborts naming the haplotype and position |
| Selection         | plan                                   | `marker_solved$report`         | one row per target; `ok` is TRUE iff a scheme designed; `attempts` has one row per (target, scheme) tried |
| Chemistry         | plan costs, allowed schemes            | chosen scheme per target       | the chosen scheme has the minimum cost among successful attempts; reversing the costs changes the choice |
| Assembly          | target + genome                        | `marker_assay` or `marker_doomed` | a designed assay passes `contract_assay_output()`; a doomed result carries a non-empty reason |
| Checks            | designed assays                        | checks table                   | every `criterion` id is in `criterion_spec()$id`; every designed assay has a row    |
| Artefact          | solved plan, directory, fixed time     | files on disk                  | with a fixed `SOURCE_DATE_EPOCH`, every text file is byte-identical on rerun; `solved.rds` is excluded by name and reason |

The LD seam has no default threshold. A default would be an uncalibrated number
in a package whose criteria registry records where every bar comes from. `NULL`
means "not applied", and that fact is written into the artefact. A caller who
passes a table or a threshold gets the prune and the record.

### 2. New arguments, not new wrappers

`marker_design_panel()` gains:

- `ld = NULL` — `NULL`, a linkage table, or a numeric `r2` threshold computed
  from `haplotypes`;
- `costs = NULL` — passed to `marker_tech()`, so chemistry ordering is a seam;
- `created = NULL` — passed to `marker_artefact()`, defaulting to
  `Sys.getenv("SOURCE_DATE_EPOCH")` and then to the current time.

The new internals are `panel_ld()`, `panel_haplotype_contract()` and
`panel_duplicate_guard()`. `marker_plan()`, `marker_solve()`,
`marker_one_pot()` and `marker_artefact()` are reused as they are; the panel
does not wrap them twice.

### 3. Reproducibility is part of the contract

`marker_artefact(created = )` writes the `created` field from an explicit value.
When `SOURCE_DATE_EPOCH` is set, that value is used, so the same input and the
same epoch produce the same text files. The gate compares `artefact.json`,
`decisions.csv`, `checks.csv`, `why.csv` and `targets.csv` byte-for-byte across
two directories. `solved.rds` is excluded, with the reason written into the
test: it serialises an igraph pointer and a random id, and a byte comparison
would test the serializer, not the design.

### 4. The fixture carries the blocker

`tests/testthat/fixtures/panel/` keeps `genome.fa` and `haplotypes.csv` (three
rows: H1 twice and H2 once). It gains:

- `haplotypes_ld.csv` — a pair of targets in perfect LD, to prove the LD seam
  collapses them when a threshold is given and leaves both when it is not;
- `haplotypes_unphased.csv` — a haplotype with two alleles at one position, to
  prove the phasing contract aborts rather than picking one or silently
  imputing a reference.

The tests use `skip_if_no_oligoarrayaux()` for the design path, as the existing
file does, but the fixture assertions themselves are not skipped: they do not
need the engine. No `panel_ld <- function(x) x` stub is added. A stub would
make the test green without making the seam real, which is the failure mode
`ZY1` was closed to avoid.

## Gates

1. **Every target accounted for.** `test-panel.R` keeps the one-row-per-target
   assertion; `test-panel-artefact.R` reads `targets.csv` and asserts the same
   count and labels as `panel$targets`.
2. **LD.** With a perfect-LD pair and `ld = 0.8`, one target remains; with
   `ld = NULL`, both remain and `artefact.json` says `ld: "not_applied"`.
3. **Phasing.** An unphased (two-allele) row aborts with the haplotype and
   position in the message; a frequency sum above one aborts per chromosome; an
   absent target allele is counted as reference, not silently invented.
4. **Chemistry.** With reversed costs, the chosen scheme changes and is the
   minimum-cost successful attempt; with no successful chemistry, the abort
   lists every target and every cause.
5. **Assembly and checks.** Every designed assay passes the existing contract,
   and every check id is a member of the registry.
6. **Artefact.** Two runs with the same `SOURCE_DATE_EPOCH` are byte-identical
   for the text files; `solved.rds` is excluded by name in the test. A second
   run into the same directory does not duplicate rows or leave stale files.

## Consequences

- The panel becomes auditable at each seam, not only at the final verdict.
- LD is opt-in and recorded; no uncalibrated threshold enters the defaults.
- `marker_artefact()` gains an explicit `created`, and `SOURCE_DATE_EPOCH` is
  honoured. This changes the manifest only when the variable is set.
- Duplicate targets and unphased haplotypes become errors. This is a breaking
  change and goes to `NEWS.md`.
- `solved.rds` stays non-byte-stable. The contract says so; a future replay
  format would be a separate decision.
- The fixture grows by two small CSVs, not by a mock engine.

## Alternatives

- **Leave the wrapper and document the seams in prose.** Rejected: the existing
  test already shows a target can be designed or explained, but nothing checks
  LD, phasing, chemistry order or the artefact bytes.
- **Default the LD threshold to 0.8.** Rejected: there is no calibration in the
  registry for a population LD bar, and a default would hide that fact.
- **Add a `marker_pipeline` class.** Rejected: `marker_plan` plus `marker_solved`
  is the object graph; a third class would duplicate it.
- **Add a `panel_ld()` stub so the test can be written now.** Rejected: the
  fixture is the blocker, and the fixture is committed. A stub would pass a test
  that proves nothing.
- **Compare `solved.rds` byte-for-byte.** Rejected: it would test serialisation
  of an igraph pointer, not the design. The text artefact is the reproducible
  product.

## Revisit trigger

A real GWAS panel with more than 100 targets, or a calibrated population LD
threshold from a named source. At that point the LD default, the per-target
chunking and the `solved.rds` replay format become worth a new decision.

## References

- [`R/marker_panel.R`](../../R/marker_panel.R) — the orchestrator.
- [`R/marker_plan.R`](../../R/marker_plan.R) — targets, context and tech layers.
- [`R/marker_solve.R`](../../R/marker_solve.R) — cheapest-first attempts.
- [`R/marker_haplotypes.R`](../../R/marker_haplotypes.R) — the current normaliser.
- [`R/marker_artefact.R`](../../R/marker_artefact.R) — the manifest and files.
- [`tests/testthat/test-panel.R`](../../tests/testthat/test-panel.R) — the existing gate.
- [`tests/testthat/fixtures/panel/`](../../tests/testthat/fixtures/panel) — the fixture.
- [`ROADMAP.md`](../../ROADMAP.md) — ZY1, ZY2, ZT T4-T5.

## Evidence at acceptance (2026-11-28)

`Rscript dev/test.R panel` ran the panel, panel-artefact and panel-LD suites:
376 assertions, 0 failures, 0 skips
(`dev/logs/test-filter-panel-20261128-065434.log`). That covers gates 1–6 as
written: target accounting, LD applied and not applied, phasing aborts, scheme
choice under reversed costs, contract membership of every check id, and the
byte-identical text set under a fixed `SOURCE_DATE_EPOCH`.

# 0009 — The panel orchestrator: GWAS, LD, haplotypes, chemistry, artefact

- Status: proposed
- Date: 2026-11-27
- Supersedes: none
- Builds on: [ADR 0006](0006-panel-pipeline.md) (the seams)

## Context

ADR 0006 made `marker_design_panel()` the one public panel entry point and
named eight seams. `ZT T4` and `ZT T5` then landed the phased-haplotype
fixture, `panel_ld()`, `panel_duplicate_guard()`, `panel_haplotype_contract()`,
`costs =`, `created =` and a byte-identical text artefact. The panel is no
longer a wrapper without checks.

The 2026-11-27 inventory still finds the object terminal: **`panel` has one
producer and zero typed consumers** (`docs/user-api-surface.md` §2, defect 1).
The only way to use the result is `panel$solved`, `panel$artefact` or
`panel$one_pot`. That is the same failure ZU1 names.

Three seams are also weaker than the ADR 0006 table claims:

- **GWAS input is a variants table, not a signal table.** `marker_targets()`
  accepts `data =` and passes extra columns through, but no contract says that
  `label`, `p`, `beta`, `se` and `maf` survive into `targets.csv`. A caller
  with summary statistics has to guess which columns are legal.
- **LD is applied after the full haplotype contract.** `panel_haplotype_contract()`
  validates every haplotype row, then `panel_ld()` prunes. An unphased row at a
  target that LD would drop aborts the panel before the prune can happen.
- **The LD record is not uniform.** A pair table records `list(table = ld)` and
  does not name the dropped targets; a numeric threshold records
  `list(threshold = ld, dropped = ...)`. The artefact therefore cannot say
  which signals were removed and why in one shape.

The current order is duplicate guard → haplotype contract → LD → plan → solve.
The declared scenario J3 order is GWAS → LD → haplotypes → chemistry →
artefact. This ADR fixes the order and the record, not the bricks.

## Decision

### 1. The pipeline is one ordered contract

`marker_design_panel()` keeps its signature and becomes the executor of five
named stages. Each stage either returns the object the next stage reads or
aborts with a classed condition.

| Stage | Input | Output | Contract |
| ----- | ----- | ------ | -------- |
| **S1 GWAS** | `targets` table or `region`; optional `variants` | `plan$targets` | one row per distinct `(chrom, pos, ref, alt)`; `label` unique; a duplicate aborts with both labels and both rows; optional `signal_id`, `p`, `beta`, `se`, `maf` survive to `targets.csv`; no p-value threshold is invented |
| **S2 LD** | target table; `ld = NULL`, a pair table, or an `r2` threshold; haplotypes as the dosage source | pruned target table + `plan$ld` record | `NULL` = `not_applied`; a pair table names both targets and their labels; a threshold computes r2 from haplotype dosage at the target positions only; every dropped target is named with `reason = "ld"`; a table naming an absent target aborts |
| **S3 Haplotypes** | surviving targets; haplotypes; optional variants | normalised haplotypes + counts | the full ADR 0006 contract applies to the surviving targets: no duplicate `(haplotype, chrom, pos)`, target `ref`/`alt` match, `sum(freq) <= 1` per chromosome, absence = reference and is counted; rows outside the variants are counted as `outside_variants`, not silently accepted |
| **S4 Chemistry** | plan + `costs` + allowed schemes | chosen scheme per target | cheapest successful scheme wins; reversing `costs` changes the choice; a target with no success carries the full cause chain |
| **S5 Artefact** | solved plan + one-pot + panel record | `artefact.json` + CSVs | fixed `created` / `SOURCE_DATE_EPOCH` gives byte-identical text files; `solved.rds` is excluded by name; the panel block records all five stages |

S2 needs a dosage matrix, so it reads haplotypes at the target positions. That
is a local phase check, not the full S3 contract: one allele per position for
the targets being tested. S3 then validates the surviving set in full. The
practical change is that an unphased row at an LD-dropped target no longer
aborts the panel, and the artefact says which target was dropped.

### 2. `marker_panel` is an intermediate object

`as_marker_solved(x)` is the one bridge:

```r
as_marker_solved(panel)   # marker_solved, with attr(x, "panel") = the record
```

`marker_pareto()`, `marker_artefact()`, `marker_export()` and `marker_replay()`
accept a `marker_panel` through it. They do not grow a second argument and do
not copy the panel logic. `marker_panel$solved` stays for compatibility, but no
scenario uses it.

`marker_panel$experimental` stays `TRUE` until J3 is green with a committed
fixture and byte-identical artefact. Only then does it become `FALSE`; that is
the graduation condition, not a version number.

### 3. The artefact records the pipeline, not only the result

`artefact.json` gains a `panel` block, written by the panel, never by the
solver:

```json
"panel": {
  "order": ["gwas", "ld", "haplotypes", "chemistry", "artefact"],
  "gwas": {"n_signals": 12, "columns": ["label", "p", "beta"]},
  "ld": {"method": "threshold", "threshold": 0.8,
          "dropped": ["sig-3", "sig-7"], "reason": "ld"},
  "haplotypes": {"n": 48, "reference_calls": 6, "outside_variants": 0},
  "chemistry": {"costs": ["snp", "hrm"], "chosen": {"sig-1": "snp"}},
  "one_pot": {"groups": 3, "unassigned": 1}
}
```

The block is part of the artefact identity because LD pruning changes the
target set. A run with `ld = NULL` and a run with `ld = 0.8` must not share an
`artefact_id`. This is compatible with the current `artefact_id()` core
(`inputs`, `plan`, `parameters`, `overrides`): the panel block enters `plan`.

### 4. No new GWAS reader, no default threshold

Summary statistics enter through `marker_targets(data = )`; the required
columns are `chrom`, `pos`, `ref`, `alt`, and `label` if the caller wants
stable ids. `p`, `beta`, `se`, `maf` and `signal_id` are optional and recorded.
A p-value cut is the caller's decision and is written into the input name or
the `gwas` record; the package does not choose one.

The LD threshold has no default for the same reason as ADR 0006: there is no
calibrated population bar in the criterion registry. `NULL` means "not
applied" and is recorded as such.

## Gates

1. **J3.** `test-scenario-j3.R` runs the five stages from a committed fixture
   and asserts the class at each step: `targets` → `plan` → `marker_solved` →
   `marker_panel`, then `marker_pareto(panel)` and `marker_artefact(panel)`.
2. **Order.** An unphased row at a target that `ld = 0.8` drops does not abort;
   the same row on a surviving target aborts with the haplotype and position.
3. **LD record.** With a pair table, both target labels are present and the
   dropped label is named; with a numeric threshold, the same shape is
   produced; with `ld = NULL`, `method = "not_applied"` and no drop.
4. **Chemistry.** Reversing `costs` changes the chosen scheme on at least one
   target; the chosen scheme is the minimum-cost successful attempt.
5. **Artefact.** Two runs with the same `SOURCE_DATE_EPOCH` are byte-identical
   for `artefact.json`, `decisions.csv`, `checks.csv`, `why.csv`,
   `targets.csv` and `signals.csv`; the `panel` block is present and lists the
   dropped targets; `solved.rds` is excluded by name.
6. **Surface.** `Rscript dev/check-user-api.R` shows typed edges
   `panel → marker_pareto`, `panel → marker_artefact`, `panel → marker_export`
   and `report → marker_solved`; no scenario uses `$solved`.

## Consequences

- The panel becomes a typed intermediate object. `ZU1` closes, and the
  "panel → Pareto → patent" path stops requiring `$solved`.
- Stage order changes behaviour: an unphased row at an LD-dropped target is no
  longer an error. This is a user-visible change and goes to `NEWS.md`.
- `artefact.json` grows a `panel` block. The text files remain byte-identical
  under a fixed epoch; the schema string is not bumped by this ADR alone.
- `marker_panel$experimental` becomes `FALSE` only after J3 is green.
- The panel still calls `marker_plan()`, `marker_solve()`, `marker_one_pot()`
  and `marker_artefact()`; it does not wrap them twice.
- No GWAS reader and no default LD threshold enter the package in this cycle.

## Alternatives

- **Make `marker_design_panel()` return `marker_solved` directly.** Rejected:
  the one-pot grouping and the panel record are real outputs; flattening them
  loses the artefact.
- **Add `marker_panel()` as a separate class and a second entry point.**
  Rejected: `marker_plan + marker_solved` is already the object graph; a third
  class duplicates it.
- **Keep the haplotype contract before LD.** Rejected: it rejects a panel for a
  target that is about to be pruned, and the error names a row the user did not
  need.
- **Default `ld = 0.8`.** Rejected: no calibration source; the registry would
  have to invent one.
- **Read GWAS summary statistics in a new `read_gwas()`.** Rejected until a
  second scenario needs a format `read_variants_vcf()` cannot express. The
  current contract is a table with optional signal columns.

## Revisit trigger

A real GWAS panel with more than 100 targets, a calibrated LD threshold from a
named source, or a second summary-statistics format. At that point the input
adapter and the LD default become a new decision.

## References

- [`R/marker_panel.R`](../../R/marker_panel.R) — the orchestrator and the three
  current seams.
- [`R/marker_plan.R`](../../R/marker_plan.R) — targets, context and tech layers.
- [`R/marker_solve.R`](../../R/marker_solve.R) — cheapest-first attempts and
  `marker_export()`.
- [`R/marker_artefact.R`](../../R/marker_artefact.R) — the manifest and text
  files.
- [`tests/testthat/test-panel.R`](../../tests/testthat/test-panel.R) — the
  existing seam tests J3 extends.
- [`tests/testthat/fixtures/panel/`](../../tests/testthat/fixtures/panel) — the
  phased and unphased fixtures.
- [`docs/user-api-surface.md`](../user-api-surface.md) — defect 1 (`panel`) and
  the class table.
- [`ROADMAP.md`](../../ROADMAP.md) — §ZU (ZU1) and §ZV (ZV6).

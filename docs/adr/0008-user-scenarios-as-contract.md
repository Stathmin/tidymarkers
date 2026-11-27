# 0008 — User scenarios J1–J7 are the executable contract

- Status: accepted
- Date: 2026-11-27
- Supersedes: none

## Context

The 2026-11-27 inventory measured the public surface: **207 exports**, **48 S3
methods**, 255 CSV rows (137 user-facing, 118 service), **934 typed edges**, 58
entry points, 116 leaves and 29 isolates. The graph answers one question — does
the output class of `A` fit an input class of `B`? — and answers it honestly.
It does not answer the question a new user asks: *which of the 207 functions do
I call, in what order, to get a marker?*

The tests are per function. The README has two quick starts and the vignettes
are organised by chemistry. Nothing in the repository fails when a named user
goal stops being reachable. `docs/user-api-surface.md` §8 lists the consequence:
`panel`, `outcomes` and `report` are terminal objects, `primers` is two
contracts under one class, and the registries have no bridge into design. Those
are not independent defects; they are points where a scenario would have to
unpack an object by hand.

A scenario is therefore not documentation. It is the acceptance surface: a
committed fixture, an ordered call path, the class expected at every step, and
a test that runs offline. If a scenario needs `$solved`, `$report` or
`tibble$scheme`, the contract is missing and the scenario fails.

## Decision

### 1. Seven scenarios, one contract each

| id | Actor and goal | Input | Public path | Output | Runtime gate |
| -- | -------------- | ----- | ----------- | ------ | ------------ |
| **J1** | Breeder: one SNP to one ordered assay | committed FASTA + one variant row | `read_genome_fasta()` → `assemble_kasp_assay()` → `marker_design_checks()` → `marker_explain()` → `write_marker_results()` / `export_primers()` | assay, checks, order sheets | `test-scenario-j1.R` |
| **J2** | Breeder: a GWAS table to a batch report you can continue | VCF fixture + genome | `read_variants_vcf()` → `design_kasp_markers()` → `as_marker_solved(report)` → `marker_design_checks()` / `marker_export()` | `marker_solved`, per-locus verdicts | `test-scenario-j2.R` |
| **J3** | Panel designer: GWAS signals to a panel with LD, haplotypes, chemistry and an artefact | targets + haplotypes + variant table | `marker_design_panel()` → `marker_pareto(panel)` / `marker_artefact(panel)` / `marker_export(panel)` | `marker_panel`, artefact | `test-scenario-j3.R` |
| **J4** | Wheat programme: a pangenome, not a genome (> 21.5 GB) | three small build fixtures A/B/D | `genome_index(strategy = "bowtie", chroms = )` → `specificity_of()` → `design_kasp_markers(BPPARAM = )` | per-shard index, verdicts, memory/time record | `test-scenario-j4.R` |
| **J5** | Lab lead: the lab answers back and the next design uses it | artefact + `outcomes.csv` | `marker_artefact_read()` → `marker_outcomes()` → `marker_priors(outcomes = )` → `marker_tech()` → `marker_solve()` | next `marker_solved`, ordering provenance | `test-scenario-j5.R` |
| **J6** | IP / release: a marker you can defend and replay | `marker_solved` | `marker_artefact(created = , scenario = "J6")` → `marker_export_rocrate()` → `marker_patent_export()` → `marker_patent_validate()` → `marker_replay()` | artefact, RO-Crate, ST.26, replay verdicts | `test-scenario-j6.R` |
| **J7** | Agent: drive the package over MCP | MCP request fixtures | `tools/list` → `tools/call` (read-only and mutating) | JSON-RPC results | `test-scenario-j7.R` |

J3, J4, J6 and J7 are the scenarios that force the four decisions in ADR 0009,
0010 and 0011. J1, J2 and J5 are the regression floor: if the quick start,
the batch or the outcome loop breaks, the rest is not a product.

### 2. Executable means a test, not a diagram

Each scenario is one `tests/testthat/test-scenario-jN.R` file with:

1. a committed fixture under `tests/testthat/fixtures/scenarios/` or an
   existing fixture (`fixtures/panel/`, `fixtures/tair10/`,
   `fixtures/scenarios/`); no network and no download;
2. the exact call chain above, with `expect_s3_class()` at every typed step;
3. no `$`, `[[` or column pull on a typed object between steps — an accessor
   such as `as_marker_solved()` or `as_primer_set()` is the only allowed
   bridge, and it must be exported and listed in the CSV;
4. a declared degradation: if an engine (OligoArrayAux, bowtie, BLAST) is
   absent, the test may skip the design call but must still assert the
   contract step that does not need the engine. A pass that is only a skip is
   not a scenario pass.

The tests are the runtime proof. A static proof lives in
`dev/check-user-api.R`: every step after the first must be a typed edge already
present in `docs/user-api-surface.csv`, and every function named in a scenario
must be an export or an S3 method. This catches a scenario that compiles only
because a new hidden bridge was added without a CSV row.

Every scenario run also writes the metrics contract in
[ADR 0011](0011-release-artefact-and-mcp-contract.md) §5: one JSONL line per
stage, numbers only. A scenario without a metrics file is not a run, and the
daily report cannot compute p95 or cost per route without it. The baseline for
J1–J7 is captured on `main` before the first code change, not at the end of the
programme.

### 3. A scenario does not add API by default

The default is to reach the goal with the current 207 exports. A new export is
allowed only when a scenario is impossible without it, and then it is named in
the ADR or `ROADMAP.md` §ZV item that owns the seam and in
`docs/user-api-surface.csv` in the same change. The scenario paths as written
need two accessors:

- `as_marker_solved()` — for J2 (`report`) and J3 (`panel`);
- `as_primer_set()` — for J1/J6 when a marker's primers are passed to a
  thermodynamic or export function that expects a named set.

The ZU items scheduled in the same programme may add the input readers and
registry adapters they name (`read_variant_genotypes()`, `as_variant_row()`,
`as_scheme_ids()`/`as_limits()`); those are owned by their ZV item, not by this
ADR. No `read_gwas()`, no `marker_pipeline()` class and no `run_scenario()`
helper are added. J2's VCF reader and J3's panel input stay the existing
`read_variants_vcf()` and `marker_design_panel()`.

### 4. Failure is a named condition, not a broken example

A scenario fails when any of these is true:

- a public step needs `$`, `[[` or a column pull to reach the next step;
- the expected class at a step is absent from `output_class` in the CSV;
- the test needs the network, an uncommitted file, or an interactive session;
- the only green path is a `skip_if_no_*()` with no asserted fallback;
- the verdict of the scenario changes on the golden set while the input does
  not (the `D1`–`D5` rule);
- the scenario's artefact is not reproducible from its own fixture;
- the scenario run has no JSONL line, or a line carries text/content instead of
  the numeric fields of ADR 0011 §5.

The test name and the failing step are the error message. "The example does not
run" is not an acceptable report.

## Gates

1. `Rscript dev/test.R scenario` runs J1–J7 and exits non-zero on any failure.
2. `Rscript dev/check-user-api.R` exits non-zero when a scenario names a
   function that is not in the CSV, or a step that is not a typed edge.
3. `Rscript dev/check.R --quick --no-render` stays green; a new export without
   a CSV row fails the user-API gate, not a review comment.
4. `Rscript dev/check-docs.R` stays green: the scenario table in this ADR and
   the plan link from `docs/INDEX.md` resolve.
5. `Rscript dev/metrics.R --check <date>` accepts that day's JSONL and
   `--report <date>` writes `data/metrics/metrics-<date>.md`; a scenario run
   without a metrics file fails the daily gate.

## Consequences

- The 207 exports keep their number unless a scenario proves a bridge is
  needed. This decision authorises two accessors; the ZU items scheduled
  alongside it own any input reader or registry adapter they add.
- `docs/user-api-surface.md` becomes the input to the scenarios, not only a
  report about them; a new export now has a scenario cost.
- J1–J2 are available on day 1 of the three-day programme, J3/J5 on day 2,
  J4/J6/J7 on day 2–3 (see `docs/plans/2026-11-27-integration-3day.md`).
- A scenario that cannot be made green within its day is a scope cut in
  `ROADMAP.md` §ZV, not a silently skipped test.

## Alternatives

- **Write a vignette per scenario.** Rejected: a vignette is not executed on
  every merge, and the existing quickstart already drifts only because
  `dev/check-docs.R` compares one chunk. A test is the contract.
- **Generate the scenarios from the graph.** Rejected: the graph has 934 edges
  and 116 leaves; a path through it is not a user goal. The seven goals are
  chosen, then checked against the graph.
- **Add a `run_scenario(id)` dispatcher.** Rejected: it would hide the call
  chain behind a string and make the test assert the dispatcher, not the API.
- **Treat a skipped engine test as a pass.** Rejected: the package already has
  honest `skip_if_no_*()` helpers; a scenario that only skips proves nothing.

## Revisit trigger

A new user-facing goal that cannot be expressed as one of J1–J7, or a scenario
that needs a third bridge. At that point the scenario set changes in a new ADR,
not by widening an existing test.

## References

- [`docs/user-api-surface.md`](../user-api-surface.md) — the inventory, the
  graph rules and the nine surface defects.
- [`dev/check-user-api.R`](../../dev/check-user-api.R) — the static gate.
- [`README.md`](../../README.md) — the two quick starts J1 mirrors.
- [`tests/testthat/test-panel.R`](../../tests/testthat/test-panel.R) — the
  existing panel test J3 extends.
- [`tests/testthat/fixtures/panel/`](../../tests/testthat/fixtures/panel) — the
  committed phased-haplotype fixture.
- [`ROADMAP.md`](../../ROADMAP.md) — §ZU (defects) and §ZV (the three-day queue).

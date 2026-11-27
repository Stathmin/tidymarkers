# tidymarkers — documentation index

Canonical entry point for every document in the repository. `dev/check-docs.R`
verifies that each document is linked from this page and that every relative
link resolves; `dev/check-artifacts.R` verifies the artifact graph.

## Project

- [README](../README.md) — what the package is, install, quick start, decision trees.
- [AGENTS](../AGENTS.md) — agent instructions and the environment to run in.
- [CONTRIBUTING](../CONTRIBUTING.md) — conventions, scope and the merge gate.
- [CREDITS](../CREDITS.md) — contributor record and how work is ranked here.
- [ROADMAP](../ROADMAP.md) — the open queue (ZA8 is this docs/artifacts layer).
- [NEWS](../NEWS.md) — changelog.
- [VALIDATION](../VALIDATION.md) — what is checked against external evidence,
  what is not, and where the gaps are.

## Vignettes

Marker chemistries:

- [KASP and allele-specific designs](../vignettes/markers-kasp.Rmd) — `snp`,
  `indel`, `multi`, `duplex`, `lax`, `asq`, `primer`.
- [Gel and melt designs](../vignettes/markers-gel.Rmd) — `caps`, `arms`, `hrm`,
  `ssr`, `indel_size`.
- [Probe designs](../vignettes/markers-probe.Rmd) — `taqman`, `dpcr`,
  `snapshot`, `castpcr`.
- [Isothermal designs](../vignettes/markers-isothermal.Rmd) — `lamp`, `rpa`.
- [CRISPR and Sanger sequencing](../vignettes/markers-crispr-sequencing.Rmd) —
  `crispr`, `sanger`.

Other vignettes:

- [Vignette index](../vignettes/index.Rmd) — one page into all of them.
- [KASP introduction](../vignettes/kasp-intro.Rmd)
- [Custom genome](../vignettes/custom-genome.Rmd)
- [Reading results](../vignettes/reading-results.Rmd)

## Developer docs

- [dev/README](../dev/README.md) — the internal checks, renderers and benchmarks.
- [dev/dag-rebuild-plan](../dev/dag-rebuild-plan.md) — the decision-DAG rebuild:
  plan, invariants (C1–C14), agent scopes.
- [docs/plans/2026-11-27-integration-3day](plans/2026-11-27-integration-3day.md)
  — the three-day integration programme: J1–J7, parallel streams, gates and
  failure criteria.
- [docs/plans/day-1-gate](plans/day-1-gate.md) — the Day-1 state report and
  merge-marker evidence.
- [docs/plans/self-repair-log](plans/self-repair-log.md) — the append-only
  blocker and repair log for the programme.
- [data/metrics/metrics-2026-11-27](../data/metrics/metrics-2026-11-27.md) —
  the generated run report: stage latencies, J1–J7 cost, failure share and the
  delta to the baseline.
- [dev/lang_purity](../dev/lang_purity.md) — the prose-tic audit: measured
  tiers, per-tier triage, and the verbosity meter roadmap.
- [docs/user-api-surface](user-api-surface.md) — the full export inventory,
  data classes, composition rules, the connectivity graph and the ZU defects.
- [docs/index-perf](index-perf.md) — one locus, memory vs per-chromosome bowtie:
  build and query times, and what the index object does not hold.

## Architecture decisions

- [adr/README](adr/README.md) — the decision records, and the full triage table
  for the 2026-11-26 tool list.
- [ADR 0001](adr/0001-llm-mt-and-agent-tooling-out-of-scope.md) — LLM, machine
  translation and agent tooling is out of scope (18 tools rejected).
- [ADR 0002](adr/0002-targets-stays-dev-only.md) — `targets` stays a dev-only
  freshness holder, not a package dependency.
- [ADR 0003](adr/0003-duckdb-deferred.md) — DuckDB/`duckplyr` is deferred until
  an aggregation is measured.
- [ADR 0004](adr/0004-mcp-ci-and-cross-project-tooling.md) — the MCP server
  stays one R file; no Python client, no published JSON Schema.
- [ADR 0005](adr/0005-index-backed-genome-access.md) — genome access is an
  index-backed strategy with declared capabilities.
- [ADR 0006](adr/0006-panel-pipeline.md) — the panel pipeline is an
  orchestrator with falsifiable seams.
- [ADR 0007](adr/0007-module-seams-and-mcp-contract.md) — module seams
  (contracts, registries, splits) and the MCP runtime contract.
- [ADR 0008](adr/0008-user-scenarios-as-contract.md) — user scenarios J1–J7
  are the executable contract.
- [ADR 0009](adr/0009-panel-orchestrator-contract.md) — the panel orchestrator:
  GWAS, LD, haplotypes, chemistry, artefact.
- [ADR 0010](adr/0010-index-strategy-at-scale.md) — the index strategy at
  scale: shards, builds and a memory budget.
- [ADR 0011](adr/0011-release-artefact-and-mcp-contract.md) — the release
  artefact and the MCP contract are the product.
- [ADR 0012](adr/0012-autonomy-and-self-repair.md) — autonomy and self-repair:
  a blocker never stops the programme.
- [ADR 0013](adr/0013-metrics-gate-is-not-a-green-delta.md) — the metrics
  report is not a green delta; a p95 gate compares comparable samples only.

## Operator decisions

- [decisions/README](decisions/README.md) — dated operator decisions that close
  a queue item with a measurement or a policy.
- [Decision 0001](decisions/0001-16gb-fixpur.md) — a 16 Gb fixpur is a disk job
  on per-chromosome shards (measured peak RSS), not a 21.5 GB RAM job.
- [Decision 0002](decisions/0002-ci-strict-gate.md) — CI fails when a linter, a
  test runner or the test files are absent; the release leg runs the local gate.

## Docs and artifacts layer

- [docs/INDEX](INDEX.md) — this page.
- [docs/decision-audit](decision-audit.md) — where each threshold's number comes
  from (`source` and the `CRITERION_EVIDENCE` rung), and where a cross-scheme
  invariant deliberately does not hold (ZW0.1, Z11.1).
- [docs/artifacts.yml](artifacts.yml) — the artifact manifest (source of truth).
  Committed artifacts only: the logo, the two README decision trees, the 20
  per-chemistry decision SVGs in `vignettes/figures` and the 20 vignette gallery
  SVGs. The `.txt` twins of the decision SVGs are committed too; only the
  `vignettes/figures/*.png` rasters are generated locally and ignored.
- [docs/generated/artifacts.md](generated/artifacts.md) — the generated artifact
  graph; one heading per artifact with links to the artifact, its generator,
  its code and the docs that consume it.
- [docs/generated/graphify-health.md](generated/graphify-health.md) — the code
  graph's health and what it does not prove; regenerated by
  `dev/graphify-health.R --write` from the local `graphify-out/graph.json`.

Refresh and check the manifest after changing a figure, its generator or the
engine it reads:

```sh
Rscript dev/artifacts.R --write    # hashes + docs/generated/artifacts.md
Rscript dev/check-artifacts.R      # stale/orphan gate (exit 1)
Rscript dev/check-docs.R           # orphans, links, include_graphics, backlinks
Rscript dev/graphify-health.R --check   # code-graph invariants (skips without a graph)
```

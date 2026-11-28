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
- [docs/plans/day-2-gate](plans/day-2-gate.md) — the Day-2 panel, loop and
  scale gate: runs, counts and the two open items.
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
- [docs/ci-bench-baseline](ci-bench-baseline.md) — the runbook for capturing
  the `github-ubuntu-latest` timing baseline on CI and committing it.

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
- [ADR 0014](adr/0014-benchmark-baseline-is-keyed-by-runner-class.md) — a
  timing baseline is keyed by runner class, not by hostname.
- [ADR 0015](adr/0015-a-record-resolves-only-its-own-directory.md) — a
  persisted record names only the files inside its own directory.
- [ADR 0016](adr/0016-one-vcf-admission-rule.md) — the two VCF readers share
  one admission rule, so a symbolic allele cannot become a fake `N`.
- [ADR 0017](adr/0017-a-raw-identifier-is-admitted-before-normalisation.md) —
  a raw identifier is admitted before `ath_chrom()` can rename it, so an
  empty chromosome cannot become a plausible contig.
- [ADR 0018](adr/0018-one-allele-and-position-rule-for-three-vcf-readers.md) —
  all three VCF readers admit the same alleles and refuse the same positions,
  so a depth table and a variant table join on one key.
- [ADR 0019](adr/0019-an-absent-allele-is-not-a-base.md) — a missing or empty
  REF/ALT is the absent state, dropped by the same rule before
  normalisation, so an empty field cannot become a base-R error or an empty
  allele.
- [ADR 0020](adr/0020-a-gt-that-names-a-missing-allele-is-not-a-call.md) — a
  `GT` naming an allele the record does not have is `NA`; the `AD` fallback
  answers only for an absent or no-call `GT`.
- [ADR 0021](adr/0021-one-path-and-one-header-rule-for-three-vcf-readers.md) —
  all three VCF readers name a bad path the same way and find the `#CHROM`
  line by scanning, not by a 1000-line cap.
- [ADR 0022](adr/0022-a-fasta-base-outside-the-alphabet-is-refused-by-every-reader.md) —
  a FASTA base outside the DNA/IUPAC alphabet is a named refusal in the
  transcript reader too, not a silently shortened sequence.
- [ADR 0023](adr/0023-a-refusal-that-a-retry-cannot-fix-carries-a-class.md) —
  a refusal a retry cannot fix (bad alphabet) carries a condition class, so a
  reader that refetches a broken cache does not retry the unfixable case.
- [ADR 0024](adr/0024-a-replay-compares-only-the-records-own-verdicts.md) —
  a replay compares only the record's own verdicts; a record with no
  `replay.csv` is refused when replayed with a trusted `solved`, instead of
  passing on the caller's own solve.
- [ADR 0025](adr/0025-a-record-member-is-resolved-not-only-named.md) — every
  member a record names (sidecar, outcomes, targets, crate descriptor, copied
  file) is resolved under the record's real root before it is opened; a
  symlink out is a named refusal, not a read of the neighbour's file.
- [ADR 0026](adr/0026-an-artefacts-id-is-the-hash-of-what-it-stores.md) — an
  artefact's id is the hash of the JSON it stores and is recomputed on read; a
  manifest edited after it was written is refused instead of replayed.
- [ADR 0027](adr/0027-a-record-directory-holds-one-design.md) — writing a
  design into a directory replaces the whole record, including `solved.rds`,
  so a rewritten manifest never replays the design it replaced.
- [ADR 0028](adr/0028-a-vcf-writer-emits-only-records-its-readers-admit.md) —
  `write_variants_vcf()` refuses a row its three readers would drop or
  reinterpret (symbolic/absent allele, POS outside the integer range, a tab in
  a field, a non-ACGTN allele, a non-representable depth), so a written VCF
  reads back as written.
- [ADR 0029](adr/0029-an-output-writer-creates-its-files-not-writes-through-a-link.md)
  — `write_marker_results()` and `export_primers()` declare every file they
  create and refuse to write through a link, the same rule and message as the
  record writers (ADR 0025 point 5), so an output set is never written outside
  its directory.
- [ADR 0030](adr/0030-a-derived-path-component-is-a-declared-name.md) —
  every path component a writer derives is a declared name: the outcomes
  handover refuses a link at `crates` or at a crate's pseudonym before the first
  crate exists, not only at the files inside it.
- [ADR 0031](adr/0031-a-named-root-is-normalised-before-a-sibling-is-derived.md)
  — the handover normalises `dir` before deriving `<dir>-identity.csv` from it;
  a trailing separator used to move the identity map inside the handover, where
  it would travel with the pseudonymised crates.
- [ADR 0032](adr/0032-the-corpus-exclusion-check-reads-the-resolved-directory.md)
  — the handover decides "outside the scanned corpus root" on `fs::path_real()`
  after creating `dir`, so a link at `dir` whose target is inside the corpus is
  refused instead of writing crates the next scan reads back.
- [ADR 0033](adr/0033-a-cache-and-a-rendered-file-are-declared-names.md) — a
  cache entry and a rendered medium are declared names too: `ensembl_gene()`,
  `ensembl_variants()`, `tair10_download()`, `marker_ensdb()` and
  `assemble_gif()` refuse a link at the name they are about to write.
- [ADR 0034](adr/0034-an-index-build-declares-every-file-it-writes.md) — a genome
  index is a declared set of names too: the bowtie, BLAST and FM builders refuse
  a link at any file the engine creates, and a linked `.1.ebwt`/`.nsq`/`.files`
  is not counted as a finished index.
- [ADR 0035](adr/0035-a-guard-refuses-a-name-a-hit-test-only-reports-one.md) —
  the writer gate's guard list names only calls that refuse a name; a cache hit
  test reports and is not a guard, and the write scanner matches the imported
  call form as well as the namespaced one.
- [ADR 0036](adr/0036-a-nested-closure-is-a-writer-of-its-own.md) — a closure
  inside a function is its own writer: the guard must be in the closure that
  writes, and a primitive built as text is not a write.
- [ADR 0037](adr/0037-the-write-vocabulary-names-every-primitive.md) — the write
  vocabulary names every primitive that creates a file: `writeXStringSet()` is
  a writer, and the self-test proves the list cannot silently shrink.
- [ADR 0038](adr/0038-a-depth-the-writer-refuses-is-not-a-depth-the-reader-reports.md) —
  the `AD` fallback admits only non-negative whole counts within integer range;
  a depth the writer refuses is `NA`, not a dosage outside `0`/`1`/`2`/`NA`.
- [ADR 0039](adr/0039-a-negative-depth-carries-no-confidence.md) — the depth
  term of `variant_confidence()` is floored at zero, so `weight` stays in
  `[0, 1]` and a negative `DP` weighs nothing rather than less than nothing.
- [ADR 0040](adr/0040-a-validator-enforces-the-unit-the-contract-names.md) —
  the metrics validator enforces the unit the contract names: integer fields
  are whole numbers, and `gpu_util` is a fraction in `[0, 1]`, not a
  percentage.
- [ADR 0041](adr/0041-a-genotype-is-its-allele-multiset.md) — a biallelic
  genotype is classified by its allele multiset, so `1/0` and `1|0` get the
  same `het` class, weight and depth-gate verdict as `0/1`.
- [ADR 0042](adr/0042-a-threshold-is-a-count-not-a-number-to-truncate.md) — a
  sample-size threshold is validated as a whole number before the cast, so
  `min_lines = 4.9` no longer passes a four-line panel.
- [ADR 0043](adr/0043-a-panel-gives-one-call-per-line-and-target.md) — a panel
  carries one call per line and target; two different calls are refused
  instead of resolved by file order.
- [ADR 0044](adr/0044-a-count-has-one-definition.md) — one
  `is_whole_number()` predicate backs `min_lines`, `max_markers`,
  `spacer_len` and `min_evidence`, so `Inf` reaches each gate's own message.
- [ADR 0045](adr/0045-a-position-is-a-whole-number-before-the-cast.md) — a
  `pos` is checked before `as.integer()`, so a panel row at 250.9 no longer
  scores the target at 250.
- [ADR 0046](adr/0046-a-design-key-names-every-input-the-design-reads.md) — a
  design key names `alt2` and `strand`, so two rows at one locus are no longer
  answered by one design.
- [ADR 0047](adr/0047-a-search-radius-is-a-whole-number-before-the-cast.md) — a
  `max_mismatch` or `k` is checked before `as.integer()`, so `1.9` no longer
  scans radius 1 and `Inf` no longer reaches a base R error.
- [ADR 0048](adr/0048-a-base-pair-threshold-is-a-whole-number-before-the-cast.md) — a
  `min_dist`, `min_gap` or `flank` is checked before `as.integer()`, so `5.9`
  no longer lowers a resolvable-gap floor to 5.

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

# 0010 — The index strategy at scale: shards, builds and a memory budget

- Status: proposed
- Date: 2026-11-27
- Supersedes: none
- Builds on: [ADR 0005](0005-index-backed-genome-access.md) (the strategy)

## Context

ADR 0005 made genome access a declared strategy and `ZT T1`/`T2` implemented
it: `genome_index()` carries `capabilities` and a 12-column `manifest`;
`specificity_of()` dispatches on the strategy; `scan_question()` refuses to mix
a 3'-seed question with a whole-primer question; a chromosome is a shard; a
shard is opened lazily; a valid shard is reused without a new mtime; and a
missing shard names its chromosome.

The measured numbers are in `docs/index-perf.md`:

- one locus, memory vs per-chromosome bowtie: median 18.64 s vs 18.32 s,
  index/memory 0.98×; maximum 21.69 s (1.18× the median), inside the 2× budget;
- `dev/bench-index-memory.R` (`VmHWM`): baseline 455.4 MB, TAIR10 in memory
  569.7 MB, one 30.4 Mb shard 563.2 MB, a synthetic 250 Mb chromosome 885.5 MB.
  That is 1.7 bytes per base over baseline; the conservative planning figure is
  5 bytes per base, so a 0.8 Gb chromosome costs about 4.5 GB.

The threshold is `BOWTIE_MONOLITH_LIMIT = 21.5 * 1024^3` (1.34 bytes per base
for a 16 Gb genome). Above it, a bowtie build is sharded even when `chroms =
NULL`. `docs/decisions/0001-16gb-fixpur.md` records that a 16 Gb fixpur is a
disk job on per-chromosome shards on a 16 GB host, not a monolithic RAM job.

Three gaps remain, and they are the reason this ADR is not already closed by
`T2`:

- **A pangenome has no explicit build identity in the query.** `marker_genome()`
  has no `build =` argument, so `genome$build` is normally `NULL` and
  `genome_build` is `NA` in every manifest row. A/B/D are a named list of
  indexes, but nothing stops a query from combining their hits.
- **The memory budget is a measurement, not a gate.** `bench-index-memory.R`
  prints `VmHWM`; no gate fails when a shard crosses 5 bytes per base or when a
  query builds a monolith.
- **The pangenome path is exercised by a three-build test, not by J4.** The
  test builds six shards and sums them; it does not assert the selection rule
  (one build per query), the missing-shard error with build and chromosome, or
  the budget on a fixture that stands in for the 16 Gb case.

## Decision

### 1. Build identity is explicit and required above one build

`marker_genome()` and `read_genome_fasta()` gain `build = NULL`. When the
caller supplies a pangenome, `build` is a short id (`"A"`, `"B"`, `"D"`, or a
release name) and it is written into every manifest row. `genome_index()` also
accepts `build =` and refuses to guess when the genome has none and the caller
asks for a pangenome.

A pangenome is a **named list of `genome_index` objects**, one per build. It is
not a new class and not one index over the union. The shard key is
`(genome_build, chrom)`; the path contains the build, so two builds with the
same chromosome name cannot collide.

### 2. One query, one build

A specificity query names one build. `specificity_of()` and
`specificity_verdict()` carry the build in the result; a hit table that mixes
builds aborts with `tidymarkers_mixed_build_error` naming both builds. There is
no implicit union: a hit on B is not evidence about A.

This is the same rule as `scan_question()`: the answer is only meaningful for
the question that was asked. A pangenome scan is N queries with N verdicts, or
a helper that loops and returns one row per build; it is never one pooled
verdict.

### 3. Lazy opening is a contract, not an optimisation

Opening a `genome_index` does not read a FASTA. The manifest is the only thing
read until a shard is queried. A shard that exists and verifies is reused
without changing its mtime. A missing shard aborts with the build, the
chromosome, the expected path and the command that rebuilds it. `index_verify()`
is the integrity check; a hash mismatch is an error, not a rebuild.

### 4. The budget is enforced

The planning budget is **5 bytes per base over baseline** for a shard and
**2× the memory-path query time** for a query on the golden locus. The measured
1.7 bytes per base is the evidence for the 5; the margin covers the query and
the R object overhead.

`dev/perf-gate.R` (the `ZT T11` gate) reads the budget. It fails when:

- a shard's `VmHWM` exceeds `baseline + 5 * bases`;
- a query on the golden locus exceeds 2× the memory-path median;
- a build above `BOWTIE_MONOLITH_LIMIT` is attempted without `chroms`;
- the budget cannot be measured and the gate is reported as `not_run` without a
  reason.

`not_run` with a reason is honest; a silent pass is not. The 16 Gb operator run
stays an operator action recorded in `docs/decisions/0001-16gb-fixpur.md`; the
unit test uses the three small builds and the synthetic 250 Mb chromosome.

### 5. J4 is the acceptance scenario

J4 runs the three-build fixture through the public path:

```r
a <- read_genome_fasta(fixture_a, build = "A")
b <- read_genome_fasta(fixture_b, build = "B")
d <- read_genome_fasta(fixture_d, build = "D")
ia <- genome_index(a, strategy = "bowtie", chroms = genome_chroms(a))
# one query per build; one verdict per build
```

The test asserts: six shards, one row per `(build, chrom)`; a query for A does
not see B; a missing B shard names `B` and the chromosome; a second object over
the same directory verifies without rebuilding; the memory record is attached
to the result, not printed.

## Gates

1. `tests/testthat/test-genome-index.R` keeps the shard, lazy-open, reuse and
   missing-shard assertions; it gains the build column and the mixed-build
   abort.
2. `tests/testthat/test-scenario-j4.R` runs J4 on the committed three-build
   fixture and asserts one verdict per build.
3. `Rscript dev/bench-index-memory.R` reports `VmHWM`; `dev/perf-gate.R` fails
   the budget with a table, not a message.
4. `Rscript dev/check.R --quick --no-render` stays green; `docs/index-perf.md`
   is regenerated when a number changes.
5. The real 16 Gb run is recorded in `docs/decisions/0001-16gb-fixpur.md`; it is
   not claimed as a unit test.

## Consequences

- A pangenome query is explicit and cannot silently pool A/B/D. The cost is one
  query per build and one verdict per build.
- `marker_genome(build = )` and `read_genome_fasta(build = )` are additive; a
  single-genome caller keeps `NA` and sees no change.
- The manifest's `genome_build` stops being dead weight; a pangenome artefact
  records which build produced which target.
- A shard over 5 bytes per base fails the gate instead of being discovered on
  the 16 Gb host.
- The 16 Gb fixpur remains an operator action; this ADR does not add a fixture
  that the host cannot hold.

## Alternatives

- **One index over the union of A/B/D.** Rejected: the union is ~64 GB and a
  hit cannot be attributed to a build. It also defeats the 21.5 GB threshold.
- **Default the build to the file name.** Rejected: a file name is not a
  release identity, and two builds can share a basename.
- **Pool hits and report the best build.** Rejected: a pooled verdict hides a
  build-specific failure and makes the answer unreproducible.
- **Use RSS from `/proc` instead of `VmHWM`.** Rejected: `VmHWM` is the peak
  and is already measured by `bench-index-memory.R`; a sampled RSS misses the
  peak.
- **Treat the 16 Gb run as a skipped test.** Rejected: the decision record is
  the honest form; a skip would look like coverage.

## Revisit trigger

A host with more than 64 GB RAM, a pangenome format that carries build identity
in the FASTA header, or a query that must compare builds in one statistical
model. At that point the named-list contract is re-opened.

## References

- [`R/marker_genome.R`](../../R/marker_genome.R) — `genome_index()` and the
  strategy dispatcher.
- [`R/marker_index.R`](../../R/marker_index.R) — the manifest, shards and
  `BOWTIE_MONOLITH_LIMIT`.
- [`R/marker_bowtie.R`](../../R/marker_bowtie.R) — the per-chromosome build.
- [`dev/bench-index-memory.R`](../../dev/bench-index-memory.R) — the `VmHWM`
  measurement.
- [`docs/index-perf.md`](../index-perf.md) — the measured query numbers.
- [`docs/decisions/0001-16gb-fixpur.md`](../decisions/0001-16gb-fixpur.md) — the
  16 Gb operator decision.
- [`tests/testthat/test-genome-index.R`](../../tests/testthat/test-genome-index.R)
  — the existing shard tests.
- [`ROADMAP.md`](../../ROADMAP.md) — `ZT T1`/`T2`, `T11`, `ZX2`, §ZV (ZV8).

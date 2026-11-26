# 0005 — Genome access is an index-backed strategy, not an in-memory object

- Status: proposed
- Date: 2026-11-26
- Supersedes: none

## Context

`marker_genome()` (`R/marker_genome.R`) already refuses to hold a genome in RAM:
it wraps an Rsamtools-indexed FASTA, reads the design window on demand, and
caches a whole chromosome only when a genome-wide step asks for it. Whole-genome
verification then dispatches inside `specificity_of()` (`R/marker_scan.R`):
`specificity_seed()` for an in-memory `DNAStringSet`, `specificity_blast()` for
a BLAST database, `specificity_bowtie()` when the genome carries a bowtie index,
and a named error when it carries only an FM index.

The dispatch is a set of optional slots (`blast_db`, `bowtie_index`, `fm_index`)
on one list, not a strategy. Today:

- the caller cannot ask what a genome can answer. The answer is inferred from
  which slot is non-NULL, and `specificity_of()` silently picks a path;
- `k` and `max_mismatch` are meaningful on the seed path and ignored by bowtie;
  the code warns and then answers a different question (D2). The BLAST path
  also ignores `max_mismatch`: it counts full-length exact seed hits only;
- `marker_bowtie_index(chroms = )` builds one index per chromosome and
  `specificity_bowtie()` sums the shards, but the set is a character vector with
  a class and a `chroms` attribute. It carries no build, version or hash record,
  and it is read eagerly: `bowtie_index_set()` reads the whole FASTA to slice it
  (the shortcut comment in `R/marker_bowtie.R`);
- a monolithic wheat index is about 21.5 GB at the measured 1.34 bytes per base.
  A pan-genome A/B/D at 16 Gb each is about 64 GB, which no host here can build
  or hold;
- `bp_lapply()` (`R/marker_thermo.R`) is a thin pass-through to
  `BiocParallel::bplapply()` plus a source carrier and a progress bar. Locus
  parallelism in `design_kasp_markers()` works (measured 40.5 s -> 8.0 s on six
  cores), but the index-backed path has no contract for workers: a shared
  read-only index, no per-worker rebuild, deterministic order and error
  isolation are all conventions, not guarantees.

`attr(x, "scan")` is the one contract that survived D1-D3: every specificity
table names the question it answers, `"3'-seed"` or `"whole-primer"`, and
`assemble_kasp_assay()` marks its `blast` field `"whole-primer"` explicitly. The
attribute exists on every path, but it remains an annotation, not a type: no
caller is forced to read it before comparing a number with a seed threshold.

## Decision

### 1. One genome object, one declared strategy

`genome_index()` becomes the constructor:

```r
genome_index(genome, strategy = c("memory", "bowtie", "blast", "fm"),
             chroms = NULL, dir = NULL, threads = 1L)
```

It returns a `genome_index` with:

- `strategy` — one of the four, never inferred at query time;
- `capabilities` — named logicals `seed_scan`, `whole_primer_scan`,
  `mismatch_tolerant`, plus `max_mismatch` and `engine`;
- `manifest` — one row per index shard: `shard`, `chrom`, `strategy`, `path`,
  `bytes`, `sha256`, `tool`, `tool_version`, `build_id`, `genome_build`,
  `source_sha256`, `created`;
- `source` — the FASTA path or `BSgenome` key, never the sequence.

`marker_genome(index = )` attaches one. The legacy `bowtie_index=`, `blast_db=`
and `fm_index=` arguments keep working for one release and are translated into a
`genome_index` with a deprecation warning.

The capability matrix is data, not prose:

| strategy | `seed_scan` | `whole_primer_scan` | `mismatch_tolerant` | engine      |
| -------- | ----------- | ------------------- | ------------------- | ----------- |
| memory   | yes         | no                  | yes, any radius     | Biostrings  |
| bowtie   | yes (seed)  | yes                 | yes, at most 3      | Rbowtie     |
| blast    | yes, exact  | yes, exact (`k=NA`) | no, exact only      | BLAST+      |
| fm       | no          | yes                 | fixed at 2          | Rsubread    |

`specificity_of()` asks the capability, not the slot. A request whose `k` or
`max_mismatch` the strategy cannot honour aborts with a classed error naming the
strategy and the missing capability. It does not warn and then answer a
different question. FM stays out of the design scan for the reason recorded in
D5: 31.7 bytes per base is about 490 GB on wheat, and `specificity_fm()` is a
read mapper.

### 2. The scan contract is typed

The wire format stays `attr(x, "scan")`, because existing tables and tests carry
it. The reader becomes explicit:

- every specificity result is classed `specificity_scan`;
- `scan_question(x)` is the only sanctioned accessor and returns one of
  `"3'-seed"` or `"whole-primer"`; an absent or unknown value is an error, never
  `NULL`;
- `specificity_verdict()` refuses a mixed table (a named genome list may carry
  one question per genome, but not two questions for one genome);
- `assemble_kasp_assay()` keeps the `blast` field labelled `"whole-primer"` and
  the aligned `specificity` field labelled `"3'-seed"`; the two fields answer
  different questions, so they are not deduplicated (D3).

This makes D1-D3 a testable contract rather than a comment: a number cannot
reach a seed-threshold comparison without its question.

### 3. Per-chromosome shards and lazy loading

A monolithic index remains the default only while the estimate fits the named
limit: 21.5 GB at 1.34 bytes per base, i.e. a 16 Gb genome. Above that limit,
`strategy = "bowtie"` and `strategy = "fm"` build one shard per chromosome, and
`strategy = "blast"` builds one database per chromosome. `chroms` selects a
subset; `NULL` means all.

The object stores only manifest paths and opens no handle. A query opens each
shard for the duration of the query and closes it; the existing chromosome cache
remains an explicit opt-out (`tidymarkers.genome_cache`), not part of the index
object. The design window keeps using `genome_seq()` and `scanFa()` for region
streaming.
A pan-genome is a named list of `genome_index` objects (one per assembly), and a
shard is keyed by `(genome_build, chrom)`, so an A/B/D query never has to build
a 64 GB monolith.

### 4. Parallelism by locus, for real

`bp_lapply()` gains `chunksize`, `on_error = c("stop", "collect")`, a
deterministic result order and the backend's `RNGseed`. `design_kasp_markers()`
uses `on_error = "collect"`: one failed locus becomes a row with
`verdict = "error"` and the condition message in `note`, instead of aborting the
whole batch. Workers open the index read-only, never build it, and inherit
`options(tidymarkers.genome_cache = FALSE)` on the index path. Nested
parallelism is forced to one core per worker, as it already is in
`assemble_kasp_assay()`.

## Gates

1. **Gold verdicts, not gold counts.** A committed gold set under
   `tests/testthat/fixtures/gold/` holds TAIR10 loci with the expected verdict
   and failing criterion ids for each. `test-genome-index.R` runs the memory and
   bowtie strategies and compares verdicts, criterion ids and `scan_question()`
   — not raw hit counts. D1 already established that an index answers a
   different question; the test asserts the label instead of pretending the
   numbers are equal.
2. **Per-locus budget.** On the TAIR10 gold set, the index path must stay within
   2x of the measured baseline: median 5.7 s per locus on one core and 1.1 s on
   six, with no single locus above 2x the median. `dev/bench.R --gate` records
   the table; a regression fails the gate rather than a later wave.
3. **Manifest integrity.** A rebuilt index reuses an existing shard only when
   `sha256` matches; a mismatch aborts naming the shard. A missing shard names
   the chromosome. Rebuilding does not change the mtime of a valid shard.
4. **Lazy and pan-genome.** A test builds three small assemblies as per-chromosome
   shards, queries them as a named list and sums the hits; no monolithic build
   runs, and the manifest lists every shard.
5. **Parallel equivalence.** Serial and 2/4-worker runs return identical
   verdicts in identical order; a failing locus is one row, not a batch abort.

## Consequences

- `specificity_of()` becomes stricter: a non-default `k` on bowtie, or
  `max_mismatch > 0` on BLAST, is an error. This is a breaking change and goes
  to `NEWS.md`.
- A wheat or pan-genome run no longer needs a 21.5-64 GB monolithic index; it
  needs per-chromosome shards and a manifest. The default location stays
  `tempdir()`; a large build is an explicit `dir=` under `marker_cache_dir()`.
- The index manifest becomes the place where build, tool version and source
  hash are recorded, so a stale shard is detectable without re-reading the
  FASTA.
- The capability table is data, so a future MCP or cross-project consumer can
  read it instead of guessing which index slot is populated. The current MCP
  registry does not expose `genome_index()`.

## Alternatives

- **Add `capabilities` to `marker_genome` and keep the slots.** Rejected: the
  slots stay the source of truth, and the query still infers the question.
- **One monolithic index with memory mapping.** Rejected: bowtie's build itself
  needs the genome and the result is 21.5 GB for wheat; mmap does not remove the
  build cost or the pan-genome product.
- **Make FM the default.** Rejected: 31.7 bytes per base (~490 GB on wheat) and
  a read mapper, not a seed counter (D5).
- **Keep the warning and answer whole-primer.** Rejected: D2 showed the number
  is then read as a seed count. A warning is not a contract.
- **Build indexes eagerly in `marker_genome()`.** Rejected: a 21.5 GB write
  must be an operator decision, not a side effect of opening a file.

## Revisit trigger

A host with enough RAM and disk for a single 64 GB pan-genome index, or an index
format that beats 1.34 bytes per base without giving up the seed question. Also
revisit if a second consumer needs a capability not in the matrix, or if
`Rbowtie`/`Rsubread` change their process contract.

## References

- [`R/marker_genome.R`](../../R/marker_genome.R) — the indexed FASTA source.
- [`R/marker_scan.R`](../../R/marker_scan.R) — `specificity_of()` dispatch.
- [`R/marker_bowtie.R`](../../R/marker_bowtie.R) — shards and whole-primer scan.
- [`R/marker_fm.R`](../../R/marker_fm.R) — the FM read mapper.
- [`R/marker_assemble.R`](../../R/marker_assemble.R) — the `scan` labels in an assay.
- [`R/marker_thermo.R`](../../R/marker_thermo.R) — `bp_lapply()` and `marker_bpparam()`.
- [`ROADMAP.md`](../../ROADMAP.md) — D1-D3, ZX1-ZX3, ZT T1-T3 and T11.

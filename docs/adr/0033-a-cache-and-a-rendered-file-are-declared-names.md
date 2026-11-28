# 0033 — A cache entry and a rendered file are declared names too

- Status: accepted (amended)
- Date: 2026-11-28
- Supersedes: none
- Builds on: [ADR 0029](0029-an-output-writer-creates-its-files-not-writes-through-a-link.md)
  (a writer creates its files, point 1: all names before the first write) and
  [ADR 0025](0025-a-record-member-is-resolved-not-only-named.md) (point 5)

## Context

ADR 0029 gave one rule to the writers of user-visible output, and ZV29–ZV39
applied it to every record, result, VCF, trace, reprex and registry writer. The
roster at the end of ZV39 said no unguarded writer path was left, but the audit
that produced it stopped at files the user asks for. Two other kinds of named
file had no guard:

- **cache entries** — `ensembl_gene()` / `ensembl_variants()` write the Ensembl
  REST body to `<cache>/<species>_<id>.json`, `tair10_download()` writes the
  download half to `<fa>.part` before moving it, and `marker_ensdb()` writes the
  GTF and the EnsDb sqlite into the cache;
- **a rendered medium** — `assemble_gif()` passes the GIF path to ImageMagick.

These are the same shape as ZV39: the function computes a name, a link at that
name redirects the write, and the outside file is replaced. Measured 2026-11-28
on the branch, with one dangling link per call and an outside file holding
`keep`:

| link at | call | result before the fix |
| ------- | ---- | --------------------- |
| `<cache>/arabidopsis_thaliana_AT4G00650.json` | `ensembl_gene()` | the response replaced the outside file (119 bytes of JSON); call succeeded |
| `<cache>/…_var_<region>.json` | `ensembl_variants()` | same shape: body written through the link |
| `<fa>.part` | `ath_genome()` download | `utils::download.file()` writes through the link |
| `<cache>/….gtf.gz` | `marker_ensdb()` | `utils::download.file()` writes through the link |
| `<cache>/….EnsDb.sqlite` | `marker_ensdb()` | `ensDbFromGtf()` writes through the link |
| `path` | `render_design_search()` | ImageMagick writes the GIF through the link |

The first row is measured end to end; the other rows are the same call shape
(`utils::download.file()`, `ensDbFromGtf(outfile=)`, an external tool with an
output argument) and are not measured here because they need the network, the
`ensembldb` stack or ImageMagick.

A cache is not a record: ZV39 let the GenBank cache be *skipped* because the
input was still readable without it. Here the write is not optional on a cache
miss — the response is the result, and `marker_ensdb()` builds from the file it
just wrote — so the call refuses instead of continuing without a cache.

## Decision

1. **A named file is guarded whoever asked for it.** Every writer that creates a
   named file calls `declared_path()` (a one-line wrapper over `refuse_links()`)
   before the first write, including a cache entry and a rendered medium.
2. **One rule, one message.** All sites use the same `refuse_links()` and its
   message ("which is a link"), so a caller reads the same sentence whichever
   writer it called.
3. **The guard covers the name the call actually writes.** `tair10_download()`
   guards the `.part` file it downloads into; the final move with
   `fs::file_move()` replaces a link at the destination rather than following
   it, so the destination needs no second guard.
4. **A reader does not take a link at a cache name as a cache hit.** The guard
   stays on creation for a path the caller supplies (ADR 0025 point 5): a
   reader follows the record or directory it was handed. A cache entry is
   different because the package named it, not the caller, and a link there
   means the file was put by someone else. `cache_present()` answers "a file
   this package created" (`file_exists()` and no link) and every hit test uses
   it. When the answer is no, points 1–3 apply: the call declares the name
   before the request, so a link is refused instead of followed, and a missing
   file is created.

## Consequences

- A link at any of the six names stops the call with a named error before the
  outside file is touched. Tests: `test-target.R`, "the Ensembl cache refuses to
  write through a link" (both endpoints, mocked HTTP); `test-output.R`,
  "declared_path refuses a link and keeps an ordinary path"; `test-search-anim.R`,
  "assemble_gif refuses to write through a link". Before the fix the Ensembl
  test created the outside file; the other two refused for the wrong reason
  (missing tool) or did not run.
- The download and EnsDb sites are wired through the same tested helper; their
  own end-to-end proof is deferred with the cache-owner question below.
- The roster is now "every writer that creates a named file", not "every writer
  of a user file". A new writer that skips `declared_path()` is a defect of the
  same family, and the audit question is one grep for the write call, not a
  judgement about who consumes the file.

## Amendment 2026-11-28

The roster was tested again the same day, and it was still short by three
groups: a set of files created by one build, a derived FASTA read back as a
cache, and an operator install.

- **An index build** writes a set of names, one of which later marks the build
  as finished. `bowtie_index_build()`, `index_write_chrom_fasta()`,
  `marker_fm_index()` and `marker_blast_index()` wrote their FASTA through a
  link, and a linked `.1.ebwt` / `.nsq` / `.files` was reported as a finished
  index. Closed by `declared_index()` and `cache_present()`; the decision is
  [ADR 0034](0034-an-index-build-declares-every-file-it-writes.md), ZV41.
- **`fasta_index()`** decompresses `<x>.gz` into `<x>` and builds `<x>.fai`.
  A link at `<x>` was accepted as the decompressed copy and read as the FASTA;
  a link at `<x>.fai` was accepted as a ready index, so random access used a
  foreign index. Both names are now declared; ZV42.
- **`marker_install_tools()`** downloads an archive, unpacks a tree and
  installs `bin`/`share` into `dest`. None of the four names was guarded, so a
  link redirected the download, the unpack or the build. All four are declared
  before the download; ZV43.

With these three, the roster sentence above holds again for every writer the
audit can find: one grep for the write call (`writeXStringSet`, `writeLines`,
`write_csv`, `file_copy`, `download.file`, `untar`, an engine output argument).

## Amendment 2026-11-28 (a cache hit)

The roster above closed the writers; the hit tests were not audited. Every one
was `fs::file_exists()`, which follows a link, so a link at a cache name was
read as the file the package would have downloaded or built:

| name | call | result before the fix |
| ---- | ---- | --------------------- |
| `<cache>/arabidopsis_thaliana_<id>.json` | `ensembl_gene()` | the link's JSON was returned as the Ensembl record |
| `<cache>/…_var_<region>.json` | `ensembl_variants()` | the same |
| `<cache>/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz` | `ath_genome()` | the linked FASTA was read as the genome |
| `<cache>/….gtf.gz`, `<cache>/….EnsDb.sqlite` | `marker_ensdb()` | the linked GTF was used to build, or the linked sqlite was opened as the annotation |

Measured 2026-11-28: the Ensembl and FASTA rows were run end to end with a
valid file behind the link, and each call returned the outside content; the
EnsDb row has the same shape (`ensDbFromGtf()`, `EnsDb(out)`) and is covered by
the refusal, which does not need `ensembldb`.

- `cache_present()` is renamed `cache_present()`: it was already the predicate
  "a file this package created", and the index builders and the caches now call
  it by that name.
- Each cache declares its name with `declared_path()` before the request or the
  build: `ensembl_gene()`/`ensembl_variants()` before the HTTP call,
  `ath_genome()` before `tair10_download()`, `marker_ensdb()` before
  `check_installed("ensembldb")`. The declaration turns a link at a hit into
  the shared refusal; the hit test alone only reports a miss.
- Tests: `test-target.R`, "the Ensembl cache does not read a hit through a
  link"; `test-ath-genome-cache.R`, "a link at the cache name is refused, not
  read as the genome"; `test-annotate.R`, "marker_ensdb refuses a link at
  either cache name". Each keeps the outside bytes and asserts "which is a
  link"; before the fix the first two returned the outside content. ZV44.

The roster now covers the read side too: a hit test that asks `file_exists()`
alone is a defect of the same family.

## Amendment 2026-11-28 (the roster is a checked ledger)

The roster was closed by hand, one audit at a time: ZV29, ZV33–ZV44. A rule
enforced by memory fails the first time nobody rereads it, so the audit is now a
gate.

- `dev/check-writers.R` parses every top-level function in `R/` and fails when a
  body calls a write primitive (`write_csv`, `write_lines`, `writeLines`,
  `saveRDS`, `save`, `ggsave`, `write_xml`, `download.file`, `file_copy`,
  `file_move`, `file.copy`, `file.rename`, `writeBin`, `write.csv`,
  `write.table`) without `refuse_links()`, `declared_path()`,
  `declared_index()` or `cache_present()` in the same body.
- `dev/writers.yml` is the ledger of the exceptions, each with a reason that
  names the declaring caller or the private temp directory. A registered symbol
  that stops writing unguarded is stale and fails too, so the ledger cannot
  quietly outlive its excuse.
- The scan is a heuristic over source text, not a proof: `reprex_from_table()`
  is registered because it builds script text that contains `writeLines()`, not
  because it writes. The gate's own `--self-test` runs one unguarded and one
  guarded synthetic body through the checker, so a scanner that stopped
  matching fails instead of passing.
- Step `2b2` in `dev/check.R`. Measured 2026-11-28: 24 writer bodies, 17 declare
  their names inline, 7 are registered (three write into `withr::local_tempdir()`,
  two are helpers whose caller declares the name, one is the script-text false
  positive, one is `fold_amplicon()`'s private FASTA). ZV45.

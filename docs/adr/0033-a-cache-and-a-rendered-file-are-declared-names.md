# 0033 — A cache entry and a rendered file are declared names too

- Status: accepted
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
4. **A reader still follows the path it is given.** The guard is on creation
   only; a cache hit reads through whatever path the caller supplied, as ADR
   0025 point 5 already says for readers.

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

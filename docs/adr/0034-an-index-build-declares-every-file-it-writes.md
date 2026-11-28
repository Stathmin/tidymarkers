# 0034 — A genome index declares every file its build writes

- Status: accepted
- Date: 2026-11-28
- Supersedes: none
- Builds on: [ADR 0029](0029-an-output-writer-creates-its-files-not-writes-through-a-link.md)
  (a writer creates its files, point 1: all names before the first write) and
  [ADR 0033](0033-a-cache-and-a-rendered-file-are-declared-names.md)
  (a cache entry is a declared name too)

## Context

ADR 0029 gave one rule to the writers of user-visible output, ZV29–ZV39 applied
it to every record, result, VCF, trace, reprex and registry writer, and ADR 0033
extended it to a cache entry and a rendered medium. The roster those two closed
still stopped at files a single call writes with one R function. A genome index
is a different shape: an external engine (`bowtie-build`, `makeblastdb`,
`Rsubread::buildindex`) creates *a set* of files next to a basename, and the
package also treats one of those files as the "build is finished" marker.

Measured 2026-11-28 on the branch, one link per call and an outside file
holding `keep`:

| link at | call | result before the fix |
| ------- | ---- | --------------------- |
| `<base>.fa` | `bowtie_index_build()` | `writeXStringSet()` replaced the outside file with the FASTA |
| `<base>.1.ebwt` | `bowtie_index_build()` | the link was followed by `file_exists()` and the foreign file was reported as a finished index ("Reusing…") |
| `<path>.fa` | `index_write_chrom_fasta()` | the shard FASTA replaced the outside file |
| `<base>.fa` | `marker_fm_index()` | the FASTA replaced the outside file; `buildindex()` then ran against it |
| `<path>.nsq` | `marker_blast_index()` | the link was followed and a foreign database was returned as the genome's index |

The first, second, third and fifth rows are cheap to reproduce and are the ones
the regression tests pin (`test-genome-index.R`).

Two properties separate this case from ADR 0033:

1. **The write is not one call.** Six bowtie files or eight BLAST files are
   created by the engine after the guard has to run. A guard on the FASTA alone
   would miss `<base>.1.ebwt` and the rest.
2. **A linked file is also a read hazard.** `fs::file_exists()` follows a link,
   so a link at the completion marker silently promotes a foreign file to "this
   genome's index". Refusing the write is not enough; the marker must be
   treated as absent first, and then the same guard refuses the rebuild.

## Decision

1. **An index is a declared set.** `declared_index(prefix, suffixes)` refuses a
   link at every name `paste0(prefix, suffixes)`; each builder calls it before
   the engine runs. The suffix sets are the engine's documented output:
   six `.ebwt` files for bowtie, eight BLAST database files, `.files`/`.log`
   for the Subread FM index.
2. **A link is not a finished index.** `cache_present(path)` is
   `file_exists()` **and** no link. The reuse checks use it, so a linked
   `.1.ebwt`, `.nsq` or `.files` does not short-circuit the build; the builder
   then refuses at `declared_index()` with the shared message.
3. **One message, one helper.** All sites use `refuse_links()` through
   `declared_path()`/`declared_index()` and the same sentence ("which is a
   link"), as ADR 0029 and ADR 0033 require.
4. **The caller's input is still the caller's.** `marker_blast_index()` on a
   `marker_genome` indexes the genome's own FASTA in place: only the index
   files beside it are declared, never the FASTA the caller passed, which may
   legitimately be a link. The FM builder does the same for `genome$path`.
5. **A reader follows its input.** `specificity_bowtie()`, `specificity_blast()`
   and `specificity_fm()` read whatever index path they are given; the guard is
   on creation, not on the query.

## Consequences

- A dangling or live link at any index file is a named error before the build,
  and the outside file is untouched. The old failure was silent: a successful
  call that returned a foreign index.
- The guard runs before the engine, so a caller without `Rbowtie`, `rBLAST` or
  BLAST+ still gets the link refusal from the writer that runs first
  (`bowtie_index_build()`, `index_write_chrom_fasta()`), which is why the
  regression tests need no engine.
- An index directory that a caller populated with real files is still reused;
  no mtime or content check is added. Reusing a foreign *regular* file remains
  possible by design — the rule is about links, not provenance.

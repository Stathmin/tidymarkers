# 0037 — The write vocabulary names every primitive that creates a file

- Status: accepted
- Date: 2026-11-28
- Supersedes: none
- Builds on: [ADR 0029](0029-an-output-writer-creates-its-files-not-writes-through-a-link.md)
  (a writer creates its files) and
  [ADR 0036](0036-a-nested-closure-is-a-writer-of-its-own.md) (the unit and the
  call rule)

## Context

The writer gate (ZV45–ZV47) is only as good as its list of write primitives.
That list was built by hand and covered the base-R, readr, fs and xml2
spellings, but nothing made it complete: a primitive that was never added is
invisible, and its unit passes because the gate sees no write at all.

An audit enumerated every call in `R/` whose name matches
`write|save|dump|export|upload|knit|create|copy|move|rename|download|fetch`
and compared the result with `write_primitive`. One file-creating call was
missing: `Biostrings::writeXStringSet()`. It appears eight times in `R/` and
writes a FASTA/sequence file wherever it is given a path:

```r
Biostrings::writeXStringSet(gs, p)   # bowtie, FM index, genome draw, fasta_index
```

With the list as it was, three units wrote a file and were not writers as far
as the gate knew: `bowtie_primer_hits()`, `specificity_fm()` and
`genome_fasta_path()`. All three write into a temp directory, so no production
path was ever unguarded; the gate simply could not have caught it if one had
been.

## Decision

1. **`writeXStringSet` is a write primitive.** It creates a sequence file at
   the path it is given, so it belongs in `write_primitive` beside
   `writeLines()` and `saveRDS()`.
2. **The completeness of the list is testable.** The self-test gains a body
   that calls `Biostrings::writeXStringSet(x, p)` and must be reported as an
   undeclared writer, so deleting the primitive from the list reddens the gate
   instead of silently shrinking what it sees.
3. **The three tempdir writers are registered, not excused.** Each gets a
   ledger entry naming its directory (`withr::local_tempdir()` or
   `fs::path_temp()`), so the fact that the name is not a caller's path is
   written down and is checked for staleness like every other entry.

## Evidence

- Before: `23 writer(s), 17 guard in place, 6 registered`. After:
  `31 writer(s), 22 guard in place, 9 registered`. The eight new hits are the
  `writeXStringSet` call sites; five units already declared their names
  (`bowtie_index_build`, `marker_fm_index`, `index_write_chrom_fasta`,
  `fasta_index`, `marker_blast_index`), three are the tempdir writers above.
- The self-test body `function(p) Biostrings::writeXStringSet(x, p)` is
  reported; with the old list the same body produced no hit.
- No production write changed. `genome_fasta_path()` writes under
  `fs::path_temp()` at a sequence-keyed name; the other two write under
  `withr::local_tempdir()`. None receives a caller's directory.

## Amendment (2026-11-28, night shift)

The audit recipe above enumerated names matching
`write|save|dump|export|upload|knit|create|copy|move|rename|download|fetch`.
That pattern cannot see an archive extractor: `utils::untar()` and
`utils::unzip()` create every file in the archive, and `tar()`/`zip()` create
an archive. None was in `write_primitive`, so a unit that unpacked a caller's
archive without `refuse_links()` would have passed the gate. `R/marker_tools.R`
already guards its extraction (`marker_install_tools()` declares the archive,
the tree, `bin` and `share`), so no live writer was unguarded; the gate could
not have caught a new one.

The four names are added, and the self-test gains
`utils::untar(a, exdir = d)` and `utils::unzip(a, exdir = d)` bodies that must
be reported. Live scan stays `31 writer(s), 23 guard in place, 8 registered`;
the two new bodies exist only in the self-test. The audit recipe is widened to
`untar|unzip|tar|zip` as well, so the next completeness check can see them.

## Consequences

- A new primitive still has to be added by hand, but the omission now shows up
  as a test failure when someone removes one, and the audit recipe (enumerate
  calls by name, diff against the list) is repeatable.
- The scan remains a source heuristic: it proves the unit contains a refusing
  call, not that the refusal runs before the write. A unit that builds an index
  through an engine (`Rbowtie::bowtie_build`, `Rsubread::buildindex`) is
  covered only because `declared_index()` sits in the same unit.
- The three ledger entries are about directories, not files: if one of those
  functions later accepts a caller's path, the entry must be removed and the
  name declared.

## References

- [`dev/check-writers.R`](../../dev/check-writers.R) — the vocabulary and the
  self-test.
- [`dev/writers.yml`](../../dev/writers.yml) — the three tempdir entries.
- [`ROADMAP.md`](../../ROADMAP.md) — `ZV48`.

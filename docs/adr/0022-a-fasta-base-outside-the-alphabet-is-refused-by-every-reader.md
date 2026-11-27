# 0022 — A FASTA base outside the alphabet is refused by every reader

- Status: accepted
- Date: 2026-11-28
- Supersedes: none
- Builds on: [ADR 0021](0021-one-path-and-one-header-rule-for-three-vcf-readers.md)
  (a rule shared by a reader family lives in one helper)

## Context

Two FASTA readers exist for two different jobs: `read_genome_fasta()` reads a
genome, `read_transcripts()` reads transcript models for alignment and
haplotypes. Both use `Biostrings::readDNAStringSet()`.

An earlier repair found that Biostrings does not refuse a letter outside its
alphabet. It warns — `ignored 1 invalid one-letter sequence codes` — and
returns the sequence without that base. Every base after it moves one
position left, so a 12-base transcript comes back as 11 bases and any
coordinate computed on it is off by one. `read_genome_fasta()` was fixed then:
`read_fasta_strict()` turns that warning into a named refusal.

`read_transcripts()` was not fixed. Measured on `>t1 / ACGTXCGTACGT`:

| reader | result |
| --- | --- |
| `read_genome_fasta()` | `rlang_error`: "has bases outside the DNA/IUPAC alphabet; refusing to shift coordinates silently" |
| `read_transcripts()` | 11-base `DNAStringSet`, warning only |

The silent branch is the transcript reader — the one whose output feeds
alignment and haplotype coordinates, where a shifted base is not visible in the
sequence itself.

## Decision

A FASTA file whose base alphabet is not the DNA/IUPAC alphabet is refused by
name, in every reader, through `read_fasta_strict()`.

- `read_genome_fasta()` and `read_transcripts()` both wrap the Biostrings call.
- A reader does not call `readDNAStringSet()` directly on a user path.
- IUPAC ambiguity codes (`R`, `Y`, `S`, `W`, `K`, `M`, `B`, `D`, `H`, `V`) stay
  legal: they are part of the alphabet, not an unknown letter.
- A `.gz` path is still read transparently; the rule is about the letters, not
  the container.
- `fasta_index()` materialises a `.gz` FASTA as a plain file before indexing;
  that copy is written from the strict read, so an invalid base is refused
  before it can be baked into the `.fai`, not only before one read.

## Consequences

- Both Biostrings readers now agree on the same file: the genome reader and the
  transcript reader refuse it with the same sentence.
- A transcript file with a stray `X`/`-`/digit now fails before alignment
  instead of producing a shortened sequence.
- Tests: `test-io.R` — "read_transcripts refuses a base outside the DNA
  alphabet" (`ACGTNXYZ` errors, `ACGTRYSWKM` round-trips) and "fasta_index
  refuses an invalid base before writing the plain copy".
- `ath_genome()` reads its cached TAIR10 download through the same helper, so a
  cached genome with a bad base is refused by name rather than shortened. The
  refusal carries class `tidymarkers_fasta_alphabet`, which lets
  `ath_genome()` tell a corrupt cache (refetch once) from a bad alphabet
  (report the file, never mask it as a download failure). Test:
  `test-ath-genome-cache.R` with `ACGTXCGT` and a mocked download that must not
  run.

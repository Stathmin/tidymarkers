# 0017 — A raw identifier is admitted before a normaliser may rename it

- Status: accepted
- Date: 2026-11-27
- Supersedes: none
- Builds on: [ADR 0016](0016-one-vcf-admission-rule.md) (one VCF admission
  rule), [ADR 0015](0015-a-record-resolves-only-its-own-directory.md) (a
  reader refuses what it cannot represent)

## Context

Several readers normalise an identifier before the rest of the package sees it.
`ath_chrom()` maps `1`, `chr1` and `Chr1` to one spelling and, for anything it
does not recognise, returns `paste0("Chr", key)`. That is the right rule for a
real name: it makes `chr1` and `1` one chromosome. It is the wrong rule for a
missing value.

A VCF record with an empty `CHROM` field therefore became `"ChrNA"` (or `"Chr"`)
— a plausible contig name, not an error. A target could then land on a sequence
that does not exist or, worse, collide with a real contig named `ChrNA`. The
same shape appeared in `as_variant_row()`, where the raw `chrom` was mapped
before admission, and in the index, where a chromosome name was used as a file
name before any check.

The defect is not in `ath_chrom()`. It is the order: normalisation ran before
admission, so the value that reached the check was already invented.

## Decision

Every reader validates the raw identifier first, then normalises it:

- `read_variants_vcf()` — both backends refuse an empty or `NA` `CHROM` before
  `ath_chrom()` (vroom: raw `CHROM` before `transmute`; VariantAnnotation: raw
  `seqnames` before the tibble).
- `read_variant_genotypes()` — the streaming loop refuses the raw first field
  before `ath_chrom()` and names the record. Found in the review of this ADR:
  the reader had kept the pre-fix order and returned the invented contig `Chr`
  for the same file the variant reader refused.
- `as_variant_row()` — the selected row must carry a chromosome and `pos >= 1`
  before `ath_chrom()` and before the `variant` class is set.
- `genome_index()` / `index_shards()` — a build or chromosome name containing a
  path separator or naming `.`/`..` is refused before a path is built.

`ath_chrom()` stays a pure, vectorised mapper. It is used on whole columns and
on already-admitted values; adding an abort inside it would turn an optional
mapping into a global failure and break callers that legitimately pass a
vector.

## Alternatives rejected

- **Make `ath_chrom()` refuse `NA`/empty.** It runs on whole columns and in
  fallback paths; one missing value would abort a read that could simply
  report the record. The mapper is also used for lookups where `NA` is a
  legitimate "no match".
- **Check after normalisation.** `"ChrNA"` is a well-formed string, so a later
  check cannot tell it from a real name. The information is gone.
- **Drop the bad record silently.** A dropped record shifts nothing but hides
  a malformed file; the readers already count and report dropped ALT alleles,
  and a missing chromosome is a harder error than a symbolic allele. The
  genotype reader therefore aborts, as the variant reader does, instead of
  counting a fourth drop class.

## Consequences

- A malformed identifier names the record and stops the read.
- A normalised name can be trusted: it came from a non-empty raw value.
- The rule is repeated in each reader on purpose; a shared helper cannot see
  the raw column after a backend has mapped it.

## Verification

- `tests/testthat/test-io.R`: `read_variants_vcf()` and
  `read_variant_genotypes()` both refuse a VCF record with no chromosome
  (commit `6474452`); the variant row refuses an empty chromosome and
  `pos < 1` (`d1fcd18`); the VCF reader refuses `POS` outside integer range
  (`e5c2c87`).
- `tests/testthat/test-genome-index.R`: build and chromosome names with a path
  separator are refused before a shard path exists.
- `tests/testthat/test-io.R`: a FASTA that names a sequence twice is refused
  (`7ea07f8`) — the same admission-before-rename rule on the sequence side.
- `Rscript dev/test.R "^io$"` — 146 pass, 0 fail; `Rscript dev/check.R --quick
  --no-render` and `Rscript dev/check-docs.R`.

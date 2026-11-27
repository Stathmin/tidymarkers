# 0018 — One allele and position rule for all three VCF readers

- Status: accepted
- Date: 2026-11-27
- Supersedes: none
- Builds on: [ADR 0016](0016-one-vcf-admission-rule.md) (one VCF admission
  rule), [ADR 0017](0017-a-raw-identifier-is-admitted-before-normalisation.md)
  (admit the raw value before a normaliser renames it)

## Context

ADR 0016 introduced one admission rule for the VCF readers and named two of
them. The package has three: `read_variants_vcf()` (the variant table),
`read_variant_genotypes()` (the streaming dosage table) and
`read_variant_depth()` (DP/AD/GT out of one sample, or `INFO/DP` for a
sites-only file). Only the first two were swept for allele admission.

The same file therefore produced different tables depending on which reader
opened it:

- A symbolic or absent ALT (`<DEL>`, `*`, `.`) was dropped and counted by the
  variant and genotype readers, and kept verbatim by both depth paths. The
  depth table carried a `<DEL>` allele that no other table had, and a join on
  `(chrom, pos, ref, alt)` silently dropped those rows.
- The depth reader never normalised `ref`/`alt`. A lower-case `a`/`g` record,
  or an IUPAC code such as `R`, stayed as written, while the other readers
  produced `A`/`G`/`N`. The same locus had two spellings.
- A `POS` that is not a positive integer was refused by the variant and depth
  readers (commit `e5c2c87`), and dropped with a count by the genotype reader.
  The genotype reader's drop was defended in its test on the grounds that a
  row with `pos = NA` would poison comparisons. That argument is satisfied by
  refusal, which the other two readers already chose.

ADR 0017's rule applies unchanged: the raw identifier is admitted before
normalisation. A symbolic allele is not a base string that `normalise_allele()`
could rename; it must be dropped before that call turns it into `N`.

## Decision

One file, one admitted set, in all three readers:

- A symbolic or absent REF or ALT is dropped and counted, with the same
  message class in every reader. The depth reader shares this rule between its
  sample path and its `INFO/DP` path (`depth_admit_alleles()` in
  `R/marker_depth.R`).
- Kept `ref`/`alt` are normalised with `normalise_allele()` (upper case, IUPAC
  folded to `N`), so the join key means one thing.
- A `POS` outside the positive integer range stops the read and names the
  record, in all three readers. No reader returns a partial table for it.
- A one-row selector refuses rather than drops. `as_variant_row()` asks for one
  variant, so a symbolic or absent allele is an error that names the allele:
  there is no second record to fall back to. Before this decision the row
  passed through `normalise_allele()` and became `NNNNN`, the fake `N` the
  readers already refuse to produce. Found by sweeping the non-streaming
  entrance to the same `variant` class.
- Differences remain only where the table shape forces them, as ADR 0016
  states: the genotype reader drops a multi-allelic record (a dosage against
  "the ALT" is ambiguous), a line shorter than the header (nothing to read)
  and a record with neither `GT` nor `AD`; the variant reader splits a
  multi-allelic ALT into one row per allele; the depth reader splits it into
  one row per allele with the matching `AD` index.

## Alternatives rejected

- **Make the depth reader abort on a symbolic allele.** Structural alleles are
  legitimate VCF content (a deletion record in a real file), not a malformed
  file. The other two readers drop and count them; refusing would make the
  depth reader unusable on the same file the variant reader accepts.
- **Keep a symbolic allele as a row and exclude it downstream.** The exclusion
  would have to be repeated at every join, and the depth table would still
  disagree with the variant table on row count.
- **Leave the genotype reader's POS drop.** It was the only reader that
  returned `0` rows with a message where the others stopped, so a caller
  switching readers to save memory silently changed a hard error into an
  empty result.

## Consequences

- A file opened by any reader yields the same `(chrom, pos, ref, alt)` set.
- A caller that joins depth onto variants no longer loses symbolic rows
  silently in one table and keeps them in the other.
- A malformed `POS` is one error class across the package.

## Verification

- `tests/testthat/test-io.R`: one temporary VCF is read by all three readers;
  a `POS` of `3000000000` is refused by each, the genotype reader naming the
  record.
- `tests/testthat/test-depth.R`: a file with `<DEL>`, `*`, lower-case `g` and
  `R` yields a depth table whose `(chrom, pos, ref, alt)` equals the variant
  table's (symbolic dropped, `G`/`N` normalised), on the sample path and the
  `INFO/DP` path.
- `tests/testthat/test-io.R`: `as_variant_row()` refuses `<DEL>`, `*`, `.` and
  `NA` with "concrete REF and ALT" instead of returning `N`.
- `Rscript dev/test.R "^(io|depth)$"` — 184 pass, 0 fail;
  `Rscript dev/test.R "scenario|plan"` — 316 pass, 0 fail.

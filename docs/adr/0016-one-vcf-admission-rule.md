# 0016 — One VCF admission rule, shared by both readers

- Status: accepted
- Date: 2026-11-27
- Supersedes: none
- Builds on: [ADR 0015](0015-a-record-resolves-only-its-own-directory.md) (a
  reader refuses what it cannot represent), commit `6a7b835`

## Context

One VCF has two readers. `read_variants_vcf()` answers *what varies at this
site*; `read_variant_genotypes()` answers *who carries the ALT*. For a join
between the two tables to mean anything, they must agree on which records they
admit.

They did not. The variant reader dropped a symbolic or absent allele
(`<DEL>`, `*`, `.`) and said how many. The genotype reader kept the record and
passed the allele to `normalise_allele()`, which rewrites every non-ACGTN
character to `N`: `<DEL>` became `NNNNN` and a missing ALT became `N`. The same
file therefore produced a genotype row for a site with no variant row, keyed by
a fake allele. A downstream `marker_informativeness()` would treat `N` as a
real allele.

This is not two independent bugs. It is one admission rule written twice, and
only one copy was updated when the rule was learned.

## Decision

1. **One predicate names the rule.** `vcf_symbolic_allele()` in
   `R/marker_io.R` decides whether a REF or ALT field is symbolic (`<DEL>`) or
   absent (`*`, `.`). Both readers call it; neither repeats the pattern.
2. **Drop the record before normalising any allele.** A symbolic or absent
   allele is unrepresentable for a short oligo, and normalisation would turn it
   into a fake `N`. The record is skipped and counted, and the count is named in
   a message.
3. **Readers may differ only where the table shape forces it.** The variant
   reader splits a multi-allelic ALT into one row per allele; the genotype
   reader drops a multi-allelic record because ALT dosage is ambiguous. Both
   differences are documented in the roxygen block of the reader that makes
   them. Everything else about admission is common.
4. **A shared-rule change is tested across both readers.** The test hands one
   VCF to both and expects the same admitted `pos` set, not just the right rows
   from one function.

## Alternatives

- **Keep the duplicated check.** Rejected: the duplication is the defect. The
  next rule learned in one reader would drift again, and no test compares the
  two readers.
- **Keep symbolic rows in the genotype table with `alt = NA`.** Rejected: `alt`
  is part of every join key and of the dosage definition. A row whose ALT
  cannot be named is a row the caller cannot use.
- **Share the whole filter pipeline.** Rejected: the multi-allelic difference
  is intentional, and factoring the pipeline would either move that difference
  into an option or force the variant reader to drop sites it can represent.

## Consequences

- The genotype table can no longer carry `N`/`NNNNN` from a VCF, and its `pos`
  set matches the variant table for every record both can represent. Test:
  `test-io.R`, "read_variant_genotypes drops a symbolic or absent ALT like the
  variant reader".
- A sites-only VCF (no sample columns) still yields an empty genotype table;
  that is an absence of samples, not a dropped record, and stays as it was.
- Residual limit: a record whose REF/ALT normalise to the same base (e.g. an
  InDel written as `A`/`A`) is still admitted by both readers. Neither reader
  validates biological plausibility, only representability.

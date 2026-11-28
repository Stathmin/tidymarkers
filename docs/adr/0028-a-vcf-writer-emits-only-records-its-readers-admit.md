# 0028 — A VCF writer emits only records its readers admit

- Status: accepted
- Date: 2026-11-28
- Supersedes: none
- Builds on: [ADR 0018](0018-one-allele-and-position-rule-for-three-vcf-readers.md)
  (one admission rule for the three readers), [ADR 0019](0019-an-absent-allele-is-not-a-base.md)
  (an absent allele is absent before it is normalised)

## Context

The three VCF readers share one admission rule (ADR 0018–0021): a symbolic or
absent allele is dropped, `POS` outside the positive integer range refuses the
file, an allele outside the DNA alphabet is refused. `write_variants_vcf()` —
the only writer in the family — checked only for `NA` and then formatted the
row. Measured 2026-11-28 on the branch, one call per row:

| row written | file body | what the readers did |
| ----------- | --------- | -------------------- |
| `alt = ""` | empty ALT field | `read_variants_vcf()` returned 0 rows |
| `alt = "<DEL>"` | `<DEL>` | returned 0 rows |
| `pos = 0L` | `POS = 0` | refused the file: `POS outside the positive integer range` |
| `pos = 3e9` | `POS = NA` plus a coercion warning | refused the file |
| `ref = "A\tB"` | `A` and `B` in two columns | read one row with `ref = "A"`, `alt = "N"` — a different allele |
| `depth = 1e10` | `DP=NA` plus a warning | read `depth = NA` |

So the writer could produce a file that its own reader refuses, or a file that
reads back as a different variant. The tab case is the worst: `alt = "G"` became
`alt = "N"` while the `G` landed in `QUAL`, with no warning from the writer.

## Decision

1. **One rule at both ends.** `write_variants_vcf()` applies the readers'
   admission predicate before it formats a row and refuses the first offending
   record by number:
   - symbolic or absent allele (`vcf_symbolic_allele()`) — the readers drop it;
   - a chromosome or allele containing a tab, CR or LF — it would become
     another column;
   - an allele outside `ACGTN` after upper-casing — the readers fold it to `N`,
     which names a base the caller did not write;
   - `POS` outside the positive integer range;
   - a non-missing `depth` outside the non-negative integer range (that is
     what `as.integer(1e10)` turns into `DP=NA`).
2. **One allele per row.** A comma ALT is refused rather than written: the
   readers expand a multi-allelic record into one row per ALT, so a compressed
   row could not read back as one row.
3. **Canonical case on write.** Alleles are upper-cased (`toupper()`), so a
   lowercase table writes the same file the readers return; this is
   normalisation, not a silent reinterpretation of a base.

## Consequences

- A variant table that cannot round-trip fails at `write_variants_vcf()` with
  the record number and the field, instead of leaving a file that drops rows or
  changes an allele. Test: `test-io.R`, "write_variants_vcf refuses a row its
  own readers drop or reinterpret".
- `write_variants_vcf()` now refuses rows it used to write. No caller in the
  package writes an empty, symbolic, lowercase-alias or multi-allelic ALT, so
  no existing design changes; a caller that did was writing a file its own
  readers would not return.
- The writers outside this family (`marker_artefact()` sidecars, RO-Crate) are
  unchanged: their admission rules are the record rules of ADR 0015/0025, not
  the VCF allele rule.
- The test and the writer name the same predicate; a future change to
  `vcf_symbolic_allele()` reaches the writer for free.

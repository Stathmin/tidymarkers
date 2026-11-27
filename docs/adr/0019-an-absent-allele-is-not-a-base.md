# 0019 — An absent allele is absent before it is normalised

- Status: accepted
- Date: 2026-11-28
- Supersedes: none
- Builds on: [ADR 0016](0016-one-vcf-admission-rule.md) (one VCF admission
  rule), [ADR 0017](0017-a-raw-identifier-is-admitted-before-normalisation.md)
  (admit the raw value before a normaliser renames it),
  [ADR 0018](0018-one-allele-and-position-rule-for-three-vcf-readers.md) (one
  allele and `POS` rule for all three readers)

## Context

ADR 0018 closed the divergence between the three VCF readers for the literal
forms of an unrepresentable allele: `<DEL>`, `*` and `.`. The predicate it
introduced, `vcf_symbolic_allele()`, was
`stringr::str_detect(x, "^<|[*.]")`, which is `FALSE` for an empty string and
`NA` for a missing value. An empty `REF` or `ALT` field is neither of the
three literal forms, and the readers handled it three different ways:

- `read_variants_vcf()` computed `drop_sym` as a logical with `NA`, then
  evaluated `if (sum(drop_sym) > 0L)`. That is `if (NA > 0L)` — a base-R
  error, `missing value where TRUE/FALSE needed`, with no naming of the
  record. The exported, canonical reader therefore failed on a file the other
  two opened.
- `read_variant_genotypes()` counted neither an empty `REF` nor an empty
  `ALT` as symbolic and kept the row, producing a dosage table with
  `alt = ""` — a record no oligo can be designed against.
- `read_variant_depth()` dropped the record with the counted message, but
  only because it carried its own patch, `drop[is.na(drop)] <- TRUE`, which
  no other reader had.

The class is the same one ZV15, ZV17 and ZV18 closed: one admission rule
stated in prose, three implementations, and one input shape the prose does
not reach. The third state of an allele field is "there is no base here",
and it arrives as `NA` from `vroom` and as `""` from the line-oriented
reader.

## Decision

The admission predicate is total over the three states of an allele field:

```r
vcf_symbolic_allele <- function(x) {
  is.na(x) | !nzchar(x) | stringr::str_detect(x, "^<|[*.]")
}
```

- A missing (`NA`) or empty (`""`) allele is **absent**, the same class as
  `.` and `*`, and is dropped before `normalise_allele()` runs. It is never
  read as "keep", because a logical `NA` is not a keeper in an `if ()`.
- The counted message is one helper, `inform_absent_alleles(n)`, used by the
  variant reader, the genotype reader and both depth paths. A count of zero
  says nothing; a count above zero names the number and the reason.
- Local patches that duplicated the rule are removed. `read_variant_depth()`
  no longer needs `drop[is.na(drop)] <- TRUE`; the predicate supplies it.
  `as_variant_row()` no longer needs a separate `is.na()` test.

The readers keep the table-shape differences ADR 0016 named (a multi-allelic
record is split by the variant and depth readers and dropped by the genotype
reader, which cannot carry an ambiguous dosage). Only admission is shared, and
it is now shared for every state of the field.

## Alternatives rejected

- **Name the missing allele in the error and stop.** A real VCF can carry a
  record the package cannot design against; the other two readers already
  drop such records with a count, so refusing in one reader reproduces the
  divergence ADR 0018 removed.
- **Keep the depth reader's local NA patch.** Two implementations of one rule
  is the defect. The patch also read as if `NA` were a depth-reader special
  case, when it is the general absent state.
- **Treat `""` as a value and let `normalise_allele()` produce `N`.** That is
  the fake `N` ADR 0017 refuses; an empty allele has no base to fold.
- **Use `is.na()` only.** The genotype reader splits lines on tabs, so an
  empty field is `""`, not `NA`, and would still be kept.

## Consequences

- One file opened by any reader yields one admitted `(chrom, pos, ref, alt)`
  set, for every state of the field: concrete, symbolic, absent and missing.
- A malformed or truncated allele field is a counted drop, not a raw base-R
  error that names no record and no reader.
- One message text; a future reader inherits the rule by calling the
  predicate and the helper.

## Verification

- `tests/testthat/test-io.R`, "an absent REF or ALT is dropped by every
  reader, not a raw failure": one temporary VCF with an empty ALT, an empty
  REF and one valid record is read by all three readers; each keeps only
  `pos = 300`, and the variant table carries no `NA` or empty allele.
  Before the change the test failed at `if (n_dropped > 0L)` with the base-R
  error, and the depth reader returned the record the other two disagreed on.
- `Rscript dev/test.R "^io$"` — 157 pass, 0 fail.
  `Rscript dev/test.R "^depth$"` — 38 pass, 0 fail.
- The same probe outside testthat: an empty ALT and an empty REF each give
  0 rows from `read_variants_vcf()`, `read_variant_genotypes()` and
  `read_variant_depth()`, with the same counted message.

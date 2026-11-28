# 0041 — A genotype is its allele multiset, not its spelling

- Status: accepted
- Date: 2026-11-28
- Supersedes: none
- Builds on: [ADR 0018](0018-one-allele-and-position-rule-for-three-vcf-readers.md)
  (one allele and `POS` rule for the three readers),
  [ADR 0039](0039-a-negative-depth-carries-no-confidence.md) (the weight domain)

## Context

`variant_confidence()` classifies a genotype so that an imbalanced heterozygote
gets the balance penalty. The classifier matched the literal strings `"0/1"`
and `"0|1"` only. A VCF may write the same heterozygous call as `"1/0"` or
`"1|0"`; the two allele orders are the same multiset `{0, 1}`, and
`genotype_call()` already treats them alike (both give dosage `1`).

The other spellings fell through to `"unknown"`, which skips the balance term.
Measured on the branch with `depth = 10`, `alt_depth = 1`:

```r
variant_confidence(tibble::tibble(
  chrom = "Chr1", pos = 1:2, depth = 10L, alt_depth = 1L,
  genotype = c("0/1", "1/0")
))
#> 0/1: genotype_class "het",     weight 0.193, confident FALSE
#> 1/0: genotype_class "unknown", weight 0.964, confident TRUE
```

The same call, with the same read counts, passed or failed the depth gate
depending on allele order. `read_variant_depth()` passes `GT` through
unchanged, so a real VCF reaches this path.

## Decision

Classify the allele multiset of a biallelic call:

- two alleles, both in `{0, 1}`: sum `0` -> `"hom_ref"`, `1` -> `"het"`,
  `2` -> `"hom_alt"`;
- any other GT (missing allele, more than two alleles, an allele other than
  `0`/`1`) stays `"unknown"`, as before;
- order and phase separator (`/` or `|`) do not matter.

The classifier reads the GT string, not the reader's dosage: GT present but
not a valid `0`/`1` call must stay `"unknown"` rather than borrow `AD`'s
answer (ADR 0020). The reader is unchanged.

## Alternatives rejected

- **Normalise `GT` in `read_variant_depth()`.** That would rewrite an
  annotation this package did not produce; the classifier is the consumer of
  the spelling, so the fix belongs there.
- **Derive the class from `alt_fraction` whenever GT is unusual.** That would
  let `AD` answer for a present but invalid GT, which ADR 0020 already
  refuses for dosage.
- **List every spelling in a lookup table.** `1/0`, `1|0`, and future
  separators make the list grow; the multiset rule is the actual contract.

## Consequences

- `0/1`, `1/0`, `0|1` and `1|0` now get one class, one weight and one gate
  verdict.
- `0/2`, `2/2`, `1` and phased multi-allelic calls remain `"unknown"` and keep
  their previous weight.
- No reader, writer or artefact field changes; replay verdicts for records
  written as `0/1` are unchanged.

## Verification

- `tests/testthat/test-depth.R`, "a genotype's class does not depend on allele
  order or phase": four spellings are `"het"` with one weight and none is
  confident at `alt_depth = 1`; `0/0`/`1/1` and the invalid calls keep their
  classes.
- `tests/testthat/test-depth.R`, "read_variant_depth() keeps a 1/0 call as a
  het": the same input through a VCF gives `"het"` and `confident = FALSE`.
- Before the fix the first test produced `"unknown"` and `confident = TRUE`
  for the `1/0` row; `Rscript dev/test.R "^depth$"` — 53 pass, 0 fail.

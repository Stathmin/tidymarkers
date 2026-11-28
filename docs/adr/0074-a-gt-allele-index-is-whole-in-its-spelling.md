# ADR 0074 — A GT allele index is whole in its spelling

- Status: accepted
- Date: 2026-11-28
- Builds on: ADR 0020 (a GT naming a nonexistent allele is not a call), ADR 0038
  (a read count is a count), ADR 0068 (a position is whole in its own spelling)

## Context

`biallelic_gt_dosage()` splits a VCF `GT` value on `/` or `|` and casts each
token with `as.integer()`, then keeps a call only when both tokens are `0` or
`1`:

```r
alleles <- lapply(parts, function(z) suppressWarnings(as.integer(z)))
good <- vapply(alleles, function(a) {
  length(a) == 2L && !anyNA(a) && all(a >= 0L & a <= 1L)
}, logical(1))
```

`as.integer()` truncates and re-parses the text, so a token that names no allele
becomes one this biallelic record does have:

| `GT` | before the fix |
|---|---|
| `0.9/1` | dosage 1 (het) |
| `1.5/1` | dosage 2 (hom alt) |
| `1e0/1` | dosage 1 (het) |

The dosage is the ALT count used by [variant_confidence()] and by the
informativeness panel, so a malformed call became a confident genotype. A
negative token was already rejected by the `a >= 0L` bound (ADR 0020
amendment).

## Decision

- An allele token is admitted only when it is one or more digits
  (`^[0-9]+$`). A token with a sign, a decimal point, an exponent, whitespace or
  any other character makes the whole call `NA`, and `AD` does not answer in its
  place (a present `GT` is authoritative, ADR 0020).
- The existing two-allele, `0`/`1` and overflow bounds stay; the digit test runs
  before the cast, on the spelling the file carried.

## Alternatives rejected

- **Judge the cast value.** `as.integer("1.5")` is `1`, a legal allele index;
  only the text shows that the file named something else.
- **Round the token.** Rounding invents a call the file did not make.
- **Fall back to `AD`.** AD is the evidence for an absent or no-call GT, not for
  a GT the file spelled as a non-allele (ADR 0020).

## Consequences

- `read_variant_genotypes()` reports `NA` genotype for a malformed GT while a
  valid `AD`/`DP` still fills `depth`/`alt_depth`.
- `0/1`, `1|0`, `0/0` and `1/1` behave exactly as before.

## Verification

- `tests/testthat/test-io.R`: a VCF with `0.9/1`, `1.5/1`, `1e0/1` and `-1/0`
  reads genotype `NA` for all four, while `0/1` gives `1` and `1/1` gives `2`,
  each with depth 10.
- `Rscript dev/test.R "^io$"` — 245 pass / 0 fail / 3 skip.

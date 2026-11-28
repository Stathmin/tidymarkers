# 0039 — A negative depth carries no confidence, not negative confidence

- Status: accepted
- Date: 2026-11-28
- Supersedes: none
- Builds on: [ADR 0018](0018-one-allele-and-position-rule-for-three-vcf-readers.md)
  (one allele and `POS` rule for the three readers)

## Context

`variant_confidence()` documents `weight` as a number in `[0, 1]`: the depth
term `1 - exp(-depth / min_depth)` rises with depth, and the balance term only
down-weights an imbalanced heterozygote. The depth term has no floor, so a
negative depth makes it negative:

```r
variant_confidence(tibble::tibble(chrom = "Chr1", pos = 1L, depth = -5))
#> weight = -4.29, confident = FALSE
```

Measured on the branch, `depth = -5` with `min_depth = 3` gave `-4.294`.
A heterozygote row happened to land on `0` because the balance term multiplies
by `pmax(0, ...)`, which is why the defect was invisible in the VCF probe: the
clamp was on the balance axis, not on the depth axis. A negative `DP` is
admissible input to the function (it is a public numeric column), so the
documented domain was not closed.

## Decision

The depth term is floored at zero:

```r
depth_weight <- pmax(0, 1 - exp(-d / min_depth))
```

- A negative depth means "no reads behind this call", which is exactly weight
  `0`, not a penalty smaller than no confidence.
- `confident` already requires `depth > 0`, so the floor changes no gate
  decision; it closes the returned value's domain.
- `depth = 0` stays `0`; positive depths are unchanged.

The reader stays raw: `read_variant_depth()` reports the `DP` value it was
given, and `NA` for an absent field. This decision is about the derived weight,
not about rewriting the annotation (contrast ADR 0038, where the reader's `AD`
fallback computes a dosage and therefore had to admit only writer-valid
counts).

## Alternatives rejected

- **Refuse a negative depth with an error.** The function's contract is to
  weight a table, not to be a validator; a negative count in real input should
  make the row unconfident, not abort a whole pipeline.
- **Set `weight <- NA` for a negative depth.** That loses the distinction
  between "no depth recorded" and "a recorded count that cannot be true"; the
  gate treats both as unconfident, and `0` is the truthful weight.

## Consequences

- `weight` is in `[0, 1]` for every finite numeric `depth`, as documented.
- `confident` is unchanged for every input: negative and zero depths remain
  `FALSE`.
- No reader or writer changes; `read_variant_depth()` still reports the raw
  `DP` and `AD` values.

## Verification

- `tests/testthat/test-depth.R`, "a negative depth weighs nothing, not less
  than nothing": depths `-5` (homozygote and heterozygote) and `0` give weight
  `0` and `confident = FALSE`; depth `3` gives `1 - exp(-1)` and is confident.
- Before the fix the same homozygote row returned `-4.294`; the new assertion
  fails on it.
- `Rscript dev/test.R "^depth$"` — 46 pass, 0 fail.

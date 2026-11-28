# ADR 0076 — A panel variant table position is whole in its spelling

- Status: accepted
- Date: 2026-11-28
- Builds on: ADR 0045 (a position is a whole number before the cast), ADR 0068
  (a position is whole in its own spelling), ADR 0075 (a fixed-input position is
  whole before the cast)

## Context

`marker_panel()` accepts a variant table as well as a VCF path. `panel_variants()`
returned a caller's table unchanged, and the haplotype contract keyed it with

```r
vk <- paste(v$chrom, as.integer(v$pos), toupper(v$ref), toupper(v$alt), sep = "\r")
```

`as.integer()` truncates, so a table row at `20.9` produced the key of `20`:

| variant table row | haplotype row | before the fix |
|---|---|---|
| `Chr1:20.9 A>G` | `Chr1:20 A>G` | "present in the table" |

The contract exists to prove every non-target haplotype variant is in the
supplied table; a shifted key let a variant the table does not carry pass the
check, and the panel then designed with it.

## Decision

- `panel_variants()` validates a table input once, at the seam, before any
  consumer: it must be a data frame with `chrom`/`pos`/`ref`/`alt`, and every
  `pos` must pass `bad_position()`. A fractional, `NA`, `Inf`, non-numeric or
  non-positive position is refused with "`variants` needs a positive whole
  `pos` on every row."
- A character path still goes through `read_variants_vcf()`, whose POS reader
  already refuses a non-whole spelling (ZV94).

## Alternatives rejected

- **Check the key after the cast.** `as.integer(20.9)` is `20`, a legal base;
  only the spelling shows the table named another one.
- **Normalise by rounding.** A rounded key matches a variant the table does not
  carry, which is the defect.
- **Validate only inside `panel_haplotype_contract()`.** The table is also read
  when no haplotypes are given; the seam covers every consumer.

## Consequences

- `marker_panel(variants = )` and `marker_design_panel(variants = )` use the
  table the caller wrote, or refuse it before LD, haplotypes or design.
- A VCF path keeps the reader's own refusal; integer tables behave as before.

## Verification

- `tests/testthat/test-panel.R`: a table at `100.9` and at `"100.9"` is refused
  with "positive whole `pos`"; the integer table is returned unchanged.
- `Rscript dev/test.R "^panel$"` — 98 pass / 0 fail.

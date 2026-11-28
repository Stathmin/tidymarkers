# ADR 0075 — A fixed-input position is whole before the cast

- Status: accepted
- Date: 2026-11-28
- Builds on: ADR 0045 (a position is a whole number before the cast), ADR 0068
  (a position is whole in its own spelling), ADR 0071 (a region bound is a
  coordinate before the cast)

## Context

Two normalisers read a one-row or table input, cast `pos` with `as.integer()`
and let the result through:

```r
# normalise_fixed_snp()
pos = as.integer(snp$pos[[1L]])
...
if (is.na(out$pos) || !nzchar(out$chrom)) abort_msg(...)

# as_variant_table()
pos = as.integer(x$pos)
```

Neither checks the spelling or the sign, so a fractional or negative position
became a plausible base:

| input | before the fix |
|---|---|
| `assemble_primer_assay(snp = <pos 10.9>)` | SNP at 10 |
| `assemble_primer_assay(snp = <pos -5>)` | SNP at -5 |
| `marker_chunk(<pos 10.9>, reference = )` | chunk variant at 10 |

`normalise_fixed_snp()` is reached by `assemble_primer_assay()`;
`as_variant_table()` by `marker_chunk.data.frame()`. The design that follows
reads the genome at the truncated base and can build an assay for a locus the
caller never named, while `-5` survives until a later window check names a
different argument.

## Decision

- `normalise_fixed_snp()` checks `bad_position(snp$pos[[1L]])` before the cast
  and refuses a fractional, `NA`, `Inf`, non-numeric or non-positive position
  with its own SNP message.
- `as_variant_table()` checks `bad_position(x$pos)` on the whole column before
  the cast and refuses any bad row with a chunk-table message.
- Both keep their existing chromosome/allele checks and their accepted integer
  inputs unchanged.

## Alternatives rejected

- **Check the cast value.** `as.integer(10.9)` is `10`, a legal base; only the
  spelling shows the caller named another one.
- **Round the position.** A rounded SNP or variant is a different locus.
- **Rely on the later window check.** `-5` and `10.9` were caught only by an
  engine message naming another argument, or not at all.

## Consequences

- `assemble_primer_assay(snp = )` and `marker_chunk(<table>, reference = )`
  read the locus the caller wrote, or refuse it before any engine or sequence
  lookup.
- Integer SNP and variant tables behave exactly as before.

## Verification

- `tests/testthat/test-marker-assemble.R`: a SNP at `pos + 0.5` and at `-5` is
  refused with "positive whole-number position".
- `tests/testthat/test-alignment.R`: a chunk table at `10.9` and `-5` is refused
  with "positive whole `pos`"; the integer table still returns `c(10L, 20L)`.
- `Rscript dev/test.R "marker-assemble|alignment"` — 72 pass / 0 fail.

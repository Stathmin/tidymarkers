# ADR 0070 — A selector is resolved outside the data mask

- Status: accepted
- Date: 2026-11-28
- Builds on: ADR 0045 (a position is a whole number before the cast), ADR 0047
  (a radius is checked before the cast), ADR 0050 (a coordinate is checked
  before the read), ADR 0068 (a selector keeps a position whole)

## Context

`read_variants_vcf(path, chrom =, start =, end =)` filtered an already-read
table:

```r
if (!is.null(chrom)) v <- dplyr::filter(v, .data$chrom == ath_chrom(chrom))
if (!is.null(start)) v <- dplyr::filter(v, .data$pos >= start)
if (!is.null(end))   v <- dplyr::filter(v, .data$pos <= end)
```

`dplyr::filter()` evaluates its second argument in a data mask, so the bare name
`chrom` in `ath_chrom(chrom)` resolved to the **column**, not to the function
argument. `ath_chrom()` mapped every column value to its own canonical name, the
comparison was true for every row, and the chromosome filter did nothing:

| input | before the fix |
|---|---|
| `read_variants_vcf(vcf, chrom = "Chr9")` on a Chr1/Chr2 file | both rows |
| `read_variants_vcf(vcf, chrom = "Chr1")` | both rows (coincidence) |
| `read_variants_vcf(vcf, chrom = NA)` | zero rows, no message |

The existing region test passed `chrom`, `start` and `end` together, so the
position window hid the no-op: it asserted only the one row inside the window.
The same call also accepted a selector that is not a position — `start = NA`,
`start = "abc"` and `end = 50.5` returned zero rows or a window nobody named, and
`start > end` returned an empty table instead of naming the inverted window.
ADR 0053 fixed the sibling case of the same class: a sample selector names a
row, it does not round.

## Decision

- Resolve `chrom` to `chrom_key <- ath_chrom(as.character(chrom))` **before**
  the filter, and compare the column with that ordinary value. The data mask
  never sees the argument name.
- Refuse a selector that cannot name one region: `chrom` must be one non-empty,
  non-`NA` name; `start` and `end` go through `check_coord()` (ADR 0050), so a
  fraction, `NA`, `Inf`, zero or an out-of-range value is refused instead of
  silently read; `start > end` is refused.
- Both backends keep the same post-filter. `VariantAnnotation`'s
  `ScanVcfParam(which=)` remains a prefilter only when all three selectors are
  given; it is a performance path, not the contract.

## Alternatives rejected

- **`!!chrom` or `.env$chrom` inside the filter.** It works but depends on the
  reader knowing the masking rule at the call site; the resolved key is the
  same value the prefilter uses and cannot be re-masked.
- **Rename the argument to `chromosome`.** It changes the public surface to
  avoid an evaluation rule; the argument name stays.
- **Validate only `start`/`end`.** The no-op filter is the user-visible wrong
  answer; leaving it would keep `chrom =` documented but ineffective.

## Consequences

- `read_variants_vcf(chrom = )` now selects the named chromosome on both
  backends; an unknown name returns zero rows because it is a valid name that
  matches nothing, while `NA`, `""` and a vector are refused.
- A window bound is a whole positive position or a refusal; an inverted window
  names both bounds in the message.
- No accepted input changes its rows: a call that passed `chrom` together with
  `start`/`end` already got the window's rows, and the chromosome predicate is
  now the one the caller asked for.

## Verification

- `tests/testthat/test-io.R`: "read_variants_vcf filters by chromosome alone" —
  a Chr1/Chr2 file returns one row for `chrom = "Chr2"`, zero for `"Chr9"`, and
  refuses `chrom = NA`; "read_variants_vcf refuses a window bound that is not a
  position" — `start = NA`, `start = 20.5`, `end = "abc"` and `start > end` are
  refused, `start = 20L` keeps the one row.
- Before the fix the seven assertions failed (`^io$`: 236 pass / 7 fail); after
  it `Rscript dev/test.R '^io$'` is 243 pass / 0 fail.

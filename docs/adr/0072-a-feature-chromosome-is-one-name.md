# ADR 0072 — A feature chromosome is one name

- Status: accepted
- Date: 2026-11-28
- Builds on: ADR 0070 (a selector is resolved outside the data mask), ADR 0071
  (a region bound is a coordinate before the cast)

## Context

`marker_features()` filters the feature table with a base comparison:

```r
if (!is.null(chrom)) {
  out <- out[!is.na(out$chrom) & out$chrom == as.character(chrom), , drop = FALSE]
}
```

`chrom` was never checked to be one name. The comparison is vectorised, so the
caller's vector was recycled against the table instead of selecting with it:

| input | before the fix |
|---|---|
| `chrom = NA` | all `NA` in the mask, zero rows, no message |
| `chrom = c("Chr4", "Chr1")` on a two-row table | row 1 kept, row 2 dropped |
| `chrom = ""` | zero rows, no message |

The caller asked for two chromosomes (or named none) and got a table that
silently answered a different question. The reader cannot tell an empty locus
from a bad selector.

## Decision

- A non-`NULL` `chrom` must be one non-empty, non-`NA` name, checked before the
  comparison, with the same message as `read_variants_vcf()`:
  `` `chrom` must be one chromosome name. ``
- The `chrom_map` translation and the comparison itself are unchanged; a name
  absent from the file still returns zero rows, because that is a real answer.

## Alternatives rejected

- **Filter with `%in%`.** It would accept a vector selector, which the rest of
  the input family does not: `read_variants_vcf(chrom = )` takes one name.
- **Treat `NA` as "no filter".** A missing chromosome is a mistyped selector,
  not a request for every chromosome; returning the whole file hides the bug.
- **Warn and drop the `NA` rows.** The table returned would still not answer
  the question the caller wrote.

## Consequences

- `marker_features(chrom = )` reads one chromosome or refuses the selector.
- `chrom = NULL` and a name absent from the file behave exactly as before.

## Verification

- `tests/testthat/test-annotate.R`: `chrom = NA`, `chrom = c("Chr4", "Chr1")`
  and `chrom = ""` are refused with "one chromosome name"; the existing
  coordinate and default-filter tests are unchanged.
- `Rscript dev/test.R '^annotate$'` is 13 pass / 0 fail.

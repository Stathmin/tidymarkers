# ADR 0084 — Annotation is asked about the rows and bases the caller named

- Status: accepted
- Date: 2026-11-28
- Builds on: ADR 0045 (a position is a whole number before the cast),
  ADR 0038 (a count the writer refuses is reported as absent)

## Context

`annotate_variants()` hands a variant table to a transcript model through one
`GRanges`, and documents its `gene` argument as "a gene filter (`$chrom`,
`$start`, `$end`) to restrict the annotation to one locus". Two failures were
measured in that function:

1. **The cast was delegated to `IRanges`.** The ranges were built as
   `IRanges::IRanges(variants$pos, width = nchar(refs))`, with `pos` never
   judged. `IRanges(269745.5, width = 1)` is `269745-269745` and emits no
   warning, so the model read the base at 269745 while the returned table still
   said `269745.5`: the consequence and `gene_id` belonged to a locus the table
   did not name. The package's whole ZV family exists to refuse exactly this,
   and this call bypassed the guard by handing the cast to Bioconductor.
2. **The documented filter was not applied.** The code compared
   `gene$chrom` only. A caller passing `start`/`end` got every variant on the
   chromosome annotated, with `gene_id` values from other genes — silently more
   than the documentation promised, and more work for `locateVariants()`.

## Decision

1. `annotate_variants()` judges `variants$pos` with `bad_position()` as its
   first act, before the `VariantAnnotation` loader and the `TxDb` load, and
   refuses with the row index and the value. The refusal therefore needs no
   annotation model at all.
2. The cast is written as `as.integer(variants$pos)` so the cast gate sees the
   guard it names, instead of leaving the conversion inside `IRanges()`.
3. The `gene` filter applies both bounds: `chrom` selects the chromosome and
   `pos >= gene$start`, `pos <= gene$end` select the locus. A bound that is
   given must be a positive whole number; a missing bound is open (1 and `Inf`).
4. The filter and the empty-table return run before the loader, so a table with
   no row inside the locus returns an empty table instead of requiring the
   annotation packages.

## Rejected alternatives

- Judge `pos` after the loader: the check is free and its message is about the
  input, so it must not depend on a Suggested package being installed.
- Filter after `locateVariants()`: the model would still be asked about rows the
  caller excluded, and the cost of a whole-chromosome annotation is the reason
  the filter exists.
- Fix the documentation instead of the filter: the argument is named `gene` and
  documents a locus; a chromosome-wide annotation is not a smaller version of
  that promise, it is a different result. The doc was the specification.

## Consequences

- A fractional `pos` is refused before any annotation work, with the row index.
- `gene = list(chrom = ..., start = ..., end = ...)` narrows the result to the
  locus; `gene = list(chrom = ...)` still narrows to the chromosome only.
- A table whose rows all fall outside the locus returns unchanged apart from
  the filter (empty), and does not load `VariantAnnotation`.
- The full annotation tests still run only with `dev/test.R --full`; the fast
  suite covers the refusals and the filter, which are the parts that do not need
  a `TxDb`.

## Verification

`Rscript dev/test.R '^annotate$'` — 17 pass, 0 fail, 5 skip (the `TxDb` tests),
including three new cases: a fractional `pos`, a textual fractional `pos`, and
a locus filter that keeps no row plus a fractional `gene$start`.
`Rscript dev/check-casts.R` — `53 position/threshold cast(s), 22 guarded in
place, 31 registered with a reason; self-test passed`.

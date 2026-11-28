# 0053 — A sample selector names a row, it does not round

- Status: accepted
- Date: 2026-11-28
- Supersedes: none
- Builds on: [ADR 0044](0044-a-count-has-one-definition.md) (one
  `is_whole_number()` predicate backs every count), [ADR 0045](0045-a-position-is-a-whole-number-before-the-cast.md)
  (a position is whole before the cast), [ADR 0051](0051-a-nucleotide-length-is-whole-before-the-design-grid.md)

## Context

`read_variant_depth(path, sample = )` takes a column name or an index into the
VCF's sample columns. `resolve_depth_sample()` cast the numeric form with
`as.integer()` and only then checked the range:

```r
if (is.numeric(sample)) {
  i <- as.integer(sample)                       # 1.9 -> 1
  if (length(i) != 1L || is.na(i) || i < 1L || i > length(sample_names)) {
    abort_msg("`sample` must be a sample index (1..N) or a column name in the VCF.")
  }
  return(sample_names[[i]])
}
```

So `sample = 1.9` on a two-sample VCF returned the depth of **S1** — a
different individual from the one the caller asked for, with no warning. The
range check could not catch it because truncation happens before the check.

## Decision

The numeric branch admits a whole index through `is_whole_number()` before the
cast; the existing range check and message stay as they are. `2` and `"S2"`
still resolve to the same column; `1.9`, `2.5`, `0`, `-1`, `Inf` and a
two-element index are refused with "`sample` must be a sample index (1..N) or a
column name in the VCF."

## Alternatives rejected

- **Round.** A rounded row is still a row the caller did not name.
- **Let vroom select by name only.** The index form is the documented default
  (`1L`) and is what a caller with no sample names in hand uses.
- **Refuse a numeric index entirely.** That removes a documented, working form
  to fix a fractional one.

## Consequences

- A fractional, zero, negative, `Inf` or vector index is refused before any
  column is read.
- A whole index and a column name are unchanged.

## Verification

- `tests/testthat/test-depth.R`, "read_variant_depth refuses a sample index it
  would have truncated": on a two-sample VCF, `1.9`, `2.5`, `0`, `-1`, `Inf`
  and `c(1, 2)` are refused; `sample = 2L` and `sample = "S2"` both read S2's
  depth (20).
- `Rscript dev/test.R "^depth$"` — 61 pass, 0 fail.

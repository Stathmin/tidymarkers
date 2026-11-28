# 0048 — A base-pair threshold is a whole number before the cast

- Status: accepted
- Date: 2026-11-28
- Supersedes: none
- Builds on: [ADR 0042](0042-a-threshold-is-a-count-not-a-number-to-truncate.md)
  (a threshold is a count, not a number to truncate), [ADR 0044](0044-a-count-has-one-definition.md)
  (one `is_whole_number()` predicate backs every count), [ADR 0047](0047-a-search-radius-is-a-whole-number-before-the-cast.md)
  (the same defect in a search radius)

## Context

Three design paths cast a caller-supplied number of base pairs with
`as.integer()` before any check, and a fourth does the same before passing it
to `snp_template()`:

| path | cast before the check | what moves |
|---|---|---|
| `design_sanger_marker()` | `min_dist` in the `m_dist` margin | 50.9 becomes 50, while `sanger_seq_primers()` compares the raw 50.9 against the same read-start bar |
| `design_indel_size_marker()` | `min_gap` in the gap floor and in `m_gap` | 5.9 becomes 5: the resolvable gap floor drops below the bar the caller named |
| `snp_template()` | `flank` on both window ends | 100.9 becomes a 100 bp half-window |
| `assemble_primer_assay()` | `flank` before `snp_template()` widens it | 150.9 becomes 150, and the template the assay locates its primers in is 0.9 bp narrower |

Measured on the branch before the fix:

```r
as.integer(50.9)   # 50
as.integer(5.9)    # 5
as.integer(100.9)  # 100
```

`min_dist` was read twice by two rules in the same call: the selection margin
used the truncated 50, while `sanger_seq_primers()` used the raw 50.9. A
candidate at distance 50.5 therefore scored as if it cleared 50 and then
failed the usable-window check — the design and its report disagreed about one
argument. `min_gap = 5.9` lowered the `gate_indel_size_gap` floor from 5.9 to
5, so an indel the caller called unresolvable was accepted as resolvable. A
fractional `flank` silently narrowed the template on which every designer
searches.

## Decision

`check_base_pairs(x, arg)` in `R/marker_validate.R` is the one gate:

```r
check_base_pairs <- function(x, arg) {
  if (!is_whole_number(x, min = 1L)) {
    abort_msg("`{arg}` must be one positive whole number of base pairs.")
  }
  as.integer(x)
}
```

- `design_sanger_marker()` checks `min_dist` before it builds the grid.
- `design_indel_size_marker()` checks `min_gap` before it reads the alleles.
- `snp_template()` checks `flank` before computing the window; every designer
  that reaches the template through it is covered by the same call.
- `assemble_primer_assay()` checks `flank` before widening it for the amplicon,
  so the fixed-primer constructor reaches `snp_template()` with a whole value.

The gate returns the integer, so no path keeps a second cast.

## Alternatives rejected

- **Round the threshold instead of refusing it.** `5.9` would test a bar the
  caller never named; ADR 0042 rejects the same move for a count.
- **Reuse `check_scan_radius()`.** Its name and its `min = 0L` default describe
  a mismatch radius; a base-pair threshold is positive and needs its own
  message. The shared predicate is still `is_whole_number()`.
- **Validate each designer only.** `snp_template()` is exported and is the
  single window shared by every scheme; leaving its cast unguarded would keep
  the defect on the path all designers take.

## Consequences

- A fractional, zero, negative, `Inf`, vector or character base-pair threshold
  is refused with the package message before any grid or sequence work.
- Every accepted whole threshold behaves as before; no accepted design changes.
- `design_sanger_marker()` now uses one value for the margin and the usable
  window, so its selection and its report cannot disagree.

## Verification

- `tests/testthat/test-sanger.R`: `min_dist = 50.9` and `flank = 100.5` abort
  with "positive whole number of base pairs".
- `tests/testthat/test-indel-size.R`: `min_gap = 5.9` and `flank = 250.5` abort
  with the same phrase.
- `tests/testthat/test-genome.R`: `snp_template(flank = 10.9)` aborts, `10L`
  still returns a 21 bp window.
- `tests/testthat/test-marker-assemble.R`: `assemble_primer_assay(flank =
  150.5)` aborts before the template is built.
- `Rscript dev/test.R "sanger|indel-size|^genome$"` — 243 pass, 0 fail, 1 skip;
  `Rscript dev/test.R "^marker-assemble$"` — 8 pass, 0 fail.

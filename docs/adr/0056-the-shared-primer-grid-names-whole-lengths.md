# 0056 — The shared primer grid names whole lengths

- Status: accepted
- Date: 2026-11-28
- Supersedes: none
- Builds on: [ADR 0055](0055-a-common-primer-search-names-whole-lengths-gaps-and-anchor.md)
  (a common-primer search names whole lengths, gaps and anchor), [ADR 0051](0051-a-nucleotide-length-is-whole-before-the-design-grid.md)
  (a nucleotide length is whole before the design grid)

## Context

`design_hrm_marker()` exposes `lengths` and hands it to `hrm_flanking_pair()`,
which builds both flanking sets with the internal `primer_candidates()` grid.
That grid did not check the value:

```r
d <- tidyr::crossing(start = starts, len = lengths) |>
  dplyr::mutate(
    primer = stringr::str_sub(s, .data$start, .data$start + .data$len - 1L),
    end = .data$start + .data$len - 1L
  )
```

`stringr::str_sub()` truncates a fractional end, but the row keeps the
fraction. Measured on the branch before the fix (`snp_fixture()`, 600 bp,
seed 2, SNP at 300, `flank = 150L`):

```r
design_hrm_marker(f$genome, f$snp, flank = 150L, lengths = 20.9)
#> forward/reverse are 20 nt, but the pair reports product 86.9 and end 200.9
```

The integer call returns the same 20 nt pair with `product` 86 and `end` 199.
The row therefore described an amplicon one base longer than the one the
primers form — the same class as ZV78 (`design_as_primers()`) and ZV79
(`design_common_primers()`), one level lower in the stack.

`primer_candidates()` has many internal callers (`marker_hrm.R`,
`marker_snapshot.R`, `marker_arms.R`, `marker_sanger.R`, `marker_dpcr.R`,
`marker_primer_pairs.R`, `marker_castpcr.R`). Each passes lengths that are
whole by construction or that come from an exported designer's `lengths`
argument (`design_hrm_marker()`, `design_sanger_marker()`,
`design_castpcr_*()`, `design_dpcr_marker()`, the snapshot and primer-pair
designers). `design_lamp_marker()` reaches the same grid for both its
`outer_lengths` and its `arm_lengths` (the latter through `lamp_arm_sites()`
→ `safe_candidates()`), but `safe_candidates()` turns the grid's refusal into
an empty set, which is why ADR 0057 also checks both at that entry.

## Decision

`primer_candidates()` checks `lengths` before it builds the grid:

```r
lengths <- check_nucleotide_length(lengths, "lengths", multi = TRUE)
```

Every exported designer that passes a user length argument into this grid
inherits one refusal instead of a per-designer copy. The message names the
grid's parameter (`lengths`), which is also the public argument name for
every such designer except `design_lamp_marker()`'s `outer_lengths`.

## Alternatives rejected

- **Cast with `as.integer()` inside the grid.** The table would still be
  built from a length the caller did not name, and the fraction would vanish
  silently.
- **Check `lengths` in each exported designer.** That duplicates the same
  guard at five or more entries and leaves the next internal caller free to
  reintroduce the truncation.
- **Check `starts` too.** Every current `starts` is a `seq_len()`/`seq.int()`
  result or a whole position computed by the caller; the measured defect is
  the length. A start guard belongs with the caller that owns the coordinate
  (ADR 0050), not with the length grid.

## Consequences

- A fractional, zero, `Inf`, vector or non-numeric `lengths` is refused with
  the package message before any primer is built.
- Whole values keep the same candidates; `len`, `end` and every product
  computed from `end` now match the molecule.
- `design_lamp_marker()`'s `arm_lengths` and `outer_lengths` reach this grid,
  but `safe_candidates()` swallows the refusal; ADR 0057 checks both names at
  the entry so the caller sees the message.

## Verification

- `tests/testthat/test-hrm.R`, "a fractional primer length is refused before
  the primer grid": `20.9` and `c(20L, Inf)` are refused with "positive whole
  numbers of nucleotides"; the integer call still returns `product` 86 and
  `end` 199.
- `dev/test.R "hrm|sanger|castpcr|dpcr|primer-pairs|snapshot|arms|lamp|taqman|rpa"`
  — 2008 pass, 0 fail, 5 skip.

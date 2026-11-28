# ADR 0079 — The cast gate judges an interval, not only a point

- Status: accepted
- Date: 2026-11-28
- Builds on: ADR 0045 (a position is a whole number before the cast), ADR 0068
  (a position is whole in its own spelling), ADR 0078 (a plotted feature
  interval is whole in its spelling)

## Context

`dev/check-casts.R` requires every position-like `as.integer()` cast in `R/` to
be guarded in the same body or registered in `dev/casts.yml` with a reason. The
scanner named a position with `pos`, `POS`, `pos_raw`, `positions` and `stab_at`
only. `start` and `end` were left out on the argument that "a range is built
from a position that its unit already judged".

ZV106 falsified that argument: `normalise_features()` accepted a caller's
`features` table and cast both bounds raw, so `start = 10.9` drew the gene from
10. The defect was invisible to the gate for exactly the excluded name, and it
was found by a manual audit, not by the check.

The same exclusion also hid the `start`/`end` casts in `normalise_region()`,
`marker_from_granges()` and the two Ensembl readers. Those turned out to be
benign — each interval is either a `GRanges` bound or validated upstream — but
"benign" was asserted in prose, not checked by the gate.

## Decision

- `position_name` includes `start` and `end`. An interval bound comes from a
  file, a caller's table or a pointer exactly like `pos`; the gate judges it the
  same way.
- Every existing unguarded `start`/`end` cast is registered in `dev/casts.yml`
  with a reason that names where its integer-ness comes from (a `GRanges`
  bound, `parse_location()`, `check_coord()` at the exported seam, a designer's
  geometry, or an external engine's output). The ledger is checked, not
  remembered: a stale entry still fails the gate.
- The gate's self-test gains an unguarded interval cast and a guarded interval
  cast, so a scanner that quietly stopped matching `start`/`end` is caught.

## Rejected alternatives

- Keep `start`/`end` out and rely on audits: ZV106 is the measured cost of that
  policy — one real defect missed by the gate.
- Flag only `start`/`end` casts whose receiver is a data-frame column: the
  scanner is syntactic and cannot know a receiver's class; a name-based rule
  with a justified ledger is the same trade the other casts already use.
- Guard every flagged site in code: a `GRanges` bound is already whole, and a
  second check there would claim a refusal the type makes impossible.

## Consequences

- The gate now covers 46 position casts (was 33), 16 guarded in place and 30
  registered; the 12 new entries each name the upstream guarantee.
- A future unguarded interval cast fails `Rscript dev/check.R`, not a manual
  audit.
- No runtime behaviour changes: the registered casts were already safe.

## Verification

`Rscript dev/check-casts.R` — `46 position cast(s), 16 guarded in place, 30
registered with a reason; self-test passed`, exit 0. Removing one ledger entry
fails the gate with `unguarded position cast ...: <symbol>`, and the new
self-test cases fail if `start`/`end` stop matching.

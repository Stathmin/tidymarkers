# ADR 0081 — A threshold is whole in its spelling before the cast

- Status: accepted
- Date: 2026-11-28
- Builds on: ADR 0048 (a base-pair threshold is checked before the cast),
  ADR 0079 (the cast gate judges an interval too), ADR 0068 (a position is
  whole in its spelling)

## Context

ADR 0048 fixed three design gates (`min_dist`, `min_gap`, `flank`) that cast a
base-pair threshold with `as.integer()` before any check, so `50.9` narrowed
the gate the caller named. The cast gate (ADR 0079) was built to make that
roster checked: a top-level function that casts a position-like expression must
guard it or be registered with a reason.

Both the roster and the gate had the same blind spot. `design_arms_marker()`
takes `min_band_gap`, a base-pair threshold for the tetra-primer ARMS gel, and
wrote `limits$band_gap <- as.numeric(min_band_gap)`. Measured on an 800 bp
target with the package's example seed:

- `min_band_gap = 20.9` accepted the call and set the arms floor to
  `10 - 20.9` bp, so a 20 bp band gap the documented default accepts was
  reported unresolved;
- `min_band_gap = Inf` was accepted and made the floor `-Inf`;
- `min_band_gap = NA` reached an `if` and failed with base R's
  `missing value where TRUE/FALSE needed`, naming no argument.

`dev/check-casts.R` did not see it for two reasons: `band_gap` is not in
`position_name`, and the cast was `as.numeric()`, not `as.integer()`. A gate
whose rule is "only a variable named `pos` cast with `as.integer()`" cannot
enforce ADR 0048; it enforces one spelling of the rule.

## Decision

1. `design_arms_marker()` passes `min_band_gap` through `check_base_pairs()`
   before it reaches the floor, the same guard `snp_template()` uses for
   `flank` and the design paths use for `min_dist`/`min_gap` (ADR 0048). A
   fractional, zero, negative, `NA`, `Inf`, vector or non-numeric value is
   refused by argument name before any primer is designed.
2. `dev/check-casts.R` treats a **direct** cast of an argument named
   `min_*`/`max_*` as a cast of the same kind, for `as.integer()`,
   `as.numeric()` and `as.double()`. `check_base_pairs()` and
   `check_worker_count()` count as guards.
3. Only a direct cast is reported. `as.numeric(x$start <= pos)` casts the
   comparison's logical result; a threshold mentioned inside it is not a
   threshold being narrowed. `marker_generic_checks()` has four such casts and
   must not need a registry entry.
4. Two callers whose guard lives in another call are registered in
   `dev/casts.yml` with the reason: `marker_ingest_literature()` (validated by
   `ingest_check_criteria()`) and the nested `has_palindrome()` scanner (the
   enclosing function checks all three counts first).

## Rejected alternatives

- Fix only `design_arms_marker()`: the next threshold cast with a different
  name or `as.numeric()` would be invisible again, exactly as `min_band_gap`
  was invisible before.
- Treat every `as.numeric()` on an expression mentioning `pos` as a position
  cast: `as.numeric(pr$start <= marker$snp$pos)` is a logical-to-numeric
  conversion, not a coordinate; it produced four false positives in
  `marker_generic_checks()` alone.
- Widen the rule to every numeric cast: a Tm, a fraction or a cost is
  legitimately continuous, and a gate that needs a registry entry for each one
  stops being read.
- Accept a fractional threshold and compare with it: `20.9` is not a number of
  base pairs, and the arms floor then differs from the documented integer by a
  fraction no caller asked for.

## Consequences

- `min_band_gap = 20L` keeps the documented floor; every non-whole spelling is
  refused by name before the design starts.
- The cast gate reports 52 casts (20 guarded, 32 registered) instead of 46
  (16/30); the two new registrations name where the guard actually lives.
- A future `min_*`/`max_*` cast with `as.numeric()` is a red gate until it is
  guarded or registered, including a threshold in a module the position roster
  does not mention.
- The gate still does not judge a threshold passed as a positional argument
  with a name outside `min_*`/`max_*` (for example `gap = 20.9`); that remains a
  review question, not a gate verdict.

## Verification

`Rscript dev/check-casts.R` — `52 position/threshold cast(s), 20 guarded in
place, 32 registered with a reason; self-test passed`. The self-test contains
an unguarded `as.numeric(x$min_gap)` and a guarded one, so a scanner that
stopped matching the threshold shape fails its own run. Removing the
`check_base_pairs()` line from `design_arms_marker()` and restoring
`as.numeric(min_band_gap)` makes the gate name `design_arms_marker`.
`Rscript dev/test.R '^arms$'` — 78 pass, 0 fail; four bad `min_band_gap`
spellings are refused by name.

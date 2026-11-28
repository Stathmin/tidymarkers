# ADR 0082 — A guard covers the cast it names, not the whole function

- Status: accepted
- Date: 2026-11-28
- Builds on: ADR 0045 (a position is a whole number before the cast),
  ADR 0068 (a position column is whole in its spelling), ADR 0079 (the cast
  gate judges an interval too)

## Context

`dev/check-casts.R` answered one question per function: did some guard appear
before the first cast? A guard for `hits$start`, written in the same function,
therefore also blessed `as.integer(hits$end)`. The unit, not the value, was
judged guarded.

`offtarget_profile()` shows the cost. It checked `hits$start` by hand:

```r
start <- hits$start
if (!is.numeric(start) || anyNA(start) || any(start < 1) || any(start != floor(start))) {
  abort_msg("`hits$start` must be a positive integer.")
}
```

and then cast both bounds: `width <- as.integer(hits$end) - as.integer(hits$start) + 1L`.
Measured on a 21-mer hit at 11–31 in a 41 bp genome, with `hits$end <- 31.5`:
the width was `21`, so the width check passed, the window was read over 11–31,
and the profile reported one hit with `n_mismatch = 0` — a perfect alignment at
a coordinate pair the caller never named. A fractional end is not a number of
bases, and the profile of a hit that does not exist is a wrong answer, not a
refusal.

## Decision

1. `dev/check-casts.R` matches guards to casts by value. Every guard argument
   text seen in the unit is kept; a cast is covered when a guard text equals
   the cast expression, prefixes it (`hits$start` covers
   `as.integer(hits$start[[i]])`), or appears in the message of the nearest
   preceding `abort_msg()`. A guard for one column no longer blesses another.
2. The gate self-test gained the two shapes that used to pass: a guarded
   `start` beside an unguarded `end` in one function, and the same pair with the
   `end` named in the abort message. Both must fail the scanner.
3. `offtarget_profile()` judges `hits$start` and `hits$end` with
   `bad_position()` before the cast, refusing a fractional, `NA`, `Inf` or
   non-positive bound by column name. The hand-written check it replaces
   accepted only numeric spellings; the shared guard also accepts a whole
   character spelling, as everywhere else (ADR 0068).
4. `as_variant_row()` is registered in `dev/casts.yml`, not fixed: its filter
   cast compares the selector with the table's `pos` before `pos_raw` is judged.
   A fractional row can only reach that abort, never a returned position.

## Rejected alternatives

- Add the self-test cases without changing the rule: the gate would then be
  honest about two shapes and blind to a guard for any third column in the same
  function, which is the defect itself.
- Require an exact expression match only: `bad_position(hits$start)` guards
  `as.integer(hits$start[[i]])` legitimately, and an indexed read is still the
  same column. Rejecting it would push correct code into the registry.
- Judge guards only before the first cast: a guard written after an unrelated
  cast still runs before the one it protects, and the message-label match
  already anchors it to a refusal.

## Consequences

- A function with one guarded and one unguarded coordinate is a red gate. The
  scanner reports 52 casts, 19 guarded in place, 33 registered; the guarded
  count fell by one because `offtarget_profile` no longer counts as guarded by
  its sibling guard alone — it now has two guards of its own, and the old
  single-guard reading had been the whole problem.
- The gate remains spelling-based: a guard in another function, or on a value
  copied into a differently named local, still needs a registry entry with a
  reason naming where the check lives.
- Behaviour changes only where a caller passed a non-whole or non-positive
  `hits$end`; every whole table behaves as before.

## Verification

`Rscript dev/check-casts.R` — `52 position/threshold cast(s), 19 guarded in
place, 33 registered with a reason; self-test passed`. Reverting the
`hits$end` guard makes the gate name `offtarget_profile`, and the self-test's
`start`-guarded/`end`-cast shape fails its own run. `Rscript dev/test.R
'^specificity-3p$'` — 52 pass, 0 fail; a hit with `end = 31.5` and one with
`start = 10.5` are both refused by column name.

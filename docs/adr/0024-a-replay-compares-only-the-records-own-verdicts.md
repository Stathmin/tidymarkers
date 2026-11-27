# 0024 — A replay compares only the record's own verdicts

- Status: accepted
- Date: 2026-11-28
- Supersedes: none
- Builds on: [ADR 0015](0015-a-record-resolves-only-its-own-directory.md)
  (an artefact's verdicts are a file, not an RDS)

## Context

ZV13.3 made a stranger's artefact checkable without opening its `solved.rds`:
`marker_artefact()` writes the verdicts as `replay.csv`, and
`marker_replay(artefact, solved = )` compares that table with a re-solve of the
trusted plan. A record written before the sidecar existed has no such table.

The fallback for a missing `replay.csv` was
`recorded <- artefact_replay_table(solved)`. When the caller supplied `solved`,
that reads the caller's own object, so the "recorded" side was the trusted solve
and the replay compared it with itself. Measured before the fix: an artefact
with `replay.csv` deleted and `solved =` passed returned `match = TRUE` with
`recorded_verdict = replayed_verdict = "doomed"` — a record that carries no
verdicts passed as if it had been checked.

ADR 0015 already said a record without `replay.csv` is replayable only where
its RDS is trusted. The code silently widened that to "or wherever the caller
has a solve of their own".

## Decision

The recorded side of a replay always comes from the artefact.

- A manifest that names `replay.csv` is compared against that file
  (unchanged).
- A manifest that names none, replayed without `solved`, reads the artefact's
  own `solved.rds` and derives the recorded table from it: the pre-existing
  artefact@1 path, local to a trusted RDS.
- A manifest that names none, replayed **with** `solved`, is refused with a
  sentence that names the limit and the two exits.
- The comparison covers the union of target labels: a verdict the record
  dropped, and a target the re-solve no longer produces, are both mismatches.
  A table with no rows claims nothing and does not pass.
- A table is refused, not reconciled, when it names no target on a row or names
  one target twice: two verdicts do not say which one is the record's claim.
- `decisions.csv` is not used to rebuild the table: a target whose every
  attempt was rejected has no `chosen` row, so the sidecar cannot reproduce the
  one-row-per-target report exactly.

## Consequences

- A record with no verdicts can no longer pass on the caller's own solve; the
  silent pass is a refusal.
- Old records stay replayable exactly where ADR 0015 allows it: from their own
  RDS, never from someone else's `solved`.
- `marker_replay()`'s `@param solved` and the ZV13.3 note now say which records
  `solved` applies to.
- Tests: `test-plan.R` — "a trusted solve cannot replay a record that recorded
  no verdicts" (refusal, then the own-RDS path still matches) and "a replay
  checks the record's row set, not just the rows it kept" (a dropped verdict,
  an empty table and the complete table) and "a replay refuses a table that
  does not name one verdict per target" (a duplicated and an empty label).

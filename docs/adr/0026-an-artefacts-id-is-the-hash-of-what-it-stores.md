# 0026 — An artefact's id is the hash of what it stores

- Status: accepted
- Date: 2026-11-28
- Supersedes: none
- Builds on: [ADR 0025](0025-a-record-member-is-resolved-not-only-named.md)
  (the manifest is the record's claim about itself)

## Context

`artefact_id()` names a design: it hashes the manifest core so two runs of the
same design share an id and a replay can be tied to the record it replayed.
The id was computed from the in-memory manifest, not from the JSON the writer
stores. Measured 2026-11-28 on the branch, on a record written by the package:

- the stored id was `77233810…f63a2a3`, and recomputing it from the record's
  own `artefact.json` gave `4efe2682…1c9c884`;
- the difference is one field. `jsonlite::toJSON()` with no `na` argument
  renders a character `NA` as the string `"NA"`; the writer used
  `na = "null"`, so the same field is `null` on disk. An absent engine's
  version is exactly such a value.

An id that cannot be re-derived from the record is a label, not a check. It
also went unverified: `marker_artefact_read()` took `manifest$id` as given, so
a manifest edited after it was written still read, replayed and accepted
outcomes under the old id.

## Decision

1. **Hash the bytes the writer emits.** `artefact_id()` serialises the core
   with the writer's own rules (`auto_unbox = TRUE`, `null = "null"`,
   `na = "null"`, `digits = NA`), so the id is recomputable from the record.
2. **Verify it on read.** `artefact_manifest()`, the single reader every
   artefact path goes through, recomputes the id for an `artefact@2` record
   and aborts when it differs:
   `The id in "…/artefact.json" does not match its own manifest: "…".` The
   hint names the recomputed hash and says the manifest was edited after it was
   written.
3. **Legacy and future records are not verified.** `artefact@1` predates the
   `contract` field, so its core cannot be re-derived and the check would
   report every old record as edited; a schema newer than this reader is
   refused for its version alone, which is the honest reason.

## Consequences

- A tampered record fails at its first read instead of feeding edited values to
  replay or to outcome recording. Test: `test-trace-schema.R`, "an artefact
  reader refuses a manifest edited after it was written"; the untouched record
  still reads.
- The frozen-hash test in `test-trace-schema.R` keeps its hash: its only `NA`
  is logical, which both spellings render as `null`.
- An `artefact@2` record written by 0.4.1-dev before this change stores an id
  that no longer verifies, so it is refused. No such record was published;
  0.4.0 wrote `artefact@1`. This is the cost of making the id checkable, and it
  ends at this branch.
- Tests that edit a record to simulate something else (a moved engine, an
  escaping member path) must keep the core consistent with the stored id; the
  alternative is to exercise the tamper refusal, not the behaviour they target.
- The rule does not make a record authentic: whoever can write the manifest can
  recompute its id. It makes the id a claim the record can be held to.

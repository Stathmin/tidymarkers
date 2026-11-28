# 0027 — A record directory holds one design

- Status: accepted
- Date: 2026-11-28
- Supersedes: none
- Builds on: [ADR 0026](0026-an-artefacts-id-is-the-hash-of-what-it-stores.md)
  (the manifest is the record's claim about itself)

## Context

`marker_artefact()` writes a design into a directory. Its writer has kept an
existing `solved.rds` since its first commit
(`if (!fs::file_exists(...)) saveRDS(...)`), while every other member is
overwritten. Measured 2026-11-28: calling it twice into one directory with two
different solves leaves the second manifest and the first solve. A replay of
that directory re-solves the first design and compares it with the second
design's `replay.csv`, so the record's own file contradicts its own manifest.

Nothing needed the old file: the manifest is rewritten either way, and the
solved object passed to the call is always in hand. The guard only made a
rewritten directory a mix of two records.

## Decision

1. **Write the solve the manifest describes.** `marker_artefact()` always calls
   `saveRDS()` for `solved.rds`, so a rewrite replaces it with the design in
   `artefact.json`.
2. **The manifest stays the list of members.** A sidecar the new design does
   not produce (for example `checks.csv` when the second solve has no checks)
   is still named as absent, and every reader takes its file list from the
   manifest; a residual file on disk is not a member of the record.

## Consequences

- Rewriting a directory is a complete replacement of the design, not a patch.
  Test: `test-trace-schema.R`, "a second write replaces solved.rds instead of
  another design's"; it fails at the previous commit.
- A replay of a rewritten directory now re-solves the manifest's design.
- The second decision is a recorded limit, not a cleanup: a stale sidecar
  stays on disk until the directory is emptied. It is inert because no member
  is read by name except through the manifest; deleting it is not worth a
  writer that may not own every file in the directory.

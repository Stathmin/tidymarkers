# 0015 — A record resolves only the paths inside its own directory

- Status: accepted
- Date: 2026-11-27
- Supersedes: none
- Builds on: [ADR 0011](0011-release-artefact-and-mcp-contract.md) (the release
  artefact) and commit `b431859` (the manifest sidecar name)

## Context

An artefact is a directory: `artefact.json` names its sidecars, and an RO-Crate
is the same record plus a metadata descriptor that travels. Both files are
inputs a stranger can hand us — a collaborator's crate, a fetched archive, a CI
cache — so every name inside them is untrusted.

Two escapes were found, both by the same question: *what does this record let
the reader or writer touch?*

1. `artefact.json` with `files$decisions = "../escape.csv"` made
   `marker_export_rocrate()` copy a file from outside the artefact directory,
   and `marker_artefact_read()` read it. Fixed in `b431859`: every sidecar name
   is one path component, refused otherwise.
2. `ro-crate-metadata.json` with a graph entity
   `{"@id": "../secret.txt", "sha256": "<the outside file's hash>"}` imported
   silently: `marker_import_rocrate()` joined the id to the crate directory,
   hashed the file it found there, saw the declared hash match, and returned the
   artefact. The crate vouched for a file it does not carry.

Both are the same defect at two ends of one contract: a record is trusted to
describe itself and not the filesystem around it.

## Decision

1. **A persisted record names only its own contents.** An artefact sidecar name
   is exactly one path component (no `/`, no `\`, not `.` or `..`). A crate
   graph id resolves inside the crate directory; absolute paths and `..`
   segments are refused with the id in the message.
2. **Refusal happens before use.** The name is checked before a file is opened,
   copied or hashed, so a rejected record leaves no partial output and reads
   nothing it did not already own.
3. **A reference outside the record is data, not a path.** The one legitimate
   outside reference — an artefact input that came from a file — is recorded as
   `contentUrl` + `sha256` (RO-Crate `File` entity with a PROV type) and is
   never opened by the import. The crate says where the input came from; it does
   not claim to carry it.
4. **Each new record reader answers the question in its test.** A reader of a
   persisted record gets a test that hands it a record naming `../something` and
   expects an error, not a successful read of the neighbour.

## Alternatives

- **Trust records written by this package.** Rejected: the reader cannot tell
  its own output from a hand-edited copy, and `artefact.json` is exactly what a
  collaborator sends when a design is being reviewed.
- **Normalise instead of refuse** (`../a/../b.csv` → `b.csv`). Rejected: silent
  normalisation hides that the record was malformed; the caller learns nothing
  and a real escape attempt looks like success.
- **Resolve symlinks and compare real paths.** Rejected as the primary check:
  `fs::path_abs()` already gives the contract we want (lexical containment), and
  a symlink inside the directory is the operator's own filesystem, not the
  record's claim. Noted below as a limit, not as a gate.

## Consequences

- A crate that escapes is a named error (`The crate declares "../x", which is
  outside ...`), not a silent hash of a neighbour. Test:
  `test-rocrate.R`, "a crate cannot declare a file outside its own directory".
- The rule is lexical first, then physical: an `@id` is checked as a path, and
  a member that exists is also resolved with `fs::path_real()` and compared with
  the crate's real root, so a symlink (or a file reached through a linked
  directory) out of the crate is refused instead of hashed. A dangling link is
  refused as a missing member. Residual limit: hard links, which are
  indistinguishable from ordinary files by path.
- Any new persisted record (trace, replay, outcome corpus) inherits item 4: a
  reader without an escape test is an open item, not a finished reader.
- Execution is a path out of the record too. An artefact's `solved.rds` is an R
  object, and `readRDS()` runs what the file says; the verdicts are therefore
  also written as `replay.csv` (three character columns), and
  `marker_replay(artefact, solved = <trusted>)` compares against that table
  without opening the RDS. A record without `replay.csv` can still be replayed
  only where its RDS is trusted (`ZV13.3`).

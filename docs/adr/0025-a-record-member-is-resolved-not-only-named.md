# 0025 — A record's member is resolved, not only named

- Status: accepted
- Date: 2026-11-28
- Supersedes: none
- Amends: [ADR 0015](0015-a-record-resolves-only-its-own-directory.md) (the
  lexical check is now the first of two)

## Context

ADR 0015 decided two things about a persisted record: a sidecar name is one
path component, and the name is refused before a file is opened. Its
"Consequences" also recorded a residual limit for artefacts — a symlink inside
the directory was called "the operator's own filesystem, not the record's
claim" — while the same section added a physical check for RO-Crate members
(`fs::path_real()` compared with the crate's real root, ZV13.2).

The two answers drifted apart. Measured 2026-11-28 on the branch: an artefact
directory whose `decisions.csv` is a symlink to a file beside it passes
`artefact_safe_name()` (the name is one component), and
`marker_artefact_read()` reads the target. `marker_export_rocrate()` copies the
target's bytes into the crate under the record's own name, and
`marker_replay()` opens a linked `solved.rds` with `readRDS()` — the one path
ADR 0015 already names as execution, not just reading. A record fetched from a
collaborator or a CI cache therefore still reaches outside its directory.

The crate check is the rule the family wants; the artefact reader was simply
not brought along. A name that is one component is necessary, not sufficient:
the file the name resolves to is what the reader trusts.

## Decision

1. **Lexical containment first, physical containment second.** A record member
   is still refused if its name leaves the directory (`artefact_safe_name()`).
   If the path exists, it is also resolved with `fs::path_real()` and must stay
   under the record's real root; a link out is refused with both paths named.
2. **Every member read goes through the same resolver.** `artefact.json`,
   every manifest sidecar, `replay.csv` and `solved.rds` are resolved by
   `artefact_member_path()` before they are opened. The crate writer uses it
   before copying, so a record cannot smuggle a neighbour's bytes into a crate.
3. **A dangling link is a missing member, not a read.** `fs::file_exists()` is
   false for a dangling link, so the reader returns its declared empty table or
   names the missing file — it does not follow the link.
4. **The test is part of the contract.** A record reader that accepts a name
   from the record gets a test with a link out; the artefact reader and the
   crate writer share one test in `test-rocrate.R`.

## Alternatives

- **Keep the lexical-only artefact rule.** Rejected: the crate and the
  artefact are the same record at two ends, and the reason given in ADR 0015
  ("the operator's own filesystem") does not hold for a record received from
  someone else, which is exactly what `marker_artefact_read()` is for.
- **Resolve and normalise the link instead of refusing.** Rejected for the same
  reason as in ADR 0015: a record that points outside is malformed, and
  silently reading the target hides the escape.
- **Check only the crate import.** Rejected: the reader and the writer are
  separate trust boundaries; a linked `solved.rds` is executed by
  `marker_replay()` before any crate exists.

## Consequences

- A linked member is a named error:
  `The artefact names "decisions.csv", which resolves to "/…/outside.csv" outside "/…/artefact"`.
  Test: `test-rocrate.R`, "an artefact sidecar cannot link out of the record";
  the reader and `marker_export_rocrate()` both refuse it.
- The rule now matches ADR 0015's crate consequence: lexical first, physical
  for an existing member. Hard links remain indistinguishable by path and stay
  a recorded limit.
- A locator that the package derives itself (a cache index, a path passed by
  the caller to `marker_trace_read()`) is still not a record claim and is not
  covered by this rule.

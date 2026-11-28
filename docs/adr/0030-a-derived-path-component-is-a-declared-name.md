# 0030 — A derived path component is a declared name too

- Status: accepted
- Date: 2026-11-28
- Supersedes: none
- Builds on: [ADR 0029](0029-an-output-writer-creates-its-files-not-writes-through-a-link.md)
  (an output writer creates its files, not writes through a link; point 4: a
  file beside the directory is still a declared name)

## Context

ADR 0029 point 4 closed `<dir>-identity.csv` because the name was derived from
the destination rather than declared inside it. The same writer,
`marker_export_outcomes()`, creates two more derived paths before it writes a
crate:

- `<dir>/crates/` — the directory that holds every crate;
- `<dir>/crates/<pseudonym>/` — one crate per artefact.

Neither was checked. `refuse_links()` ran on `outcomes.csv` inside a crate, and
the identity map was guarded, but a link at `crates` or at one pseudonym was
followed by `fs::dir_create()` and every file of the crate landed in its target.

Measured 2026-11-28 on the branch, with the corpus fixture of
`test-outcomes-crate.R`:

| link at | result before the fix |
| ------- | --------------------- |
| `<dir>/crates` | both crates created in the linked directory (8 files outside `dir`); the export reported success |
| `<dir>/crates/site-…` | that crate written in the linked directory (8 files outside `dir`) |

The identity map was still written, so the record claimed a crate under `dir`
that was not there. This is the ADR 0025/0029 rule one component above the
names it already covered: the guard named files, not the directories the writer
derives on the way to them.

## Decision

1. **Every path component the writer derives is a declared name.** Before the
   first crate, `marker_export_outcomes()` guards `crates` inside `dir` and each
   pseudonym inside `crates` with the same `refuse_links()`.
2. **The guard runs before any crate exists.** Both checks precede the loop, so
   a link at the second pseudonym does not leave the first crate on disk; the
   handover is all-or-nothing for the set it is about to write.
3. **One rule, one message.** The refusal is the existing "which is a link"
   sentence; no new error class.

## Consequences

- A link at `crates` or at a crate's pseudonym stops the call before a single
  file is written, and the outside directory keeps its contents. Test:
  `test-outcomes-crate.R`, "the handover refuses to write its crates through a
  link" — both link positions, the outside directory stays empty, and no
  identity map is written for the first case; before the fix both cases wrote 8
  files outside `dir` (37 pass / 0 fail after, 5 fail before).
- The rule is now stated for the whole path the writer creates: a file name
  inside a directory (0029 point 1), a file beside it (0029 point 4), and a
  directory between the two (this ADR).
- The caller-named root `dir` is not guarded: the caller chose that path, and
  following a link there is the ordinary meaning of naming it. This ADR covers
  only components the package derives itself.

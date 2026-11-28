# 0029 — An output writer creates its files, not writes through a link

- Status: accepted
- Date: 2026-11-28
- Supersedes: none
- Builds on: [ADR 0025](0025-a-record-member-is-resolved-not-only-named.md)
  (a record member is resolved, not only named; point 5: every writer of the
  catalogue refuses to write through a link)

## Context

ADR 0025 point 5 made one rule for the writers that fill a record directory:
`marker_artefact()`, `marker_export()`, the RO-Crate export and the outcomes
writers each call `refuse_links()` on every name they are about to create, before
the first write. ZV25.7 closed `marker_export()` for exactly this reason: a link
at a member name sent the write to a file outside the directory, and `keep`
became a CSV.

The two result writers were never brought into the family:

- `write_marker_results()` writes `marker_design.csv`,
  `marker_primer_qc.csv`, `marker_dimer_dg.csv` (and `marker_blast_hits.csv`
  when the assay has hits), then calls `export_primers()`;
- `export_primers()` is exported on its own and writes `{base}_idt.csv`,
  `{base}_lgc.csv`, `{base}_tidy.csv`.

Measured 2026-11-28 on the branch, one link per call, an outside file holding
`keep`:

| link at | call | result before the fix |
| ------- | ---- | --------------------- |
| `marker_design.csv` | `write_marker_results()` | wrote the design table over the outside file; `keep` lost |
| `marker_FRI_Chr4_269090_tidy.csv` | `write_marker_results()` | wrote the order table over the outside file; `keep` lost |

The outcomes handover has the same shape one level up: `marker_export_outcomes()`
writes one crate per artefact and then the identity map *beside* the export
(`<dir>-identity.csv`, deliberately outside every crate). That name was derived,
not declared, so no guard saw it: a link at `<dir>-identity.csv` was followed
after every crate had already been written, and the map replaced the outside
file.

`safe_name_part()` already keeps the *name* inside `dir`; it says nothing about
an existing link at that name. The output writers therefore disagreed with the
record writers on the same question, as the VCF family did before ADR 0028.

## Decision

1. **Every writer that creates named files in a directory calls
   `refuse_links()` on all names before the first write.** `write_marker_results()`
   declares its own three (or four, with BLAST hits) names plus the three order
   sheets it delegates, so a refusal happens before any file in the set is
   written. `export_primers()` declares its own three names, because it is
   exported and can be called without `write_marker_results()`.
2. **One rule, one message.** Both call the same `refuse_links()` the record
   writers use, with its message ("which is a link"). A caller that meets the
   refusal reads the same sentence whichever writer it called.
3. **A conditional output is guarded when it is written.** `marker_blast_hits.csv`
   is only declared when `assay$blast` is not `NULL`, so a link at a name the
   call would not write is not refused.
4. **A file beside the directory is still a declared name.**
   `marker_export_outcomes()` computes `<dir>-identity.csv` before anything is
   created, guards it, and only then writes the crates; the guard covers a name
   derived from the destination, not only a name inside it.

## Consequences

- A link at any output name stops the call with a named error and leaves the
  outside file untouched. Test: `test-output.R`, "the output writers refuse to
  write through a link" — all six names of `write_marker_results()` and the
  direct `export_primers()` path; before the fix the outside file lost `keep`.
- The writer is all-or-nothing for its own name set: the guard runs before
  `marker_design.csv`, so a link at an order sheet does not leave three files of
  a half-written result set behind. The handover is the same: a link at
  `<dir>-identity.csv` refuses before the first crate exists. Test:
  `test-outcomes-crate.R`, "the handover refuses to write its identity map
  through a link".
- The rule now names no writer by file: it is "all names before the first
  write", the same sentence ADR 0025 uses for records and ADR 0028 uses for a
  VCF row.

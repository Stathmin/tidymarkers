# Self-repair log

Append-only log of blockers met by the three-day programme. The protocol is in
[ADR 0012](../adr/0012-autonomy-and-self-repair.md); the programme is in
[the plan](2026-11-27-integration-3day.md).

Every detection, attempt, bypass and quarantine gets one line. No edits: a
correction is a new line. No secrets, prompts, sequences or file content — ids
and actions only.

| time (UTC) | blocker | class | action | result | commit | attempt | status | next task |
| ---------- | ------- | ----- | ------ | ------ | ------ | ------- | ------ | --------- |
| 2026-11-27 16:24 | `R CMD check` is ~7 min and the release gate needs it | transient | run in the background, take the result later; no second run on the same tree | bypassed: gate ran once | `c74ef0b` | 1 | bypassed | day-2 gate |
| 2026-11-27 16:27 | `check-docs.R`: ADR index and ADR record disagree on status | logic | green gate prints the same status on both sides; stale inventory rows fixed | fixed: `check-docs.R` ok | `c2938e1` | 1 | fixed | ZV11 metrics gate |
| 2026-11-27 16:38 | `dev/bench.R` compares against a hostname, so CI is `not_run` forever | environment | host key from `TIDYMARKERS_BENCH_HOST`, stable runner-class name, baseline committed per key | bypassed: local `not_run` is honest; baseline on the runner is open (ZV12) | `8d2b35a` | 1 | quarantined | CI step reads the key |
| 2026-11-27 16:46 | a criteria threshold is parsed as a call, not a number | logic | parse to numeric at the boundary, refuse a non-numeric threshold | fixed: `test-criteria.R` | `dd57511` | 1 | fixed | hardening notes |
| 2026-11-27 17:22 | a VCF POS outside integer range becomes `NA` and the record disappears | logic | refuse the file, name count and first record | fixed: `test-io.R` | `e5c2c87` | 1 | fixed | FASTA path audit |
| 2026-11-27 17:37 | an unknown FASTA base is dropped by Biostrings, shifting every later coordinate | logic | both readers abort on a base outside the alphabet, name the file | fixed: `test-io.R` | `4ced5e1` | 1 | fixed | VCF diagnostics |
| 2026-11-27 17:40 | a missing VCF column errors without naming the column | logic | the error names the column and the file | fixed: `test-io.R` | `1857f69` | 1 | fixed | artefact record audit |
| 2026-11-27 17:42 | `artefact.json` sidecar name `../escape.csv` is read and packaged | logic | sidecar name is one path component, refused before any open or copy | fixed: `test-rocrate.R` | `b431859` | 1 | fixed | crate import audit |
| 2026-11-27 17:47 | RO-Crate `@id` `../outside-the-crate.txt` is hashed and accepted | logic | graph id resolves inside the crate directory or the import refuses it | fixed: `test-rocrate.R`; rule in ADR 0015 | `3944ca2` | 1 | fixed | ZV13.1 other readers |

The first nine rows are a backfill written 2026-11-27 17:50 UTC from the branch's
own commit history: the placeholder row ("no blockers recorded yet") was in place
while the fixes above were already committed. Times are commit times converted to
UTC; the commit column is the evidence. Nothing here was re-timed or re-worded
after the fact — a later correction would be a new row.

Status values: `open`, `fixed`, `bypassed`, `quarantined`, `watchdog`.
Class values: `transient`, `environment`, `data`, `access`, `logic`, `spec`,
`external`.

Daily metrics derived from this log and from `stage = "repair"` JSONL:
blockers/hour, MTTR, auto-resolved share, tasks per blocker. They are written to
`data/metrics/metrics-<date>.md`.

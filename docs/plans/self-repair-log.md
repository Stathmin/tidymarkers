# Self-repair log

Append-only log of blockers met by the three-day programme. The protocol is in
[ADR 0012](../adr/0012-autonomy-and-self-repair.md); the programme is in
[the plan](2026-11-27-integration-3day.md).

Every detection, attempt, bypass and quarantine gets one line. No edits: a
correction is a new line. No secrets, prompts, sequences or file content — ids
and actions only.

| time (UTC) | blocker | class | action | result | commit | attempt | status | next task |
| ---------- | ------- | ----- | ------ | ------ | ------ | ------- | ------ | --------- |
| — | — | — | no blockers recorded yet | — | — | — | — | — |

Status values: `open`, `fixed`, `bypassed`, `quarantined`, `watchdog`.
Class values: `transient`, `environment`, `data`, `access`, `logic`, `spec`,
`external`.

Daily metrics derived from this log and from `stage = "repair"` JSONL:
blockers/hour, MTTR, auto-resolved share, tasks per blocker. They are written to
`data/metrics/metrics-<date>.md`.

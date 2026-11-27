# 0012 — Autonomy and self-repair: a blocker never stops the programme

- Status: proposed
- Date: 2026-11-27
- Supersedes: none
- Builds on: [ADR 0008](0008-user-scenarios-as-contract.md) (scenarios),
  [ADR 0011](0011-release-artefact-and-mcp-contract.md) (metrics)

## Context

The three-day programme is deliberately parallel: seven scenarios, ten queue
items, three release gates and a metrics baseline. Its known failure mode is
not a wrong design; it is **idleness**. One missing Bioconductor package, one
engine that crashes on a fixture, one missing full text, one slow
`R CMD check` can stop the lane and leave the other nine items waiting. The
earlier plan rule — “a day ends green, a red gate is not carried over” — makes
that worse: it turns one local blocker into a programme-wide stop.

The repository already contains the raw material for recovery:
`dev/install-deps.R` installs the hard dependencies and the package;
`internal_duplex_dg()` / `internal_hairpin_dg()` are pure-R fallbacks;
`specificity_of()` dispatches memory/bowtie/blast/FM; `docs/decisions/` records
operator decisions; `tests/testthat/fixtures/` holds the committed data; and
`docs/plans/2026-11-27-integration-3day.md` names the failure criteria. What is
missing is a rule that turns a blocker into a *repair task* while the main lane
continues.

A blocker is information. It is never a reason to stop the programme.

## Decision

### 1. One rule

**A blocker never stops the programme.** The blocked task leaves the main lane,
a repair task is created, and the lane immediately takes the next task. The
programme stops only when the operator says so, when a prohibition would be
violated, or when the final acceptance criteria are known to be unreachable and
the operator accepts that result.

Nothing in this ADR weakens a gate. It changes what happens *while* a gate is
red: repair or quarantine, not idle.

### 2. Detection and classification

Every blocker is classified before any fix. The class selects the default
action and the budget.

| Class | Signals | Default action |
| ----- | ------- | -------------- |
| **transient** | HTTP 429, timeout, rate limit, connection reset, a single flaky process | retry with backoff; no code change; record the retry count |
| **environment** | missing package, missing binary, wrong R library, engine absent, `R_LIBS` order | install or point the environment (`dev/install-deps.R`, `marker_install_tools()`, `R_LIBS_SITE`); never skip silently |
| **data** | missing full text, fixture, file, column, or a table that does not parse | obtain from a named source; if impossible, generate a fixture from measured data with the assumption recorded; never invent a value |
| **access** | missing key, expired credential, no permission, read-only mount | use the read-only fallback; ask the human; do not exfiltrate or guess a secret |
| **logic** | reproducible bug, failing test, wrong verdict, exception in code | reproduce, minimal fix at the shared root, add the regression test |
| **spec** | ambiguous requirement, two readings of a threshold, no named source | choose the conservative default, write a revisit note, ask the human asynchronously; do not invent a policy |
| **external** | upstream tool crashes, service down, incompatible release, XSD missing | fallback or replacement, pin the version, quarantine the claim if no fallback exists |

A blocker whose class is unclear is `spec` until evidence moves it. A blocker
that changes class is re-classified and the log records both classes.

### 3. Repair task, budget and priority

A repair task is a first-class queue item:

```
R-<n> | blocker: <id> | class: <class> | priority: P0|P1|P2
      | lane: repair | budget: 30 min/attempt, max 3 attempts
      | owner: <role> | status: open|fixed|bypassed|quarantined
      | next: <task id>
```

- **P0** blocks J1, J6, J7 or a release gate; repair starts immediately.
- **P1** blocks a day item; repair runs beside the main lane.
- **P2** is polish or a non-blocking warning; repair waits for a free slot.

The budget is **30 minutes per attempt and at most 3 attempts** (90 minutes
total) per blocker. A fix that needs a fourth attempt is a new blocker with a
new decision, not an unlimited retry.

### 4. Fix, bypass, quarantine

The repair lane tries, in order:

1. **Fix** the root cause and add the regression test.
2. **Bypass** with a declared fallback: pure-R thermodynamics instead of
   OligoArrayAux, memory instead of an index, a cached fixture instead of a
   download, a background `R CMD check` instead of a blocking one. The bypass
   is recorded with its limit.
3. **Quarantine** when neither works. Quarantine writes the risk, the impact,
   the conservative default (if any), the question to the human and a revisit
   date. The item stays open; it is not marked done.

After quarantine the lane **immediately takes the next task**. A quarantined
item does not hold the queue.

### 5. The self-repair log

Every detection, attempt, bypass and quarantine is appended to
[`docs/plans/self-repair-log.md`](../plans/self-repair-log.md). The columns are:
time (UTC), blocker id, class, action, result, commit, attempt, status, next
task. The log is append-only; a correction is a new line, not an edit. It
contains no secrets and no data content, only ids and actions.

The log is the source for the self-repair metrics in §8. An entry is written
within 10 minutes of detection; a blocker with no line is an unlogged stop and
fails the daily gate.

### 6. Anti-idle rules

- **Three next tasks.** The queue always holds at least three ready tasks. If
  fewer remain, the owner creates the next tasks before the lane drains.
- **Time-box.** Every task has a start time and a box; the box is shorter than
  the 30-minute repair attempt only when the task is smaller.
- **Watchdog.** No commit and no log progress for **45 minutes** during the
  working window means the lane is stuck: switch to the next task and write a
  watchdog line. The watchdog is a signal, not a punishment.
- **Repair lane does not block the main lane.** Repair runs in its own lane
  with its own budget; the main lane never waits for it. If the same agent
  owns both, the main task is paused only for the length of one attempt and
  then resumed or handed off.

### 7. Prohibitions

These are never a valid repair:

1. **Do not invent data.** A missing value is `NA` with a reason, or a fixture
   generated from a named source with the assumption recorded. A synthetic
   scale fixture is labelled synthetic and is never reported as measured.
2. **Do not weaken tests.** No skipping a failing assertion, no loosening a
   threshold, no deleting a golden case. A test that cannot run is quarantined,
   not edited to pass.
3. **Do not push to `main` without the gate.** Repair commits go to the repair
   branch; `dev/check.R`, `check-user-api.R`, `check-mcp-parity.R` and the
   metrics gates run before any merge.
4. **Do not take an irreversible action without a human.** Deleting data,
   rewriting history, changing a published contract, spending money or
   rotating a credential waits for the operator.

### 8. Metrics

The daily report adds a self-repair panel:

- **blockers/hour** — new blocker ids divided by elapsed programme hours;
- **MTTR** — mean time from detection to `fixed` or `bypassed`; quarantined
  blockers are reported separately and are not folded into MTTR;
- **auto-resolved share** — `(fixed + bypassed) / all blockers`;
- **tasks per blocker** — repair tasks divided by distinct blocker ids.

Repair events also write `stage = "repair"` lines to the JSONL contract of
[ADR 0011](0011-release-artefact-and-mcp-contract.md) §5, so duration, retries,
status and queue depth are measured by the same collector. The class lives in
the log and in `object_id` as `blocker:<id>`; it is not free text in the JSONL.

A day with a higher blockers/hour than the previous day is reported, not
hidden. A quarantine share above 30% is a process signal, not a code failure.

### 9. Known blockers and their bypasses

These are the blockers this repository has already met. Each one has a
declared route; none is allowed to stop the lane.

| Blocker | Class | Route |
| ------- | ----- | ----- |
| No Bioconductor package (`GenomicRanges`, `Rsamtools`, `DECIPHER`, …) | environment | `Rscript dev/install-deps.R`; `BiocManager::install()` with the pinned release; verify with `requireNamespace()`. If the host cannot install, quarantine the external suite and keep the pure-R tests green. |
| External binary crashes or is absent (OligoArrayAux, BLAST, bowtie, ImageMagick) | environment / external | use the pure-R fallback (`internal_duplex_dg()`, `internal_hairpin_dg()`), the memory specificity path, or another engine (`specificity_of()` dispatches bowtie/blast/FM). Pin the version; if no fallback exists, quarantine the claim and say so. |
| No fixture or no full text | data | generate a fixture from a named source with the assumption recorded in `docs/decisions/` or `inst/extdata/`; a synthetic scale fixture is labelled synthetic. Never invent a measured number. |
| A human decision is needed (threshold, policy, release) | spec / access | take the conservative default already in the registry (no LD default, no p-value cut, `calibrated = FALSE`), write a revisit note, ask the human asynchronously, continue on the next task. |
| `R CMD check` is slow (about 7 minutes) | transient / environment | run it in the background, continue the lane, collect the result, never run it twice on the same tree. |
| No 16 Gb fixpur on the host | data / environment | use the three-build fixture and the synthetic 250 Mb chromosome; record the real run as an operator action in `docs/decisions/0001-16gb-fixpur.md`; never claim a measured 16 Gb result. |
| Network / Ensembl Plants is down | external | use the committed FASTA/transcript fixture or the cache; tests never require a live download. |
| Engine absent for a scenario | environment | the contract step that does not need the engine still asserts; the design step is a declared skip, not a pass, and the scenario is marked degraded in the report. |
| A refactor changes a verdict | logic | stop merging that item, reproduce on the golden set, fix the root cause; never relax the golden verdict. |
| Missing schema/XSD for ST.26 validation | data / external | validate the structure with the committed schema if present; otherwise quarantine the validation with the XSD question, keep the export. |

### 10. Gates

1. Every blocker has a log line within 10 minutes; the daily report lists the
   count by class.
2. The queue holds at least three ready tasks at every checkpoint.
3. No 45-minute gap without a commit or a log line during the working window;
   the watchdog line names the task switched to.
4. The repair lane never consumes the only main-lane slot.
5. `blockers/hour`, `MTTR`, `auto-resolved share` and `tasks per blocker` are in
   `data/metrics/metrics-<date>.md`.
6. No repair commit reaches `main` without `dev/check.R`,
   `check-user-api.R`, `check-mcp-parity.R` and the metrics gates.
7. No test was weakened, skipped without a quarantine line, or deleted to make
   a repair green.

## Consequences

- The programme keeps moving through blockers. Some items end quarantined; that
  is visible and honest, not hidden.
- Every blocker costs a log line and a bounded repair attempt. The cost is
  small compared with a stopped lane.
- Repair commits are isolated on a repair branch and still pass the release
  gate; autonomy does not mean unchecked pushes.
- A quarantine is not a closure. The item stays in `ROADMAP.md` §ZV with the
  human question and the revisit date.
- The self-repair metrics are part of the same JSONL and the same daily report,
  so “we fixed it” is a number, not a claim.

## Alternatives

- **Stop the line on the first red gate.** Rejected: it is the failure mode
  this ADR exists to remove; nine independent items wait for one missing
  package.
- **Retry without a cap.** Rejected: a 30-minute cap and three attempts force
  a quarantine decision instead of an infinite loop.
- **Silent workaround.** Rejected: an unlogged bypass is indistinguishable from
  a wrong result. The log is the contract.
- **Weaken the test to keep the gate green.** Rejected: it converts a visible
  blocker into an invisible defect.
- **Wait for the human on every ambiguity.** Rejected: the conservative default
  plus a revisit note keeps the lane moving; the human answers asynchronously.

## Revisit trigger

More than 30% of blockers quarantined, an MTTR above 90 minutes, or a repair
that cannot be bounded by the 30-minute/3-attempt budget. At that point the
classification or the budget is wrong, not the rule that the programme
continues.

## References

- [`docs/plans/2026-11-27-integration-3day.md`](../plans/2026-11-27-integration-3day.md)
  — the programme and its gates.
- [`docs/plans/self-repair-log.md`](../plans/self-repair-log.md) — the
  append-only repair log.
- [`dev/install-deps.R`](../../dev/install-deps.R) — the environment route.
- [`R/marker_exec.R`](../../R/marker_exec.R) — timeouts and external processes.
- [`R/marker_dg_internal.R`](../../R/marker_dg_internal.R) — the pure-R
  thermodynamic fallback (`internal_duplex_dg()`, `internal_hairpin_dg()`).
- [`docs/decisions/0001-16gb-fixpur.md`](../decisions/0001-16gb-fixpur.md) —
  the 16 Gb operator decision.
- [`ROADMAP.md`](../../ROADMAP.md) — §ZV, the queue the repair lane serves.

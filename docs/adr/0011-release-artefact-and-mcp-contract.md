# 0011 — The release artefact and the MCP contract are the product

- Status: proposed
- Date: 2026-11-27
- Supersedes: none
- Builds on: [ADR 0004](0004-mcp-ci-and-cross-project-tooling.md) and
  [ADR 0007](0007-module-seams-and-mcp-contract.md)

## Context

The package already produces two machine-readable objects that a user does not
read as R output:

- **The artefact.** `marker_artefact()` writes `artefact.json`,
  `decisions.csv`, `checks.csv`, `why.csv`, `targets.csv`, optional
  `signals.csv` and `solved.rds`. The manifest records the package version, git
  sha, `renv.lock` hash, platform, engine stamp, hashed inputs, plan layers,
  costs, LD, parameters, overrides, counts and file names. `created` comes from
  `SOURCE_DATE_EPOCH` when set, and `artefact_id()` deliberately excludes
  `created`, so the text files are byte-identical on a rerun. `marker_replay()`
  reads `solved.rds`, re-solves the plan and compares the verdict per target.
- **The MCP contract.** `inst/mcp/tools.R` holds **38 tools**; each has an
  `inputSchema`, a handler and annotations. `marker_technology` is the only
  mutating tool. `mcp_validate_args()` validates a call against the schema and
  returns JSON-RPC `-32602` with `tool`, `path`, `expected` and `got`.
  `dev/mcp-schema.R` publishes `$defs.inputSchemas` and `annotations` into
  `inst/mcp/schema.json`; `tests/testthat/test-mcp.R` fails on drift.

Those are two products, but the repository does not treat them as one release
surface:

- **The artefact has no scenario and no contract version.** A reader cannot
  tell whether a JSON file came from J1, J3 or J6, or which contract the file
  obeys. `schema = "tidymarkers/artefact@1"` is the only version marker.
- **MCP parity is tested in pieces.** `test-mcp.R` compares the registry with
  `schema.json`; `test-mcp-validation.R` exercises 38 tools; no gate compares
  the tool set, the handler arguments, the annotations and the user-facing
  inventory at once. A tool can exist in `tools.R` and be absent from
  `docs/user-api-surface.csv` without a failure.
- **The release gate is not named.** `dev/check.R` runs the local checks and
  `--full` adds `R CMD check`; `check-user-api.R` is a separate command. A
  release is a sequence an operator remembers, not one gate.
- **No run metrics.** Nothing records p50/p95 per stage, throughput, cost per
  J1–J7 route, cache hits or the failure share. A regression in any of them is
  invisible until a reviewer notices it, and there is no baseline to compare
  against. Collecting this at the end of the programme would measure the
  finished tree, not the three days that produced it.
- **`solved.rds` is non-byte-stable by design.** It carries an igraph pointer
  and a random id; the current test excludes it by name. That exclusion is a
  contract fact and must stay visible in the artefact, not only in a test.

## Decision

### 1. Two products, one release gate

The release product is:

1. a **self-documenting panel artefact** — JSON plus CSVs, reproducible from
   the recorded inputs under a fixed `SOURCE_DATE_EPOCH`, with an optional
   RO-Crate and an optional WIPO ST.26 export;
2. the **MCP `tools/list` contract** — `schema.json` plus annotations, the
   machine interface an agent reads before it calls anything.

A release is not “`R CMD check` is green”. It is all of:

```sh
Rscript dev/check.R                 # local gate, including docs/artifacts
Rscript dev/check-user-api.R        # 213 exports / 52 S3 / CSV / DOT
Rscript dev/check-mcp-parity.R      # tools / schema / handlers / inventory
Rscript dev/test.R --full           # the full suite
Rscript dev/check.R --full          # adds R CMD check --as-cran
```

Any red gate blocks the release. `check.R --full` already includes the test
suite and `R CMD check`; the other three are explicit because they prove
different contracts.

### 2. Artefact contract @2

The manifest moves to `schema = "tidymarkers/artefact@2"` and gains:

```json
"contract": {"version": "tidymarkers/contract@1", "scenario": "J6"},
"panel": { ... },
"files": {
  "solved": {"name": "solved.rds", "excluded": "igraph pointer and random id"}
}
```

Rules:

- `scenario` is one of `J1`–`J7` or `"unspecified"`. It is set by the caller
  or by the panel: `marker_artefact(scenario = "J6")` in J6, `"J3"` from the
  panel, and `"unspecified"` when a direct call does not name one. A wrong
  attribution is worse than an honest gap.
- `contract.version` is independent of the package version. It bumps when a
  field is removed or renamed, a field becomes required, or a verdict changes
  meaning; additive optional fields do not bump it.
- `inputs` hashes every input it records, including the target table. The
  current `artefact_inputs()` hashes the genome, variants and haplotypes; the
  target table is part of `plan` and is added to `inputs` with
  `role = "targets"`.
- `created` is excluded from `artefact_id()` as today. `scenario` and
  `contract.version` are included, because two scenarios are different
  products even when the solved plan is the same.
- The byte-identical set is exactly the text files:
  `artefact.json`, `decisions.csv`, `checks.csv`, `why.csv`, `targets.csv`,
  `signals.csv`. `solved.rds` is excluded by name and reason. A second writer
  must reproduce the text set, not the RDS.
- `marker_replay()` compares the per-target verdict and the report, not the RDS
  bytes. It succeeds when every target gets the same verdict under the recorded
  knobs; a different verdict is a failure with the target and both verdicts.

This is a breaking change to the manifest. It goes to `NEWS.md` and is the
first field a reader sees; `artefact@1` files remain readable by
`marker_artefact_read()` with a migration note.

### 3. MCP parity is one gate

`dev/check-mcp-parity.R` is the single source of truth for the machine
interface. It fails, with the offending name, when:

1. the tool names in `mcp_tool_specs()` differ from `schema.json`
   `$defs.inputSchemas` or from the handlers, or a published `inputSchema`
   differs from the runtime registry when both are rendered as JSON;
2. a schema `required` argument is not a formal of the handler, or a handler
   argument is not declared in the schema;
3. an `enum`, `minimum`, `maximum`, `pattern` or `additionalProperties` in the
   schema is weaker than the validator actually enforces;
4. the annotations disagree with `MCP_MUTATING_TOOLS` (only
   `marker_technology` is non-read-only);
5. a tool's input or output class is absent from
   `docs/user-api-surface.csv`, or the tool names a function that is not an
   export;
6. a J1–J7 scenario has no tool path and is not listed as `human_only` in the
   parity table.

The gate is static: it reads files, it does not start a server. Runtime
behaviour stays in `test-mcp.R` and `test-mcp-validation.R`: the former
compares all 38 published schemas with the registry, the latter sends
`tools/list` and `tools/call` for five representative tools (`primer_tm`,
`marker_schemes`, `design_markers`, `screen_snp_targets`, `marker_technology`)
and checks the wire result. A tool that is removed or renamed bumps
`contract.version` and is a breaking change.

### 4. No second interface

The following stay out of scope, consistent with ADR 0004:

- no Python client and no second schema dialect;
- no runtime read of `schema.json` by the server (the registry is
  `mcp_tool_specs()`; the file is a published contract);
- no hosted service and no network tool (`openWorldHint = FALSE`);
- no MCP tool that writes outside the artefact directory.

A new transport is a new ADR, not a new file in `inst/mcp/`.

### 5. Run metrics are collected from the first day

Metrics are a product contract, not a post-mortem. The first action of day 1
is the baseline run on `main`; every later run writes the same format. A
scenario run without a metrics line is not a run, and a daily slice without a
baseline is `not_run` with a reason, not a green gate.

**File and format.** One JSONL file per run:

```
data/metrics/<date>-<run>.jsonl
```

One JSON object per line, one line per stage event. The fields are fixed:

| Field | Type | Meaning |
| ----- | ---- | ------- |
| `run_id` | string | stable id of the run; deterministic for a reproducibility run |
| `stage` | enum | `scenario`, `ingest`, `design`, `specificity`, `thermo`, `checks`, `artefact`, `mcp`, `gate`, `repair` |
| `object_id` | string | stable id, artefact id or `scenario:JN`; never free text |
| `duration_ms` | integer | wall-clock duration of the stage |
| `status` | enum | `ok`, `warn`, `fail`, `skip`, `error` |
| `retries` | integer | attempts after the first |
| `tokens_in` | integer or null | LLM input tokens; null when no model is called |
| `tokens_out` | integer or null | LLM output tokens; null when no model is called |
| `cost_usd` | number or null | billed cost; `0` for local compute, null when unknown |
| `cache_hit` | boolean | the stage reused a shard, artefact or cached value |
| `queue_depth` | integer | tasks waiting at stage start |
| `language` | string | `en`, `ru` or `und`; the run language, not its text |

**No content.** Prompts, sequences, target labels, file paths and error text do
not enter the metrics. `object_id` is an id or a hash, not a label. The
validator rejects unknown keys, so a `prompt` or `content` field fails the
gate instead of being silently stored.

**Who writes what.** Package stages write `tokens_in = null`,
`tokens_out = null` and `cost_usd = 0`; they are local R. The agent/LLM wrapper
that drives the programme writes the same fields with real numbers. A route
that mixes both keeps one file and one schema.

**Daily slices.** `dev/metrics.R --report <date>` writes
`data/metrics/metrics-<date>.md` with:

- p50 and p95 of `duration_ms` per stage, nearest-rank; a stage with fewer than
  20 samples is reported as `insufficient_n`, not as a percentile;
- throughput: completed runs per hour and per day;
- cost per J1–J7 route: sum of `cost_usd` divided by the number of runs of that
  route;
- failure share: `(fail + error) / total` per stage and per route;
- cache-hit rate and p95 `queue_depth`.

The final report of the programme is `data/metrics/metrics-<date>.md`, with the
baseline delta for every stage and route.

**Gates.**

1. **p95 does not get worse.** For every stage with a baseline and at least 20
   samples, `p95(stage) <= p95_baseline(stage)`. A stage below 20 samples is
   `not_run`, with the count and the reason.
2. **Cost does not grow.** `cost_per_route(JN) <= cost_per_route_baseline(JN)`
   for J1–J7. When the baseline is zero, any billed cost fails the gate.
3. **Rerun is byte-identical.** The scenario output (artefact text files and
   CSVs) is byte-identical under a fixed `SOURCE_DATE_EPOCH`, as in §2. The
   raw JSONL cannot be byte-identical because `duration_ms` and `queue_depth`
   vary; the gate compares a canonical projection that drops those two fields
   and `cache_hit`, and requires the projection to be identical. The timing
   fields are compared by p50/p95, not by bytes. This is the honest form of
   "rerun is byte-identical"; a byte compare of timings would fail every
   correct rerun.
4. **Failure share is reported, not hidden.** A day with a higher failure share
   than the baseline fails the daily report even when p95 is green; the report
   names the stage and the route.

## Gates

1. **J6.** `test-scenario-j6.R` writes an artefact with a fixed epoch, reads it
   with `marker_artefact_read()`, exports an RO-Crate, validates an ST.26 file,
   and replays the artefact to the same verdicts. It asserts the `scenario` and
   `contract.version` fields.
2. **Byte identity.** Two runs with the same `SOURCE_DATE_EPOCH` produce
   identical text files; `solved.rds` is present, excluded by name and reason,
   and is not compared.
3. **Parity.** `Rscript dev/check-mcp-parity.R` exits 1 on each of the six
   conditions above; a deliberately removed tool and a deliberately added
   undeclared argument both fail the gate.
4. **Surface.** `Rscript dev/check-user-api.R` stays green; a new tool export
   without a CSV row fails the user-API gate, and the parity gate names it
   first.
5. **Release dry run.** `Rscript dev/check.R --full` plus the other three
   commands are green on a clean tree; the tarball builds and installs.
6. **Metrics.** `Rscript dev/metrics.R --check <date>` validates the day's
   JSONL against the schema and rejects content fields; `--report <date>`
   writes `data/metrics/metrics-<date>.md`. The p95, cost and rerun gates of §5
   are green, or `not_run` with the sample count and the reason.

## Consequences

- An artefact is attributable to a scenario and a contract version. A reader
  can tell a J3 panel from a J6 release without opening the RDS.
- `artefact@2` is a breaking manifest change; old files stay readable with a
  migration note.
- `dev/check-mcp-parity.R` becomes the fourth release gate and the only place
  that compares the 38 tools, the schema, the handlers and the inventory.
- A new MCP tool now costs a schema entry, a CSV row and a parity-table row.
  That is the intended cost: the machine interface is a product.
- Every run from day 1 leaves one JSONL file; the baseline is captured before
  the first code change, so p95 and cost deltas describe the programme, not the
  finished tree.
- The metrics file is numbers only. A prompt, a sequence or a target label in
  it is a schema violation, not a privacy review.
- No Python client, no hosted service and no runtime schema read are added.

## Alternatives

- **Keep `artefact@1` and add fields without a version bump.** Rejected: a
  reader cannot distinguish an old file from a new one, and `scenario` changes
  identity.
- **Put `scenario` only in the directory name.** Rejected: the manifest travels
  without its directory, and a copied JSON must still say what it is.
- **Extend `test-mcp.R` instead of a new gate.** Rejected: the test loads the
  package and runs handlers; parity must run before a package load and must
  compare the CSV, exactly as `check-user-api.R` does.
- **Make `schema.json` the runtime source.** Rejected by ADR 0004: the registry
  is R; the file is a published contract, and two sources would drift.
- **Compare `solved.rds` bytes.** Rejected: it serialises an igraph pointer and
  a random id; the verdict comparison is the real replay.
- **Collect metrics at the end of the three days.** Rejected: it measures the
  finished tree and cannot show a regression introduced on day 2. The baseline
  must be captured before the first code change.
- **Store prompts or error text in the metrics.** Rejected: the contract is
  numbers only; content in a metrics file is a privacy and size problem, and
  the validator rejects it.

## Blocker

Status stays `proposed` because gate 5 (the release dry run) has not been run
since this ADR was written. Gate 5 needs `dev/check.R --full`, which runs
`R CMD check --as-cran`; the 2026-11-28 night shift may not run the full suite.
Gate 6 was run on 2026-11-28: the day's JSONL passes
`dev/metrics.R --check 2026-11-28` (29 rows) and `--report` writes
`data/metrics/metrics-2026-11-28.md`. Its §5 comparison cells are `not_run`
with n = 0 and the reason "no current rows", because no J1–J7 route ran that
day; the day's rows are all `stage = "repair"`. That is the named-`not_run`
form gate 6 allows, not a green p95 or cost claim. Gates 1–4 are green:
`test-scenario-j6.R` passes, `dev/check-mcp-parity.R` prints `ok mcp-parity`,
and `dev/check-user-api.R` prints `ok user-api` at 213 exports / 52 S3 methods.
The ADR becomes `accepted` in the commit that runs the dry run; until then
`proposed` names a specific unmet gate, per plan §8.

## Revisit trigger

A second consumer that cannot read R and needs the artefact schema as a
standalone JSON Schema, or a tool that changes the artefact directory outside
`marker_artefact()`. At that point the contract version, not the package
version, is the thing to re-open.

## References

- [`R/marker_artefact.R`](../../R/marker_artefact.R) — the manifest, text files
  and `marker_replay()`.
- [`R/marker_export_rocrate.R`](../../R/marker_export_rocrate.R) — the RO-Crate
  writer.
- [`R/marker_patent.R`](../../R/marker_patent.R) — the ST.26 export and
  validation.
- [`inst/mcp/tools.R`](../../inst/mcp/tools.R) — the 38 tool specs and handlers.
- [`inst/mcp/schema.json`](../../inst/mcp/schema.json) — the published
  `tools/list` contract.
- [`dev/mcp-schema.R`](../../dev/mcp-schema.R) — the generator and `--check`.
- [`tests/testthat/test-mcp.R`](../../tests/testthat/test-mcp.R) — the runtime
  conformance suite.
- [`docs/user-api-surface.md`](../user-api-surface.md) — the export and class
  inventory the parity gate reads.
- [`ROADMAP.md`](../../ROADMAP.md) — `ZT T8`, `T10`, §ZV (ZV9, ZV10).

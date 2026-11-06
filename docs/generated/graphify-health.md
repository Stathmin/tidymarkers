# Graphify health (generated)

Run `Rscript dev/graphify-health.R --write` after `graphify update .`.
Do not edit by hand.

The code graph is a **limitation of orientation, not proof of correctness**: an empty graph and clean code look the same.

## Structural invariants

These must stay at zero. `dev/graphify-health.R --check` fails if one moves.

| metric | value |
| --- | ---: |
| dangling endpoints | 0 |
| exact duplicate edges | 0 |
| same-endpoint (collapsed) edges | 0 |
| import cycles | 0 |
| self-loops outside `.tools/` | 0 |

## Known, accepted

Not zero, and not gated: the cause is understood and lives outside this package.

| metric | value | status |
| --- | ---: | --- |
| self-loops | 1 | accepted: one recursive `collectEdges()` self-edge in the vendored `.tools/elk/elk-layout.js`; `--check` fails on any self-loop outside `.tools/` |

## Community labels

A label that echoes a filename is honest for a single-file community and misleading for a multi-file one. `--check` warns (never fails) when any label is `missing`, `duplicate` or `echo-stale`.

| status | communities |
| --- | ---: |
| `missing` | 0 |
| `duplicate` | 0 |
| `echo-stale` | 0 |
| `echo-honest` | 0 |
| `named` | 401 |

No missing, duplicate or multi-file file-echo labels.

### Relabel provenance

| date | model | requested | relabelled | kept | note |
| --- | --- | ---: | ---: | ---: | --- |
| 2026-11-06 | `Ternary-Bonsai-2-27B-tight` | 74 | 69 | 5 |  |
| 2026-11-06 | `Ternary-Bonsai-2-27B-tight` | 5 | 5 | 0 |  |
| 2026-11-06 | `Ternary-Bonsai-2-27B-tight` | 3 | 3 | 0 |  |
| 2026-11-06 | `Ternary-Bonsai-2-27B-tight` | 251 | 249 | 2 |  |
| 2026-11-06 | `Ternary-Bonsai-2-27B-tight` | 2 | 2 | 0 |  |
| 2026-11-06 | `manual` | 1 | 1 | 0 | cid 408: orphan external yaml node; the model kept proposing a rejected placeholder, fixed by hand |
| 2026-11-06 | `Ternary-Bonsai-2-27B-tight` | 4 | 2 | 2 | cids 210/217 re-asked with --cids after the unchanged-check fix; the model repeated the same honest single-file name, so both stay |
| 2026-11-06 | `Ternary-Bonsai-2-27B-tight` | 38 | 37 | 1 | after the .graphifyignore rebuild renumbered communities and dropped 38 label keys; 37 named here, cid 264 fixed by hand |
| 2026-11-06 | `manual` | 1 | 1 | 0 | cid 264: single-file SECURITY.md community; model proposed the member echo 'Security Policy', fixed by hand |
| 2026-11-06 | `Ternary-Bonsai-2-27B-tight` | 4 | 4 | 0 | targeted --cids=190,210,217,246; thin but honest names deepened to say what the cluster is for |

## Size

| metric | value |
| --- | ---: |
| nodes | 7 754 |
| edges | 10 361 |
| communities | 401 |
| nodes without `source_file` | 368 |

## Nodes without a `source_file`

Two kinds, both expected. `external_function` is the AST record of a call
into another package: it has no file here, and `explain` on one shows its
callers, not its definition. `package_function` is published deliberately by
`dev/graphify-handrolled.R` for a hand-roll's replacement; when the same
`pkg::fn` already exists as an AST node the injection reuses that node, so
`explain "pkg::fn"` resolves to one node and shows the `supersedes` edge.
The total is a drift signal, not a defect.

| kind | count |
| --- | ---: |
| external_function | 297 |
| package | 55 |
| concept | 14 |
| external_class | 1 |
| package_function | 1 |

## God nodes

Defined here as the top targets by **in-degree over `calls`**: how often the
rest of the code calls this. (Graphify's own ranking mixes document nodes and
is not reproduced.)

| node | callers |
| --- | ---: |
| `tibble::tibble()` | 180 |
| `cli::cli_abort()` | 105 |
| `glue::glue()` | 59 |
| `dplyr::mutate()` | 48 |
| `dplyr::filter()` | 46 |
| `purrr::list_rbind()` | 46 |
| `Biostrings::DNAStringSet()` | 41 |
| `utils::read.csv()` | 40 |
| `purrr::map()` | 38 |
| `fs::path()` | 36 |

## Known gaps

Adopted from GWASILIY's graph practice.

| gap | status | note |
| --- | --- | --- |
| hand-roll ledger | data | `dev/handrolled.yml` holds the decisions; `dev/check-handrolled.R` reads it |
| `supersedes` graph edges | adopted | `dev/graphify-handrolled.R` injects one `supersedes` edge per entry, reusing the AST `pkg::fn()` node when it exists and creating a `package::function` node otherwise; idempotent, and `dev/check.R` runs it |
| default label prompt on rebuild | guarded | `graphify update .` re-runs the weak built-in prompt on changed communities; the strict pass and the origin ledger record and repair it |
| `dev/relabel-bonsai.R` helper | shipped | repeatable strict pass; defaults to missing/duplicate/echo-stale, `--all-echo` also renames the honest single-file echoes, `--cids=` repairs one name, and a proposal equal to the current label is kept, not counted as a rename. A rebuild that renumbers communities drops label keys; the pass restores them as `missing` and names them |
| single-file echo labels (254) | skipped | a file-echo name is honest for a one-file community; the strict pass renames only missing/duplicate/echo-stale, `--all-echo` opts in |

## Configuration levers

What is wired into graphify here, and what is deliberately not. `open` is usable in the installed build but not adopted; `skipped` has no consumer.

| lever | status | note |
| --- | --- | --- |
| `graphify update .` | shipped | code-only re-extract, no LLM; first step of Rerun below |
| strict label pass + origin ledger | shipped | `dev/relabel-bonsai.R`; all communities named, every pass recorded |
| `supersedes` edges from the hand-roll ledger | shipped | see Known gaps; answers "what replaced this hand-roll" |
| local LLM semantic backend | open | `dev/graphify-semantic.sh` holds the safe flags; `extract --backend=openai` reads `OPENAI_BASE_URL`/`OPENAI_MODEL` (bonsai). Measured: the served model has an 8192-token context and graphify asks for up to 16384 output tokens, so a full 63-document pass truncates and splits for hours. Opt-in until a faster/larger-context model is available |
| `export wiki` | shipped | 411 articles in gitignored `graphify-out/wiki/`, one per community plus an index; template-only, no LLM, so its quality is the curated labels |
| git hook (`graphify hook install`) | shipped | post-commit/checkout/merge/rewrite plus the graph merge driver; it runs `update` only (no LLM), so `dev/graphify-health.R --write` stays manual |
| `.graphifyignore` | shipped | `data-raw/proposed/` (37 files, 335 nodes), `results/`, `.tmp/` are design scratch; a rebuild prunes them, and it can only ever add exclusions |
| `--mode deep` | skipped | richer INFERRED edges, but a full re-extract and it re-names changed communities with the weak built-in prompt; no consumer here |
| `query` / `path` / `explain` / `affected` | shipped | verified from the shell: `explain "Biostrings::complement"` resolves one node and shows its `supersedes` edge, `affected "normalise_allele()"` lists its caller. 0.9.92 has no MCP server, so agents use this CLI |
| `export neo4j` / `graphml` / `obsidian` / `svg` | skipped | no consumer in this repository |

## Provenance

Built at commit `fd0780c1d4eb9e0896098c5a98cc6caf20c3d315`.
The graph matches `HEAD`.

## Rerun

```sh
graphify update .
Rscript dev/graphify-handrolled.R
Rscript dev/graphify-health.R --write
graphify export wiki
graphify export html
# when a label warns: Rscript dev/relabel-bonsai.R [--all-echo]
# one wrong name: Rscript dev/relabel-bonsai.R --cids=N
# semantic pass (opt-in, local LLM): dev/graphify-semantic.sh
# refresh the analysis sidecar: graphify cluster-only . --no-label
```

Do not pass `--no-cluster` to `update`: it drops the community
assignment and duplicates some edges, so the health numbers above
mislead. `cluster-only` is deterministic on an unchanged graph but can
renumber after a rebuild, so re-check the labels and rerun
`dev/relabel-bonsai.R` if it warns.

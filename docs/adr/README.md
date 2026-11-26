# Architecture decisions

One file per decision that changes what this repository depends on. Numbered in
the order they were accepted. `docs/` is `.Rbuildignore`d, so nothing here ships
in the package.

Status values used here: `accepted`, `rejected` (the proposal was refused),
`superseded`.

| ADR                                                                    | Decision                                                                    | Status   |
| ---------------------------------------------------------------------- | --------------------------------------------------------------------------- | -------- |
| [0001](0001-llm-mt-and-agent-tooling-out-of-scope.md)                  | LLM, MT and agent tooling is out of scope                                   | accepted |
| [0002](0002-targets-stays-dev-only.md)                                 | `targets` stays a dev-only freshness holder, not a package dependency       | accepted |
| [0003](0003-duckdb-deferred.md)                                        | DuckDB/`duckplyr` is deferred until an aggregation is measured              | accepted |
| [0004](0004-mcp-ci-and-cross-project-tooling.md)                       | The MCP server stays one R file; no Python client, no published JSON Schema | accepted |

## Triage of the 2026-11-26 tool list

Twenty-five tools were proposed for adoption. The verdicts are below; "no call
site" means the repository has nothing that would run it, which is the condition
`CONTRIBUTING.md` refuses (a dependency is accepted when existing code reaches a
caller).

| #  | Tool                        | Verdict         | Reason                                                                                       |
| -- | --------------------------- | --------------- | -------------------------------------------------------------------------------------------- |
| 1  | DSPy                        | Reject          | No prompt anywhere; rubric optimisation has no object to optimise (ADR 0001).                |
| 2  | vLLM / sglang               | Reject          | No model is served; the package runs no inference (ADR 0001).                                |
| 3  | Instructor / Outlines       | Reject          | Nothing is decoded; MCP input is validated against a JSON Schema in R (ADR 0001).            |
| 4  | promptfoo                   | Reject          | No prompts to regress; the gate is `dev/check.R`.                                            |
| 5  | OpenTelemetry + Phoenix     | Reject          | `Rprof` is the profiling standard here (`AGENTS.md`); there is no request or span to trace.  |
| 6  | `targets`                   | Already in use  | `dev/targets/_targets.R` holds artifact freshness; stays dev-only (ADR 0002).                |
| 7  | Letta                       | Reject          | No agent runtime; there are no sessions to persist.                                          |
| 8  | DeepPavlov                  | Reject          | No Russian text; the package processes DNA.                                                   |
| 9  | NLLB-200                    | Reject          | No translation surface.                                                                      |
| 10 | ragas / DeepEval            | Reject          | No RAG and no LLM output to score.                                                            |
| 11 | mcp-use / mcp-hub / FastMCP | Reject          | The server is R; `tests/testthat/test-mcp.R` is its conformance suite (ADR 0004).            |
| 12 | Shared Pydantic registry    | Reject          | The canonical schemas are R: `marker_trace_schema()`, `mcp_tool_specs()` (ADR 0004).         |
| 13 | uv workspaces               | Reject          | No Python runtime, only a dev cross-check script (ADR 0004).                                 |
| 14 | nektos/act                  | Reject          | `dev/check.R` and CI already run the same gate (ADR 0004).                                   |
| 15 | lm-eval-harness             | Reject          | No model to evaluate.                                                                        |
| 16 | COMET-Kiwi                  | Reject          | No translation output.                                                                        |
| 17 | xCOMET                      | Reject          | No translation output.                                                                        |
| 18 | inspect                     | Reject          | No trajectories; there is no agent.                                                           |
| 19 | OpenHands                   | Reject          | Editor-side tool, not a repository dependency.                                                |
| 20 | aider                       | Reject          | Same: the git workflow belongs to the operator, not to the package.                           |
| 21 | AgentOps                    | Reject          | No agents to observe.                                                                         |
| 22 | duckplyr / DuckDB           | Defer           | Dev-only analytics; the trigger is not met (ADR 0003).                                        |
| 23 | Quarto                      | Already in use  | `_quarto.yml`, `dev/lecture/`; no change.                                                    |
| 24 | Apertium                    | Reject          | No translation surface.                                                                      |
| 25 | LibreTranslate              | Reject          | No translation surface.                                                                      |

The tools in rows 1–5, 7–10 and 15–21 are not wrong tools. They belong to the
consumer of this package: a pipeline that drives these MCP tools with a model
can optimise its prompts or constrain its decoding there, against the tool
schemas published here. This repository owes that consumer a stable, tested
contract — `tests/testthat/test-mcp.R` — not a copy of the consumer's stack.

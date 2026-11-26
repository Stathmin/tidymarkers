# 0004 — Keep the MCP server in R, the gate local, and the schema unpublished

- Status: accepted
- Date: 2026-11-26
- Supersedes: none

## Context

The MCP surface here is one R file of tool specs (`inst/mcp/tools.R`) and a thin
stdio loop (`inst/mcp/server.R`: one JSON-RPC line in, `jsonlite` in/out,
`mcp_dispatch()` for `initialize`, `tools/list`, `tools/call`). Its conformance
suite is `tests/testthat/test-mcp.R`: the `tools/list` shape, one call per
spec, and enum membership checked against the exported registries.

Two canonical schemas exist in this package and both are R objects:
`marker_trace_schema()` (the versioned artefact record) and `mcp_tool_specs()`
(the tool registry). The only Python is `dev/as-tm-crosscheck.py`, a dev
cross-check of melting temperatures.

The proposal offered `mcp-use` and `mcp-hub`, FastMCP, `nektos/act`, uv
workspaces, a shared Pydantic schema registry, and Quarto. The stated aim was to
compose three projects' MCP servers into one agent surface.

## Decision

- **No Python MCP client, hub or FastMCP.** The server is R. A second client
  implementation in another language would be a second definition of
  `tools/call`, and `tests/testthat/test-mcp.R` is already the contract. A cross
  project composition, if it happens, reads this server over stdio like any
  other client; it does not need a copy of it.
- **No `nektos/act`.** `dev/check.R` is the local gate (style, lint, roxygen
  sync, docs/artifacts, figure freshness, install smoke) and
  `.github/workflows/R-CMD-check.yaml` runs `R CMD check --as-cran` on push and
  pull request. `act` would run a containerised subset of the same gate and add
  a Docker requirement to a repository whose checks already run offline on a
  pinned R 4.6.1 tree.
- **No uv workspace and no shared Pydantic registry.** There is no Python
  workspace to unify. The overlapping entities such a registry would hold (DOI,
  Person, Work, Species, Gene) are not what this package stores; the entities it
  does own are R-side and versioned (`marker_trace_schema()`,
  `"tidymarkers/trace@1"`).
- **No published JSON Schema yet.** `test-mcp.R` holds the shape without a file
  on disk. `inst/mcp/schema.json` becomes worth emitting only when a consumer
  must validate tool specs without loading R.
- **Quarto: no change.** It is already the docs/lecture stack (`_quarto.yml`,
  `dev/lecture/`).

## Consequences

- The MCP surface stays offline, stdio-only and testable with the package's own
  suite; no new runtime or container dependency.
- A cross-project consumer gets the tool specs through `tools/list`, not through
  a generated schema file; if it needs the latter, see the trigger.
- CI does not gain a second, slower path that can disagree with `dev/check.R`.

## Alternatives

- Emit `inst/mcp/schema.json` now, "for future consumers". Rejected: a generated
  file with no reader is an artifact the gate must then keep in sync, and
  `test-mcp.R` already fails on the same drift.
- Port the server to Python behind FastMCP. Rejected: it would replace a working
  R implementation and its test with a second language in a package that has no
  Python runtime.

## Revisit trigger

A second consumer of the MCP contract that cannot load R. Then emit
`inst/mcp/schema.json` from `mcp_tool_specs()` and add a round-trip test; the
server itself stays R.

## References

- [`inst/mcp/tools.R`](../../inst/mcp/tools.R) — the tool registry.
- [`inst/mcp/server.R`](../../inst/mcp/server.R) — the stdio loop.
- [`tests/testthat/test-mcp.R`](../../tests/testthat/test-mcp.R) — the contract.
- [`dev/check.R`](../../dev/check.R) — the local gate.
- [`_quarto.yml`](../../_quarto.yml) — the docs/lecture stack.

# 0001 — LLM, machine-translation and agent tooling is out of scope

- Status: accepted
- Date: 2026-11-26
- Supersedes: none

## Context

`tidymarkers` designs, checks and audits plant SNP markers. Its runtime is
exported R functions over Bioconductor: sequence handling (`Biostrings`), ranges
(`GenomicRanges`), nearest-neighbour thermodynamics (`TmCalculator`,
OligoArrayAux) and exact/approximate matching for specificity. Every verdict is
deterministic: the same input and the same criteria table produce the same
artefact, and `marker_trace_schema()` (`"tidymarkers/trace@1"`) records the
evidence behind it.

The package exposes that API over MCP as an offline stdio server
(`inst/mcp/tools.R`, `inst/mcp/server.R`): a registry of tool specs, each with a
JSON-Schema input and an R handler. It answers JSON-RPC 2.0 for `initialize`,
`tools/list` and `tools/call`. It calls no model, no network service and no GPU.

Python appears once, as a dev cross-check of melting temperatures against
`Bio.SeqUtils.MeltingTemp` (`dev/as-tm-crosscheck.py`, driven by
`dev/as-tm-crosscheck.R`). It is not a runtime.

The 2026-11-26 proposal offered eighteen tools aimed at LLM application work:
DSPy, vLLM/sglang, Instructor/Outlines, promptfoo, OpenTelemetry+Phoenix, Letta,
DeepPavlov, NLLB-200, ragas/DeepEval, lm-eval-harness, COMET-Kiwi, xCOMET,
inspect, OpenHands, aider, AgentOps, Apertium, LibreTranslate.

## Decision

None is adopted. Each would be a dependency with no call site, which is the
condition `CONTRIBUTING.md` refuses under "Wiring": a dependency is accepted
when existing code reaches a caller.

The reasons are specific to this repository:

1. Prompt optimisation (DSPy) and constrained decoding (Instructor, Outlines)
   need a prompt or a decoder. There is neither. MCP input is validated against
   a schema in R and the handler returns an R object; no token sampling happens
   in this package.
2. Serving (vLLM, sglang) and output scoring (lm-eval-harness, ragas, DeepEval,
   COMET-Kiwi, xCOMET, inspect) need model output. There is none.
3. Translation (NLLB-200, Apertium, LibreTranslate) and Russian NLP (DeepPavlov)
   need text. The package processes DNA.
4. Agent memory and observation (Letta, AgentOps) and agent IDEs (OpenHands,
   aider) assume an agent runtime, which the package does not have.
5. Tracing (OpenTelemetry + Phoenix) duplicates a standard already fixed in
   `AGENTS.md`: `Rprof`, then read the innermost frame together with its parent.
   A package that builds S4 objects has profiles, not spans.

Most of these belong to the consumer. A pipeline that drives these MCP tools
with a model — retrieval, extraction, translation — can use DSPy or Outlines
there, against the tool schemas this package publishes. That decision is that
repository's, and this package's obligation is a stable, tested tool contract,
which `tests/testthat/test-mcp.R` already enforces.

## Consequences

- The dependency surface stays what `DESCRIPTION` declares. Nothing is installed.
- The MCP server keeps working offline and stays one R file plus its test.
- If a model call ever lands inside the package, this ADR has to be superseded,
  not quietly contradicted.

## Alternatives

- Adopt a subset "for later". Rejected: an unused dependency is not a
  contribution here, and it ages against the R 4.6.1 / Bioconductor 3.23 lock.
- Adopt only the evaluation frameworks, to score the decision graph. Rejected:
  the graph's output is a verdict plus a cause chain, held by `testthat`
  assertions, not sampled text.

## Revisit trigger

A merged feature in which the package itself calls a model, serves inference or
publishes generation prompts. Until then none of the eighteen has a caller here.

## References

- [`inst/mcp/tools.R`](../../inst/mcp/tools.R) — the tool registry.
- [`tests/testthat/test-mcp.R`](../../tests/testthat/test-mcp.R) — the contract.
- [`CONTRIBUTING.md`](../../CONTRIBUTING.md) — the "Wiring" rule.
- [`R/marker_trace_schema.R`](../../R/marker_trace_schema.R) — the versioned record.
- [`AGENTS.md`](../../AGENTS.md) — the profiling standard.

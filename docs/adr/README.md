# Architecture decisions

One file per decision that changes what this repository depends on, or how a
module is bounded. Numbered in the order they were accepted. `docs/` is
`.Rbuildignore`d, so nothing here ships in the package.

Status values used here: `proposed` (the decision is written, the code gate is
not met), `accepted`, `rejected` (the proposal was refused), `superseded`,
`accepted (amended)` (the decision stands and a dated amendment at the end of
the file records the measurement or event that re-tested it).

| ADR                                                                    | Decision                                                                    | Status   |
| ---------------------------------------------------------------------- | --------------------------------------------------------------------------- | -------- |
| [0001](0001-llm-mt-and-agent-tooling-out-of-scope.md)                  | LLM, MT and agent tooling is out of scope                                   | accepted |
| [0002](0002-targets-stays-dev-only.md)                                 | `targets` stays a dev-only freshness holder, not a package dependency       | accepted (amended) |
| [0003](0003-duckdb-deferred.md)                                        | DuckDB/`duckplyr` is deferred until an aggregation is measured              | accepted (amended) |
| [0004](0004-mcp-ci-and-cross-project-tooling.md)                       | The MCP server stays one R file; no Python client, no published JSON Schema | accepted (amended) |
| [0005](0005-index-backed-genome-access.md)                             | Genome access is an index-backed strategy with declared capabilities       | accepted |
| [0006](0006-panel-pipeline.md)                                         | The panel pipeline is an orchestrator with falsifiable seams                | accepted |
| [0007](0007-module-seams-and-mcp-contract.md)                          | Module seams (contracts, registries, splits) and the MCP runtime contract   | proposed |
| [0008](0008-user-scenarios-as-contract.md)                             | User scenarios J1–J7 are the executable contract                            | accepted |
| [0009](0009-panel-orchestrator-contract.md)                            | The panel orchestrator: GWAS, LD, haplotypes, chemistry, artefact           | accepted |
| [0010](0010-index-strategy-at-scale.md)                                | The index strategy at scale: shards, builds and a memory budget            | accepted |
| [0011](0011-release-artefact-and-mcp-contract.md)                      | The release artefact and the MCP contract are the product                   | proposed |
| [0012](0012-autonomy-and-self-repair.md)                               | Autonomy and self-repair: a blocker never stops the programme               | proposed |
| [0013](0013-metrics-gate-is-not-a-green-delta.md)                      | The metrics report is not a gate; a mixed-sample p95 is not a regression     | accepted |
| [0014](0014-benchmark-baseline-is-keyed-by-runner-class.md)            | A timing baseline is keyed by runner class, not by hostname                  | accepted |
| [0015](0015-a-record-resolves-only-its-own-directory.md)               | A record resolves only the paths inside its own directory                    | accepted |
| [0016](0016-one-vcf-admission-rule.md)                                 | One VCF admission rule, shared by both readers                               | accepted |
| [0017](0017-a-raw-identifier-is-admitted-before-normalisation.md)       | A raw identifier is admitted before a normaliser may rename it               | accepted |
| [0018](0018-one-allele-and-position-rule-for-three-vcf-readers.md)      | One allele and `POS` rule for all three VCF readers                          | accepted |
| [0019](0019-an-absent-allele-is-not-a-base.md)                         | An absent allele is absent before it is normalised                           | accepted |
| [0020](0020-a-gt-that-names-a-missing-allele-is-not-a-call.md)          | A GT that names a missing allele is not a call                               | accepted |
| [0021](0021-one-path-and-one-header-rule-for-three-vcf-readers.md)      | One path and one header rule for three VCF readers                           | accepted |
| [0022](0022-a-fasta-base-outside-the-alphabet-is-refused-by-every-reader.md) | A FASTA base outside the alphabet is refused by every reader                 | accepted |
| [0023](0023-a-refusal-that-a-retry-cannot-fix-carries-a-class.md)       | A refusal that a retry cannot fix carries a class                            | accepted |
| [0024](0024-a-replay-compares-only-the-records-own-verdicts.md)        | A replay compares only the record's own verdicts                            | accepted |
| [0025](0025-a-record-member-is-resolved-not-only-named.md)             | A record's member is resolved, not only named                               | accepted |
| [0026](0026-an-artefacts-id-is-the-hash-of-what-it-stores.md)          | An artefact's id is the hash of what it stores                              | accepted |
| [0027](0027-a-record-directory-holds-one-design.md)                    | A record directory holds one design                                          | accepted |
| [0028](0028-a-vcf-writer-emits-only-records-its-readers-admit.md)      | A VCF writer emits only records its readers admit                            | accepted |
| [0029](0029-an-output-writer-creates-its-files-not-writes-through-a-link.md) | An output writer creates its files, not writes through a link          | accepted |
| [0030](0030-a-derived-path-component-is-a-declared-name.md)             | A derived path component is a declared name                                 | accepted |
| [0031](0031-a-named-root-is-normalised-before-a-sibling-is-derived.md)  | A named root is normalised before a sibling is derived from it              | accepted |
| [0032](0032-the-corpus-exclusion-check-reads-the-resolved-directory.md) | The corpus-exclusion check reads the resolved directory                     | accepted |
| [0033](0033-a-cache-and-a-rendered-file-are-declared-names.md) | A cache entry and a rendered file are declared names too | accepted (amended) |
| [0034](0034-an-index-build-declares-every-file-it-writes.md) | A genome index declares every file its build writes | accepted |
| [0035](0035-a-guard-refuses-a-name-a-hit-test-only-reports-one.md) | A guard refuses a name; a hit test only reports one | accepted (amended) |
| [0036](0036-a-nested-closure-is-a-writer-of-its-own.md) | A nested closure is a writer of its own | accepted |
| [0037](0037-the-write-vocabulary-names-every-primitive.md) | The write vocabulary names every primitive that creates a file | accepted |
| [0038](0038-a-depth-the-writer-refuses-is-not-a-depth-the-reader-reports.md) | A depth the writer refuses is not a depth the reader reports | accepted |
| [0039](0039-a-negative-depth-carries-no-confidence.md) | A negative depth carries no confidence, not negative confidence | accepted |
| [0040](0040-a-validator-enforces-the-unit-the-contract-names.md) | A validator enforces the unit the contract names | accepted |
| [0041](0041-a-genotype-is-its-allele-multiset.md) | A genotype is its allele multiset, not its spelling | accepted |
| [0042](0042-a-threshold-is-a-count-not-a-number-to-truncate.md) | A threshold is a count, not a number to truncate | accepted |
| [0043](0043-a-panel-gives-one-call-per-line-and-target.md) | A panel gives one call per line and target | accepted |
| [0044](0044-a-count-has-one-definition.md) | A count has one definition | accepted |
| [0045](0045-a-position-is-a-whole-number-before-the-cast.md) | A position is a whole number before the cast | accepted |
| [0046](0046-a-design-key-names-every-input-the-design-reads.md) | A design key names every input the design reads | accepted |
| [0047](0047-a-search-radius-is-a-whole-number-before-the-cast.md) | A search radius is a whole number before the cast | accepted |
| [0048](0048-a-base-pair-threshold-is-a-whole-number-before-the-cast.md) | A base-pair threshold is a whole number before the cast | accepted |
| [0049](0049-a-worker-count-is-a-whole-number-before-the-cast.md) | A worker count is a positive whole number before the cast | accepted |
| [0050](0050-a-genomic-interval-is-whole-before-the-read.md) | A genomic interval is whole before the read | accepted |
| [0051](0051-a-nucleotide-length-is-whole-before-the-design-grid.md) | A nucleotide length is whole before the design grid | accepted |
| [0052](0052-a-palindrome-count-is-whole-before-the-scan.md) | A palindrome's counts are whole before `findPalindromes()` | accepted |
| [0053](0053-a-sample-selector-names-a-row-it-does-not-round.md) | A sample selector names a row, it does not round | accepted |
| [0054](0054-a-stabilising-mismatch-is-a-whole-distance-inside-the-primer.md) | A stabilising mismatch is a whole distance inside the primer | accepted |
| [0055](0055-a-common-primer-search-names-whole-lengths-gaps-and-anchor.md) | A common-primer search names whole lengths, gaps and anchor | accepted |
| [0056](0056-the-shared-primer-grid-names-whole-lengths.md) | The shared primer grid names whole lengths | accepted |
| [0057](0057-a-lamp-length-is-checked-before-the-swallowed-grid.md) | A LAMP length is checked before the swallowed grid | accepted |
| [0058](0058-a-repeat-scan-names-whole-units-and-copies.md) | A repeat scan names whole units and copies | accepted |
| [0059](0059-an-as-primer-anchor-is-whole-before-the-window.md) | An AS-primer anchor and linked offset are whole before the window | accepted |
| [0060](0060-a-target-window-names-its-whole-locus.md) | A target window names its whole locus before the cut | accepted |

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
| 22 | duckplyr / DuckDB           | Defer           | Dev-only analytics; measured 2026-11-26, `dplyr` wins (ADR 0003).                                        |
| 23 | Quarto                      | Already in use  | `_quarto.yml`, `dev/lecture/`; no change.                                                    |
| 24 | Apertium                    | Reject          | No translation surface.                                                                      |
| 25 | LibreTranslate              | Reject          | No translation surface.                                                                      |

The tools in rows 1–5, 7–10 and 15–21 are not wrong tools. They belong to the
consumer of this package: a pipeline that drives these MCP tools with a model
can optimise its prompts or constrain its decoding there, against the tool
schemas published here. This repository owes that consumer a stable, tested
contract — `tests/testthat/test-mcp.R` — not a copy of the consumer's stack.

## Amendments (2026-11-26)

The three deferred follow-ups from the triage were executed; none changed its
decision, so each ADR keeps `accepted` with a dated amendment at the end:

- [0002](0002-targets-stays-dev-only.md#amendment-2026-11-26-the-trigger-was-measured-and-it-was-not-met) —
  step cost measured; no step nears ten minutes, so `targets` is not widened.
- [0003](0003-duckdb-deferred.md#amendment-2026-11-26-the-ladder-was-run-dplyr-wins-at-every-scale) —
  aggregation ladder run; `dplyr` wins at every scale.
- [0004](0004-mcp-ci-and-cross-project-tooling.md#amendment-2026-11-26-the-schema-is-published-by-operator-decision) —
  `inst/mcp/schema.json` published by operator decision plus a drift test.

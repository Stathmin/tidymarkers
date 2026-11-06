#!/bin/sh
# =============================================================================
#  dev/graphify-semantic.sh  --  opt-in semantic pass with the local LLM
# =============================================================================
#  Graphify's AST pass is local and deterministic; this is the only step that
#  calls a model.  It adds document "concept" nodes and INFERRED edges.  Run
#  `graphify update .` first and re-run this script, not `update`, for prose
#  changes.
#
#  Measured on this machine (2026-11-06, Ternary-Bonsai-2-27B-tight,
#  llama-server -c 65536 -np 2, i.e. 32768 tokens per slot, 2 slots):
#    * ~32 output tokens/second; a ~3.3 KB document prompt takes ~90 s and
#      ~3000 output tokens, finish_reason=stop, no reasoning tokens;
#    * --token-budget 6000 fits one slot without truncation;
#    * --max-concurrency 2 uses both slots.
#  A full 54-document pass is still hours; scope it to a subtree by passing a
#  path, e.g. `docs`.
#
#  NOTE: graphify indexes .md/.txt/.yaml/etc but NOT .Rmd or .Rd, so
#  vignettes/ and man/ are invisible to this pass.  Convert or export first if
#  those must enter the graph.
#
#  Run:  dev/graphify-semantic.sh [PATH]   (default: .)
# =============================================================================
set -eu

export OPENAI_BASE_URL="${OPENAI_BASE_URL:-http://127.0.0.1:8000/v1}"
export OPENAI_MODEL="${OPENAI_MODEL:-Ternary-Bonsai-2-27B-tight}"
export OPENAI_API_KEY="${OPENAI_API_KEY:-local}"
export GRAPHIFY_MAX_OUTPUT_TOKENS="${GRAPHIFY_MAX_OUTPUT_TOKENS:-4000}"
export GRAPHIFY_API_TIMEOUT="${GRAPHIFY_API_TIMEOUT:-600}"

exec graphify extract "${1:-.}" \
  --backend=openai \
  --max-concurrency 2 \
  --token-budget 6000 \
  --out .

#!/bin/sh
# =============================================================================
#  dev/graphify-semantic.sh  --  opt-in semantic pass with the local LLM
# =============================================================================
#  Graphify's AST pass is local and deterministic; this is the only step that
#  calls a model.  It adds document "concept" nodes and INFERRED edges.  Run
#  `graphify update .` first and re-run this script, not `update`, for prose
#  changes.
#
#  Measured on this machine (2026-11-10, Ternary-Bonsai-2-27B-tight,
#  llama-server -c 65536 -np 2, i.e. 32768 tokens per slot, 2 slots):
#    * the server ignores graphify's GRAPHIFY_DISABLE_THINKING shape
#      ({"thinking":{"type":"disabled"}} still emitted 105 tokens of reasoning
#      prose); `chat_template_kwargs.enable_thinking=false` and
#      `reasoning_effort: low` both returned clean JSON in 9-14 output tokens;
#    * a large doc (18-28 KB) needs more than 4000 output tokens: at that cap
#      graphify logs "truncated at max_completion_tokens", splits the chunk and
#      retries, and keeps partial results uncached, so the pass never converges.
#      16000 fits one slot (input budget 12000 + output);
#    * --max-concurrency 2 uses both slots.
#  The thinking-off shape is a custom provider, registered once in the trusted
#  global ~/.graphify/providers.json as `bonsai` (extra_body + reasoning_effort).
#  Pass a path to scope the pass, e.g. `docs`; a path-scoped pass stamps only
#  that subtree in the manifest.
#
#  NOTE: `vignettes/*.Rmd` **are** indexed -- `.graphifyrc` maps `.Rmd` to
#  markdown.  `man/*.Rd` is not: use `dev/graphify-docs.R` to export plain `.md`
#  stubs first if the manuals must enter this pass.
#
#  Run:  dev/graphify-semantic.sh [PATH]   (default: .)
# =============================================================================
set -eu

# `bonsai` comes from ~/.graphify/providers.json; it pins the local base_url,
# the model, and the thinking-off extra_body. Register it once with:
#   graphify provider add bonsai --base-url http://127.0.0.1:8000/v1 \
#     --default-model Ternary-Bonsai-2-27B-tight --env-key BONSAI_API_KEY
# then add the "extra_body" and "reasoning_effort" keys to that entry.
export BONSAI_API_KEY="${BONSAI_API_KEY:-local}"
# Override to use the built-in openai backend instead (env OPENAI_BASE_URL/OPENAI_MODEL).
export GRAPHIFY_LLM_BACKEND="${GRAPHIFY_LLM_BACKEND:-bonsai}"
export GRAPHIFY_MAX_OUTPUT_TOKENS="${GRAPHIFY_MAX_OUTPUT_TOKENS:-16000}"
export GRAPHIFY_API_TIMEOUT="${GRAPHIFY_API_TIMEOUT:-1200}"

exec graphify extract "${1:-.}" \
  --backend="$GRAPHIFY_LLM_BACKEND" \
  --max-concurrency 2 \
  --token-budget 12000 \
  --out .

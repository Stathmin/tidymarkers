#!/bin/sh
# Render the lecture deck with OligoArrayAux on PATH (as run.sh does), so the
# design chunks can score PCR efficiencies with DECIPHER.
set -e
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
export PATH="$ROOT/.tools/bin:$PATH"
export UNAFOLDDAT="$ROOT/.tools/share/oligoarrayaux"
cd "$ROOT"
exec quarto render lecture/kasp-slides.qmd "$@"

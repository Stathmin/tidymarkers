#!/bin/sh
# Render the lecture deck with OligoArrayAux on PATH (as run.sh does), so the
# design chunks can score PCR efficiencies with DECIPHER.
#
# Windows: run this from Git-Bash (or the Rtools/MSYS2 shell) - `sh` and the
# build tools come with Rtools; cairo for ggtools/ggmsa is part of Rtools too.
set -e
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
export PATH="$ROOT/.tools/bin:$PATH"
export UNAFOLDDAT="$ROOT/.tools/share/oligoarrayaux"
cd "$ROOT"
exec quarto render dev/lecture/kasp-slides.qmd "$@"

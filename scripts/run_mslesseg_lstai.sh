#!/usr/bin/env bash
# Same subjects, but mscard runs LST-AI itself (docker), so the segmentation is graded too.
# Usage: scripts/run_mslesseg_lstai.sh P2 P3 ...
set -uo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
MANIFEST="${MANIFEST:-$HOME/repos/lesiontrack/work/manifest_mslesseg.tsv}"
OUT="$ROOT/derivatives/mslesseg-lstai"
LOG="$ROOT/work/logs/lstai-$(date -u +%Y%m%d-%H%M%S).log"
mkdir -p "$OUT" "$ROOT/work/logs"
ARGS=()
for s in "$@"; do ARGS+=(--subject "$s"); done
nohup "$HOME/.venvs/mscard/bin/mscard" run "$MANIFEST" --out "$OUT" --threads 4 --segmenter lst-ai "${ARGS[@]}" > "$LOG" 2>&1 &
echo "pid $! log $LOG"

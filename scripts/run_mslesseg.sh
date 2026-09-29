#!/usr/bin/env bash
# Run mscard over the longitudinal MSLesSeg subjects (expert masks, --segmenter given), detached.
# Usage: scripts/run_mslesseg.sh [subject ...]     (default: every subject in the manifest)
set -uo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
MANIFEST="${MANIFEST:-$HOME/repos/lesiontrack/work/manifest_mslesseg.tsv}"
OUT="$ROOT/derivatives/mslesseg"
LOG="$ROOT/work/logs/cohort-$(date -u +%Y%m%d-%H%M%S).log"
mkdir -p "$OUT" "$ROOT/work/logs"
ARGS=()
for s in "$@"; do ARGS+=(--subject "$s"); done
nohup "$HOME/.venvs/mscard/bin/mscard" run "$MANIFEST" --out "$OUT" --threads 4 "${ARGS[@]}" > "$LOG" 2>&1 &
echo "pid $! log $LOG"

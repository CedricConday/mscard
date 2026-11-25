#!/usr/bin/env bash
# Both MSLesSeg cohorts on lesiontrack 0.3.2 (lesion masks resliced without the similarity scale).
# Registrations are cached; only the masks are resliced again.
# Usage: scripts/run_cohorts_lt032.sh     log: work/logs/lt032.log, ends "lt032 done"
set -uo pipefail
cd "$(dirname "$0")/.."
MS=~/.venvs/mscard/bin
MANIFEST=~/repos/lesiontrack/work/manifest_mslesseg.tsv
LOG=work/logs/lt032.log
say() { echo "$(date "+%F %T") $*" | tee -a "$LOG"; }
say "start $(git rev-parse --short HEAD), $($MS/pip list 2>/dev/null | grep -i "^lesiontrack " | tr -s " ")"
for prof in default smooth; do
  out=derivatives/mslesseg; [ $prof = smooth ] && out=derivatives/mslesseg-smooth
  start=$(date +%s)
  $MS/mscard run $MANIFEST --out $out --threads 4 --reg-profile $prof > work/logs/lt032_$prof.log 2>&1
  say "$prof exit=$? seconds=$(( $(date +%s) - start ))"
  $MS/python scripts/results_table.py $out > work/results_lt032_$prof.md 2>> "$LOG"
done
say "lt032 done"

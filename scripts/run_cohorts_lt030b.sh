#!/usr/bin/env bash
# Both MSLesSeg cohorts again on lesiontrack 0.3.0 after the brain-volume fix (similarity
# volume factor). Registrations are cached from lt030. Waits for the RepairGauge backtests.
# Usage: scripts/run_cohorts_lt030b.sh     log: work/logs/lt030b.log, ends "lt030b done"
set -uo pipefail
cd "$(dirname "$0")/.."
MS=~/.venvs/mscard/bin
MANIFEST=~/repos/lesiontrack/work/manifest_mslesseg.tsv
LOG=work/logs/lt030b.log
say() { echo "$(date "+%F %T") $*" | tee -a "$LOG"; }
while pgrep -f "[r]un_reported.sh" > /dev/null; do sleep 60; done
say "start $(git rev-parse --short HEAD), $($MS/pip list 2>/dev/null | grep -i "^lesiontrack " | tr -s " ")"
for prof in default smooth; do
  out=derivatives/mslesseg; [ $prof = smooth ] && out=derivatives/mslesseg-smooth
  start=$(date +%s)
  $MS/mscard run $MANIFEST --out $out --threads 4 --reg-profile $prof > work/logs/lt030b_$prof.log 2>&1
  say "$prof exit=$? seconds=$(( $(date +%s) - start ))"
  $MS/python scripts/results_table.py $out > work/results_lt030b_$prof.md 2>> "$LOG"
done
say "lt030b done"

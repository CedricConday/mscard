#!/usr/bin/env bash
# Both MSLesSeg cohorts (default and --reg-profile smooth) again on lesiontrack 0.3.0, whose
# defaults differ from the 0.2.0 the README numbers were measured with. Waits for
# lesiontrack's ANTs chain so the two do not share the cores.
# Usage: scripts/run_cohorts_lt030.sh     log: work/logs/lt030.log, ends "lt030 done"
set -uo pipefail
cd "$(dirname "$0")/.."
MS=~/.venvs/mscard/bin
MANIFEST=~/repos/lesiontrack/work/manifest_mslesseg.tsv
mkdir -p work/logs
LOG=work/logs/lt030.log
say() { echo "$(date '+%F %T') $*" | tee -a "$LOG"; }
while pgrep -f "[r]un_ants030.sh" > /dev/null; do sleep 60; done
say "start, $($MS/pip list 2>/dev/null | grep -i '^lesiontrack ' | tr -s ' ')"
start=$(date +%s)
$MS/mscard run $MANIFEST --out derivatives/mslesseg --threads 4 > work/logs/lt030_default.log 2>&1
say "default exit=$? seconds=$(( $(date +%s) - start ))"
start=$(date +%s)
$MS/mscard run $MANIFEST --out derivatives/mslesseg-smooth --threads 4 --reg-profile smooth > work/logs/lt030_smooth.log 2>&1
say "smooth exit=$? seconds=$(( $(date +%s) - start ))"
$MS/python scripts/results_table.py derivatives/mslesseg > work/results_lt030_default.md 2>> "$LOG"
$MS/python scripts/results_table.py derivatives/mslesseg-smooth > work/results_lt030_smooth.md 2>> "$LOG"
say "lt030 done"

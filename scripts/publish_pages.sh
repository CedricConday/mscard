#!/usr/bin/env bash
# Copy the gallery and every self-contained report into docs/ for GitHub Pages.
# Usage: scripts/publish_pages.sh [derivatives-dir]   (default derivatives/mslesseg)
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
SRC="${1:-$ROOT/derivatives/mslesseg}"
DOCS="$ROOT/docs"
mkdir -p "$DOCS"
cp "$SRC/index.html" "$DOCS/index.html"
n=0
for r in "$SRC"/*/report.html; do
  s=$(basename "$(dirname "$r")")
  mkdir -p "$DOCS/$s"
  cp "$r" "$DOCS/$s/report.html"
  n=$((n+1))
done
touch "$DOCS/.nojekyll"
echo "docs/: gallery + $n reports"

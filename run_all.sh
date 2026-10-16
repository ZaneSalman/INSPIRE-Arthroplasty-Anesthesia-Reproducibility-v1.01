#!/usr/bin/env bash
# Runs the v1.02 analysis in order. Requires INSPIRE_DATA_DIR and INSPIRE_WORK_DIR (see README).
set -euo pipefail
: "${INSPIRE_DATA_DIR:?set INSPIRE_DATA_DIR to the folder holding the INSPIRE v1.4.2 csv.gz files}"
: "${INSPIRE_WORK_DIR:?set INSPIRE_WORK_DIR to a private local folder for intermediates}"
cd "$(dirname "$0")/code"
for s in 0*.py 1*.py; do
  echo "=== $s"
  python3 "$s"
done

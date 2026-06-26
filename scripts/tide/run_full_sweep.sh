#!/bin/bash
# run_full_sweep.sh — chain sweep_process → collect → sweep_cv → collect
#
# Usage (from the repo root):
#   bash scripts/tide/run_full_sweep.sh
#
# Sequence:
#   1. sweep_process: upload + submit + wait + fetch
#   2. sweep_process collect.py → _data/processed/sweep/sweep_all.csv
#   3. sweep_cv:     upload + submit + wait + fetch
#   4. sweep_cv collect.R → sweep_results / sweep_summary / sweep_optima
#
# Each sub-sweep is independently re-runnable via its own submit.sh.
# This orchestrator is just a sequencer for the common end-to-end path.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

echo "=== run_full_sweep: phase 1/4 — radius sweep on Tide ==="
bash "${SCRIPT_DIR}/sweep_process/submit.sh" all

echo ""
echo "=== run_full_sweep: phase 2/4 — assembling sweep_all.csv locally ==="
python "${SCRIPT_DIR}/sweep_process/collect.py"

echo ""
echo "=== run_full_sweep: phase 3/4 — CV sweep on Tide ==="
bash "${SCRIPT_DIR}/sweep_cv/submit.sh" all

echo ""
echo "=== run_full_sweep: phase 4/4 — assembling CV summary + optima locally ==="
Rscript "${SCRIPT_DIR}/sweep_cv/collect.R"

echo ""
echo "=== Done ==="
echo "Optima written to _data/processed/sweep_cv/sweep_optima.csv"

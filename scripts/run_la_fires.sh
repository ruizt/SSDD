#!/usr/bin/env bash
# Stage 1 SSDD computation for the LA fires (Palisades + Eaton).
#
# Reproduces the runs that produced the files currently in
# _data/processed/{palisades,eaton}/.
#
# Inputs: LARIAC6 building footprints + DINS public-view points already staged
#         under _data/raw/{buildings,dins}/.
# Outputs: per-fire {raw_metrics.csv, buildings.gpkg, compute_log.txt}.
#
# Usage:
#   conda activate ssdd
#   ./scripts/run_la_fires.sh                # both fires, default params
#   ./scripts/run_la_fires.sh palisades      # one fire
#   ./scripts/run_la_fires.sh eaton --dins-only --r-d 150   # extra flags pass through
#
# All paths are resolved relative to the repo root, regardless of where the
# script is invoked from.

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
RAW_BLD="${REPO_ROOT}/_data/raw/buildings"
RAW_DINS="${REPO_ROOT}/_data/raw/dins"
OUT_ROOT="${REPO_ROOT}/_data/processed"
COMPUTE="${REPO_ROOT}/src/ssdd_compute.py"

run_fire() {
    local fire="$1"; shift
    local fire_cap
    fire_cap="$(tr '[:lower:]' '[:upper:]' <<< "${fire:0:1}")${fire:1}"

    local buildings="${RAW_BLD}/LARIAC6_Buildings_2020_${fire}.shp"
    local dins="${RAW_DINS}/DINS_2025_${fire_cap}_Public_View.geojson"
    local output="${OUT_ROOT}/${fire}"

    echo "=== Stage 1 :: ${fire} ==="
    echo "  buildings: ${buildings}"
    echo "  dins     : ${dins}"
    echo "  output   : ${output}"
    mkdir -p "${output}"

    python "${COMPUTE}" \
        --buildings "${buildings}" \
        --dins      "${dins}" \
        --output    "${output}" \
        --run-name  "${fire}" \
        "$@"
}

# Pick fires from arg list; bare flags (--*) pass through to ssdd_compute.py.
fires=()
extra_args=()
for arg in "$@"; do
    case "${arg}" in
        palisades|eaton) fires+=("${arg}") ;;
        *)               extra_args+=("${arg}") ;;
    esac
done
if [[ ${#fires[@]} -eq 0 ]]; then
    fires=(palisades eaton)
fi

for fire in "${fires[@]}"; do
    # ${arr[@]+"${arr[@]}"} expands to the array if set, to nothing if empty —
    # avoids the "unbound variable" error from set -u.
    run_fire "${fire}" ${extra_args[@]+"${extra_args[@]}"}
    echo
done

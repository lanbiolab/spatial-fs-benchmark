#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

RESULTS_ROOT="results/biological_subsets_rebuild_v1"
LOG_DIR="$RESULTS_ROOT/queue_logs"
DRIVER_LOG="$LOG_DIR/feature_precompute_driver.log"
SENTINEL="$RESULTS_ROOT/feature_precompute.done"
mkdir -p "$LOG_DIR"
rm -f "$SENTINEL"

run_one() {
    local family="$1"
    local subset="$2"
    local config="configs/rebuild_v1/biological_subsets/feature_only/${family}/${subset}.yaml"
    local log="$LOG_DIR/feature_${subset}_${family}.log"
    echo "[$(date --iso-8601=seconds)] START ${subset}/${family}" | tee -a "$log" "$DRIVER_LOG"
    OMP_NUM_THREADS=6 OPENBLAS_NUM_THREADS=6 MKL_NUM_THREADS=6 \
        R_FUTURE_PLAN=sequential \
        .conda-env/bin/spatial-fs-run --config "$config" >>"$log" 2>&1
    echo "[$(date --iso-8601=seconds)] DONE ${subset}/${family}" | tee -a "$log" "$DRIVER_LOG"
}

echo "[$(date --iso-8601=seconds)] FEATURE PRECOMPUTE START" >"$DRIVER_LOG"
for family in non_spatial svg; do
    run_one "$family" stomics_0212_immune_subset
    run_one "$family" stomics_0212_epithelial_subset
done
touch "$SENTINEL"
echo "[$(date --iso-8601=seconds)] FEATURE PRECOMPUTE DONE" | tee -a "$DRIVER_LOG"

#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

RESULTS_ROOT="results/biological_subsets_rebuild_v1"
LOG_DIR="$RESULTS_ROOT/queue_logs"
LOG="$LOG_DIR/somde_repair.log"
DONE="$RESULTS_ROOT/somde_repair.done"
FAILED="$RESULTS_ROOT/somde_repair.failed"
mkdir -p "$LOG_DIR"
rm -f "$DONE" "$FAILED"
trap 'touch "$FAILED"' ERR

echo "[$(date --iso-8601=seconds)] SOMDE REPAIR START" >"$LOG"
for subset in stomics_0212_immune_subset stomics_0212_epithelial_subset; do
    echo "[$(date --iso-8601=seconds)] START ${subset}" | tee -a "$LOG"
    OMP_NUM_THREADS=6 OPENBLAS_NUM_THREADS=6 MKL_NUM_THREADS=6 \
        R_FUTURE_PLAN=sequential \
        .conda-env/bin/spatial-fs-run \
        --config "configs/rebuild_v1/biological_subsets/feature_only/svg/${subset}.yaml" \
        >>"$LOG" 2>&1
    echo "[$(date --iso-8601=seconds)] DONE ${subset}" | tee -a "$LOG"
done

touch "$DONE"
trap - ERR
echo "[$(date --iso-8601=seconds)] SOMDE REPAIR DONE" | tee -a "$LOG"

tmux new-session -d -s bio_subset_immune \
    "cd '$ROOT' && scripts/run_biological_subset_queue.sh stomics_0212_immune_subset"
tmux new-session -d -s bio_subset_epithelial \
    "cd '$ROOT' && scripts/run_biological_subset_queue.sh stomics_0212_epithelial_subset"

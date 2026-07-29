#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

RESULTS_ROOT="results/biological_subsets_rebuild_v1"
LOG="$RESULTS_ROOT/queue_logs/task_metadata_refresh.log"
DONE="$RESULTS_ROOT/task_metadata_refreshed.done"
FAILED="$RESULTS_ROOT/task_metadata_refreshed.failed"
rm -f "$DONE" "$FAILED"
trap 'touch "$FAILED"' ERR
echo "[$(date --iso-8601=seconds)] TASK METADATA REFRESH START" >"$LOG"

for subset in stomics_0212_immune_subset stomics_0212_epithelial_subset; do
    for family in non_spatial svg; do
        echo "[$(date --iso-8601=seconds)] START ${subset}/${family}" | tee -a "$LOG"
        CUDA_VISIBLE_DEVICES="" \
            OMP_NUM_THREADS=6 OPENBLAS_NUM_THREADS=6 MKL_NUM_THREADS=6 \
            .conda-env/bin/spatial-fs-run \
            --config "configs/rebuild_v1/biological_subsets/${family}/${subset}.yaml" \
            >>"$LOG" 2>&1
        echo "[$(date --iso-8601=seconds)] DONE ${subset}/${family}" | tee -a "$LOG"
    done
done

.conda-env/bin/python scripts/rebuild_results_from_records.py \
    --results-root "$RESULTS_ROOT" \
    --include-glob 'stomics_0212_*_subset' \
    --merged-output "$RESULTS_ROOT/merged_results.csv" >>"$LOG" 2>&1
.conda-env/bin/python scripts/audit_svg_downstream_results.py \
    --results-root "$RESULTS_ROOT" \
    --config-root configs/rebuild_v1/biological_subsets \
    --output-dir "$RESULTS_ROOT/audit" >>"$LOG" 2>&1

touch "$DONE"
trap - ERR
echo "[$(date --iso-8601=seconds)] TASK METADATA REFRESH DONE" | tee -a "$LOG"

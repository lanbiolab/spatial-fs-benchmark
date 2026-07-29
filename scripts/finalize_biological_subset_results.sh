#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

RESULTS_ROOT="results/biological_subsets_rebuild_v1"
LOG_DIR="$RESULTS_ROOT/queue_logs"
LOG="$LOG_DIR/finalize.log"
SUBSETS=(stomics_0212_immune_subset stomics_0212_epithelial_subset)
mkdir -p "$LOG_DIR"

echo "[$(date --iso-8601=seconds)] FINALIZER START" >"$LOG"
while true; do
    ready=1
    for subset in "${SUBSETS[@]}"; do
        if [[ -f "$RESULTS_ROOT/${subset}.failed" ]]; then
            echo "[$(date --iso-8601=seconds)] ABORT failed subset=${subset}" | tee -a "$LOG"
            exit 1
        fi
        [[ -f "$RESULTS_ROOT/${subset}.done" ]] || ready=0
    done
    [[ "$ready" -eq 1 ]] && break
    echo "[$(date --iso-8601=seconds)] WAIT subset queues" >>"$LOG"
    sleep 120
done

.conda-env/bin/python scripts/rebuild_results_from_records.py \
    --results-root "$RESULTS_ROOT" \
    --include-glob 'stomics_0212_*_subset' \
    --merged-output "$RESULTS_ROOT/merged_results.csv" >>"$LOG" 2>&1

.conda-env/bin/python scripts/audit_svg_downstream_results.py \
    --results-root "$RESULTS_ROOT" \
    --config-root configs/rebuild_v1/biological_subsets \
    --output-dir "$RESULTS_ROOT/audit" >>"$LOG" 2>&1

touch "$RESULTS_ROOT/finalized.done"
echo "[$(date --iso-8601=seconds)] FINALIZER DONE" | tee -a "$LOG"

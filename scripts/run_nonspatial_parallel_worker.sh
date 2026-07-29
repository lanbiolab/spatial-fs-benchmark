#!/usr/bin/env bash
set -euo pipefail

if [[ "$#" -lt 3 ]]; then
    echo "Usage: $0 GPU_ID WORKER_ID CONFIG [CONFIG ...]" >&2
    exit 2
fi

GPU_ID="$1"
WORKER_ID="$2"
shift 2

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

LOG_DIR="results/spatial_svg_rebuild_v1/non_spatial_parallel_logs"
DRIVER_LOG="results/spatial_svg_rebuild_v1/non_spatial_parallel_${WORKER_ID}.log"
mkdir -p "$LOG_DIR"

echo "[$(date --iso-8601=seconds)] WORKER START id=${WORKER_ID} gpu=${GPU_ID}" >>"$DRIVER_LOG"
for config in "$@"; do
    phase="$(basename "$(dirname "$config")")"
    name="$(basename "$config" .yaml)"
    log="$LOG_DIR/${WORKER_ID}_${phase}_${name}.log"
    echo "[$(date --iso-8601=seconds)] START ${phase}/${name}" | tee "$log" | tee -a "$DRIVER_LOG"
    if CUDA_VISIBLE_DEVICES="$GPU_ID" \
        OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8 \
        R_FUTURE_PLAN=sequential \
        .conda-env/bin/spatial-fs-run --config "$config" >>"$log" 2>&1; then
        echo "[$(date --iso-8601=seconds)] DONE ${phase}/${name}" | tee -a "$log" | tee -a "$DRIVER_LOG"
    else
        status=$?
        echo "[$(date --iso-8601=seconds)] FAILED(${status}) ${phase}/${name}" | tee -a "$log" | tee -a "$DRIVER_LOG"
        exit "$status"
    fi
done
echo "[$(date --iso-8601=seconds)] WORKER DONE id=${WORKER_ID}" >>"$DRIVER_LOG"

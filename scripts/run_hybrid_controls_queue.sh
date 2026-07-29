#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

LOG_DIR="results/hybrid_controls_v1/logs"
DRIVER_LOG="results/hybrid_controls_v1/queue.log"
MAX_EXISTING_MEMORY_MIB="${HYBRID_MAX_EXISTING_MEMORY_MIB:-2500}"
mkdir -p "$LOG_DIR"

choose_idle_gpu() {
    nvidia-smi --query-gpu=index,memory.used,utilization.gpu --format=csv,noheader,nounits |
        awk -F',' -v limit="$MAX_EXISTING_MEMORY_MIB" '
            {
                gsub(/ /, "", $1); gsub(/ /, "", $2); gsub(/ /, "", $3)
                if ($2 < limit && $3 < 15) { print $1; exit }
            }
        '
}

wait_for_gpu() {
    local gpu=""
    while [[ -z "$gpu" ]]; do
        gpu="$(choose_idle_gpu || true)"
        if [[ -z "$gpu" ]]; then
            echo "[$(date --iso-8601=seconds)] WAIT no GPU below ${MAX_EXISTING_MEMORY_MIB} MiB" >>"$DRIVER_LOG"
            sleep 120
        fi
    done
    printf '%s' "$gpu"
}

for dataset in mouse_brain_serial_sections dlpfc stomics_0212 stomics_0218 stomics_0224 e8p5_embryo e9p5_embryo; do
    config="configs/rebuild_v1/hybrid_controls/${dataset}.yaml"
    log="$LOG_DIR/${dataset}.log"
    gpu="$(wait_for_gpu)"
    echo "[$(date --iso-8601=seconds)] START ${dataset} gpu=${gpu}" | tee -a "$DRIVER_LOG" "$log"
    CUDA_VISIBLE_DEVICES="$gpu" \
        OMP_NUM_THREADS=6 OPENBLAS_NUM_THREADS=6 MKL_NUM_THREADS=6 \
        .conda-env/bin/spatial-fs-run --config "$config" >>"$log" 2>&1
    echo "[$(date --iso-8601=seconds)] DONE ${dataset} gpu=${gpu}" | tee -a "$DRIVER_LOG" "$log"
done

echo "[$(date --iso-8601=seconds)] QUEUE DONE" | tee -a "$DRIVER_LOG"

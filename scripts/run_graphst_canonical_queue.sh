#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

LOG_DIR="results/graphst_canonical_v1/logs"
mkdir -p "$LOG_DIR"

run_queue() {
    local gpu="$1"
    shift
    for dataset in "$@"; do
        local config="configs/rebuild_v1/downstream/graphst_canonical/${dataset}.yaml"
        local log="$LOG_DIR/${dataset}.log"
        echo "[$(date --iso-8601=seconds)] START dataset=${dataset} gpu=${gpu}" | tee -a "$log"
        CUDA_VISIBLE_DEVICES="$gpu" \
            OMP_NUM_THREADS=6 OPENBLAS_NUM_THREADS=6 MKL_NUM_THREADS=6 \
            .conda-env/bin/spatial-fs-run --config "$config" >>"$log" 2>&1
        echo "[$(date --iso-8601=seconds)] DONE dataset=${dataset} gpu=${gpu}" | tee -a "$log"
    done
}

run_queue 0 mouse_brain_serial_sections stomics_0218 &
pid0=$!
run_queue 1 stomics_0212 e8p5_embryo dlpfc &
pid1=$!
run_queue 2 stomics_0224 e9p5_embryo &
pid2=$!

wait "$pid0" "$pid1" "$pid2"
echo "[$(date --iso-8601=seconds)] ALL GRAPHST QUEUES COMPLETE" | tee -a "$LOG_DIR/queue.log"

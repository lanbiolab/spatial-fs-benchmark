#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

LOG_DIR="results/spatial_svg_rebuild_v1/downstream_queue_logs"
mkdir -p "$LOG_DIR"
SVG_DOWNSTREAM_GPU="${SVG_DOWNSTREAM_GPU:-2}"

run_one() {
    local config="$1"
    local phase="$2"
    local name
    name="$(basename "$config" .yaml)"
    local log="$LOG_DIR/${phase}_${name}.log"
    echo "[$(date --iso-8601=seconds)] START ${phase}/${name}" | tee "$log"
    if CUDA_VISIBLE_DEVICES="$SVG_DOWNSTREAM_GPU" \
        OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8 \
        .conda-env/bin/spatial-fs-run --config "$config" >>"$log" 2>&1; then
        echo "[$(date --iso-8601=seconds)] DONE ${phase}/${name}" | tee -a "$log"
    else
        status=$?
        echo "[$(date --iso-8601=seconds)] FAILED(${status}) ${phase}/${name}" | tee -a "$log"
        return "$status"
    fi
}

DATASET_ORDER=(
    mouse_brain_serial_sections
    dlpfc
    stomics_0212
    stomics_0218
    stomics_0224
    e8p5_embryo
    e9p5_embryo
)

for phase in canonical feature_number; do
    for dataset in "${DATASET_ORDER[@]}"; do
        run_one "configs/rebuild_v1/downstream/${phase}/${dataset}.yaml" "$phase"
    done
done

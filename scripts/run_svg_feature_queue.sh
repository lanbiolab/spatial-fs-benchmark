#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

mkdir -p results/spatial_svg_rebuild_v1/queue_logs

SVG_QUEUE_JOBS="${SVG_QUEUE_JOBS:-3}"
SVG_METHOD_THREADS="${SVG_METHOD_THREADS:-4}"

run_one() {
    local config="$1"
    local name
    name="$(basename "$config" .yaml)"
    local log="results/spatial_svg_rebuild_v1/queue_logs/${name}.log"
    echo "[$(date --iso-8601=seconds)] START ${name}" | tee "$log"
    if OMP_NUM_THREADS="$SVG_METHOD_THREADS" \
        OPENBLAS_NUM_THREADS="$SVG_METHOD_THREADS" \
        MKL_NUM_THREADS="$SVG_METHOD_THREADS" \
        .conda-env/bin/spatial-fs-run --config "$config" >>"$log" 2>&1; then
        echo "[$(date --iso-8601=seconds)] DONE ${name}" | tee -a "$log"
    else
        status=$?
        echo "[$(date --iso-8601=seconds)] FAILED(${status}) ${name}" | tee -a "$log"
        return "$status"
    fi
}

export -f run_one
export ROOT
export SVG_METHOD_THREADS

find configs/rebuild_v1 -maxdepth 1 -name '*_svg_features.yaml' -print0 \
    | sort -z \
    | xargs -0 -n1 -P"$SVG_QUEUE_JOBS" bash -c 'run_one "$1"' _

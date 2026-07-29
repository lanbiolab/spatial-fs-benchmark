#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

if [[ $# -ne 1 ]]; then
    echo "Usage: $0 {stomics_0212_immune_subset|stomics_0212_epithelial_subset}" >&2
    exit 2
fi

SUBSET="$1"
case "$SUBSET" in
    stomics_0212_immune_subset|stomics_0212_epithelial_subset) ;;
    *) echo "Unknown subset: $SUBSET" >&2; exit 2 ;;
esac

RESULTS_ROOT="results/biological_subsets_rebuild_v1"
LOG_DIR="$RESULTS_ROOT/queue_logs"
DRIVER_LOG="$LOG_DIR/${SUBSET}_driver.log"
FEATURE_SENTINEL="$RESULTS_ROOT/feature_precompute.done"
DONE_SENTINEL="$RESULTS_ROOT/${SUBSET}.done"
FAILED_SENTINEL="$RESULTS_ROOT/${SUBSET}.failed"
mkdir -p "$LOG_DIR"
rm -f "$DONE_SENTINEL" "$FAILED_SENTINEL"
trap 'touch "$FAILED_SENTINEL"' ERR

gpu_has_compute_process() {
    local gpu="$1"
    [[ -n "$(nvidia-smi -i "$gpu" --query-compute-apps=pid --format=csv,noheader,nounits 2>/dev/null | tr -d '[:space:]')" ]]
}

acquire_free_gpu() {
    while true; do
        for gpu in 0 1 2 3; do
            exec {lock_fd}>"/tmp/spatial_fs_benchmark_gpu_${gpu}.lock"
            if flock -n "$lock_fd"; then
                if ! gpu_has_compute_process "$gpu"; then
                    GPU="$gpu"
                    GPU_LOCK_FD="$lock_fd"
                    export GPU GPU_LOCK_FD
                    return 0
                fi
                flock -u "$lock_fd"
            fi
            eval "exec ${lock_fd}>&-"
        done
        echo "[$(date --iso-8601=seconds)] WAIT no unoccupied GPU" | tee -a "$DRIVER_LOG"
        sleep 120
    done
}

acquire_shared_gpu() {
    local gpu="${BIOLOGICAL_SUBSET_GPU:?BIOLOGICAL_SUBSET_GPU is required}"
    local max_used_mib="${BIOLOGICAL_SUBSET_MAX_EXISTING_MIB:-12000}"
    if [[ ! "$gpu" =~ ^[0-3]$ ]]; then
        echo "Invalid BIOLOGICAL_SUBSET_GPU: $gpu" >&2
        return 2
    fi
    exec {lock_fd}>"/tmp/spatial_fs_benchmark_gpu_${gpu}.lock"
    flock "$lock_fd"
    while true; do
        used_mib="$(nvidia-smi -i "$gpu" --query-gpu=memory.used --format=csv,noheader,nounits)"
        if (( used_mib <= max_used_mib )); then
            GPU="$gpu"
            GPU_LOCK_FD="$lock_fd"
            export GPU GPU_LOCK_FD
            return 0
        fi
        echo "[$(date --iso-8601=seconds)] WAIT gpu=${gpu} memory=${used_mib}MiB" | tee -a "$DRIVER_LOG"
        sleep 120
    done
}

run_one() {
    local family="$1"
    local config="configs/rebuild_v1/biological_subsets/${family}/${SUBSET}.yaml"
    local log="$LOG_DIR/${SUBSET}_${family}.log"
    echo "[$(date --iso-8601=seconds)] START ${family} gpu=${GPU}" | tee -a "$log" "$DRIVER_LOG"
    CUDA_VISIBLE_DEVICES="$GPU" \
        OMP_NUM_THREADS=6 OPENBLAS_NUM_THREADS=6 MKL_NUM_THREADS=6 \
        R_FUTURE_PLAN=sequential \
        .conda-env/bin/spatial-fs-run --config "$config" >>"$log" 2>&1
    echo "[$(date --iso-8601=seconds)] DONE ${family} gpu=${GPU}" | tee -a "$log" "$DRIVER_LOG"
}

echo "[$(date --iso-8601=seconds)] QUEUE START subset=${SUBSET}" >"$DRIVER_LOG"
while [[ ! -f "$FEATURE_SENTINEL" ]]; do
    echo "[$(date --iso-8601=seconds)] WAIT feature precompute" | tee -a "$DRIVER_LOG"
    sleep 60
done
if [[ -n "${BIOLOGICAL_SUBSET_GPU:-}" ]]; then
    acquire_shared_gpu
else
    acquire_free_gpu
fi
echo "[$(date --iso-8601=seconds)] ACQUIRED gpu=${GPU}" | tee -a "$DRIVER_LOG"
run_one non_spatial
run_one svg
touch "$DONE_SENTINEL"
trap - ERR
echo "[$(date --iso-8601=seconds)] QUEUE DONE subset=${SUBSET}" | tee -a "$DRIVER_LOG"

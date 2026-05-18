#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT_DIR"

CONFIGS=(
  "configs/benchmark/stomics_0212_spatial_main_native_fullfs_repair.yaml"
  "configs/benchmark/stomics_0218_spatial_main_native_fullfs_repair.yaml"
  "configs/benchmark/stomics_0224_spatial_main_native_fullfs_repair.yaml"
)

gpu_count="${STOMICS_FULLFS_GPU_COUNT:-3}"
max_jobs="${STOMICS_FULLFS_MAX_JOBS:-3}"
launched=0

for config in "${CONFIGS[@]}"; do
  if (( launched >= max_jobs )); then
    break
  fi
  session="bench_$(basename "$config" .yaml)"
  gpu=$(( launched % gpu_count ))
  if tmux has-session -t "$session" 2>/dev/null; then
    echo "tmux session already exists: $session"
    launched=$(( launched + 1 ))
    continue
  fi
  tmux new-session -d -s "$session" "cd '$ROOT_DIR' && CUDA_VISIBLE_DEVICES=$gpu ./.conda-env/bin/python scripts/run_benchmark.py --config '$config'"
  echo "started $session on GPU $gpu"
  launched=$(( launched + 1 ))
done

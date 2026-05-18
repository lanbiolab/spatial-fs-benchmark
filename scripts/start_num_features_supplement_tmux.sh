#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT_DIR"

./.conda-env/bin/python scripts/generate_num_features_supplement_configs.py

CONFIGS=(
  "configs/benchmark/dlpfc_spatial_main_native_num_features_extra.yaml"
  "configs/benchmark/mouse_brain_serial_sections_spatial_main_native_num_features_extra.yaml"
  "configs/benchmark/stomics_visium_5samples_spatial_main_native_num_features_extra.yaml"
  "configs/benchmark/e8p5_embryo_spatial_main_native_num_features_extra.yaml"
  "configs/benchmark/e9p5_embryo_spatial_main_native_num_features_extra.yaml"
  "configs/benchmark/stomics_0218_wilcoxon_spatial_main_native_num_features_extra.yaml"
  "configs/benchmark/stomics_0224_wilcoxon_spatial_main_native_num_features_extra.yaml"
  "configs/benchmark/stomics_0212_wilcoxon_spatial_main_native_num_features_extra.yaml"
  "configs/benchmark/dlpfc_spatial_main_native_all_features.yaml"
  "configs/benchmark/mouse_brain_serial_sections_spatial_main_native_all_features.yaml"
  "configs/benchmark/stomics_visium_5samples_spatial_main_native_all_features.yaml"
  "configs/benchmark/e8p5_embryo_spatial_main_native_all_features.yaml"
  "configs/benchmark/e9p5_embryo_spatial_main_native_all_features.yaml"
  "configs/benchmark/stomics_0218_wilcoxon_spatial_main_native_all_features.yaml"
  "configs/benchmark/stomics_0224_wilcoxon_spatial_main_native_all_features.yaml"
  "configs/benchmark/stomics_0212_wilcoxon_spatial_main_native_all_features.yaml"
)

gpu_count="${NUM_FEATURES_SUPPLEMENT_GPU_COUNT:-4}"
max_jobs="${NUM_FEATURES_SUPPLEMENT_MAX_JOBS:-$gpu_count}"
start_offset="${NUM_FEATURES_SUPPLEMENT_START_OFFSET:-0}"
index=0
launched=0

for config in "${CONFIGS[@]}"; do
  if (( index < start_offset )); then
    index=$(( index + 1 ))
    continue
  fi
  if (( launched >= max_jobs )); then
    break
  fi
  session="bench_$(basename "$config" .yaml)"
  gpu=$(( launched % gpu_count ))
  if tmux has-session -t "$session" 2>/dev/null; then
    echo "tmux session already exists: $session"
    index=$(( index + 1 ))
    continue
  fi
  tmux new-session -d -s "$session" "cd '$ROOT_DIR' && CUDA_VISIBLE_DEVICES=$gpu ./.conda-env/bin/python scripts/run_benchmark.py --config '$config'"
  echo "started $session on GPU $gpu"
  index=$(( index + 1 ))
  launched=$(( launched + 1 ))
done

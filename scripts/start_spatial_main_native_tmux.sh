#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT_DIR"

CONFIGS=(
  "configs/benchmark/dlpfc_spatial_main_native.yaml"
  "configs/benchmark/mouse_brain_serial_sections_spatial_main_native.yaml"
  "configs/benchmark/stomics_visium_5samples_spatial_main_native.yaml"
  "configs/benchmark/e8p5_embryo_spatial_main_native.yaml"
  "configs/benchmark/e9p5_embryo_spatial_main_native.yaml"
)

for config in "${CONFIGS[@]}"; do
  session="bench_$(basename "$config" .yaml)"
  if tmux has-session -t "$session" 2>/dev/null; then
    echo "tmux session already exists: $session"
    continue
  fi
  tmux new-session -d -s "$session" "cd '$ROOT_DIR' && ./.conda-env/bin/python scripts/run_benchmark.py --config '$config'"
  echo "started $session"
done

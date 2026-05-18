#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT_DIR"

export OMP_NUM_THREADS=8
export OPENBLAS_NUM_THREADS=8
export MKL_NUM_THREADS=8
export NUMEXPR_NUM_THREADS=8
export VECLIB_MAXIMUM_THREADS=8
export BLIS_NUM_THREADS=8

declare -A CONFIGS=(
  [bench_dlpfc_atlas_fs]="configs/benchmark/dlpfc_atlas_fs.yaml"
  [bench_e8p5_embryo_atlas_fs]="configs/benchmark/e8p5_embryo_atlas_fs.yaml"
  [bench_e9p5_embryo_atlas_fs]="configs/benchmark/e9p5_embryo_atlas_fs.yaml"
  [bench_mouse_atlas_fs]="configs/benchmark/mouse_brain_serial_sections_atlas_fs.yaml"
  [bench_sagittal_atlas_fs]="configs/benchmark/sagittal_region_atlas_fs.yaml"
  [bench_stomics_visium_5samples_atlas_fs]="configs/benchmark/stomics_visium_5samples_atlas_fs.yaml"
)

for session in "${!CONFIGS[@]}"; do
  if tmux has-session -t "$session" 2>/dev/null; then
    echo "tmux session already exists: $session"
    continue
  fi
  config="${CONFIGS[$session]}"
  tmux new-session -d -s "$session" "cd '$ROOT_DIR' && ./.conda-env/bin/python scripts/run_benchmark.py --config '$config'"
  echo "started $session -> $config"
done

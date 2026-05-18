#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

mkdir -p results/lineage_subsets/logs

./.conda-env/bin/python scripts/create_stomics0212_lineage_subsets.py

tmux kill-session -t bench_stomics_0212_epithelial_subset_scvi 2>/dev/null || true
tmux kill-session -t bench_stomics_0212_immune_subset_scvi 2>/dev/null || true

tmux new-session -d -s bench_stomics_0212_epithelial_subset_scvi \
  "cd '$ROOT' && CUDA_VISIBLE_DEVICES=0 ./.conda-env/bin/python scripts/run_benchmark.py --config configs/benchmark/stomics_0212_epithelial_subset_scvi.yaml 2>&1 | tee results/lineage_subsets/logs/stomics_0212_epithelial_subset_scvi.log"

tmux new-session -d -s bench_stomics_0212_immune_subset_scvi \
  "cd '$ROOT' && CUDA_VISIBLE_DEVICES=1 ./.conda-env/bin/python scripts/run_benchmark.py --config configs/benchmark/stomics_0212_immune_subset_scvi.yaml 2>&1 | tee results/lineage_subsets/logs/stomics_0212_immune_subset_scvi.log"

echo "Started bench_stomics_0212_epithelial_subset_scvi on GPU 0"
echo "Started bench_stomics_0212_immune_subset_scvi on GPU 1"

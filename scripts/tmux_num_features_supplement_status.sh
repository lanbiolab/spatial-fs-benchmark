#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT_DIR"

echo "tmux sessions:"
tmux ls 2>/dev/null | grep 'bench_.*num_features_extra\|bench_.*all_features' || true

python - <<'PY'
from pathlib import Path
import yaml

targets = [
    Path("configs/benchmark/dlpfc_spatial_main_native_num_features_extra.yaml"),
    Path("configs/benchmark/mouse_brain_serial_sections_spatial_main_native_num_features_extra.yaml"),
    Path("configs/benchmark/stomics_visium_5samples_spatial_main_native_num_features_extra.yaml"),
    Path("configs/benchmark/e8p5_embryo_spatial_main_native_num_features_extra.yaml"),
    Path("configs/benchmark/e9p5_embryo_spatial_main_native_num_features_extra.yaml"),
    Path("configs/benchmark/stomics_0218_wilcoxon_spatial_main_native_num_features_extra.yaml"),
    Path("configs/benchmark/stomics_0224_wilcoxon_spatial_main_native_num_features_extra.yaml"),
    Path("configs/benchmark/stomics_0212_wilcoxon_spatial_main_native_num_features_extra.yaml"),
    Path("configs/benchmark/dlpfc_spatial_main_native_all_features.yaml"),
    Path("configs/benchmark/mouse_brain_serial_sections_spatial_main_native_all_features.yaml"),
    Path("configs/benchmark/stomics_visium_5samples_spatial_main_native_all_features.yaml"),
    Path("configs/benchmark/e8p5_embryo_spatial_main_native_all_features.yaml"),
    Path("configs/benchmark/e9p5_embryo_spatial_main_native_all_features.yaml"),
    Path("configs/benchmark/stomics_0218_wilcoxon_spatial_main_native_all_features.yaml"),
    Path("configs/benchmark/stomics_0224_wilcoxon_spatial_main_native_all_features.yaml"),
    Path("configs/benchmark/stomics_0212_wilcoxon_spatial_main_native_all_features.yaml"),
]
for config_path in targets:
    if not config_path.exists():
        continue
    cfg = yaml.safe_load(config_path.read_text())
    if not cfg:
        continue
    root = Path(cfg["output_dir"])
    total = (
        len(cfg["feature_selection_methods"])
        * len(cfg["integration_methods"])
        * len(cfg["tasks"])
        * len(cfg["n_features"])
        * len(cfg["seeds"])
    )
    done = sum(1 for _ in root.rglob("*_records.json")) if root.exists() else 0
    pct = 100.0 * done / total if total else 0.0
    print(f"{config_path.stem}: {done}/{total} ({pct:.1f}%)")
PY

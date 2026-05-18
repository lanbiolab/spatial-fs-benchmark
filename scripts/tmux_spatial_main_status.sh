#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT_DIR"

echo "tmux sessions:"
tmux ls 2>/dev/null | grep 'bench_.*_spatial_main_rerun' || true

python - <<'PY'
from pathlib import Path
import yaml

targets = [
    ("DLPFC", Path("results/spatial_main_rerun/dlpfc_spatial_main"), Path("configs/benchmark/dlpfc_spatial_main_rerun.yaml")),
    ("Mouse", Path("results/spatial_main_rerun/mouse_brain_serial_sections_spatial_main"), Path("configs/benchmark/mouse_brain_serial_sections_spatial_main_rerun.yaml")),
    ("Sagittal", Path("results/spatial_main_rerun/sagittal_region_spatial_main"), Path("configs/benchmark/sagittal_region_spatial_main_rerun.yaml")),
    ("STOmicsVisium5", Path("results/spatial_main_rerun/stomics_visium_5samples_spatial_main"), Path("configs/benchmark/stomics_visium_5samples_spatial_main_rerun.yaml")),
    ("E8p5Embryo", Path("results/spatial_main_rerun/e8p5_embryo_region_spatial_main"), Path("configs/benchmark/e8p5_embryo_spatial_main_rerun.yaml")),
    ("E9p5Embryo", Path("results/spatial_main_rerun/e9p5_embryo_region_spatial_main"), Path("configs/benchmark/e9p5_embryo_spatial_main_rerun.yaml")),
]
for name, root, config_path in targets:
    cfg = yaml.safe_load(config_path.read_text())
    total = (
        len(cfg["feature_selection_methods"])
        * len(cfg["integration_methods"])
        * len(cfg["tasks"])
        * len(cfg["n_features"])
        * len(cfg["seeds"])
    )
    done = sum(1 for _ in root.rglob("*_records.json")) if root.exists() else 0
    pct = 100.0 * done / total if total else 0.0
    print(f"{name}: {done}/{total} ({pct:.1f}%)")
PY

#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT_DIR"

echo "tmux sessions:"
tmux ls 2>/dev/null | grep 'bench_.*_atlas_fs' || true

python - <<'PY'
from pathlib import Path
import yaml

targets = [
    ("DLPFC", Path("results/dlpfc_atlas_fs"), Path("configs/benchmark/dlpfc_atlas_fs.yaml")),
    ("E8p5Embryo", Path("results/e8p5_embryo_region_atlas_fs"), Path("configs/benchmark/e8p5_embryo_atlas_fs.yaml")),
    ("E9p5Embryo", Path("results/e9p5_embryo_region_atlas_fs"), Path("configs/benchmark/e9p5_embryo_atlas_fs.yaml")),
    ("Mouse", Path("results/mouse_brain_serial_sections_atlas_fs"), Path("configs/benchmark/mouse_brain_serial_sections_atlas_fs.yaml")),
    ("Sagittal", Path("results/sagittal_region_atlas_fs"), Path("configs/benchmark/sagittal_region_atlas_fs.yaml")),
    ("STOmicsVisium5", Path("results/stomics_visium_5samples_atlas_fs"), Path("configs/benchmark/stomics_visium_5samples_atlas_fs.yaml")),
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
    done = sum(1 for _ in root.rglob("*_records.json"))
    pct = 100.0 * done / total if total else 0.0
    print(f"{name}: {done}/{total} ({pct:.1f}%)")
PY

from __future__ import annotations

from pathlib import Path

import yaml


ROOT = Path(__file__).resolve().parents[1]
OUTPUT_DIR = ROOT / "configs" / "rebuild_v1"

DATASETS = {
    "dlpfc": "configs/datasets/dlpfc.yaml",
    "mouse_brain_serial_sections": "configs/datasets/mouse_brain_serial_sections.yaml",
    "stomics_0212": "configs/datasets/stomics_0212_wilcoxon.yaml",
    "stomics_0218": "configs/datasets/stomics_0218_wilcoxon.yaml",
    "stomics_0224": "configs/datasets/stomics_0224_wilcoxon.yaml",
    "stomics_visium_5samples": "configs/datasets/stomics_visium_5samples.yaml",
    "e8p5_embryo": "configs/datasets/e8p5_embryo.yaml",
    "e9p5_embryo": "configs/datasets/e9p5_embryo.yaml",
}

ALL_N = [100, 200, 500, 1000, 2000, 5000, 10000]
SVG_METHOD_THREADS = 4


def spatial_methods(has_integer_counts: bool) -> list[dict]:
    methods = [{"name": "morans_i", "params": {"n_neighbors": 8}, "n_features": ALL_N}]
    if has_integer_counts:
        methods.extend(
            [
                {"name": "nnsvg", "params": {"n_threads": SVG_METHOD_THREADS}, "n_features": ALL_N},
                {
                    "name": "spatialde",
                    "params": {"max_cells_per_slice": 800},
                    "n_features": ALL_N,
                },
                {"name": "sparkx", "params": {"n_threads": SVG_METHOD_THREADS}, "n_features": ALL_N},
                {"name": "somde", "params": {"spots_per_node": 20}, "n_features": ALL_N},
            ]
        )
    return methods


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    for key, dataset_path in DATASETS.items():
        has_integer_counts = key != "stomics_visium_5samples"
        config = {
            "name": f"{key}_svg_feature_selection_rebuild_v1",
            "output_dir": f"results/spatial_svg_rebuild_v1/{key}",
            "datasets": [dataset_path],
            "feature_selection_methods": spatial_methods(has_integer_counts),
            "integration_methods": [],
            "tasks": [],
            "n_features": ALL_N,
            "seeds": [0, 1, 2],
            "save_embeddings": False,
        }
        destination = OUTPUT_DIR / f"{key}_svg_features.yaml"
        destination.write_text(yaml.safe_dump(config, sort_keys=False), encoding="utf-8")


if __name__ == "__main__":
    main()

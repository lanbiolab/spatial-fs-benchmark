from __future__ import annotations

from pathlib import Path

import yaml


ROOT = Path(__file__).resolve().parents[1]
CONFIG_DIR = ROOT / "configs" / "benchmark"

SOURCE_CONFIGS = [
    "dlpfc_atlas_fs.yaml",
    "mouse_brain_serial_sections_atlas_fs.yaml",
    "sagittal_region_atlas_fs.yaml",
    "stomics_visium_5samples_atlas_fs.yaml",
    "e8p5_embryo_atlas_fs.yaml",
    "e9p5_embryo_atlas_fs.yaml",
]

SPATIAL_MAIN_INTEGRATORS = [
    {"name": "scvi", "params": {"n_latent": 30, "max_epochs": 100, "batch_size": 2048}},
    {"name": "cellcharter", "params": {"n_latent": 30, "nhood_layers": 4, "spatial_neighbors": 6, "max_epochs": 100}},
    {"name": "gpsa", "params": {"n_input_dims": 20, "m_x_per_view": 50, "m_g": 50, "max_epochs": 200, "learning_rate": 0.01, "fixed_view_idx": 0}},
    {"name": "staligner", "params": {"hidden_dims": [128, 30], "n_epochs": 600, "lr": 0.001, "knn_neigh": 50, "radius_cutoff": 50.0}},
]


def main() -> None:
    for source_name in SOURCE_CONFIGS:
        source_path = CONFIG_DIR / source_name
        payload = yaml.safe_load(source_path.read_text(encoding="utf-8"))
        payload["integration_methods"] = SPATIAL_MAIN_INTEGRATORS
        payload["tasks"] = [
            task for task in payload["tasks"] if str(task.get("name")) != "slice_representation_eval"
        ]
        payload["name"] = str(payload["name"]).replace("_atlas_fs", "_spatial_main")
        payload["output_dir"] = str(payload["output_dir"]).replace("_atlas_fs", "_spatial_main")
        destination_name = source_name.replace("_atlas_fs", "_spatial_main")
        destination_path = CONFIG_DIR / destination_name
        destination_path.write_text(yaml.safe_dump(payload, sort_keys=False), encoding="utf-8")
        print(destination_path.relative_to(ROOT))


if __name__ == "__main__":
    main()

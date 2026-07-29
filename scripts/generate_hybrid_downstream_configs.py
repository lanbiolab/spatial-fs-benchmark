from __future__ import annotations

from pathlib import Path

import yaml


ROOT = Path(__file__).resolve().parents[1]
OUTPUT_DIR = ROOT / "configs/rebuild_v1/hybrid_controls"
RESULTS_ROOT = "results/hybrid_controls_v1"
DATASETS = {
    "dlpfc": ("configs/datasets/dlpfc.yaml", 7, True),
    "mouse_brain_serial_sections": ("configs/datasets/mouse_brain_serial_sections.yaml", 9, True),
    "stomics_0212": ("configs/datasets/stomics_0212_wilcoxon.yaml", 7, False),
    "stomics_0218": ("configs/datasets/stomics_0218_wilcoxon.yaml", 9, False),
    "stomics_0224": ("configs/datasets/stomics_0224_wilcoxon.yaml", 8, False),
    "e8p5_embryo": ("configs/datasets/e8p5_embryo.yaml", 12, False),
    "e9p5_embryo": ("configs/datasets/e9p5_embryo.yaml", 12, False),
}


def task_config(n_clusters: int, alignment: bool) -> list[dict]:
    tasks = [
        {"name": "integration_eval", "params": {"n_neighbors": 15, "max_eval_cells": 12000}},
        {
            "name": "clustering_eval",
            "params": {
                "n_clusters": n_clusters,
                "n_neighbors": 8,
                "silhouette_sample_size": 10000,
                "max_eval_cells": 12000,
            },
        },
    ]
    if alignment:
        tasks.append(
            {
                "name": "alignment_eval",
                "params": {"n_neighbors": 10, "max_points_per_slice": 2500},
            }
        )
    return tasks


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    for key, (dataset, n_clusters, alignment) in DATASETS.items():
        payload = {
            "name": f"{key}_hybrid_controls_v1",
            "output_dir": f"{RESULTS_ROOT}/{key}",
            "datasets": [dataset],
            "feature_selection_methods": [
                {"name": "hybrid_hvg_svg_union", "n_features": [2000]},
                {"name": "hybrid_hvg_svg_intersection", "n_features": [2000]},
            ],
            "integration_methods": [
                {
                    "name": "scvi",
                    "params": {"n_latent": 30, "max_epochs": 100, "batch_size": 2048},
                },
                {
                    "name": "cellcharter",
                    "params": {
                        "n_latent": 30,
                        "nhood_layers": 4,
                        "spatial_neighbors": 6,
                        "max_epochs": 100,
                        "batch_size": 2048,
                    },
                },
            ],
            "tasks": task_config(n_clusters, alignment),
            "n_features": [2000],
            "seeds": [0, 1, 2],
            "save_embeddings": True,
        }
        (OUTPUT_DIR / f"{key}.yaml").write_text(yaml.safe_dump(payload, sort_keys=False))


if __name__ == "__main__":
    main()

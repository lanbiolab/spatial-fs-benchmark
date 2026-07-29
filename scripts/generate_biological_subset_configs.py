from __future__ import annotations

from pathlib import Path

import yaml


ROOT = Path(__file__).resolve().parents[1]
CONFIG_ROOT = ROOT / "configs" / "rebuild_v1" / "biological_subsets"
RESULTS_ROOT = "results/biological_subsets_rebuild_v1"

DATASETS = {
    "stomics_0212_immune_subset": (
        "configs/datasets/stomics_0212_immune_subset.yaml",
        2,
    ),
    "stomics_0212_epithelial_subset": (
        "configs/datasets/stomics_0212_epithelial_subset.yaml",
        4,
    ),
}


def non_spatial_methods() -> list[dict]:
    methods = [
        {"name": "all_features", "n_features": [1]},
        {"name": "random", "n_features": [500]},
        {"name": "TFs", "n_features": [2000]},
        {"name": "scanpy_seurat"},
        {"name": "scanpy_seurat_batch"},
        {"name": "scanpy_seurat_v3"},
        {"name": "scanpy_seurat_v3_batch"},
        {"name": "scanpy_cell_ranger"},
        {"name": "scanpy_cell_ranger_batch"},
        {"name": "scanpy_pearson"},
        {"name": "scanpy_pearson_batch"},
        {"name": "seurat_vst"},
        {"name": "seurat_mvp"},
        {"name": "seurat_disp"},
        {"name": "seurat_sct", "params": {"max_cells": 3000}},
        {"name": "scsegindex", "n_features": [500]},
        {"name": "dubstepr", "params": {"max_cells": 5000}},
        {"name": "nbumi", "params": {"max_cells": 5000}},
        {"name": "osca"},
        {"name": "scry"},
        {"name": "singleCellHaystack", "params": {"max_cells": 5000}},
        {"name": "Brennecke"},
        {"name": "scPNMF", "params": {"max_cells": 3000}},
        {"name": "triku", "params": {"max_cells": 5000}},
        {"name": "hotspot", "params": {"max_cells": 3000}},
        {"name": "anticor", "params": {"max_cells": 3000}},
        {"name": "statistic_mean"},
        {"name": "statistic_variance"},
        {"name": "wilcoxon"},
    ]
    for method in methods:
        method.setdefault("n_features", [2000])
    return methods


def spatial_methods() -> list[dict]:
    return [
        {"name": "morans_i", "params": {"n_neighbors": 8}, "n_features": [2000]},
        {"name": "nnsvg", "params": {"n_threads": 4}, "n_features": [2000]},
        {
            "name": "spatialde",
            "params": {"max_cells_per_slice": 800},
            "n_features": [2000],
        },
        {"name": "sparkx", "params": {"n_threads": 4}, "n_features": [2000]},
        {"name": "somde", "params": {"spots_per_node": 20}, "n_features": [2000]},
    ]


def payload(key: str, dataset_path: str, n_clusters: int, family: str) -> dict:
    methods = spatial_methods() if family == "svg" else non_spatial_methods()
    return {
        "name": f"{key}_{family}_canonical_rebuild_v1",
        "output_dir": f"{RESULTS_ROOT}/{key}",
        "datasets": [dataset_path],
        "feature_selection_methods": methods,
        "integration_methods": [
            {
                "name": "scvi",
                "params": {"n_latent": 30, "max_epochs": 100, "batch_size": 2048},
            }
        ],
        "tasks": [
            {
                "name": "integration_eval",
                "params": {"n_neighbors": 15, "max_eval_cells": 12000},
            },
            {
                "name": "clustering_eval",
                "params": {
                    "n_clusters": n_clusters,
                    "n_neighbors": 8,
                    "silhouette_sample_size": 10000,
                    "max_eval_cells": 12000,
                },
            },
        ],
        "n_features": [2000],
        "seeds": [0, 1, 2],
        "save_embeddings": True,
    }


def main() -> None:
    for family in ("non_spatial", "svg"):
        family_dir = CONFIG_ROOT / family
        feature_dir = CONFIG_ROOT / "feature_only" / family
        family_dir.mkdir(parents=True, exist_ok=True)
        feature_dir.mkdir(parents=True, exist_ok=True)
        for key, (dataset_path, n_clusters) in DATASETS.items():
            config_payload = payload(key, dataset_path, n_clusters, family)
            destination = family_dir / f"{key}.yaml"
            destination.write_text(
                yaml.safe_dump(config_payload, sort_keys=False),
                encoding="utf-8",
            )
            feature_payload = dict(config_payload)
            feature_payload["name"] = f"{key}_{family}_feature_only_rebuild_v1"
            feature_payload["integration_methods"] = []
            feature_payload["tasks"] = []
            feature_payload["save_embeddings"] = False
            (feature_dir / f"{key}.yaml").write_text(
                yaml.safe_dump(feature_payload, sort_keys=False),
                encoding="utf-8",
            )


if __name__ == "__main__":
    main()

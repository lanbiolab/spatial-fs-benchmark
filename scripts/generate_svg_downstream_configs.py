from __future__ import annotations

from pathlib import Path

import yaml


ROOT = Path(__file__).resolve().parents[1]
OUTPUT_DIR = ROOT / "configs" / "rebuild_v1" / "downstream"
RESULTS_ROOT = "results/spatial_svg_rebuild_v1"

DATASETS = {
    "dlpfc": ("configs/datasets/dlpfc.yaml", 7, True),
    "mouse_brain_serial_sections": ("configs/datasets/mouse_brain_serial_sections.yaml", 9, True),
    "stomics_0212": ("configs/datasets/stomics_0212_wilcoxon.yaml", 7, False),
    "stomics_0218": ("configs/datasets/stomics_0218_wilcoxon.yaml", 9, False),
    "stomics_0224": ("configs/datasets/stomics_0224_wilcoxon.yaml", 8, False),
    "e8p5_embryo": ("configs/datasets/e8p5_embryo.yaml", 12, False),
    "e9p5_embryo": ("configs/datasets/e9p5_embryo.yaml", 12, False),
}

CANONICAL_N = [2000]
FEATURE_NUMBER_N = [100, 200, 500, 1000, 5000, 10000]


def spatial_methods(n_features: list[int]) -> list[dict]:
    return [
        {"name": "morans_i", "params": {"n_neighbors": 8}, "n_features": n_features},
        {"name": "nnsvg", "params": {"n_threads": 4}, "n_features": n_features},
        {
            "name": "spatialde",
            "params": {"max_cells_per_slice": 800},
            "n_features": n_features,
        },
        {"name": "sparkx", "params": {"n_threads": 4}, "n_features": n_features},
        {"name": "somde", "params": {"spots_per_node": 20}, "n_features": n_features},
    ]


def integration_methods() -> list[dict]:
    return [
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
    ]


def tasks(n_clusters: int, include_alignment: bool) -> list[dict]:
    output = [
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
    if include_alignment:
        output.append(
            {
                "name": "alignment_eval",
                "params": {"n_neighbors": 10, "max_points_per_slice": 2500},
            }
        )
    return output


def write_config(
    key: str,
    dataset_path: str,
    n_clusters: int,
    include_alignment: bool,
    phase: str,
    n_features: list[int],
    seeds: list[int],
) -> None:
    payload = {
        "name": f"{key}_svg_{phase}_downstream",
        "output_dir": f"{RESULTS_ROOT}/{key}",
        "datasets": [dataset_path],
        "feature_selection_methods": spatial_methods(n_features),
        "integration_methods": integration_methods(),
        "tasks": tasks(n_clusters, include_alignment),
        "n_features": n_features,
        "seeds": seeds,
        "save_embeddings": True,
    }
    phase_dir = OUTPUT_DIR / phase
    phase_dir.mkdir(parents=True, exist_ok=True)
    destination = phase_dir / f"{key}.yaml"
    destination.write_text(yaml.safe_dump(payload, sort_keys=False), encoding="utf-8")


def main() -> None:
    for key, (dataset_path, n_clusters, include_alignment) in DATASETS.items():
        write_config(
            key,
            dataset_path,
            n_clusters,
            include_alignment,
            "canonical",
            CANONICAL_N,
            [0, 1, 2],
        )
        write_config(
            key,
            dataset_path,
            n_clusters,
            include_alignment,
            "feature_number",
            FEATURE_NUMBER_N,
            [0],
        )


if __name__ == "__main__":
    main()

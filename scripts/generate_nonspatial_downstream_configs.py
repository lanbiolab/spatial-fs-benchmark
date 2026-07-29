from __future__ import annotations

from pathlib import Path

import yaml


ROOT = Path(__file__).resolve().parents[1]
OUTPUT_DIR = ROOT / "configs" / "rebuild_v1" / "non_spatial_downstream"
RESULTS_ROOT = "results/spatial_svg_rebuild_v1"

DATASETS = {
    "dlpfc": ("configs/datasets/dlpfc.yaml", 7, True, True),
    "mouse_brain_serial_sections": (
        "configs/datasets/mouse_brain_serial_sections.yaml",
        9,
        True,
        True,
    ),
    "stomics_0212": ("configs/datasets/stomics_0212_wilcoxon.yaml", 7, False, True),
    "stomics_0218": ("configs/datasets/stomics_0218_wilcoxon.yaml", 9, False, True),
    "stomics_0224": ("configs/datasets/stomics_0224_wilcoxon.yaml", 8, False, True),
    "e8p5_embryo": ("configs/datasets/e8p5_embryo.yaml", 12, False, False),
    "e9p5_embryo": ("configs/datasets/e9p5_embryo.yaml", 12, False, False),
}

FEATURE_NUMBER_N = [100, 200, 500, 1000, 5000, 10000]


def benchmark_methods(n_features: list[int], include_wilcoxon: bool) -> list[dict]:
    methods = [
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
        {"name": "scsegindex"},
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
    ]
    if include_wilcoxon:
        methods.append({"name": "wilcoxon"})
    for method in methods:
        method["n_features"] = n_features
    return methods


def canonical_methods(include_wilcoxon: bool) -> list[dict]:
    methods = benchmark_methods([2000], include_wilcoxon)
    for method in methods:
        if method["name"] == "scsegindex":
            method["n_features"] = [500]
    return [
        {"name": "all_features", "n_features": [1]},
        {"name": "random", "n_features": [500]},
        {"name": "TFs", "n_features": [2000]},
        *methods,
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
    include_wilcoxon: bool,
    phase: str,
) -> None:
    if phase == "canonical":
        methods = canonical_methods(include_wilcoxon)
        seeds = [0, 1, 2]
        default_n = [2000]
    else:
        methods = benchmark_methods(FEATURE_NUMBER_N, include_wilcoxon)
        seeds = [0]
        default_n = FEATURE_NUMBER_N
    payload = {
        "name": f"{key}_non_spatial_{phase}_rebuild_v1",
        "output_dir": f"{RESULTS_ROOT}/{key}",
        "datasets": [dataset_path],
        "feature_selection_methods": methods,
        "integration_methods": integration_methods(),
        "tasks": tasks(n_clusters, include_alignment),
        "n_features": default_n,
        "seeds": seeds,
        "save_embeddings": True,
    }
    phase_dir = OUTPUT_DIR / phase
    phase_dir.mkdir(parents=True, exist_ok=True)
    destination = phase_dir / f"{key}.yaml"
    destination.write_text(yaml.safe_dump(payload, sort_keys=False), encoding="utf-8")


def main() -> None:
    for key, (dataset_path, n_clusters, include_alignment, include_wilcoxon) in DATASETS.items():
        write_config(
            key,
            dataset_path,
            n_clusters,
            include_alignment,
            include_wilcoxon,
            "canonical",
        )
        write_config(
            key,
            dataset_path,
            n_clusters,
            include_alignment,
            include_wilcoxon,
            "feature_number",
        )


if __name__ == "__main__":
    main()

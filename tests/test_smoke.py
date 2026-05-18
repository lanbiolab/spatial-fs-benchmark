from __future__ import annotations

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import yaml

from spatial_fs_benchmark.benchmark.runner import BenchmarkRunner
from spatial_fs_benchmark.config import load_benchmark_config


def make_synthetic_dataset(path: Path) -> None:
    rng = np.random.default_rng(7)
    n_slices = 3
    cells_per_slice = 18
    genes = 40
    labels = np.array(["A", "B", "C"])
    matrices = []
    obs_frames = []
    coords = []
    for slice_idx in range(n_slices):
        slice_name = f"slice_{slice_idx}"
        slice_labels = np.repeat(labels, cells_per_slice // len(labels))
        base = rng.poisson(lam=1.0, size=(cells_per_slice, genes)).astype(float)
        for label_idx, label in enumerate(labels):
            mask = slice_labels == label
            base[mask, label_idx * 6 : (label_idx + 1) * 6] += 5.0
        base += slice_idx * 0.1
        matrices.append(base)
        obs_frames.append(
            pd.DataFrame(
                {
                    "slices": slice_name,
                    "original_domain": slice_labels,
                },
                index=[f"{slice_name}_cell_{i}" for i in range(cells_per_slice)],
            )
        )
        xy = rng.normal(size=(cells_per_slice, 2)) + slice_idx * 2.0
        coords.append(xy)
    adata = ad.AnnData(
        X=np.vstack(matrices),
        obs=pd.concat(obs_frames),
        var=pd.DataFrame(index=[f"gene_{i}" for i in range(genes)]),
    )
    adata.obsm["spatial"] = np.vstack(coords)
    adata.write_h5ad(path)


def test_smoke_benchmark(tmp_path: Path) -> None:
    dataset_path = tmp_path / "synthetic.h5ad"
    make_synthetic_dataset(dataset_path)
    dataset_config = {
        "name": "Synthetic",
        "path": str(dataset_path),
        "platform": "test",
        "slice_key": "slices",
        "label_key": "original_domain",
        "coord_key": "spatial",
        "preprocess": {
            "normalize_total": False,
            "log1p": False,
            "scale": False,
        },
    }
    dataset_config_path = tmp_path / "dataset.yaml"
    dataset_config_path.write_text(yaml.safe_dump(dataset_config), encoding="utf-8")
    benchmark_config = {
        "name": "synthetic_smoke",
        "output_dir": str(tmp_path / "results"),
        "datasets": [str(dataset_config_path)],
        "feature_selection_methods": [{"name": "highly_expressed", "params": {}}],
        "integration_methods": [{"name": "pca", "params": {"n_components": 8}}],
        "tasks": [
            {"name": "integration_eval", "params": {"n_neighbors": 5}},
            {"name": "clustering_eval", "params": {"n_clusters": 3, "n_neighbors": 4}},
            {"name": "alignment_eval", "params": {"n_neighbors": 3}},
            {"name": "slice_representation_eval", "params": {"n_clusters": 3}},
        ],
        "n_features": [12],
        "seeds": [0],
        "save_embeddings": True,
    }
    benchmark_config_path = tmp_path / "benchmark.yaml"
    benchmark_config_path.write_text(yaml.safe_dump(benchmark_config), encoding="utf-8")
    runner = BenchmarkRunner(load_benchmark_config(benchmark_config_path))
    runner.run()
    results_path = tmp_path / "results" / "results.csv"
    assert results_path.exists()
    results = pd.read_csv(results_path)
    assert not results.empty
    assert {"integration_eval", "clustering_eval", "alignment_eval", "slice_representation_eval"} <= set(results["task"])

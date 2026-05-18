#!/usr/bin/env python
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from spatial_fs_benchmark.config import load_dataset_config
from spatial_fs_benchmark.data.io import load_dataset
from spatial_fs_benchmark.feature_selection import build_feature_selector


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Smoke test a feature selector on a small dataset subset.")
    parser.add_argument("--dataset-config", required=True, help="Path to dataset YAML config.")
    parser.add_argument("--method", required=True, help="Feature selector name.")
    parser.add_argument("--n-features", type=int, default=20, help="Requested number of selected features.")
    parser.add_argument("--max-obs", type=int, default=120, help="Maximum number of observations to keep.")
    parser.add_argument("--max-vars", type=int, default=1500, help="Maximum number of genes to keep.")
    parser.add_argument("--seed", type=int, default=0, help="Random seed for downsampling.")
    parser.add_argument("--params", default="{}", help="JSON object with selector keyword arguments.")
    args = parser.parse_args()
    args.dataset_config = str(Path(args.dataset_config).resolve())
    return args


def subset_dataset(dataset_config: str, max_obs: int, max_vars: int, seed: int):
    dataset = load_dataset(load_dataset_config(dataset_config))
    adata = dataset.adata
    rng = np.random.default_rng(seed)
    if adata.n_obs > max_obs:
        obs_idx = np.sort(rng.choice(adata.n_obs, size=max_obs, replace=False))
        adata = adata[obs_idx].copy()
    else:
        adata = adata.copy()
    if adata.n_vars > max_vars:
        counts = adata.layers["counts"] if "counts" in adata.layers else adata.X
        dense = counts.toarray() if hasattr(counts, "toarray") else np.asarray(counts)
        var_idx = np.argsort(dense.var(axis=0))[::-1][:max_vars]
        adata = adata[:, np.sort(var_idx)].copy()
    dataset.adata = adata
    return dataset


def main() -> int:
    args = parse_args()
    params = json.loads(args.params)
    dataset = subset_dataset(args.dataset_config, args.max_obs, args.max_vars, args.seed)
    selector = build_feature_selector(args.method, **params)
    result = selector.select(dataset, n_features=args.n_features, random_seed=args.seed)
    payload = {
        "dataset": dataset.name,
        "method": args.method,
        "n_selected": len(result.feature_names),
        "preview": result.feature_names[:5],
        "metadata": result.metadata,
    }
    print(json.dumps(payload, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

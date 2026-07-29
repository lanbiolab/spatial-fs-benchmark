from __future__ import annotations

import argparse
import json
import time
from hashlib import md5
from pathlib import Path

import numpy as np
import yaml

from spatial_fs_benchmark.config import load_benchmark_config
from spatial_fs_benchmark.data.io import load_dataset
from spatial_fs_benchmark.data.preprocess import preprocess_dataset
from spatial_fs_benchmark.data.spatial_object import SpatialDataset
from spatial_fs_benchmark.feature_selection import build_feature_selector


ROOT = Path(__file__).resolve().parents[1]
NONSPATIAL_CONFIG = ROOT / "configs/rebuild_v1/non_spatial_downstream/canonical/mouse_brain_serial_sections.yaml"
SPATIAL_CONFIG = ROOT / "configs/rebuild_v1/downstream/canonical/mouse_brain_serial_sections.yaml"


def method_config(method_name: str) -> tuple[dict[str, object], int]:
    candidates: dict[str, dict] = {}
    for path in (NONSPATIAL_CONFIG, SPATIAL_CONFIG):
        payload = yaml.safe_load(path.read_text(encoding="utf-8"))
        candidates.update({item["name"]: item for item in payload["feature_selection_methods"]})
    if method_name not in candidates:
        raise KeyError(f"Unknown canonical method: {method_name}")
    item = candidates[method_name]
    requested = int(item.get("n_features", [2000])[0])
    return item.get("params", {}), requested


def standardized_dataset(seed: int = 1729) -> SpatialDataset:
    benchmark = load_benchmark_config(NONSPATIAL_CONFIG)
    dataset_config = benchmark.datasets[0]
    dataset = preprocess_dataset(load_dataset(dataset_config), dataset_config.preprocess)
    rng = np.random.default_rng(seed)
    selected_obs: list[int] = []
    for slice_name in dict.fromkeys(dataset.slice_ids.tolist()):
        indices = np.flatnonzero(dataset.slice_ids == slice_name)
        selected_obs.extend(rng.choice(indices, size=min(1000, indices.size), replace=False).tolist())
    selected_obs = sorted(selected_obs)

    counts = dataset.adata.layers["counts"][selected_obs]
    detection = np.asarray((counts > 0).sum(axis=0)).ravel()
    gene_order = np.lexsort((np.asarray(dataset.gene_names), -detection))[:5000]
    adata = dataset.adata[selected_obs, gene_order].copy()
    return SpatialDataset(
        name="MouseBrainSerialSections_resource_profile",
        adata=adata,
        slice_key=dataset.slice_key,
        coord_key=dataset.coord_key,
        label_key=dataset.label_key,
        slice_class_key=dataset.slice_class_key,
        platform=dataset.platform,
        species=dataset.species,
        source_path=dataset.source_path,
        alignment_pairs=list(dataset.alignment_pairs),
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--method", required=True)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    params, requested_n = method_config(args.method)
    dataset = standardized_dataset()
    if args.method == "all_features":
        requested_n = dataset.n_vars
    selector = build_feature_selector(args.method, **params)
    start = time.perf_counter()
    result = selector.select(dataset, n_features=requested_n, random_seed=0)
    elapsed = time.perf_counter() - start
    payload = {
        "method": args.method,
        "requested_n_features": requested_n,
        "effective_n_features": len(result.feature_names),
        "n_spots": dataset.n_obs,
        "n_genes": dataset.n_vars,
        "wall_seconds_internal": elapsed,
        "selector_params": params,
        "selector_implementation_version": getattr(selector, "implementation_version", None),
        "feature_hash": md5("\n".join(result.feature_names).encode("utf-8")).hexdigest(),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")


if __name__ == "__main__":
    main()

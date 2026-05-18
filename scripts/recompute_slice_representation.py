from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd

from spatial_fs_benchmark.benchmark.experiment import metric_records_from_task_output
from spatial_fs_benchmark.benchmark.result_schema import MetricRecord
from spatial_fs_benchmark.config import load_benchmark_config
from spatial_fs_benchmark.data.io import load_dataset
from spatial_fs_benchmark.data.preprocess import preprocess_dataset
from spatial_fs_benchmark.feature_selection.base import FeatureSelectionResult
from spatial_fs_benchmark.integration.base import IntegrationResult
from spatial_fs_benchmark.tasks import build_task


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Recompute slice representation metrics from saved embeddings.")
    parser.add_argument("--config", required=True, help="Benchmark config path.")
    parser.add_argument("--results-dir", required=True, help="Benchmark results directory.")
    parser.add_argument("--dataset", required=True, help="Dataset name to update.")
    return parser


def load_selection(path: Path) -> FeatureSelectionResult:
    payload = json.loads(path.read_text(encoding="utf-8"))
    return FeatureSelectionResult(
        method_name=payload["method_name"],
        feature_names=payload["feature_names"],
        feature_indices=payload["feature_indices"],
        scores=payload["scores"],
        metadata=payload.get("metadata", {}),
    )


def main() -> None:
    args = build_parser().parse_args()
    config = load_benchmark_config(args.config)
    results_dir = Path(args.results_dir)
    results_path = results_dir / "results.csv"
    results = pd.read_csv(results_path)

    dataset_config = next(dataset for dataset in config.datasets if dataset.name == args.dataset)
    dataset = preprocess_dataset(load_dataset(dataset_config), dataset_config.preprocess)

    slice_task_config = next(task for task in config.tasks if task.name == "slice_representation_eval")
    task = build_task(slice_task_config.name, **slice_task_config.params)

    dataset_dir = results_dir / dataset.name.lower()
    new_records: list[MetricRecord] = []
    combinations = (
        results[
            (results["dataset"] == dataset.name)
            & (results["task"] != "slice_representation_eval")
        ][["fs_method", "n_features", "integration_method", "random_seed"]]
        .drop_duplicates()
        .sort_values(["fs_method", "n_features", "integration_method", "random_seed"])
    )

    for combo in combinations.itertuples(index=False):
        selection_path = (
            dataset_dir
            / "feature_selection"
            / combo.fs_method
            / f"n{combo.n_features}"
            / f"seed{combo.random_seed}"
            / "selected_features.json"
        )
        embedding_path = (
            dataset_dir
            / combo.integration_method
            / combo.fs_method
            / f"n{combo.n_features}"
            / f"seed{combo.random_seed}"
            / "embedding.npz"
        )
        selection = load_selection(selection_path)
        embedding_payload = np.load(embedding_path)
        integration_result = IntegrationResult(
            method_name=combo.integration_method,
            embedding=embedding_payload["embedding"],
        )
        start = time.perf_counter()
        task_output = task.evaluate(
            dataset,
            selection,
            integration_result,
            random_seed=int(combo.random_seed),
        )
        runtime = time.perf_counter() - start
        new_records.extend(
            metric_records_from_task_output(
                dataset=dataset,
                fs_method=str(combo.fs_method),
                n_features=int(combo.n_features),
                effective_n_features=int(selection.metadata.get("effective_n_features", len(selection.feature_names))),
                integration_method=str(combo.integration_method),
                task_name="slice_representation_eval",
                task_output=task_output,
                random_seed=int(combo.random_seed),
                runtime=runtime,
            )
        )

    updated = results[
        ~(
            (results["dataset"] == dataset.name)
            & (results["task"] == "slice_representation_eval")
        )
    ].copy()
    updated = pd.concat([updated, pd.DataFrame([record.to_dict() for record in new_records])], ignore_index=True)
    updated.sort_values(
        ["dataset", "fs_method", "n_features", "integration_method", "task", "metric_name", "random_seed"],
        inplace=True,
    )
    updated.to_csv(results_path, index=False)
    print(f"Recomputed {len(new_records)} slice representation metric rows for {dataset.name}")


if __name__ == "__main__":
    main()

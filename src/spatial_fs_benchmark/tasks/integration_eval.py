from __future__ import annotations

import numpy as np

from spatial_fs_benchmark.data.spatial_object import SpatialDataset
from spatial_fs_benchmark.feature_selection.base import FeatureSelectionResult
from spatial_fs_benchmark.integration.base import IntegrationResult
from spatial_fs_benchmark.metrics.integration_paper import (
    graph_connectivity_score,
    isolated_label_f1,
    scaled_lisi,
    scaled_silhouette_batch,
    scaled_silhouette_labels,
)
from spatial_fs_benchmark.tasks.base import BenchmarkTask, TaskOutput


def _sample_by_slice(slice_ids: np.ndarray, max_cells: int, random_seed: int) -> np.ndarray:
    if len(slice_ids) <= max_cells:
        return np.arange(len(slice_ids))
    rng = np.random.default_rng(random_seed)
    unique_slices = np.unique(slice_ids)
    per_slice = max(1, max_cells // max(1, len(unique_slices)))
    selected = []
    for slice_id in unique_slices:
        idx = np.flatnonzero(slice_ids == slice_id)
        take = min(len(idx), per_slice)
        selected.append(rng.choice(idx, size=take, replace=False))
    sampled = np.concatenate(selected)
    if len(sampled) > max_cells:
        sampled = rng.choice(sampled, size=max_cells, replace=False)
    return np.sort(sampled)


class IntegrationEvaluationTask(BenchmarkTask):
    name = "integration_eval"

    def __init__(self, n_neighbors: int = 15, max_eval_cells: int | None = None) -> None:
        self.n_neighbors = n_neighbors
        self.max_eval_cells = max_eval_cells

    def evaluate(
        self,
        dataset: SpatialDataset,
        selected_features: FeatureSelectionResult,
        integration_result: IntegrationResult,
        random_seed: int = 0,
    ) -> TaskOutput:
        indices = (
            _sample_by_slice(dataset.slice_ids, self.max_eval_cells, random_seed)
            if self.max_eval_cells is not None
            else np.arange(dataset.n_obs)
        )
        embedding = integration_result.embedding[indices]
        slice_ids = dataset.slice_ids[indices]
        labels = dataset.labels[indices] if dataset.labels is not None else None
        metrics = {
            "bASW": scaled_silhouette_batch(
                embedding,
                slice_ids,
                labels if labels is not None else slice_ids,
            ),
            "iLISI": scaled_lisi(
                embedding,
                slice_ids,
                n_neighbors=self.n_neighbors,
                mode="batch",
            ),
        }
        if labels is not None:
            metrics["dASW"] = scaled_silhouette_labels(embedding, labels)
            metrics["dLISI"] = scaled_lisi(
                embedding,
                labels,
                n_neighbors=self.n_neighbors,
                mode="label",
            )
            metrics["ILL"] = isolated_label_f1(
                embedding,
                labels,
                slice_ids,
                n_neighbors=self.n_neighbors,
            )
            metrics["GC"] = graph_connectivity_score(
                embedding,
                labels,
                n_neighbors=self.n_neighbors,
            )
        else:
            metrics["dASW"] = float("nan")
            metrics["dLISI"] = float("nan")
            metrics["ILL"] = float("nan")
            metrics["GC"] = float("nan")
        return TaskOutput(metrics=metrics)

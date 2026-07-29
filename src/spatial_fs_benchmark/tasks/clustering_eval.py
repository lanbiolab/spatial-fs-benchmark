from __future__ import annotations

import numpy as np

from spatial_fs_benchmark.data.spatial_object import SpatialDataset
from spatial_fs_benchmark.feature_selection.base import FeatureSelectionResult
from spatial_fs_benchmark.integration.base import IntegrationResult
from spatial_fs_benchmark.metrics.clustering_metrics import cluster_embedding, clustering_scores
from spatial_fs_benchmark.metrics.spatial_clustering import chaos_score, pas_score
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


class ClusteringEvaluationTask(BenchmarkTask):
    name = "clustering_eval"
    implementation_version = "v2_kmeans_slice_balanced_spatial_metrics"

    def __init__(
        self,
        n_clusters: int = 6,
        n_neighbors: int = 8,
        silhouette_sample_size: int | None = 10000,
        max_eval_cells: int | None = None,
    ) -> None:
        self.n_clusters = n_clusters
        self.n_neighbors = n_neighbors
        self.silhouette_sample_size = silhouette_sample_size
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
        predicted = cluster_embedding(embedding, self.n_clusters, random_seed=random_seed)
        true_labels = dataset.labels[indices] if dataset.labels is not None else None
        metrics = clustering_scores(
            embedding,
            predicted_labels=predicted,
            true_labels=true_labels,
            random_seed=random_seed,
            silhouette_sample_size=self.silhouette_sample_size,
        )
        per_slice_chaos = []
        per_slice_pas = []
        sampled_slices = dataset.slice_ids[indices]
        sampled_coords = dataset.coords[indices]
        for slice_id in np.unique(sampled_slices):
            mask = sampled_slices == slice_id
            per_slice_chaos.append(chaos_score(predicted[mask], sampled_coords[mask]))
            per_slice_pas.append(pas_score(predicted[mask], sampled_coords[mask], n_neighbors=self.n_neighbors))
        metrics["CHAOS"] = float(np.mean(per_slice_chaos))
        metrics["PAS"] = float(np.mean(per_slice_pas))
        return TaskOutput(metrics=metrics, artifacts={"predicted_clusters": predicted.tolist()})

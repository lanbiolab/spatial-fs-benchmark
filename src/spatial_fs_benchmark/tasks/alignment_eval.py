from __future__ import annotations

from spatial_fs_benchmark.data.spatial_object import SpatialDataset
from spatial_fs_benchmark.feature_selection.base import FeatureSelectionResult
from spatial_fs_benchmark.integration.base import IntegrationResult
from spatial_fs_benchmark.metrics.alignment_metrics import alignment_accuracy_ratio
from spatial_fs_benchmark.tasks.base import BenchmarkTask, TaskOutput


class AlignmentEvaluationTask(BenchmarkTask):
    name = "alignment_eval"
    implementation_version = "v2_explicit_slice_pairs"

    def __init__(self, n_neighbors: int = 10, max_points_per_slice: int | None = None) -> None:
        self.n_neighbors = n_neighbors
        self.max_points_per_slice = max_points_per_slice

    def evaluate(
        self,
        dataset: SpatialDataset,
        selected_features: FeatureSelectionResult,
        integration_result: IntegrationResult,
        random_seed: int = 0,
    ) -> TaskOutput:
        if dataset.labels is None:
            return TaskOutput(metrics={"Accuracy": float("nan"), "Ratio": float("nan")})
        representation = integration_result.metadata.get("aligned_representation", integration_result.embedding)
        accuracy, ratio = alignment_accuracy_ratio(
            representation,
            dataset.labels,
            dataset.slice_ids,
            alignment_pairs=dataset.alignment_pairs,
            max_points_per_slice=self.max_points_per_slice,
            random_seed=random_seed,
        )
        metrics = {"Accuracy": accuracy, "Ratio": ratio}
        return TaskOutput(
            metrics=metrics,
            artifacts={"alignment_pairs": [list(pair) for pair in dataset.alignment_pairs]},
        )

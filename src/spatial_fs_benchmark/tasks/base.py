from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

from spatial_fs_benchmark.data.spatial_object import SpatialDataset
from spatial_fs_benchmark.feature_selection.base import FeatureSelectionResult
from spatial_fs_benchmark.integration.base import IntegrationResult


@dataclass(slots=True)
class TaskOutput:
    metrics: dict[str, float]
    artifacts: dict[str, Any] = field(default_factory=dict)


class BenchmarkTask:
    name = "base"

    def evaluate(
        self,
        dataset: SpatialDataset,
        selected_features: FeatureSelectionResult,
        integration_result: IntegrationResult,
        random_seed: int = 0,
    ) -> TaskOutput:
        raise NotImplementedError

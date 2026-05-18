from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import numpy as np
import pandas as pd

from spatial_fs_benchmark.data.spatial_object import SpatialDataset
from spatial_fs_benchmark.feature_selection.base import FeatureSelectionResult


@dataclass(slots=True)
class IntegrationResult:
    method_name: str
    embedding: np.ndarray
    corrected_matrix: np.ndarray | None = None
    slice_embedding: pd.DataFrame | None = None
    metadata: dict[str, Any] = field(default_factory=dict)


class SpatialIntegrator:
    name = "base"

    def fit_transform(
        self,
        dataset: SpatialDataset,
        selected_features: FeatureSelectionResult,
        random_seed: int = 0,
    ) -> IntegrationResult:
        raise NotImplementedError

    @staticmethod
    def _build_slice_embedding(dataset: SpatialDataset, embedding: np.ndarray) -> pd.DataFrame:
        frame = pd.DataFrame(embedding, index=dataset.adata.obs_names)
        frame["slice_id"] = dataset.slice_ids
        return frame.groupby("slice_id", observed=True).mean(numeric_only=True)

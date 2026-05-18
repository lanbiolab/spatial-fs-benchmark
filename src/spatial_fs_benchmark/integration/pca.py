from __future__ import annotations

from hashlib import md5

from sklearn.decomposition import PCA

from spatial_fs_benchmark.data.spatial_object import SpatialDataset
from spatial_fs_benchmark.feature_selection.base import FeatureSelectionResult
from spatial_fs_benchmark.integration.base import IntegrationResult, SpatialIntegrator


class PCAIntegrator(SpatialIntegrator):
    name = "pca"

    def __init__(self, n_components: int = 30) -> None:
        self.n_components = n_components

    def fit_transform(
        self,
        dataset: SpatialDataset,
        selected_features: FeatureSelectionResult,
        random_seed: int = 0,
    ) -> IntegrationResult:
        subset = dataset.subset_features(selected_features.feature_names)
        matrix = subset.to_dense_matrix()
        pca = PCA(n_components=min(self.n_components, matrix.shape[1], matrix.shape[0] - 1), random_state=random_seed)
        embedding = pca.fit_transform(matrix)
        return IntegrationResult(
            method_name=self.name,
            embedding=embedding,
            slice_embedding=self._build_slice_embedding(dataset, embedding),
            metadata={
                "explained_variance_ratio": pca.explained_variance_ratio_.tolist(),
                "n_input_features": int(matrix.shape[1]),
                "input_feature_hash": md5("\n".join(selected_features.feature_names).encode("utf-8")).hexdigest(),
            },
        )

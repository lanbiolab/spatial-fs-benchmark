from __future__ import annotations

from hashlib import md5

import anndata as ad
import numpy as np
import scanorama
from scipy import sparse

from spatial_fs_benchmark.data.spatial_object import SpatialDataset
from spatial_fs_benchmark.feature_selection.base import FeatureSelectionResult
from spatial_fs_benchmark.integration.base import IntegrationResult, SpatialIntegrator


class ScanoramaIntegrator(SpatialIntegrator):
    name = "scanorama"

    def __init__(self, n_components: int = 30) -> None:
        self.n_components = n_components

    def fit_transform(
        self,
        dataset: SpatialDataset,
        selected_features: FeatureSelectionResult,
        random_seed: int = 0,
    ) -> IntegrationResult:
        subset = dataset.subset_features(selected_features.feature_names)
        adatas: list[ad.AnnData] = []
        slice_order: list[str] = []
        for slice_id in np.unique(subset.slice_ids):
            slice_mask = subset.adata.obs[subset.slice_key].astype(str) == slice_id
            slice_adata = subset.adata[slice_mask].copy()
            if sparse.issparse(slice_adata.X):
                slice_adata.X = slice_adata.X.tocsr()
            adatas.append(slice_adata)
            slice_order.append(str(slice_id))
        corrected = scanorama.correct_scanpy(adatas, return_dimred=True, dimred=self.n_components)
        per_slice_embeddings = [adata_slice.obsm["X_scanorama"] for adata_slice in corrected]
        embedding = np.vstack(per_slice_embeddings)
        obs_index = np.concatenate([adata_slice.obs_names.to_numpy() for adata_slice in corrected])
        order = {obs_name: idx for idx, obs_name in enumerate(obs_index)}
        reordered = embedding[[order[name] for name in subset.adata.obs_names]]
        return IntegrationResult(
            method_name=self.name,
            embedding=reordered,
            slice_embedding=self._build_slice_embedding(dataset, reordered),
            metadata={
                "slice_order": slice_order,
                "n_components": self.n_components,
                "n_input_features": int(subset.n_vars),
                "input_feature_hash": md5("\n".join(selected_features.feature_names).encode("utf-8")).hexdigest(),
            },
        )

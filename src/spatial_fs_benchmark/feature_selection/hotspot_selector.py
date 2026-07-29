from __future__ import annotations

import numpy as np
import scanpy as sc
from scipy import sparse

from spatial_fs_benchmark.data.spatial_object import SpatialDataset
from spatial_fs_benchmark.feature_selection.base import FeatureSelectionResult, FeatureSelector


class HotspotSelector(FeatureSelector):
    name = "hotspot"
    implementation_version = "v2_seeded_subsampling"

    def __init__(self, n_neighbors: int = 30, max_cells: int | None = None) -> None:
        self.n_neighbors = n_neighbors
        self.max_cells = max_cells
        self.stochastic_selection = max_cells is not None

    def select(self, dataset: SpatialDataset, n_features: int, random_seed: int = 0) -> FeatureSelectionResult:
        import hotspot

        gene_names = dataset.adata.var_names.to_numpy()
        score_map = {gene: 0.0 for gene in gene_names}

        adata = self._maybe_subsample_obs(self._adata_from_counts(dataset), self.max_cells, random_seed)
        sc.pp.filter_cells(adata, min_counts=1)
        sc.pp.filter_genes(adata, min_cells=1)
        adata.obs["total_counts"] = np.asarray(adata.X.sum(axis=1)).ravel()
        adata.layers["counts"] = adata.X.toarray() if sparse.issparse(adata.X) else np.asarray(adata.X.copy())

        sc.pp.normalize_total(adata)
        sc.pp.log1p(adata)
        sc.pp.scale(adata)
        sc.tl.pca(adata)

        hs = hotspot.Hotspot(
            adata,
            layer_key="counts",
            model="danb",
            latent_obsm_key="X_pca",
            umi_counts_obs_key="total_counts",
        )
        hs.create_knn_graph(weighted_graph=False, n_neighbors=self.n_neighbors)
        hs_results = hs.compute_autocorrelations()
        if "FDR" in hs_results.columns:
            hs_results = hs_results.loc[hs_results["FDR"] < 0.05].copy()
        for gene_name, score in hs_results["Z"].items():
            score_map[str(gene_name)] = float(score)

        scores = np.asarray([score_map[gene] for gene in gene_names], dtype=float)
        return self._build_result(
            method_name=self.name,
            gene_names=gene_names,
            scores=scores,
            n_features=n_features,
            metadata={"n_neighbors": self.n_neighbors, "max_cells": self.max_cells},
        )

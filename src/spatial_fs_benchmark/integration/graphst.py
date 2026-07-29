from __future__ import annotations

import sys
from hashlib import md5
from pathlib import Path

import numpy as np
import scanpy as sc
import torch
from scipy import sparse
from sklearn.decomposition import PCA
from sklearn.neighbors import NearestNeighbors

from spatial_fs_benchmark.data.spatial_object import SpatialDataset
from spatial_fs_benchmark.data.preprocess import is_nonnegative_integer_matrix
from spatial_fs_benchmark.feature_selection.base import FeatureSelectionResult
from spatial_fs_benchmark.integration.base import IntegrationResult, SpatialIntegrator


def _graphst_module_path() -> Path:
    return (
        Path(__file__).resolve().parents[3]
        / "external"
        / "iSTBench"
        / "Benchmark"
        / "RunModel"
        / "GraphST"
    )


def _build_sparse_slice_graph(
    coordinates: np.ndarray,
    slice_ids: np.ndarray,
    n_neighbors: int,
) -> tuple[sparse.csr_matrix, sparse.csr_matrix]:
    """Build directed-neighbor and symmetric-adjacency matrices without cross-slice edges."""
    neighbor_blocks: list[sparse.csr_matrix] = []
    adjacency_blocks: list[sparse.csr_matrix] = []
    for slice_name in dict.fromkeys(slice_ids.tolist()):
        local = coordinates[slice_ids == slice_name]
        n_spots = local.shape[0]
        k = min(n_neighbors, max(n_spots - 1, 0))
        if k == 0:
            directed = sparse.csr_matrix((n_spots, n_spots), dtype=np.float32)
        else:
            indices = NearestNeighbors(n_neighbors=k + 1).fit(local).kneighbors(
                local, return_distance=False
            )
            rows = np.repeat(np.arange(n_spots), k)
            cols = indices[:, 1 : k + 1].reshape(-1)
            directed = sparse.csr_matrix(
                (np.ones(rows.size, dtype=np.float32), (rows, cols)),
                shape=(n_spots, n_spots),
            )
        neighbor_blocks.append(directed)
        adjacency_blocks.append(directed.maximum(directed.T).tocsr())
    return (
        sparse.block_diag(neighbor_blocks, format="csr", dtype=np.float32),
        sparse.block_diag(adjacency_blocks, format="csr", dtype=np.float32),
    )


class GraphSTIntegrator(SpatialIntegrator):
    name = "graphst"
    implementation_version = "v2_sparse_counts_official_preprocessing_pca20"

    def __init__(
        self,
        epochs: int = 600,
        dim_output: int = 32,
        pca_components: int = 20,
        n_neighbors: int = 3,
        device: str = "cpu",
    ) -> None:
        self.epochs = epochs
        self.dim_output = dim_output
        self.pca_components = pca_components
        self.n_neighbors = n_neighbors
        self.device = device

    def fit_transform(
        self,
        dataset: SpatialDataset,
        selected_features: FeatureSelectionResult,
        random_seed: int = 0,
    ) -> IntegrationResult:
        module_root = _graphst_module_path()
        if str(module_root) not in sys.path:
            sys.path.insert(0, str(module_root))
        from GraphST import GraphST as graphst_module
        from GraphST.preprocess import add_contrastive_label, get_feature

        subset = dataset.subset_features(selected_features.feature_names)
        adata = subset.adata.copy()
        adata.obs["slices"] = adata.obs[subset.slice_key].astype(str)

        # GraphST assumes observations from each spatial graph are contiguous.
        original_slice_ids = adata.obs["slices"].to_numpy()
        slice_order = list(dict.fromkeys(original_slice_ids.tolist()))
        grouped_order = np.concatenate(
            [np.flatnonzero(original_slice_ids == slice_name) for slice_name in slice_order]
        )
        adata = adata[grouped_order].copy()
        grouped_slice_ids = adata.obs["slices"].to_numpy()

        # Follow the official GraphST expression preprocessing after feature selection.
        if "counts" in adata.layers:
            adata.X = adata.layers["counts"].copy()
            adata.uns.pop("log1p", None)
        if not is_nonnegative_integer_matrix(adata.X):
            source = adata.uns.get("spatial_fs_benchmark", {}).get("counts_source", "unknown")
            raise ValueError(
                f"GraphST requires non-negative integer counts before normalization, but "
                f"dataset '{dataset.name}' uses non-count input from '{source}'."
            )
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        sc.pp.scale(adata, zero_center=False, max_value=10)
        adata.var["highly_variable"] = True

        graph_neigh, adjacency = _build_sparse_slice_graph(
            np.asarray(adata.obsm["spatial"], dtype=np.float32),
            grouped_slice_ids,
            self.n_neighbors,
        )
        adata.obsm["adj"] = adjacency
        adata.obsm["graph_neigh"] = graph_neigh
        add_contrastive_label(adata)
        get_feature(adata)

        device = torch.device(self.device if self.device != "auto" else ("cuda" if torch.cuda.is_available() else "cpu"))
        model = graphst_module.GraphST(
            adata,
            device=device,
            epochs=self.epochs,
            dim_output=self.dim_output,
            random_seed=random_seed,
            datatype="Slide",
        )
        trained = model.train()
        reconstructed = np.asarray(trained.obsm["emb"], dtype=np.float32)
        n_components = min(self.pca_components, reconstructed.shape[0] - 1, reconstructed.shape[1])
        if n_components < 2:
            raise ValueError("GraphST requires at least three observations and two selected features")
        grouped_embedding = PCA(
            n_components=n_components,
            svd_solver="randomized",
            random_state=random_seed,
        ).fit_transform(reconstructed).astype(np.float32)
        embedding = np.empty_like(grouped_embedding)
        embedding[grouped_order] = grouped_embedding
        return IntegrationResult(
            method_name=self.name,
            embedding=embedding,
            slice_embedding=self._build_slice_embedding(dataset, embedding),
            metadata={
                "n_input_features": int(adata.n_vars),
                "input_feature_hash": md5("\n".join(selected_features.feature_names).encode("utf-8")).hexdigest(),
                "epochs": self.epochs,
                "dim_output": self.dim_output,
                "pca_components": n_components,
                "n_neighbors": self.n_neighbors,
                "device": str(device),
                "implementation_version": self.implementation_version,
                "input_assay": "counts",
                "counts_source": adata.uns.get("spatial_fs_benchmark", {}).get("counts_source", "unknown"),
            },
        )

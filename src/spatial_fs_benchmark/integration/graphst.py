from __future__ import annotations

import sys
from hashlib import md5
from pathlib import Path

import numpy as np
import torch
from scipy import linalg

from spatial_fs_benchmark.data.spatial_object import SpatialDataset
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


class GraphSTIntegrator(SpatialIntegrator):
    name = "graphst"

    def __init__(
        self,
        epochs: int = 200,
        dim_output: int = 32,
        n_neighbors: int = 3,
        device: str = "cpu",
    ) -> None:
        self.epochs = epochs
        self.dim_output = dim_output
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
        from GraphST.preprocess import add_contrastive_label, construct_interaction, get_feature

        subset = dataset.subset_features(selected_features.feature_names)
        adata = subset.adata.copy()
        adata.obs["slices"] = adata.obs[subset.slice_key].astype(str)
        adata.var["highly_variable"] = True

        # Build a block-diagonal spatial graph so slices are not spuriously connected by shared coordinates.
        adj_blocks = []
        neigh_blocks = []
        feat_blocks = []
        feat_aug_blocks = []
        label_blocks = []
        for slice_name in np.unique(subset.slice_ids):
            slice_adata = adata[adata.obs["slices"] == slice_name].copy()
            construct_interaction(slice_adata, n_neighbors=self.n_neighbors)
            add_contrastive_label(slice_adata)
            get_feature(slice_adata)
            adj_blocks.append(np.asarray(slice_adata.obsm["adj"], dtype=np.float32))
            neigh_blocks.append(np.asarray(slice_adata.obsm["graph_neigh"], dtype=np.float32))
            feat_blocks.append(np.asarray(slice_adata.obsm["feat"], dtype=np.float32))
            feat_aug_blocks.append(np.asarray(slice_adata.obsm["feat_a"], dtype=np.float32))
            label_blocks.append(np.asarray(slice_adata.obsm["label_CSL"], dtype=np.float32))

        adata.obsm["adj"] = linalg.block_diag(*adj_blocks).astype(np.float32)
        adata.obsm["graph_neigh"] = linalg.block_diag(*neigh_blocks).astype(np.float32)
        adata.obsm["feat"] = np.vstack(feat_blocks).astype(np.float32)
        adata.obsm["feat_a"] = np.vstack(feat_aug_blocks).astype(np.float32)
        adata.obsm["label_CSL"] = np.vstack(label_blocks).astype(np.float32)

        device = torch.device(self.device if self.device != "auto" else ("cuda" if torch.cuda.is_available() else "cpu"))
        model = graphst_module.GraphST(
            adata,
            device=device,
            epochs=self.epochs,
            dim_output=self.dim_output,
            random_seed=random_seed,
        )
        trained = model.train()
        embedding = np.asarray(trained.obsm["emb"], dtype=np.float32)
        return IntegrationResult(
            method_name=self.name,
            embedding=embedding,
            slice_embedding=self._build_slice_embedding(dataset, embedding),
            metadata={
                "n_input_features": int(adata.n_vars),
                "input_feature_hash": md5("\n".join(selected_features.feature_names).encode("utf-8")).hexdigest(),
                "epochs": self.epochs,
                "dim_output": self.dim_output,
                "n_neighbors": self.n_neighbors,
                "device": str(device),
            },
        )

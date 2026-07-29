from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse


@dataclass(slots=True)
class SpatialDataset:
    name: str
    adata: ad.AnnData
    slice_key: str
    coord_key: str
    label_key: str | None = None
    slice_class_key: str | None = None
    platform: str | None = None
    species: str | None = None
    source_path: str | None = None
    alignment_pairs: list[tuple[str, str]] = field(default_factory=list)

    @property
    def n_obs(self) -> int:
        return int(self.adata.n_obs)

    @property
    def n_vars(self) -> int:
        return int(self.adata.n_vars)

    @property
    def gene_names(self) -> list[str]:
        return self.adata.var_names.tolist()

    @property
    def slice_ids(self) -> np.ndarray:
        return self.adata.obs[self.slice_key].astype(str).to_numpy()

    @property
    def labels(self) -> np.ndarray | None:
        if self.label_key is None or self.label_key not in self.adata.obs:
            return None
        return self.adata.obs[self.label_key].astype(str).to_numpy()

    @property
    def slice_classes(self) -> np.ndarray | None:
        if self.slice_class_key is None or self.slice_class_key not in self.adata.obs:
            return None
        return self.adata.obs[self.slice_class_key].astype(str).to_numpy()

    @property
    def coords(self) -> np.ndarray:
        return np.asarray(self.adata.obsm[self.coord_key], dtype=float)

    def to_dense_matrix(self) -> np.ndarray:
        matrix = self.adata.X
        if sparse.issparse(matrix):
            return matrix.toarray().astype(np.float32, copy=False)
        return np.asarray(matrix, dtype=np.float32)

    def subset_features(self, feature_names: list[str]) -> "SpatialDataset":
        subset = self.adata[:, feature_names].copy()
        return SpatialDataset(
            name=self.name,
            adata=subset,
            slice_key=self.slice_key,
            coord_key=self.coord_key,
            label_key=self.label_key,
            slice_class_key=self.slice_class_key,
            platform=self.platform,
            species=self.species,
            source_path=self.source_path,
            alignment_pairs=list(self.alignment_pairs),
        )

    def copy(self) -> "SpatialDataset":
        return SpatialDataset(
            name=self.name,
            adata=self.adata.copy(),
            slice_key=self.slice_key,
            coord_key=self.coord_key,
            label_key=self.label_key,
            slice_class_key=self.slice_class_key,
            platform=self.platform,
            species=self.species,
            source_path=self.source_path,
            alignment_pairs=list(self.alignment_pairs),
        )

    def with_embedding(self, embedding: np.ndarray, key: str = "X_emb") -> "SpatialDataset":
        dataset = self.copy()
        dataset.adata.obsm[key] = embedding
        return dataset

    def slice_index(self) -> pd.Series:
        return self.adata.obs[self.slice_key].astype(str)

    def save_h5ad(self, path: str | Path) -> None:
        destination = Path(path)
        destination.parent.mkdir(parents=True, exist_ok=True)
        self.adata.write_h5ad(destination)

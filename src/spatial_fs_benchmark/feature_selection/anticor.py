from __future__ import annotations

from pathlib import Path
from tempfile import gettempdir

import numpy as np
import scanpy as sc

from spatial_fs_benchmark.data.spatial_object import SpatialDataset
from spatial_fs_benchmark.feature_selection.base import FeatureSelectionResult, FeatureSelector


class AnticorSelector(FeatureSelector):
    name = "anticor"

    def __init__(self, temp_dir: str | None = None, max_cells: int | None = None) -> None:
        self.temp_dir = temp_dir
        self.max_cells = max_cells

    def select(self, dataset: SpatialDataset, n_features: int, random_seed: int = 0) -> FeatureSelectionResult:
        from anticor_features.anticor_features import get_anti_cor_genes

        gene_names = dataset.adata.var_names.to_numpy()
        score_map = {gene: 0.0 for gene in gene_names}

        adata = self._maybe_subsample_obs(self._adata_from_counts(dataset), self.max_cells, random_seed)
        sc.pp.normalize_total(adata)
        sc.pp.log1p(adata)

        scratch_dir = Path(self.temp_dir or Path(gettempdir()) / "anticor" / dataset.name.lower())
        scratch_dir.mkdir(parents=True, exist_ok=True)
        anticor = get_anti_cor_genes(
            adata.X.T,
            adata.var_names.tolist(),
            pre_remove_pathways=[],
            scratch_dir=str(scratch_dir),
        )

        fdr = anticor.get("FDR")
        if fdr is not None:
            fdr = np.asarray(fdr, dtype=float)
            scores = -np.log10(np.clip(fdr, 1e-300, None))
        else:
            scores = anticor.get("selected", False).astype(float).to_numpy(dtype=float)
        scores = np.nan_to_num(scores, nan=0.0, posinf=0.0, neginf=0.0)
        for gene_name, score in zip(anticor["gene"], scores, strict=True):
            score_map[str(gene_name)] = float(score)

        score_array = np.asarray([score_map[gene] for gene in gene_names], dtype=float)
        return self._build_result(
            method_name=self.name,
            gene_names=gene_names,
            scores=score_array,
            n_features=n_features,
            metadata={"temp_dir": str(scratch_dir), "max_cells": self.max_cells},
        )

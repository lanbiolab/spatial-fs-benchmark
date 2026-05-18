from __future__ import annotations

import math

import numpy as np
import pandas as pd
import scanpy as sc

from spatial_fs_benchmark.data.spatial_object import SpatialDataset
from spatial_fs_benchmark.feature_selection.base import FeatureSelectionResult, FeatureSelector


class WilcoxonSelector(FeatureSelector):
    name = "wilcoxon"

    def select(self, dataset: SpatialDataset, n_features: int, random_seed: int = 0) -> FeatureSelectionResult:
        if dataset.label_key is None or dataset.labels is None:
            raise ValueError("WilcoxonSelector requires dataset labels.")

        adata = self._adata_from_counts(dataset)
        adata.obs["Label"] = pd.Categorical(dataset.labels)
        adata.obs["Label"] = adata.obs["Label"].cat.remove_unused_categories()
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)

        sc.tl.rank_genes_groups(adata, groupby="Label", method="wilcoxon", tie_correct=True)

        n_labels = len(adata.obs["Label"].cat.categories)
        per_label = max(1, math.ceil(n_features / max(1, n_labels)))
        gene_scores = pd.Series(0.0, index=adata.var_names, dtype=float)

        for label in adata.obs["Label"].cat.categories:
            frame = sc.get.rank_genes_groups_df(adata, group=str(label))
            frame = frame[frame["names"].notnull()]
            if "pvals_adj" in frame.columns:
                frame = frame[frame["pvals_adj"] <= 0.01]
            if "logfoldchanges" in frame.columns:
                frame = frame[frame["logfoldchanges"] > 0]
            if frame.empty:
                continue
            frame = frame.sort_values(by="logfoldchanges", ascending=False).head(per_label)
            scores = np.linspace(per_label, 1, num=len(frame), dtype=float)
            for gene_name, score in zip(frame["names"], scores, strict=True):
                gene_scores.loc[gene_name] = max(float(score), float(gene_scores.loc[gene_name]))

        return self._build_result(
            method_name=self.name,
            gene_names=adata.var_names.to_numpy(),
            scores=gene_scores.to_numpy(dtype=float),
            n_features=n_features,
            metadata={"selection_mode": "per_label_balanced"},
        )

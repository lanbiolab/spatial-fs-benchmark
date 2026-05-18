from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from spatial_fs_benchmark.data.spatial_object import SpatialDataset
from spatial_fs_benchmark.feature_selection.base import FeatureSelectionResult, FeatureSelector


class TFSelector(FeatureSelector):
    name = "TFs"

    def __init__(self, tfs_file: str | None = None) -> None:
        project_root = Path(__file__).resolve().parents[3]
        self.tfs_file = Path(tfs_file) if tfs_file is not None else project_root / "data" / "resources" / "human_tfs.tsv"

    def select(self, dataset: SpatialDataset, n_features: int, random_seed: int = 0) -> FeatureSelectionResult:
        if dataset.species is None:
            raise ValueError("TFSelector requires dataset.species to be set.")
        if not self.tfs_file.exists():
            raise FileNotFoundError(f"TF gene list not found: {self.tfs_file}")

        tfs = pd.read_csv(self.tfs_file, sep="\t")
        species = str(dataset.species).title()
        tfs = tfs[tfs["Species"] == species].copy()
        if tfs.empty:
            raise ValueError(f"No TF genes found for species '{species}' in {self.tfs_file}.")

        gene_names = dataset.adata.var_names.to_numpy()
        id_col = "ENSEMBL" if all(str(name).startswith("ENS") for name in gene_names[: min(50, len(gene_names))]) else "Gene"
        if id_col == "ENSEMBL":
            gene_set = set(map(str, gene_names.tolist()))
            ranked_features = [str(gene) for gene in tfs[id_col].dropna().tolist() if str(gene) in gene_set]
        else:
            canonical_to_feature = {str(gene).upper(): str(gene) for gene in gene_names.tolist()}
            ranked_features = []
            seen: set[str] = set()
            for gene in tfs[id_col].dropna().astype(str).tolist():
                matched = canonical_to_feature.get(gene.upper())
                if matched is None or matched in seen:
                    continue
                seen.add(matched)
                ranked_features.append(matched)
        if not ranked_features:
            raise ValueError(f"No TF genes from {self.tfs_file} were present in dataset '{dataset.name}'.")

        scores = np.linspace(len(ranked_features), 1, num=len(ranked_features), dtype=float).tolist()
        return self._build_ranked_result(
            method_name=self.name,
            gene_names=gene_names,
            ranked_features=ranked_features,
            ranked_scores=scores,
            n_features=n_features,
            metadata={"species": species, "tfs_file": str(self.tfs_file)},
        )

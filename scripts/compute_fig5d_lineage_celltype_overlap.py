#!/usr/bin/env python3

from __future__ import annotations

import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc


ROOT = Path(__file__).resolve().parents[1]

FULL_H5AD = ROOT / "data/raw/wilcoxon_candidate_0212/stomics_0212_multislice.h5ad"
METHODS_TSV = ROOT / "results/fig5a_spatial_lineages/figures/fig5a_lineage_methods.tsv"
OUT_SCORES = ROOT / "results/fig5d_spatial_lineages/figures/fig5d_lineage_celltype_overlap.tsv"
OUT_DIFFS = ROOT / "results/fig5d_spatial_lineages/figures/fig5d_lineage_celltype_overlap_diff.tsv"

LINEAGES = {
    "Immune": ["Langerhans Cell", "T Cell"],
    "Epithelial": ["Keratinocyte", "Merkel Cell", "Stem Cell", "Progenitor Cell"],
}

DATASET_PATHS = {
    "Full": ROOT / "results/spatial_main_native_seed0_fix3/stomics_0212_wilcoxon_spatial_main/stomics0212/feature_selection",
    "Immune": ROOT / "results/lineage_subsets/stomics_0212_immune_subset_scvi/stomics0212immune/feature_selection",
    "Epithelial": ROOT / "results/lineage_subsets/stomics_0212_epithelial_subset_scvi/stomics0212epithelial/feature_selection",
}


def load_markers() -> dict[str, set[str]]:
    adata = ad.read_h5ad(FULL_H5AD)
    if "raw_count" in adata.layers:
        adata.X = adata.layers["raw_count"].copy()

    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.tl.rank_genes_groups(
        adata,
        groupby="cell_type",
        method="wilcoxon",
        use_raw=False,
        n_genes=100,
    )

    markers: dict[str, set[str]] = {}
    cell_types = (
        list(adata.obs["cell_type"].cat.categories)
        if hasattr(adata.obs["cell_type"], "cat")
        else sorted(map(str, adata.obs["cell_type"].unique()))
    )
    for cell_type in cell_types:
        names = sc.get.rank_genes_groups_df(adata, group=cell_type)["names"].dropna().head(100)
        markers[cell_type] = set(names.astype(str))
    return markers


def feature_set(dataset_name: str, method_base: str, n_features: str) -> set[str]:
    path = DATASET_PATHS[dataset_name] / method_base / f"n{n_features}" / "seed0" / "selected_features.json"
    obj = json.loads(path.read_text())
    return set(obj["feature_names"])


def main() -> None:
    methods = pd.read_csv(METHODS_TSV, sep="\t").copy()
    methods["n_features"] = methods["Method"].str.replace(r"^.*-N", "", regex=True)
    markers = load_markers()

    lineage_by_label = {
        label: lineage
        for lineage, labels in LINEAGES.items()
        for label in labels
    }
    full_labels = LINEAGES["Immune"] + LINEAGES["Epithelial"]

    rows = []
    for _, row in methods.iterrows():
        method_name = row["Name"]
        method_base = row["MethodBase"]
        n_features = row["n_features"]

        for dataset_name, labels in {
            "Full": full_labels,
            "Immune": LINEAGES["Immune"],
            "Epithelial": LINEAGES["Epithelial"],
        }.items():
            selected = feature_set(dataset_name, method_base, n_features)
            for label in labels:
                marker_set = markers[label]
                prop = len(selected & marker_set) / max(len(marker_set), 1)
                rows.append(
                    {
                        "Method": method_name,
                        "Dataset": dataset_name,
                        "Label": label,
                        "Lineage": lineage_by_label[label],
                        "Value": float(prop),
                    }
                )

    scores = pd.DataFrame(rows)

    full_scores = scores[scores["Dataset"] == "Full"][["Method", "Label", "Value"]].rename(columns={"Value": "FullValue"})
    diffs = (
        scores[scores["Dataset"].isin(["Immune", "Epithelial"])]
        .merge(full_scores, on=["Method", "Label"], how="left")
        .assign(Difference=lambda d: d["Value"] - d["FullValue"])
        [["Method", "Dataset", "Label", "Lineage", "Difference"]]
    )

    OUT_SCORES.parent.mkdir(parents=True, exist_ok=True)
    scores.to_csv(OUT_SCORES, sep="\t", index=False)
    diffs.to_csv(OUT_DIFFS, sep="\t", index=False)


if __name__ == "__main__":
    main()

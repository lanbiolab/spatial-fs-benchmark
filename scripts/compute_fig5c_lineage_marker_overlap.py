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
OUT_TSV = ROOT / "results/fig5c_spatial_lineages/figures/fig5c_lineage_marker_overlap.tsv"

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
    for cell_type in adata.obs["cell_type"].cat.categories if hasattr(adata.obs["cell_type"], "cat") else sorted(adata.obs["cell_type"].unique()):
        names = sc.get.rank_genes_groups_df(adata, group=cell_type)["names"].dropna().head(100)
        markers[cell_type] = set(names.astype(str))
    return markers


def feature_set(dataset_name: str, method_base: str, n_features: str) -> set[str]:
    path = DATASET_PATHS[dataset_name] / method_base / f"n{n_features}" / "seed0" / "selected_features.json"
    obj = json.loads(path.read_text())
    return set(obj["feature_names"])


def main() -> None:
    methods = pd.read_csv(METHODS_TSV, sep="\t")
    methods = methods.copy()
    methods["n_features"] = methods["Method"].str.replace(r"^.*-N", "", regex=True)

    markers = load_markers()

    rows = []
    for _, row in methods.iterrows():
        method_name = row["Name"]
        method_base = row["MethodBase"]
        n_features = row["n_features"]

        for lineage, cell_types in LINEAGES.items():
            for dataset_name in ["Full", "Immune", "Epithelial"]:
                if method_base == "random":
                    mean_prop = np.nan
                    sd_prop = np.nan
                else:
                    selected = feature_set(dataset_name, method_base, n_features)
                    props = []
                    for cell_type in cell_types:
                        marker_set = markers[cell_type]
                        props.append(len(selected & marker_set) / max(len(marker_set), 1))
                    mean_prop = float(np.mean(props))
                    sd_prop = float(np.std(props, ddof=0))

                rows.append(
                    {
                        "Method": method_name,
                        "Lineage": lineage,
                        "Dataset": dataset_name,
                        "MeanProp": mean_prop,
                        "SDProp": sd_prop,
                    }
                )

    out = pd.DataFrame(rows)
    OUT_TSV.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(OUT_TSV, sep="\t", index=False)


if __name__ == "__main__":
    main()

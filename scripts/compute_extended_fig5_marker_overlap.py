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
OUT_SUMMARY = ROOT / "results/extended_data_fig5/figures/extended_fig5_markers_summary.tsv"
OUT_CELLTYPE = ROOT / "results/extended_data_fig5/figures/extended_fig5_celltype_markers.tsv"

LINEAGES = {
    "Endothelial": ["Endothelial Cell"],
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
    cats = adata.obs["cell_type"].cat.categories if hasattr(adata.obs["cell_type"], "cat") else sorted(adata.obs["cell_type"].unique())
    for cell_type in cats:
        names = sc.get.rank_genes_groups_df(adata, group=cell_type)["names"].dropna().head(100)
        markers[str(cell_type)] = set(names.astype(str))
    return markers


def feature_set(dataset_name: str, method_base: str, n_features: str) -> set[str]:
    path = DATASET_PATHS[dataset_name] / method_base / f"n{n_features}" / "seed0" / "selected_features.json"
    obj = json.loads(path.read_text())
    return set(obj["feature_names"])


def main() -> None:
    methods = pd.read_csv(METHODS_TSV, sep="\t").copy()
    methods["n_features"] = methods["Method"].str.replace(r"^.*-N", "", regex=True)
    methods = methods[methods["MethodBase"] != "random"].copy()

    markers = load_markers()

    summary_rows = []
    celltype_rows = []

    for _, row in methods.iterrows():
        method_name = row["Name"]
        method_base = row["MethodBase"]
        n_features = row["n_features"]

        for dataset_name in ["Full", "Immune", "Epithelial"]:
            selected = feature_set(dataset_name, method_base, n_features)

            for lineage, cell_types in LINEAGES.items():
                props = []
                for cell_type in cell_types:
                    marker_set = markers[cell_type]
                    prop = len(selected & marker_set) / max(len(marker_set), 1)
                    props.append(prop)
                    celltype_rows.append(
                        {
                            "Method": method_name,
                            "Dataset": dataset_name,
                            "Lineage": lineage,
                            "Label": cell_type,
                            "Value": float(prop),
                        }
                    )

                summary_rows.append(
                    {
                        "Method": method_name,
                        "Dataset": dataset_name,
                        "Lineage": lineage,
                        "MeanProp": float(np.mean(props)),
                        "SDProp": float(np.std(props, ddof=0)),
                    }
                )

    OUT_SUMMARY.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(summary_rows).to_csv(OUT_SUMMARY, sep="\t", index=False)
    pd.DataFrame(celltype_rows).to_csv(OUT_CELLTYPE, sep="\t", index=False)


if __name__ == "__main__":
    main()

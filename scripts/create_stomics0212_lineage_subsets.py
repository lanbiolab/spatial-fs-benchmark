#!/usr/bin/env python3

import anndata as ad
from pathlib import Path


SOURCE = Path("data/raw/wilcoxon_candidate_0212/stomics_0212_multislice.h5ad")
OUT_DIR = SOURCE.parent

SUBSETS = {
    "epithelial": {
        "outfile": OUT_DIR / "stomics_0212_epithelial_subset.h5ad",
        "labels": ["Keratinocyte", "Merkel Cell", "Stem Cell", "Progenitor Cell"],
    },
    "immune": {
        "outfile": OUT_DIR / "stomics_0212_immune_subset.h5ad",
        "labels": ["Langerhans Cell", "T Cell"],
    },
}


def main() -> None:
    adata = ad.read_h5ad(SOURCE)
    adata.obs["cell_type"] = adata.obs["cell_type"].astype(str)

    for subset_name, cfg in SUBSETS.items():
        labels = set(cfg["labels"])
        sub = adata[adata.obs["cell_type"].isin(labels)].copy()
        sub.write_h5ad(cfg["outfile"])
        print(
            f"{subset_name}: shape={sub.shape}, "
            f"labels={sorted(sub.obs['cell_type'].astype(str).unique().tolist())}"
        )


if __name__ == "__main__":
    main()

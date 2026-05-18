from __future__ import annotations

import argparse
from pathlib import Path

import anndata as ad
import pandas as pd


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Combine 5 STOmics Visium processed h5ad files into one multislice h5ad.")
    parser.add_argument("--root-dir", required=True, help="Directory containing ST*_10x_Visium_processed.h5ad files.")
    parser.add_argument("--output", required=True, help="Output h5ad path.")
    return parser


def main() -> None:
    args = build_parser().parse_args()
    root_dir = Path(args.root_dir).resolve()
    output = Path(args.output).resolve()

    files = sorted(root_dir.glob("ST*_10x_Visium_processed.h5ad"))
    if len(files) == 0:
        raise FileNotFoundError(f"No ST*_10x_Visium_processed.h5ad files found in {root_dir}")

    adatas: list[ad.AnnData] = []
    for path in files:
        sample_id = path.name.split("_")[0]
        adata = ad.read_h5ad(path)
        adata.var_names_make_unique()
        adata.obs_names = pd.Index([f"{sample_id}_{obs}" for obs in adata.obs_names], dtype=str)
        adata.obs["sample_id"] = pd.Categorical([sample_id] * adata.n_obs)
        adata.obs["slice"] = pd.Categorical([sample_id] * adata.n_obs)
        if "clusters" in adata.obs:
            adata.obs["clusters"] = adata.obs["clusters"].astype(str).astype("category")
        adatas.append(adata)

    combined = ad.concat(adatas, join="inner", merge="same")
    combined.var_names_make_unique()
    combined.obs_names_make_unique()
    combined.obs["sample_id"] = combined.obs["sample_id"].astype("category")
    combined.obs["slice"] = combined.obs["slice"].astype("category")

    output.parent.mkdir(parents=True, exist_ok=True)
    combined.write_h5ad(output, compression="gzip")
    print(f"Wrote {output} with shape {combined.shape}")


if __name__ == "__main__":
    main()

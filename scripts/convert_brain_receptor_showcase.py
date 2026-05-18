from __future__ import annotations

import argparse
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Convert BrainReceptorShowcase MERFISH CSVs into a unified h5ad.")
    parser.add_argument("--input-dir", default="data/raw/BrainReceptorShowcase")
    parser.add_argument("--output", default="data/raw/BrainReceptorShowcase/brain_receptor_showcase.h5ad")
    return parser


def load_replicate(replicate_dir: Path, slice_id: str, replicate_id: str) -> ad.AnnData:
    expr_path = replicate_dir / f"cell_by_gene_S{slice_id[-1]}R{replicate_id[-1]}.csv"
    meta_path = replicate_dir / f"cell_metadata_S{slice_id[-1]}R{replicate_id[-1]}.csv"

    expr = pd.read_csv(expr_path)
    meta = pd.read_csv(meta_path)

    cell_id_col = expr.columns[0]
    meta_id_col = meta.columns[0]
    expr[cell_id_col] = expr[cell_id_col].astype(str)
    meta[meta_id_col] = meta[meta_id_col].astype(str)

    merged = meta.merge(expr, left_on=meta_id_col, right_on=cell_id_col, how="inner", validate="one_to_one")
    feature_start = len(meta.columns)
    feature_frame = merged.iloc[:, feature_start:].copy()
    feature_names = feature_frame.columns.astype(str).tolist()

    # Drop blank control probes from the benchmark object.
    keep_columns = [name for name in feature_names if not name.startswith("Blank-")]
    feature_frame = feature_frame[keep_columns]

    obs = merged[meta.columns].copy()
    obs.rename(columns={meta_id_col: "cell_id"}, inplace=True)
    obs["slice"] = pd.Categorical([slice_id] * len(obs))
    obs["replicate"] = pd.Categorical([replicate_id] * len(obs))
    obs["sample_id"] = pd.Categorical([f"{slice_id}_{replicate_id}"] * len(obs))
    obs.index = [f"{slice_id}_{replicate_id}_{cell_id}" for cell_id in obs["cell_id"].astype(str)]

    matrix = sparse.csr_matrix(feature_frame.to_numpy(dtype=np.float32))
    var = pd.DataFrame(index=keep_columns)

    adata = ad.AnnData(X=matrix, obs=obs, var=var)
    adata.obsm["spatial"] = obs[["center_x", "center_y"]].to_numpy(dtype=np.float32)
    return adata


def main() -> None:
    args = build_parser().parse_args()
    input_dir = Path(args.input_dir)

    adatas: list[ad.AnnData] = []
    for slice_dir in sorted(path for path in input_dir.iterdir() if path.is_dir() and path.name.startswith("Slice")):
        for replicate_dir in sorted(path for path in slice_dir.iterdir() if path.is_dir() and path.name.startswith("Replicate")):
            slice_id = slice_dir.name
            replicate_id = replicate_dir.name
            expr_path = replicate_dir / f"cell_by_gene_S{slice_id[-1]}R{replicate_id[-1]}.csv"
            meta_path = replicate_dir / f"cell_metadata_S{slice_id[-1]}R{replicate_id[-1]}.csv"
            if not expr_path.exists() or not meta_path.exists():
                continue
            adatas.append(load_replicate(replicate_dir, slice_id, replicate_id))

    if not adatas:
        raise FileNotFoundError("No complete replicate CSV pairs found under BrainReceptorShowcase.")

    merged = ad.concat(adatas, join="outer", merge="same")
    merged.obs["sample_id"] = merged.obs["sample_id"].astype("category")
    merged.obs["slice"] = merged.obs["slice"].astype("category")
    merged.obs["replicate"] = merged.obs["replicate"].astype("category")
    merged.var_names_make_unique()

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    merged.write_h5ad(output)
    print(f"Wrote {output} with shape {merged.shape}")


if __name__ == "__main__":
    main()

from __future__ import annotations

import argparse
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Convert embryo digital_expression bundles into a unified multislice h5ad.")
    parser.add_argument("--root-dir", required=True, help="Directory containing *.digital_expression.txt.gz files.")
    parser.add_argument("--stage", required=True, help="Stage label to store in obs['stage'], e.g. E8.5 or E9.5.")
    parser.add_argument("--output", required=True, help="Output h5ad path.")
    parser.add_argument("--chunksize", type=int, default=512, help="Gene rows per read_csv chunk.")
    return parser


def _sample_id_from_expr(path: Path) -> str:
    return path.name.replace(".digital_expression.txt.gz", "")


def _coord_path_for_expr(path: Path) -> Path:
    return path.with_name(path.name.replace(".digital_expression.txt.gz", "_matched_bead_locations.txt.gz"))


def _load_expression(path: Path, chunksize: int) -> tuple[sparse.csr_matrix, list[str], list[str]]:
    matrices: list[sparse.csr_matrix] = []
    genes: list[str] = []
    barcodes: list[str] | None = None
    reader = pd.read_csv(
        path,
        sep="\t",
        compression="gzip",
        index_col=0,
        chunksize=chunksize,
    )
    for chunk in reader:
        if barcodes is None:
            barcodes = chunk.columns.astype(str).tolist()
        genes.extend(chunk.index.astype(str).tolist())
        values = chunk.to_numpy(dtype=np.int32, copy=False)
        matrices.append(sparse.csr_matrix(values.T))
    if barcodes is None:
        raise ValueError(f"No expression rows found in {path}")
    matrix = sparse.hstack(matrices, format="csr", dtype=np.int32)
    return matrix, barcodes, genes


def _load_coordinates(path: Path, expected_obs: int, sample_id: str, stage: str, barcodes: list[str]) -> pd.DataFrame:
    coords = pd.read_csv(
        path,
        sep="\t",
        compression="gzip",
        header=None,
        names=["in_tissue", "x", "y"],
    )
    if len(coords) != expected_obs:
        raise ValueError(
            f"Coordinate rows ({len(coords)}) do not match expression columns ({expected_obs}) for sample {sample_id}"
        )
    obs_index = pd.Index([f"{sample_id}:{barcode}" for barcode in barcodes], name="spot_id")
    obs = pd.DataFrame(index=obs_index)
    obs["barcode"] = pd.Categorical(barcodes)
    obs["sample_id"] = pd.Categorical([sample_id] * expected_obs)
    obs["slice"] = pd.Categorical([sample_id] * expected_obs)
    obs["stage"] = pd.Categorical([stage] * expected_obs)
    obs["in_tissue"] = coords["in_tissue"].to_numpy(dtype=np.int8)
    obs["x"] = coords["x"].to_numpy(dtype=np.float32)
    obs["y"] = coords["y"].to_numpy(dtype=np.float32)
    return obs


def _load_sample(expr_path: Path, stage: str, chunksize: int) -> ad.AnnData:
    sample_id = _sample_id_from_expr(expr_path)
    coord_path = _coord_path_for_expr(expr_path)
    if not coord_path.exists():
        raise FileNotFoundError(f"Missing coordinate file for {expr_path.name}: {coord_path}")
    matrix, barcodes, genes = _load_expression(expr_path, chunksize)
    obs = _load_coordinates(coord_path, matrix.shape[0], sample_id, stage, barcodes)
    adata = ad.AnnData(
        X=matrix,
        obs=obs,
        var=pd.DataFrame(index=pd.Index(genes, name="gene")),
    )
    adata.var_names_make_unique()
    adata.obs_names_make_unique()
    adata.obsm["spatial"] = obs[["x", "y"]].to_numpy(dtype=np.float32)
    return adata


def main() -> None:
    args = build_parser().parse_args()
    root_dir = Path(args.root_dir).resolve()
    output = Path(args.output).resolve()
    expr_paths = sorted(root_dir.glob("*.digital_expression.txt.gz"))
    if not expr_paths:
        raise FileNotFoundError(f"No *.digital_expression.txt.gz files found in {root_dir}")

    adatas = [_load_sample(expr_path, args.stage, args.chunksize) for expr_path in expr_paths]
    combined = ad.concat(adatas, join="inner", merge="same")
    combined.obs["sample_id"] = combined.obs["sample_id"].astype("category")
    combined.obs["slice"] = combined.obs["slice"].astype("category")
    combined.obs["stage"] = combined.obs["stage"].astype("category")
    output.parent.mkdir(parents=True, exist_ok=True)
    combined.write_h5ad(output, compression="gzip")
    print(f"Wrote {output} with shape {combined.shape}")


if __name__ == "__main__":
    main()

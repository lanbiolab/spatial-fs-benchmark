from __future__ import annotations

import argparse
import csv
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Convert sagittal CSV bundles into a unified h5ad.")
    parser.add_argument("--root-dir", required=True, help="Directory containing sagittal1/2/3 folders.")
    parser.add_argument(
        "--samples",
        nargs="*",
        default=None,
        help="Optional sample directory names to convert, e.g. sagittal1 sagittal2.",
    )
    parser.add_argument("--output", required=True, help="Output h5ad path.")
    return parser


def _count_gene_rows(path: Path) -> int:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return max(sum(1 for _ in handle) - 1, 0)


def _load_expression(path: Path) -> tuple[np.ndarray, list[str], list[str]]:
    n_genes = _count_gene_rows(path)
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle)
        header = next(reader)
        cell_ids = header[1:]
        matrix = np.empty((len(cell_ids), n_genes), dtype=np.float32)
        genes: list[str] = []
        for gene_idx, row in enumerate(reader):
            genes.append(row[0])
            matrix[:, gene_idx] = np.asarray(row[1:], dtype=np.float32)
    return matrix, cell_ids, genes


def _load_sample(sample_dir: Path) -> ad.AnnData:
    expr_path = next(sample_dir.glob("*processed_expression_pd.csv"))
    spatial_path = next(sample_dir.glob("*_spatial.csv"))

    matrix, cell_ids, genes = _load_expression(expr_path)
    spatial = pd.read_csv(spatial_path, skiprows=[1]).set_index("NAME")
    obs = spatial.reindex(cell_ids).copy()
    if obs.isnull().any().any():
        missing = obs.index[obs.isnull().any(axis=1)][:5].tolist()
        raise ValueError(f"Missing spatial metadata for cells: {missing}")
    obs["sample_id"] = sample_dir.name
    obs["slice"] = sample_dir.name
    obs.index = cell_ids

    adata = ad.AnnData(X=matrix, obs=obs, var=pd.DataFrame(index=pd.Index(genes, name="gene")))
    adata.var_names_make_unique()
    adata.obs_names_make_unique()
    adata.obsm["spatial"] = obs[["X", "Y"]].to_numpy(dtype=np.float32)
    return adata


def main() -> None:
    args = build_parser().parse_args()
    root_dir = Path(args.root_dir).resolve()
    output = Path(args.output).resolve()

    sample_dirs = sorted(path for path in root_dir.iterdir() if path.is_dir())
    if args.samples:
        keep = set(args.samples)
        sample_dirs = [path for path in sample_dirs if path.name in keep]
    if not sample_dirs:
        raise FileNotFoundError(f"No sagittal sample directories found in {root_dir}")

    adatas = [_load_sample(sample_dir) for sample_dir in sample_dirs]
    reference_genes = adatas[0].var_names
    for adata in adatas[1:]:
        if not reference_genes.equals(adata.var_names):
            raise ValueError("Gene orders differ across sagittal samples; current converter expects identical genes.")

    combined = ad.concat(adatas, join="inner", merge="same")
    combined.obs["sample_id"] = combined.obs["sample_id"].astype("category")
    combined.obs["slice"] = combined.obs["slice"].astype("category")
    output.parent.mkdir(parents=True, exist_ok=True)
    combined.write_h5ad(output, compression="gzip")
    print(f"Wrote {output} with shape {combined.shape}")


if __name__ == "__main__":
    main()

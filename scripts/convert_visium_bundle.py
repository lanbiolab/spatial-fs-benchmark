from __future__ import annotations

import argparse
import tarfile
import tempfile
from pathlib import Path

import anndata as ad
import pandas as pd
import scanpy as sc


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Convert a 10x Visium bundle into h5ad.")
    parser.add_argument("--dataset-dir", required=True, help="Directory containing the 10x files.")
    parser.add_argument("--matrix-h5", required=True, help="filtered_feature_bc_matrix.h5 filename.")
    parser.add_argument("--spatial-tar", required=True, help="spatial.tar.gz filename.")
    parser.add_argument("--analysis-tar", default=None, help="analysis.tar.gz filename.")
    parser.add_argument("--sample-id", required=True, help="Sample identifier to store in obs.")
    parser.add_argument("--output", required=True, help="Output h5ad path.")
    return parser


def _merge_clustering_outputs(adata: ad.AnnData, analysis_tar: Path) -> None:
    with tarfile.open(analysis_tar, "r:gz") as archive:
        cluster_members = [m for m in archive.getmembers() if m.name.endswith("/clusters.csv")]
        for member in cluster_members:
            handle = archive.extractfile(member)
            if handle is None:
                continue
            frame = pd.read_csv(handle)
            if "Barcode" not in frame.columns or "Cluster" not in frame.columns:
                continue
            column = member.name.split("/")[-2].replace("gene_expression_", "")
            series = frame.set_index("Barcode")["Cluster"].astype(str)
            adata.obs[column] = series.reindex(adata.obs_names).astype("category")


def main() -> None:
    args = build_parser().parse_args()
    dataset_dir = Path(args.dataset_dir).resolve()
    matrix_h5 = dataset_dir / args.matrix_h5
    spatial_tar = dataset_dir / args.spatial_tar
    analysis_tar = dataset_dir / args.analysis_tar if args.analysis_tar is not None else None
    output = Path(args.output).resolve()

    with tempfile.TemporaryDirectory(prefix="visium_convert_") as temp_dir:
        temp_path = Path(temp_dir)
        spatial_dir = temp_path / "spatial"
        with tarfile.open(spatial_tar, "r:gz") as archive:
            archive.extractall(temp_path)
        adata = sc.read_visium(path=temp_path, count_file=str(matrix_h5))

    adata.obs["sample_id"] = pd.Categorical([args.sample_id] * adata.n_obs)
    adata.obs["in_tissue"] = adata.obs["in_tissue"].astype(str).astype("category")
    adata.var_names_make_unique()
    if analysis_tar is not None and analysis_tar.exists():
        _merge_clustering_outputs(adata, analysis_tar)

    output.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(output)
    print(f"Wrote {output}")


if __name__ == "__main__":
    main()

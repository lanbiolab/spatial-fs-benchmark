from __future__ import annotations

import argparse
import tarfile
import tempfile
from pathlib import Path

import anndata as ad
import pandas as pd
import scanpy as sc


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Convert two serial Visium sections into one h5ad.")
    parser.add_argument("--root-dir", required=True, help="Directory containing the two section folders.")
    parser.add_argument("--output", required=True, help="Output h5ad path.")
    return parser


def _merge_clustering_outputs(adata: ad.AnnData, analysis_tar: Path) -> None:
    with tarfile.open(analysis_tar, "r:gz") as archive:
        cluster_members = [member for member in archive.getmembers() if member.name.endswith("/clusters.csv")]
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


def _load_section(
    dataset_dir: Path,
    matrix_h5: str,
    spatial_tar: str,
    analysis_tar: str,
    sample_id: str,
) -> ad.AnnData:
    with tempfile.TemporaryDirectory(prefix="visium_section_") as temp_dir:
        temp_path = Path(temp_dir)
        with tarfile.open(dataset_dir / spatial_tar, "r:gz") as archive:
            archive.extractall(temp_path)
        adata = sc.read_visium(path=temp_path, count_file=str(dataset_dir / matrix_h5))
    adata.var_names_make_unique()
    adata.obs["sample_id"] = pd.Categorical([sample_id] * adata.n_obs)
    adata.obs["in_tissue"] = adata.obs["in_tissue"].astype(str).astype("category")
    _merge_clustering_outputs(adata, dataset_dir / analysis_tar)
    return adata


def main() -> None:
    args = build_parser().parse_args()
    root_dir = Path(args.root_dir).resolve()
    output = Path(args.output).resolve()

    samples = [
        {
            "dataset_dir": root_dir / "Mouse Brain Serial Section 1 (Sagittal-Anterior)",
            "matrix_h5": "V1_Mouse_Brain_Sagittal_Anterior_filtered_feature_bc_matrix.h5",
            "spatial_tar": "V1_Mouse_Brain_Sagittal_Anterior_spatial.tar.gz",
            "analysis_tar": "V1_Mouse_Brain_Sagittal_Anterior_analysis.tar.gz",
            "sample_id": "section1",
        },
        {
            "dataset_dir": root_dir / "Mouse Brain Serial Section 2 (Sagittal-Anterior)",
            "matrix_h5": "V1_Mouse_Brain_Sagittal_Anterior_Section_2_filtered_feature_bc_matrix.h5",
            "spatial_tar": "V1_Mouse_Brain_Sagittal_Anterior_Section_2_spatial.tar.gz",
            "analysis_tar": "V1_Mouse_Brain_Sagittal_Anterior_Section_2_analysis.tar.gz",
            "sample_id": "section2",
        },
    ]

    adatas = [
        _load_section(
            dataset_dir=sample["dataset_dir"],
            matrix_h5=sample["matrix_h5"],
            spatial_tar=sample["spatial_tar"],
            analysis_tar=sample["analysis_tar"],
            sample_id=sample["sample_id"],
        )
        for sample in samples
    ]
    combined = ad.concat(adatas, join="inner", merge="same")
    combined.var_names_make_unique()
    combined.obs_names_make_unique()
    combined.obs["sample_id"] = combined.obs["sample_id"].astype("category")
    output.parent.mkdir(parents=True, exist_ok=True)
    combined.write_h5ad(output, compression="gzip")
    print(f"Wrote {output} with shape {combined.shape}")


if __name__ == "__main__":
    main()

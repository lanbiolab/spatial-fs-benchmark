from __future__ import annotations

import argparse
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import io, sparse


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Convert exported DLPFC matrix files into h5ad.")
    parser.add_argument("--input-dir", default="data/raw/dlpfc")
    parser.add_argument("--output", default="data/raw/dlpfc/dlpfc_multislice.h5ad")
    return parser


def main() -> None:
    args = build_parser().parse_args()
    input_dir = Path(args.input_dir)

    counts = io.mmread(input_dir / "counts.mtx")
    if not sparse.issparse(counts):
        counts = sparse.csr_matrix(counts)
    else:
        counts = counts.tocsr()

    obs = pd.read_csv(input_dir / "obs.csv", index_col=0)
    var = pd.read_csv(input_dir / "var.csv", index_col=0)
    spatial = pd.read_csv(input_dir / "spatial.csv", index_col=0)

    obs_names = [line.strip() for line in (input_dir / "obs_names.txt").read_text(encoding="utf-8").splitlines() if line.strip()]
    var_names = [line.strip() for line in (input_dir / "var_names.txt").read_text(encoding="utf-8").splitlines() if line.strip()]

    obs.index = obs_names
    var.index = var_names
    spatial.index = obs_names

    if "key" in obs.columns:
        unique_obs_names = obs["key"].astype(str).tolist()
        obs.index = unique_obs_names
        spatial.index = unique_obs_names

    if counts.shape != (len(obs_names), len(var_names)):
        counts = counts.T
    if counts.shape != (obs.shape[0], var.shape[0]):
        raise ValueError(f"Matrix shape {counts.shape} does not match obs/var lengths")

    adata = ad.AnnData(X=counts, obs=obs, var=var)
    adata.obsm["spatial"] = spatial.to_numpy(dtype=np.float32)
    adata.uns["selected_sample_ids"] = [
        line.strip()
        for line in (input_dir / "selected_sample_ids.txt").read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(output)
    print(f"Wrote {output}")


if __name__ == "__main__":
    main()

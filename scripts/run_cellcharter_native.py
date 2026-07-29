from __future__ import annotations

import argparse
import json
from pathlib import Path

import anndata as ad
import cellcharter as cc
import numpy as np
import scvi
import squidpy as sq


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run native CellCharter aggregation in the dedicated environment.")
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--metadata-output", required=True, type=Path)
    parser.add_argument("--slice-key", required=True)
    parser.add_argument("--coord-key", required=True)
    parser.add_argument("--n-latent", required=True, type=int)
    parser.add_argument("--nhood-layers", required=True, type=int)
    parser.add_argument("--max-epochs", required=True, type=int)
    parser.add_argument("--batch-size", required=True, type=int)
    parser.add_argument("--seed", required=True, type=int)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    adata = ad.read_h5ad(args.input)
    if "counts" in adata.layers:
        adata.X = adata.layers["counts"].copy()
        adata.uns.pop("log1p", None)
    matrix = adata.X
    matrix_min = float(matrix.min()) if hasattr(matrix, "min") else float(np.min(np.asarray(matrix)))
    if matrix_min < 0:
        raise ValueError("CellCharter requires non-negative count-like input in X/counts.")
    values = matrix.data if hasattr(matrix, "data") else np.asarray(matrix).ravel()
    if values.size and not np.all(np.abs(values - np.rint(values)) <= 1e-6):
        raise ValueError("CellCharter/scVI requires integer counts in the counts layer.")

    adata.obs[args.slice_key] = adata.obs[args.slice_key].astype("category")
    scvi.settings.seed = args.seed
    scvi.model.SCVI.setup_anndata(adata, layer="counts", batch_key=args.slice_key)
    model = scvi.model.SCVI(adata, n_latent=min(args.n_latent, max(2, adata.n_vars)))
    model.train(
        max_epochs=args.max_epochs,
        batch_size=args.batch_size,
        check_val_every_n_epoch=None,
        enable_progress_bar=False,
    )
    adata.obsm["X_scvi"] = model.get_latent_representation(adata).astype(np.float32)

    sq.gr.spatial_neighbors(
        adata,
        spatial_key=args.coord_key,
        library_key=args.slice_key,
        coord_type="generic",
        delaunay=True,
    )
    cc.gr.aggregate_neighbors(
        adata,
        n_layers=args.nhood_layers,
        use_rep="X_scvi",
        sample_key=args.slice_key,
        out_key="X_cellcharter",
    )

    embedding = np.asarray(adata.obsm["X_cellcharter"], dtype=np.float32)
    np.save(args.output, embedding)
    args.metadata_output.write_text(
        json.dumps(
            {
                "n_latent": args.n_latent,
                "nhood_layers": args.nhood_layers,
                "max_epochs": args.max_epochs,
                "batch_size": args.batch_size,
                "output_dim": int(embedding.shape[1]),
                "native_pipeline": "scvi + squidpy.gr.spatial_neighbors + cellcharter.gr.aggregate_neighbors",
            },
            indent=2,
        ),
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()

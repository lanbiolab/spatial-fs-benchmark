from __future__ import annotations

import argparse
from pathlib import Path

import anndata as ad
import numpy as np


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Crop a central spatial region from each sagittal slice.")
    parser.add_argument("--input", required=True, help="Input sagittal multislice h5ad.")
    parser.add_argument("--output", required=True, help="Output cropped h5ad.")
    parser.add_argument("--lower-quantile", type=float, default=0.35, help="Lower coordinate quantile per slice.")
    parser.add_argument("--upper-quantile", type=float, default=0.65, help="Upper coordinate quantile per slice.")
    return parser


def main() -> None:
    args = build_parser().parse_args()
    input_path = Path(args.input).resolve()
    output_path = Path(args.output).resolve()

    adata = ad.read_h5ad(input_path)
    spatial = np.asarray(adata.obsm["spatial"], dtype=np.float32)
    sample_ids = adata.obs["sample_id"].astype(str).to_numpy()

    keep = np.zeros(adata.n_obs, dtype=bool)
    summary: list[tuple[str, int, int]] = []
    for sample_id in np.unique(sample_ids):
        mask = sample_ids == sample_id
        coords = spatial[mask]
        x0, x1 = np.quantile(coords[:, 0], [args.lower_quantile, args.upper_quantile])
        y0, y1 = np.quantile(coords[:, 1], [args.lower_quantile, args.upper_quantile])
        sample_keep = mask.copy()
        sample_keep[mask] = (
            (coords[:, 0] >= x0)
            & (coords[:, 0] <= x1)
            & (coords[:, 1] >= y0)
            & (coords[:, 1] <= y1)
        )
        keep |= sample_keep
        summary.append((sample_id, int(mask.sum()), int(sample_keep.sum())))

    cropped = adata[keep].copy()
    cropped.uns["region_crop"] = {
        "lower_quantile": float(args.lower_quantile),
        "upper_quantile": float(args.upper_quantile),
        "summary": [
            {"sample_id": sample_id, "n_before": before, "n_after": after}
            for sample_id, before, after in summary
        ],
    }
    output_path.parent.mkdir(parents=True, exist_ok=True)
    cropped.write_h5ad(output_path, compression="gzip")
    print(f"Wrote {output_path} with shape {cropped.shape}")
    for sample_id, before, after in summary:
        print(f"{sample_id}: {after}/{before}")


if __name__ == "__main__":
    main()

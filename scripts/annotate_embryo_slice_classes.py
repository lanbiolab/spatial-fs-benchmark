from __future__ import annotations

import argparse
from pathlib import Path

import anndata as ad
import pandas as pd


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Annotate embryo multislice h5ad with serial slice classes.")
    parser.add_argument("--input", required=True, help="Input embryo h5ad path.")
    parser.add_argument("--output", required=True, help="Output annotated h5ad path.")
    parser.add_argument(
        "--groups",
        default="3,3",
        help="Comma-separated serial group sizes over lexicographically sorted sample_id values, e.g. 3,3.",
    )
    parser.add_argument(
        "--labels",
        default="serial_group_1,serial_group_2",
        help="Comma-separated labels for the serial groups.",
    )
    return parser


def main() -> None:
    args = build_parser().parse_args()
    input_path = Path(args.input).resolve()
    output_path = Path(args.output).resolve()

    group_sizes = [int(item) for item in args.groups.split(",") if item.strip()]
    group_labels = [item.strip() for item in args.labels.split(",") if item.strip()]
    if len(group_sizes) != len(group_labels):
        raise ValueError("groups and labels must have the same length")

    adata = ad.read_h5ad(input_path)
    sample_ids = sorted(adata.obs["sample_id"].astype(str).unique().tolist())
    if sum(group_sizes) != len(sample_ids):
        raise ValueError(
            f"Group sizes sum to {sum(group_sizes)}, but dataset has {len(sample_ids)} unique sample_id values"
        )

    assignments: dict[str, str] = {}
    cursor = 0
    for size, label in zip(group_sizes, group_labels, strict=True):
        for sample_id in sample_ids[cursor : cursor + size]:
            assignments[sample_id] = label
        cursor += size

    adata.obs["slice_class"] = pd.Categorical(adata.obs["sample_id"].astype(str).map(assignments))
    # Keep metadata h5ad-serializable: lists must stay homogeneous.
    adata.uns["slice_class_definition"] = {
        "sample_order": sample_ids,
        "group_labels": group_labels,
        "group_sizes": group_sizes,
        "rule": "sorted sample_id split into contiguous serial groups",
    }
    output_path.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(output_path, compression="gzip")
    print(f"Wrote {output_path} with slice_class labels: {group_labels}")


if __name__ == "__main__":
    main()

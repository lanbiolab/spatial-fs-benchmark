from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


KEY = [
    "dataset",
    "fs_method",
    "n_features",
    "integration_method",
    "task",
    "metric_name",
    "random_seed",
]


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Filter the audited merged result table to the two primary integrators."
    )
    parser.add_argument(
        "--input",
        default="results/spatial_svg_rebuild_v1/merged_results.csv",
        type=Path,
    )
    parser.add_argument(
        "--output",
        default="results/spatial_svg_rebuild_v1/merged_primary_results.csv",
        type=Path,
    )
    parser.add_argument("--expected-rows", type=int, default=43874)
    args = parser.parse_args()

    frame = pd.read_csv(args.input)
    primary = frame.loc[frame["integration_method"].isin(["scvi", "cellcharter"])].copy()
    duplicates = int(primary.duplicated(KEY).sum())
    if duplicates:
        raise RuntimeError(f"Primary merged table contains {duplicates} duplicate metric keys")
    if len(primary) != args.expected_rows:
        raise RuntimeError(f"Expected {args.expected_rows} primary metric rows, observed {len(primary)}")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    primary.to_csv(args.output, index=False)
    print(f"primary_rows={len(primary)} duplicates={duplicates} output={args.output}")


if __name__ == "__main__":
    main()

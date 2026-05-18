from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Summarize benchmark result tables.")
    parser.add_argument("--results-dir", required=True, help="Benchmark output directory containing results.csv.")
    parser.add_argument("--output", required=True, help="Path to the summary CSV to write.")
    return parser


def main() -> None:
    args = build_parser().parse_args()
    results_path = Path(args.results_dir) / "results.csv"
    results = pd.read_csv(results_path)
    summary = (
        results.groupby(
            ["dataset", "fs_method", "n_features", "integration_method", "task", "metric_name"],
            observed=True,
        )["metric_value"]
        .agg(["mean", "std"])
        .reset_index()
        .rename(columns={"mean": "metric_mean", "std": "metric_std"})
    )
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    summary.to_csv(output_path, index=False)
    print(f"Wrote summary to {output_path}")


if __name__ == "__main__":
    main()

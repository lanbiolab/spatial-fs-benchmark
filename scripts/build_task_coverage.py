from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build task coverage tables from benchmark results.")
    parser.add_argument("--results", required=True, help="Path to combined benchmark results.csv")
    parser.add_argument("--output-dir", required=True, help="Directory to write coverage tables")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    results = pd.read_csv(args.results)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    results["has_value"] = results["metric_value"].notna()

    dataset_task = (
        results.groupby(["dataset", "task", "metric_name"], dropna=False)["has_value"]
        .agg(available="sum", total="count")
        .reset_index()
    )
    dataset_task["coverage"] = dataset_task["available"] / dataset_task["total"]
    dataset_task.to_csv(output_dir / "dataset_task_metric_coverage.csv", index=False)

    dataset_task_summary = (
        results.groupby(["dataset", "task"], dropna=False)["has_value"]
        .agg(available="sum", total="count")
        .reset_index()
    )
    dataset_task_summary["coverage"] = dataset_task_summary["available"] / dataset_task_summary["total"]
    dataset_task_summary.to_csv(output_dir / "dataset_task_coverage.csv", index=False)

    method_task_summary = (
        results.groupby(["fs_method", "task"], dropna=False)["has_value"]
        .agg(available="sum", total="count")
        .reset_index()
    )
    method_task_summary["coverage"] = method_task_summary["available"] / method_task_summary["total"]
    method_task_summary.to_csv(output_dir / "method_task_coverage.csv", index=False)


if __name__ == "__main__":
    main()

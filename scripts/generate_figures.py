from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from spatial_fs_benchmark.plotting.summary_plots import save_performance_heatmap, save_ranking_plot


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Generate summary figures from benchmark results.")
    parser.add_argument("--results", required=True, help="Path to a benchmark results CSV or summary CSV.")
    parser.add_argument("--output-dir", required=True, help="Directory where plots will be written.")
    return parser


def main() -> None:
    args = build_parser().parse_args()
    results = pd.read_csv(args.results)
    if "metric_value" not in results.columns and "metric_mean" in results.columns:
        results = results.rename(columns={"metric_mean": "metric_value"})
    output_dir = Path(args.output_dir)
    heatmap_path = save_performance_heatmap(results, output_dir)
    ranking_path = save_ranking_plot(results, output_dir)
    print(f"Wrote {heatmap_path}")
    print(f"Wrote {ranking_path}")


if __name__ == "__main__":
    main()

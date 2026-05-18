from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

from spatial_fs_benchmark.metrics.ranking import add_within_task_ranks


def save_performance_heatmap(results: pd.DataFrame, output_dir: str | Path) -> Path:
    path = Path(output_dir)
    path.mkdir(parents=True, exist_ok=True)
    results = results.copy()
    results["column_key"] = results["integration_method"] + " | " + results["task"] + " | " + results["metric_name"]
    pivot = results.pivot_table(
        index="fs_method",
        columns="column_key",
        values="metric_value",
        aggfunc="mean",
    )
    figure = plt.figure(figsize=(max(10, pivot.shape[1] * 0.8), max(4, pivot.shape[0] * 0.6)))
    sns.heatmap(pivot, cmap="viridis", annot=False)
    plt.tight_layout()
    output_path = path / "performance_heatmap.png"
    figure.savefig(output_path, dpi=200)
    plt.close(figure)
    return output_path


def save_ranking_plot(results: pd.DataFrame, output_dir: str | Path) -> Path:
    path = Path(output_dir)
    path.mkdir(parents=True, exist_ok=True)
    ranked = add_within_task_ranks(results)
    summary = ranked.groupby("fs_method", observed=True)["rank"].mean().sort_values()
    figure = plt.figure(figsize=(8, 4))
    summary.plot(kind="bar", color="#2c7fb8")
    plt.ylabel("Average Rank")
    plt.tight_layout()
    output_path = path / "feature_selection_ranking.png"
    figure.savefig(output_path, dpi=200)
    plt.close(figure)
    return output_path

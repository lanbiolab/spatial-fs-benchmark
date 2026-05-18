from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd


INVALID_METHOD_NOTES = {
    "seurat_mvp": "Marked NA on SagittalAtlas: selector degenerates on processed values with negatives; selected scores are all NaN.",
    "seurat_disp": "Marked NA on SagittalAtlas: selector degenerates on processed values with negatives; selected scores are all NaN.",
    "statistic_mean": "Marked NA on SagittalAtlas: selector applies log1p to processed values with negatives; selected scores are all NaN.",
    "statistic_variance": "Marked NA on SagittalAtlas: selector applies log1p to processed values with negatives; selected scores are all NaN.",
}


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Curate SagittalAtlas benchmark results.")
    parser.add_argument(
        "--results-dir",
        default="results/sagittal_atlas_fs",
        help="Benchmark output directory containing results.csv and feature_selection artifacts.",
    )
    return parser


def load_effective_features(results_dir: Path) -> pd.DataFrame:
    base = results_dir / "sagittalatlas" / "feature_selection"
    rows: list[dict[str, object]] = []
    for path in base.rglob("selected_features.json"):
        seed_part = path.parent.name
        n_part = path.parent.parent.name
        method = path.parent.parent.parent.name
        payload = json.loads(path.read_text(encoding="utf-8"))
        metadata = payload.get("metadata", {})
        rows.append(
            {
                "dataset": "SagittalAtlas",
                "fs_method": method,
                "n_features": int(n_part.removeprefix("n")),
                "random_seed": int(seed_part.removeprefix("seed")),
                "effective_n_features": int(metadata.get("effective_n_features", len(payload.get("feature_names", [])))),
                "feature_name_hash": metadata.get("feature_name_hash", ""),
            }
        )
    return pd.DataFrame(rows)


def append_note(existing: pd.Series, note: str) -> pd.Series:
    existing = existing.fillna("").astype(str)
    return existing.where(existing == "", existing + " | ") + note


def main() -> None:
    args = build_parser().parse_args()
    results_dir = Path(args.results_dir)
    results_path = results_dir / "results.csv"
    raw_backup_path = results_dir / "results_raw.csv"
    curated_path = results_dir / "results_curated.csv"

    results = pd.read_csv(results_path)
    if "notes" not in results.columns:
        results["notes"] = ""
    results["notes"] = results["notes"].fillna("").astype(str)
    if "effective_n_features" in results.columns:
        results = results.drop(columns=["effective_n_features"])

    feature_meta = load_effective_features(results_dir)
    curated = results.merge(
        feature_meta[["dataset", "fs_method", "n_features", "random_seed", "effective_n_features"]],
        on=["dataset", "fs_method", "n_features", "random_seed"],
        how="left",
    )
    curated["effective_n_features"] = curated["effective_n_features"].fillna(curated["n_features"]).astype(int)

    for method, note in INVALID_METHOD_NOTES.items():
        mask = curated["fs_method"] == method
        if not mask.any():
            continue
        curated.loc[mask, "metric_value"] = np.nan
        curated.loc[mask, "notes"] = append_note(curated.loc[mask, "notes"].astype(str), note)

    if not raw_backup_path.exists():
        results.to_csv(raw_backup_path, index=False)

    curated.to_csv(curated_path, index=False)
    curated.to_csv(results_path, index=False)
    print(f"Wrote curated Sagittal results to {curated_path}")
    print(f"Updated canonical results file at {results_path}")


if __name__ == "__main__":
    main()

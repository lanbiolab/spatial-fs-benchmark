from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from build_frozen_scores import (  # noqa: E402
    aggregate_setting_scores,
    build_seed_task_scores,
    metric_directions,
    scale_metrics,
)


FROZEN = ROOT / "results/spatial_svg_rebuild_v1/frozen_scores"
HYBRID = ROOT / "results/hybrid_controls_v1"
MORAN_HELDOUT = ROOT / "results/validation_v1/heldout_slice"
NNSVG_HELDOUT = ROOT / "results/validation_v1/heldout_slice_nnsvg_reference"
OUTPUT = ROOT / "results/submission_validations_v1"
COMPETITIVE_GROUPS = {"expression_driven", "spatially_informed"}


def hybrid_scores() -> tuple[
    pd.DataFrame,
    pd.DataFrame,
    pd.DataFrame,
    pd.DataFrame,
    pd.DataFrame,
    pd.DataFrame,
    pd.DataFrame,
]:
    results = pd.read_csv(HYBRID / "merged_results.csv")
    directions = metric_directions()
    results = results.merge(
        directions[["Metric", "DirectionMultiplier"]],
        left_on="metric_name",
        right_on="Metric",
        how="left",
        validate="many_to_one",
    )
    if results["DirectionMultiplier"].isna().any():
        raise RuntimeError("Hybrid results contain unknown metric directions.")
    results["OrientedValue"] = results["metric_value"] * results["DirectionMultiplier"]
    ranges = pd.read_csv(FROZEN / "frozen_metric_ranges.tsv", sep="\t")
    scaled = scale_metrics(results, ranges, "OrientedValue")
    if scaled.loc[scaled["metric_value"].notna(), "FrozenScaleDenominator"].isna().any():
        raise RuntimeError("A finite hybrid metric lacks a frozen scaling range.")
    seed_scores = build_seed_task_scores(scaled)
    seed_scores["MethodGroup"] = "hybrid"
    setting_scores = aggregate_setting_scores(seed_scores)
    setting_scores["MethodGroup"] = "hybrid"

    existing = pd.read_csv(FROZEN / "representative_task_scores.tsv", sep="\t")
    existing = existing[existing["MethodGroup"].isin(COMPETITIVE_GROUPS)].copy()
    common = [column for column in existing.columns if column in setting_scores.columns]
    combined = pd.concat([existing[common], setting_scores[common]], ignore_index=True)
    for score in ["Integration", "Clustering", "Alignment", "CoreOverall"]:
        combined[f"{score}Rank"] = combined.groupby(
            ["dataset", "integration_method"], observed=True
        )[f"{score}Mean"].rank(method="average", ascending=False)
    global_summary = (
        combined.groupby(["fs_method", "integration_method", "MethodGroup"], observed=True)
        .agg(
            MeanCoreOverall=("CoreOverallMean", "mean"),
            SDCoreOverall=("CoreOverallMean", "std"),
            MeanDatasetRank=("CoreOverallRank", "mean"),
            NDatasets=("dataset", "nunique"),
        )
        .reset_index()
    )
    global_summary["GlobalRankWithHybrids"] = global_summary.groupby(
        "integration_method", observed=True
    )["MeanDatasetRank"].rank(method="average")

    effective_rows = []
    for path in HYBRID.rglob("selected_features.meta.json"):
        payload = json.loads(path.read_text())
        effective_rows.append(
            {
                "dataset": payload["dataset"]["dataset_name"],
                "fs_method": payload["selector_name"],
                "requested_n_features": payload["n_features"],
                "effective_n_features": payload["effective_n_features"],
            }
        )
    effective = pd.DataFrame(effective_rows).drop_duplicates()
    global_summary = global_summary.merge(
        effective.groupby("fs_method", observed=True)["effective_n_features"]
        .agg(["min", "median", "max"])
        .reset_index()
        .rename(
            columns={
                "min": "EffectiveFeaturesMin",
                "median": "EffectiveFeaturesMedian",
                "max": "EffectiveFeaturesMax",
            }
        ),
        on="fs_method",
        how="left",
    )
    comparison_methods = {
        "hybrid_hvg_svg_intersection",
        "hybrid_hvg_svg_union",
        "scanpy_seurat_v3_batch",
        "morans_i",
        "triku",
        "somde",
        "seurat_vst",
    }
    comparison = combined[combined["fs_method"].isin(comparison_methods)].copy()
    comparison = comparison.merge(
        effective[["dataset", "fs_method", "effective_n_features"]],
        on=["dataset", "fs_method"],
        how="left",
    )
    comparison["effective_n_features"] = comparison["effective_n_features"].fillna(2000)
    parents = comparison[comparison["fs_method"].isin({"scanpy_seurat_v3_batch", "morans_i"})]
    best_parent = (
        parents.groupby(["dataset", "integration_method"], observed=True)["CoreOverallMean"]
        .max()
        .rename("BestParentCoreOverall")
        .reset_index()
    )
    hybrid_differences = comparison[comparison["MethodGroup"].eq("hybrid")].merge(
        best_parent,
        on=["dataset", "integration_method"],
        how="left",
        validate="many_to_one",
    )
    hybrid_differences["DifferenceFromBestParent"] = (
        hybrid_differences["CoreOverallMean"] - hybrid_differences["BestParentCoreOverall"]
    )
    return scaled, seed_scores, setting_scores, global_summary, effective, comparison, hybrid_differences


def load_heldout(root: Path, reference: str) -> pd.DataFrame:
    frames = [pd.read_csv(path, sep="\t") for path in sorted(root.glob("*/metrics.tsv"))]
    frame = pd.concat(frames, ignore_index=True)
    frame["Reference"] = reference
    names = {
        "heldout_topn_jaccard": "Jaccard",
        "heldout_moran_percentile": "MeanReferencePercentile",
        "heldout_rank_spearman": "RankSpearman",
        "heldout_nnsvg_topn_jaccard": "Jaccard",
        "heldout_nnsvg_percentile": "MeanReferencePercentile",
        "heldout_nnsvg_rank_spearman": "RankSpearman",
        "heldout_ari": "ARI",
        "heldout_nmi": "NMI",
    }
    frame["Metric"] = frame["metric"].map(names)
    return frame[frame["Metric"].notna()].copy()


def cluster_bootstrap_correlation(
    frame: pd.DataFrame,
    x_metric: str,
    y_metric: str,
    draws: int = 5000,
    seed: int = 20260729,
) -> dict[str, float]:
    wide = frame.pivot_table(
        index=["dataset", "heldout_slice", "method"],
        columns="Metric",
        values="value",
    ).reset_index()
    wide = wide.dropna(subset=[x_metric, y_metric])
    observed = spearmanr(wide[x_metric], wide[y_metric]).statistic
    fold = wide["dataset"].astype(str) + "::" + wide["heldout_slice"].astype(str)
    folds = fold.unique()
    rng = np.random.default_rng(seed)
    estimates = []
    for _ in range(draws):
        sampled = rng.choice(folds, len(folds), replace=True)
        indices = np.concatenate([np.flatnonzero(fold.to_numpy() == item) for item in sampled])
        estimate = spearmanr(wide.iloc[indices][x_metric], wide.iloc[indices][y_metric]).statistic
        if np.isfinite(estimate):
            estimates.append(estimate)
    return {
        "NMethodFold": len(wide),
        "SpearmanRho": observed,
        "Lower95": np.quantile(estimates, 0.025),
        "Upper95": np.quantile(estimates, 0.975),
    }


def heldout_cross_reference() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    combined = pd.concat(
        [
            load_heldout(MORAN_HELDOUT, "Moran's I"),
            load_heldout(NNSVG_HELDOUT, "nnSVG"),
        ],
        ignore_index=True,
    )
    summary = (
        combined.groupby(["Reference", "method", "Metric"], observed=True)
        .agg(Mean=("value", "mean"), SD=("value", "std"), NFolds=("value", "count"))
        .reset_index()
    )
    rankings = summary[summary["Metric"].isin(["Jaccard", "MeanReferencePercentile", "RankSpearman"])].copy()
    rankings["Rank"] = rankings.groupby(["Reference", "Metric"], observed=True)["Mean"].rank(
        method="average", ascending=False
    )
    concordance_rows = []
    for metric, part in rankings.groupby("Metric", observed=True):
        wide = part.pivot(index="method", columns="Reference", values="Rank").dropna()
        concordance_rows.append(
            {
                "Metric": metric,
                "NMethods": len(wide),
                "SpearmanRho": spearmanr(wide.iloc[:, 0], wide.iloc[:, 1]).statistic,
            }
        )
    for reference, part in combined.groupby("Reference", observed=True):
        for metric in ["Jaccard", "MeanReferencePercentile", "RankSpearman"]:
            correlation = cluster_bootstrap_correlation(part, metric, "ARI")
            concordance_rows.append(
                {
                    "Metric": f"{metric} vs ARI ({reference} reference)",
                    "NMethods": correlation.pop("NMethodFold"),
                    **correlation,
                }
            )
    return combined, summary, pd.DataFrame(concordance_rows)


def main() -> None:
    OUTPUT.mkdir(parents=True, exist_ok=True)
    (
        scaled,
        seed_scores,
        setting_scores,
        global_summary,
        effective,
        comparison,
        hybrid_differences,
    ) = hybrid_scores()
    heldout, heldout_summary, heldout_concordance = heldout_cross_reference()
    tables = {
        "hybrid_scaled_seed_metrics.tsv": scaled,
        "hybrid_seed_task_scores.tsv": seed_scores,
        "hybrid_setting_task_scores.tsv": setting_scores,
        "hybrid_global_rank_summary.tsv": global_summary,
        "hybrid_effective_features.tsv": effective,
        "hybrid_comparison_task_scores.tsv": comparison,
        "hybrid_differences_from_best_parent.tsv": hybrid_differences,
        "heldout_cross_reference_values.tsv": heldout,
        "heldout_cross_reference_summary.tsv": heldout_summary,
        "heldout_cross_reference_concordance.tsv": heldout_concordance,
    }
    for name, frame in tables.items():
        frame.to_csv(OUTPUT / name, sep="\t", index=False)
    manifest = {
        "hybrid_metric_rows": len(scaled),
        "hybrid_seed_task_rows": len(seed_scores),
        "hybrid_setting_rows": len(setting_scores),
        "heldout_value_rows": len(heldout),
        "heldout_folds": int(heldout[["dataset", "heldout_slice"]].drop_duplicates().shape[0]),
        "frozen_ranges_reused": str(FROZEN / "frozen_metric_ranges.tsv"),
    }
    (OUTPUT / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")


if __name__ == "__main__":
    main()

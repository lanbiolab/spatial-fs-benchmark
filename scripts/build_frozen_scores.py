from __future__ import annotations

import argparse
import hashlib
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd


ROOT_DIR = Path(__file__).resolve().parents[1]
SRC_DIR = ROOT_DIR / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from spatial_fs_benchmark.config import load_benchmark_config


METRIC_METADATA = {
    "bASW": ("Integration", "BatchMixing", True),
    "iLISI": ("Integration", "BatchMixing", True),
    "dASW": ("Integration", "BiologicalConservation", True),
    "dLISI": ("Integration", "BiologicalConservation", True),
    "ILL": ("Integration", "BiologicalConservation", True),
    "GC": ("Integration", "BiologicalConservation", True),
    "silhouette": ("Clustering", "LabelIndependentClustering", True),
    "CHAOS": ("Clustering", "LabelIndependentClustering", False),
    "PAS": ("Clustering", "LabelIndependentClustering", False),
    "ari": ("Clustering", "LabelAgreement", True),
    "nmi": ("Clustering", "LabelAgreement", True),
    "Accuracy": ("Alignment", "Alignment", True),
    "Ratio": ("Alignment", "Alignment", False),
}

COMPONENT_METRICS = {
    "BatchMixing": ["bASW", "iLISI"],
    "BiologicalConservation": ["dASW", "dLISI", "ILL", "GC"],
    "LabelIndependentClustering": ["silhouette", "CHAOS", "PAS"],
    "LabelAgreement": ["ari", "nmi"],
    "Alignment": ["Accuracy", "Ratio"],
}

SPATIAL_METHODS = {"morans_i", "sparkx", "nnsvg", "spatialde", "somde"}
CONTROL_METHODS = {"all_features", "random", "TFs", "scsegindex"}
LABEL_INFORMED_METHODS = {"wilcoxon"}

METRIC_KEY = [
    "dataset",
    "fs_method",
    "n_features",
    "integration_method",
    "task",
    "metric_name",
    "random_seed",
]
SETTING_KEY = ["dataset", "fs_method", "n_features", "integration_method"]
SETTING_METRIC_KEY = [*SETTING_KEY, "task", "metric_name"]


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Build frozen metric ranges, task scores, and rankings from audited results."
    )
    parser.add_argument(
        "--results",
        default="results/spatial_svg_rebuild_v1/merged_results.csv",
        help="Merged benchmark result CSV.",
    )
    parser.add_argument(
        "--output-dir",
        default="results/spatial_svg_rebuild_v1/frozen_scores",
        help="Directory for frozen scoring tables and manifest.",
    )
    parser.add_argument(
        "--canonical-config-root",
        action="append",
        default=[],
        help="Canonical config root. May be supplied more than once.",
    )
    return parser


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def method_group(method: str) -> str:
    if method in SPATIAL_METHODS:
        return "spatially_informed"
    if method in CONTROL_METHODS:
        return "control"
    if method in LABEL_INFORMED_METHODS:
        return "label_informed"
    return "expression_driven"


def normalize_n_features(value: object) -> str:
    text = str(value)
    if text.lower() == "all":
        return "all"
    number = float(text)
    if not number.is_integer():
        raise ValueError(f"Non-integral feature count: {value}")
    return str(int(number))


def metric_directions() -> pd.DataFrame:
    rows = []
    for metric, (task_group, component, higher_better) in METRIC_METADATA.items():
        rows.append(
            {
                "Metric": metric,
                "TaskGroup": task_group,
                "Component": component,
                "HigherBetter": higher_better,
                "DirectionMultiplier": 1 if higher_better else -1,
            }
        )
    return pd.DataFrame(rows)


def seed_averaged_metrics(results: pd.DataFrame) -> pd.DataFrame:
    return (
        results.groupby(SETTING_METRIC_KEY, dropna=False, observed=True)
        .agg(
            ValueMeanRaw=("metric_value", "mean"),
            ValueSDRaw=("metric_value", "std"),
            OrientedMean=("OrientedValue", "mean"),
            OrientedSD=("OrientedValue", "std"),
            NSeeds=("random_seed", "nunique"),
            EffectiveFeaturesMin=("effective_n_features", "min"),
            EffectiveFeaturesMax=("effective_n_features", "max"),
            RuntimeMean=("runtime", "mean"),
        )
        .reset_index()
    )


def frozen_ranges(setting_metrics: pd.DataFrame) -> pd.DataFrame:
    ranges = (
        setting_metrics.groupby(["dataset", "task", "metric_name"], dropna=False, observed=True)
        .agg(
            Lower=("OrientedMean", "min"),
            Upper=("OrientedMean", "max"),
            NFiniteSettings=("OrientedMean", lambda values: int(np.isfinite(values).sum())),
        )
        .reset_index()
    )
    ranges["Range"] = ranges["Upper"] - ranges["Lower"]
    ranges["RangeStatus"] = np.select(
        [
            ranges["NFiniteSettings"].eq(0),
            ranges["Range"].eq(0) & ranges["Range"].notna(),
        ],
        ["unavailable", "constant"],
        default="ok",
    )
    ranges["FrozenScaleDenominator"] = ranges["Range"].where(ranges["RangeStatus"].eq("ok"))
    ranges["RangeSource"] = "seed-averaged setting-level oriented metric values"
    return ranges


def scale_metrics(
    metrics: pd.DataFrame,
    ranges: pd.DataFrame,
    value_column: str,
) -> pd.DataFrame:
    scaled = metrics.merge(
        ranges[
            [
                "dataset",
                "task",
                "metric_name",
                "Lower",
                "Upper",
                "Range",
                "RangeStatus",
                "FrozenScaleDenominator",
            ]
        ],
        on=["dataset", "task", "metric_name"],
        how="left",
        validate="many_to_one",
    )
    scaled["ScaledValue"] = (
        scaled[value_column] - scaled["Lower"]
    ) / scaled["FrozenScaleDenominator"]
    return scaled


def strict_mean(frame: pd.DataFrame, columns: list[str]) -> pd.Series:
    values = frame.reindex(columns=columns)
    complete = values.notna().all(axis=1)
    return values.mean(axis=1).where(complete)


def available_equal_mean(first: pd.Series, second: pd.Series) -> pd.Series:
    values = pd.concat([first, second], axis=1)
    return values.mean(axis=1, skipna=True).where(values.notna().any(axis=1))


def build_seed_task_scores(seed_scaled: pd.DataFrame) -> pd.DataFrame:
    index_columns = [*SETTING_KEY, "random_seed"]
    wide = seed_scaled.pivot(
        index=index_columns,
        columns="metric_name",
        values="ScaledValue",
    ).reset_index()
    wide.columns.name = None

    wide["BatchMixing"] = strict_mean(wide, COMPONENT_METRICS["BatchMixing"])
    wide["BiologicalConservation"] = strict_mean(
        wide, COMPONENT_METRICS["BiologicalConservation"]
    )
    wide["Integration"] = available_equal_mean(
        wide["BatchMixing"], wide["BiologicalConservation"]
    )
    wide["LabelIndependentClustering"] = strict_mean(
        wide, COMPONENT_METRICS["LabelIndependentClustering"]
    )
    wide["LabelAgreement"] = strict_mean(wide, COMPONENT_METRICS["LabelAgreement"])
    wide["Clustering"] = available_equal_mean(
        wide["LabelIndependentClustering"], wide["LabelAgreement"]
    )
    wide["Alignment"] = strict_mean(wide, COMPONENT_METRICS["Alignment"])

    core = wide[["Integration", "Clustering"]]
    wide["CoreOverall"] = core.mean(axis=1).where(core.notna().all(axis=1))
    alignment_eligible = wide[["Integration", "Clustering", "Alignment"]]
    wide["AlignmentEligibleOverall"] = alignment_eligible.mean(axis=1).where(
        alignment_eligible.notna().all(axis=1)
    )
    wide["MethodGroup"] = wide["fs_method"].map(method_group)
    return wide


def aggregate_setting_scores(seed_scores: pd.DataFrame) -> pd.DataFrame:
    score_columns = [
        "BatchMixing",
        "BiologicalConservation",
        "Integration",
        "LabelIndependentClustering",
        "LabelAgreement",
        "Clustering",
        "Alignment",
        "CoreOverall",
        "AlignmentEligibleOverall",
    ]
    aggregations: dict[str, tuple[str, str]] = {}
    for column in score_columns:
        aggregations[f"{column}Mean"] = (column, "mean")
        aggregations[f"{column}SD"] = (column, "std")
        aggregations[f"{column}NSeeds"] = (column, "count")
    summary = (
        seed_scores.groupby(SETTING_KEY, dropna=False, observed=True)
        .agg(**aggregations)
        .reset_index()
    )
    summary["MethodGroup"] = summary["fs_method"].map(method_group)
    return summary


def representative_settings(config_roots: list[Path]) -> pd.DataFrame:
    rows: dict[tuple[str, str], dict] = {}
    for config_root in config_roots:
        for path in sorted(config_root.glob("*.yaml")):
            config = load_benchmark_config(path)
            for dataset in config.datasets:
                for method in config.feature_selection_methods:
                    feature_counts = method.n_features or config.n_features
                    if len(feature_counts) != 1:
                        raise ValueError(
                            f"Canonical config {path} has multiple feature counts for {method.name}: "
                            f"{feature_counts}"
                        )
                    n_features = (
                        "all"
                        if method.name == "all_features"
                        else normalize_n_features(feature_counts[0])
                    )
                    key = (dataset.name, method.name)
                    candidate = {
                        "dataset": dataset.name,
                        "fs_method": method.name,
                        "n_features": n_features,
                        "MethodGroup": method_group(method.name),
                        "CanonicalConfig": str(path.relative_to(ROOT_DIR)),
                    }
                    previous = rows.get(key)
                    if previous is not None and previous["n_features"] != n_features:
                        raise ValueError(f"Conflicting canonical settings for {key}: {previous} vs {candidate}")
                    rows[key] = candidate
    return pd.DataFrame(rows.values()).sort_values(["dataset", "fs_method"]).reset_index(drop=True)


def add_dataset_ranks(scores: pd.DataFrame) -> pd.DataFrame:
    ranked = scores.copy()
    score_names = [
        "BatchMixing",
        "BiologicalConservation",
        "Integration",
        "LabelIndependentClustering",
        "LabelAgreement",
        "Clustering",
        "Alignment",
        "CoreOverall",
        "AlignmentEligibleOverall",
    ]
    for score in score_names:
        value_column = f"{score}Mean"
        rank_column = f"{score}Rank"
        ranked[rank_column] = ranked.groupby(
            ["dataset", "integration_method"], observed=True
        )[value_column].rank(method="average", ascending=False, na_option="keep")
    return ranked


def global_rank_summary(dataset_ranks: pd.DataFrame) -> pd.DataFrame:
    score_names = ["Integration", "Clustering", "Alignment", "CoreOverall", "AlignmentEligibleOverall"]
    aggregations: dict[str, tuple[str, str]] = {}
    for score in score_names:
        aggregations[f"Mean{score}Score"] = (f"{score}Mean", "mean")
        aggregations[f"SD{score}Score"] = (f"{score}Mean", "std")
        aggregations[f"Mean{score}Rank"] = (f"{score}Rank", "mean")
        aggregations[f"SD{score}Rank"] = (f"{score}Rank", "std")
        aggregations[f"N{score}Datasets"] = (f"{score}Rank", "count")
    summary = (
        dataset_ranks.groupby(
            ["fs_method", "integration_method", "MethodGroup"],
            dropna=False,
            observed=True,
        )
        .agg(**aggregations)
        .reset_index()
    )
    summary["GlobalCoreOverallRank"] = summary.groupby(
        "integration_method", observed=True
    )["MeanCoreOverallRank"].rank(method="average", ascending=True, na_option="keep")
    return summary.sort_values(
        ["integration_method", "GlobalCoreOverallRank", "fs_method"]
    ).reset_index(drop=True)


def audit_rows(
    results: pd.DataFrame,
    setting_metrics: pd.DataFrame,
    ranges: pd.DataFrame,
    setting_scaled: pd.DataFrame,
    setting_scores: pd.DataFrame,
    representatives: pd.DataFrame,
    representative_scores: pd.DataFrame,
    competitive_scores: pd.DataFrame,
) -> pd.DataFrame:
    input_duplicates = int(results.duplicated(METRIC_KEY).sum())
    unknown_metrics = sorted(set(results["metric_name"]) - set(METRIC_METADATA))
    missing_ranges = int(
        setting_scaled.loc[setting_scaled["OrientedMean"].notna(), "RangeStatus"].isna().sum()
    )
    outside = setting_scaled.loc[
        setting_scaled["ScaledValue"].notna()
        & ~setting_scaled["ScaledValue"].between(-1e-10, 1 + 1e-10)
    ]
    core_missing = int(setting_scores["CoreOverallMean"].isna().sum())
    alignment_outside = int(
        setting_scores.loc[
            ~setting_scores["dataset"].isin(["DLPFC", "MouseBrainSerialSections"]),
            "AlignmentMean",
        ].notna().sum()
    )
    representative_duplicates = int(
        representatives.duplicated(["dataset", "fs_method"]).sum()
    )
    expected_representative_scores = len(representatives) * results["integration_method"].nunique()
    noncompetitive_in_competitive = int(
        (
            ~competitive_scores["MethodGroup"].isin(
                ["expression_driven", "spatially_informed"]
            )
        ).sum()
    )

    expected_range_rows = results.groupby(["dataset", "task", "metric_name"]).ngroups
    checks = [
        ("input metric-key duplicates", 0, input_duplicates, input_duplicates == 0),
        ("metric direction coverage", 0, len(unknown_metrics), len(unknown_metrics) == 0),
        (
            "frozen range rows",
            expected_range_rows,
            len(ranges),
            len(ranges) == expected_range_rows,
        ),
        ("finite setting metrics missing ranges", 0, missing_ranges, missing_ranges == 0),
        ("setting metric means outside [0,1]", 0, len(outside), len(outside) == 0),
        ("setting scores missing CoreOverall", 0, core_missing, core_missing == 0),
        ("alignment scores outside eligible datasets", 0, alignment_outside, alignment_outside == 0),
        ("duplicate representative settings", 0, representative_duplicates, representative_duplicates == 0),
        (
            "representative score coverage",
            expected_representative_scores,
            len(representative_scores),
            len(representative_scores) == expected_representative_scores,
        ),
        (
            "controls/oracles in competitive ranking",
            0,
            noncompetitive_in_competitive,
            noncompetitive_in_competitive == 0,
        ),
    ]
    rows = [
        {
            "Check": check,
            "Expected": expected,
            "Observed": observed,
            "Status": "PASS" if passed else "FAIL",
        }
        for check, expected, observed, passed in checks
    ]
    constant_ranges = int(ranges["RangeStatus"].eq("constant").sum())
    rows.append(
        {
            "Check": "constant frozen ranges",
            "Expected": 0,
            "Observed": constant_ranges,
            "Status": "PASS" if constant_ranges == 0 else "WARN",
        }
    )
    return pd.DataFrame(rows)


def main() -> None:
    args = build_parser().parse_args()
    results_path = (ROOT_DIR / args.results).resolve()
    output_dir = (ROOT_DIR / args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    canonical_roots = args.canonical_config_root or [
        "configs/rebuild_v1/downstream/canonical",
        "configs/rebuild_v1/non_spatial_downstream/canonical",
    ]
    canonical_config_roots = [(ROOT_DIR / value).resolve() for value in canonical_roots]

    results = pd.read_csv(results_path)
    results["n_features"] = results["n_features"].map(normalize_n_features)
    results["metric_value"] = pd.to_numeric(results["metric_value"], errors="coerce")
    results["random_seed"] = pd.to_numeric(results["random_seed"], errors="raise").astype(int)
    results["DirectionMultiplier"] = results["metric_name"].map(
        {metric: 1 if metadata[2] else -1 for metric, metadata in METRIC_METADATA.items()}
    )
    if results["DirectionMultiplier"].isna().any():
        unknown = sorted(results.loc[results["DirectionMultiplier"].isna(), "metric_name"].unique())
        raise ValueError(f"Missing metric directions: {unknown}")
    results["OrientedValue"] = results["metric_value"] * results["DirectionMultiplier"]

    directions = metric_directions()
    setting_metrics = seed_averaged_metrics(results)
    ranges = frozen_ranges(setting_metrics)
    setting_scaled = scale_metrics(setting_metrics, ranges, "OrientedMean")
    seed_scaled = scale_metrics(results, ranges, "OrientedValue")
    seed_scores = build_seed_task_scores(seed_scaled)
    setting_scores = aggregate_setting_scores(seed_scores)

    representatives = representative_settings(canonical_config_roots)
    representative_scores = representatives.merge(
        setting_scores,
        on=["dataset", "fs_method", "n_features", "MethodGroup"],
        how="left",
        validate="one_to_many",
    )
    dataset_ranks = add_dataset_ranks(representative_scores)
    global_ranks = global_rank_summary(dataset_ranks)
    competitive_scores = representative_scores.loc[
        representative_scores["MethodGroup"].isin(
            ["expression_driven", "spatially_informed"]
        )
    ].copy()
    competitive_dataset_ranks = add_dataset_ranks(competitive_scores)
    competitive_global_ranks = global_rank_summary(competitive_dataset_ranks)
    audit = audit_rows(
        results,
        setting_metrics,
        ranges,
        setting_scaled,
        setting_scores,
        representatives,
        representative_scores,
        competitive_scores,
    )

    tables = {
        "metric_directions.tsv": directions,
        "setting_metric_seed_summary.tsv": setting_metrics,
        "frozen_metric_ranges.tsv": ranges,
        "scaled_setting_metrics.tsv": setting_scaled,
        "scaled_seed_metrics.tsv": seed_scaled,
        "seed_task_scores.tsv": seed_scores,
        "setting_task_scores.tsv": setting_scores,
        "representative_settings.tsv": representatives,
        "representative_task_scores.tsv": representative_scores,
        "dataset_representative_ranks.tsv": dataset_ranks,
        "global_method_ranks.tsv": global_ranks,
        "dataset_competitive_ranks.tsv": competitive_dataset_ranks,
        "global_competitive_method_ranks.tsv": competitive_global_ranks,
        "scoring_audit.tsv": audit,
    }
    for filename, frame in tables.items():
        frame.to_csv(output_dir / filename, sep="\t", index=False)

    manifest = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "input_results": str(results_path),
        "input_sha256": sha256_file(results_path),
        "canonical_config_roots": [str(path) for path in canonical_config_roots],
        "scaling_unit": "Dataset x Metric",
        "range_source": "seed-averaged setting-level oriented metric values",
        "seed_scaled_values_clipped": False,
        "core_overall": "equal mean of Integration and Clustering",
        "alignment_eligible_overall": "equal mean of Integration, Clustering and Alignment",
        "tables": {filename: len(frame) for filename, frame in tables.items()},
    }
    (output_dir / "frozen_score_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )

    print(audit.to_string(index=False))
    print(f"\nFrozen scoring tables: {output_dir}")
    if audit["Status"].eq("FAIL").any():
        raise SystemExit(1)


if __name__ == "__main__":
    main()

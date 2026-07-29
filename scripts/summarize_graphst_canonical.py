from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "scripts") not in sys.path:
    sys.path.insert(0, str(ROOT / "scripts"))
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))

from build_frozen_scores import (  # noqa: E402
    METRIC_KEY,
    METRIC_METADATA,
    add_dataset_ranks,
    aggregate_setting_scores,
    build_seed_task_scores,
    global_rank_summary,
    normalize_n_features,
    scale_metrics,
    seed_averaged_metrics,
)
from spatial_fs_benchmark.benchmark.experiment import read_metric_records  # noqa: E402
from spatial_fs_benchmark.benchmark.result_schema import records_to_frame  # noqa: E402


RESULTS_ROOT = ROOT / "results/spatial_svg_rebuild_v1"
OUTPUT = ROOT / "results/graphst_canonical_v1/summary"
FROZEN = RESULTS_ROOT / "frozen_scores"


def read_graphst_records() -> pd.DataFrame:
    records = []
    for path in RESULTS_ROOT.rglob("graphst/*/*/*/*_records.json"):
        records.extend(read_metric_records(path))
    frame = records_to_frame(records)
    frame = frame.loc[frame["integration_method"].eq("graphst")].copy()
    frame["n_features"] = frame["n_features"].map(normalize_n_features)
    frame["metric_value"] = pd.to_numeric(frame["metric_value"], errors="coerce")
    frame["random_seed"] = pd.to_numeric(frame["random_seed"], errors="raise").astype(int)
    frame["DirectionMultiplier"] = frame["metric_name"].map(
        {metric: 1 if metadata[2] else -1 for metric, metadata in METRIC_METADATA.items()}
    )
    frame["OrientedValue"] = frame["metric_value"] * frame["DirectionMultiplier"]
    return frame


def feature_hash_audit() -> pd.DataFrame:
    rows = []
    for graphst_meta in RESULTS_ROOT.rglob("graphst/*/*/*/embedding.meta.json"):
        graphst = json.loads(graphst_meta.read_text(encoding="utf-8"))
        if graphst.get("integrator_params", {}).get("epochs") != 600:
            continue
        graphst_hash = graphst.get("selection_feature_hash", "")
        relative = graphst_meta.relative_to(RESULTS_ROOT)
        parts = list(relative.parts)
        graphst_index = parts.index("graphst")
        for primary in ("scvi", "cellcharter"):
            primary_parts = parts.copy()
            primary_parts[graphst_index] = primary
            primary_meta = RESULTS_ROOT.joinpath(*primary_parts)
            if not primary_meta.exists():
                rows.append(
                    {
                        "graphst_meta": str(graphst_meta),
                        "primary_integrator": primary,
                        "primary_meta": str(primary_meta),
                        "graphst_feature_hash": graphst_hash,
                        "primary_feature_hash": "",
                        "status": "MISSING_PRIMARY_META",
                    }
                )
                continue
            primary_payload = json.loads(primary_meta.read_text(encoding="utf-8"))
            primary_hash = primary_payload.get("selection_feature_hash", "")
            rows.append(
                {
                    "graphst_meta": str(graphst_meta),
                    "primary_integrator": primary,
                    "primary_meta": str(primary_meta),
                    "graphst_feature_hash": graphst_hash,
                    "primary_feature_hash": primary_hash,
                    "status": "PASS" if graphst_hash == primary_hash and graphst_hash else "HASH_MISMATCH",
                }
            )
    return pd.DataFrame(rows)


def integrator_summaries(
    dataset_ranks: pd.DataFrame,
    global_ranks: pd.DataFrame,
    bootstrap_draws: int = 5000,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    score_columns = ["CoreOverallMean", "IntegrationMean", "ClusteringMean", "AlignmentMean"]
    mean_rows = []
    for integrator, frame in dataset_ranks.groupby("integration_method", observed=True):
        row = {"integration_method": integrator, "comparison": "absolute mean"}
        row.update({column: frame[column].mean() for column in score_columns})
        mean_rows.append(row)
        if integrator == "scvi":
            continue
        reference = dataset_ranks.loc[
            dataset_ranks["integration_method"].eq("scvi"),
            ["dataset", "fs_method", *score_columns],
        ]
        matched = frame.merge(
            reference,
            on=["dataset", "fs_method"],
            suffixes=("", "_scvi"),
            validate="one_to_one",
        )
        diff_row = {"integration_method": integrator, "comparison": "difference from scVI"}
        diff_row.update(
            {column: (matched[column] - matched[f"{column}_scvi"]).mean() for column in score_columns}
        )
        mean_rows.append(diff_row)

    wide = global_ranks.pivot(
        index="fs_method",
        columns="integration_method",
        values="MeanCoreOverallRank",
    ).dropna()
    rng = np.random.default_rng(20260729)
    correlation_rows = []
    pairs = [("scvi", "cellcharter"), ("scvi", "graphst"), ("cellcharter", "graphst")]
    for first, second in pairs:
        observed = float(spearmanr(wide[first], wide[second]).statistic)
        bootstrap = np.empty(bootstrap_draws, dtype=float)
        values = wide[[first, second]].to_numpy()
        for draw in range(bootstrap_draws):
            sampled = values[rng.integers(0, len(values), size=len(values))]
            bootstrap[draw] = spearmanr(sampled[:, 0], sampled[:, 1]).statistic
        correlation_rows.append(
            {
                "integrator_1": first,
                "integrator_2": second,
                "n_methods": len(wide),
                "spearman_rho": observed,
                "bootstrap_lower_95": float(np.nanquantile(bootstrap, 0.025)),
                "bootstrap_upper_95": float(np.nanquantile(bootstrap, 0.975)),
                "bootstrap_draws": bootstrap_draws,
            }
        )
    return pd.DataFrame(mean_rows), pd.DataFrame(correlation_rows)


def main() -> None:
    OUTPUT.mkdir(parents=True, exist_ok=True)
    results = read_graphst_records()
    ranges = pd.read_csv(FROZEN / "frozen_metric_ranges.tsv", sep="\t")
    seed_scaled = scale_metrics(results, ranges, "OrientedValue")
    seed_scores = build_seed_task_scores(seed_scaled)
    setting_scores = aggregate_setting_scores(seed_scores)
    hash_audit = feature_hash_audit()

    original = pd.read_csv(FROZEN / "representative_task_scores.tsv", sep="\t")
    original_competitive = original.loc[
        original["MethodGroup"].isin(["expression_driven", "spatially_informed"])
    ].copy()
    combined = pd.concat([original_competitive, setting_scores], ignore_index=True, sort=False)
    dataset_ranks = add_dataset_ranks(combined)
    global_ranks = global_rank_summary(dataset_ranks)
    integrator_means, integrator_correlations = integrator_summaries(dataset_ranks, global_ranks)

    embedding_count = len(list(RESULTS_ROOT.rglob("graphst/*/*/*/embedding.npz")))
    task_record_count = len(list(RESULTS_ROOT.rglob("graphst/*/*/*/*_records.json")))
    duplicates = int(results.duplicated(METRIC_KEY).sum())
    nonfinite_mask = ~np.isfinite(results["metric_value"])
    nonfinite = int(nonfinite_mask.sum())
    allowed_nonfinite = (
        results["dataset"].isin(["E8p5Embryo", "E9p5Embryo"])
        & results["metric_name"].isin(["dASW", "dLISI", "ILL", "GC"])
    )
    unexpected_nonfinite = int((nonfinite_mask & ~allowed_nonfinite).sum())
    missing_ranges = int(seed_scaled.loc[seed_scaled["OrientedValue"].notna(), "RangeStatus"].isna().sum())
    hash_failures = int(hash_audit["status"].ne("PASS").sum())
    audit = pd.DataFrame(
        [
            ("GraphST embeddings", 609, embedding_count),
            ("GraphST task records", 1392, task_record_count),
            ("GraphST metric rows", 6699, len(results)),
            ("duplicate metric keys", 0, duplicates),
            ("expected label-dependent missing values", 696, nonfinite),
            ("unexpected non-finite metric values", 0, unexpected_nonfinite),
            ("finite values missing frozen ranges", 0, missing_ranges),
            ("GraphST setting scores", 203, len(setting_scores)),
            ("GraphST-primary feature-hash comparisons", 1218, len(hash_audit)),
            ("GraphST-primary feature-hash failures", 0, hash_failures),
        ],
        columns=["Check", "Expected", "Observed"],
    )
    audit["Status"] = audit["Expected"].eq(audit["Observed"]).map({True: "PASS", False: "FAIL"})

    tables = {
        "graphst_metric_records.tsv": results,
        "graphst_scaled_seed_metrics.tsv": seed_scaled,
        "graphst_seed_task_scores.tsv": seed_scores,
        "graphst_setting_task_scores.tsv": setting_scores,
        "three_integrator_dataset_ranks.tsv": dataset_ranks,
        "three_integrator_global_ranks.tsv": global_ranks,
        "integrator_mean_scores_and_differences.tsv": integrator_means,
        "integrator_rank_correlations.tsv": integrator_correlations,
        "feature_hash_audit.tsv": hash_audit,
        "audit.tsv": audit,
    }
    for filename, frame in tables.items():
        frame.to_csv(OUTPUT / filename, sep="\t", index=False)
    (OUTPUT / "manifest.json").write_text(
        json.dumps(
            {
                "range_source": str(FROZEN / "frozen_metric_ranges.tsv"),
                "ranges_recomputed": False,
                "scaled_values_clipped": False,
                "tables": {name: len(frame) for name, frame in tables.items()},
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )
    print(audit.to_string(index=False))
    if audit["Status"].eq("FAIL").any():
        raise SystemExit(1)


if __name__ == "__main__":
    main()

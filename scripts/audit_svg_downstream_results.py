from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd


ROOT_DIR = Path(__file__).resolve().parents[1]
SRC_DIR = ROOT_DIR / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from spatial_fs_benchmark.config import load_benchmark_config


INTEGRATOR_VERSIONS = {
    "scvi": "v2_integer_counts_batch_covariate",
    "cellcharter": "v2_native_scvi_spatial_aggregation_counts",
}
TASK_METRICS = {
    "integration_eval": ("bASW", "iLISI", "dASW", "dLISI", "ILL", "GC"),
    "clustering_eval": ("silhouette", "ari", "nmi", "CHAOS", "PAS"),
    "alignment_eval": ("Accuracy", "Ratio"),
}
UNLABELLED_INTEGRATION_METRICS = {"dASW", "dLISI", "ILL", "GC"}
UNLABELLED_CLUSTERING_METRICS = {"silhouette", "CHAOS", "PAS"}
KEY_COLUMNS = [
    "dataset",
    "fs_method",
    "n_features",
    "integration_method",
    "task",
    "metric_name",
    "random_seed",
]


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Audit rebuilt downstream benchmark results against frozen configs."
    )
    parser.add_argument(
        "--results-root",
        default="results/spatial_svg_rebuild_v1",
        help="Root containing merged_results.csv and per-dataset artifacts.",
    )
    parser.add_argument(
        "--config-root",
        action="append",
        default=[],
        help=(
            "Root containing canonical and feature-number benchmark configs. "
            "May be supplied more than once."
        ),
    )
    parser.add_argument(
        "--output-dir",
        default="",
        help="Audit output directory; defaults to RESULTS_ROOT/audit.",
    )
    return parser


def read_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def expected_tables(config_roots: list[Path]) -> tuple[pd.DataFrame, pd.DataFrame]:
    run_rows: dict[tuple, dict] = {}
    metric_rows: dict[tuple, dict] = {}

    for config_root in config_roots:
        for config_path in sorted(config_root.glob("*/*.yaml")):
            config = load_benchmark_config(config_path)
            for dataset in config.datasets:
                has_labels = bool(dataset.label_key)
                for selector in config.feature_selection_methods:
                    feature_counts = selector.n_features or config.n_features
                    for n_features in feature_counts:
                        for seed in config.seeds:
                            for integrator in config.integration_methods:
                                run_key = (
                                    dataset.name,
                                    selector.name,
                                    int(n_features),
                                    integrator.name,
                                    int(seed),
                                )
                                run_rows[run_key] = {
                                    "dataset": dataset.name,
                                    "fs_method": selector.name,
                                    "n_features": int(n_features),
                                    "integration_method": integrator.name,
                                    "random_seed": int(seed),
                                    "has_labels": has_labels,
                                    "dataset_dir": str(Path(config.output_dir) / dataset.name.lower()),
                                }
                                for task in config.tasks:
                                    metrics = TASK_METRICS[task.name]
                                    if task.name == "clustering_eval" and not has_labels:
                                        metrics = tuple(UNLABELLED_CLUSTERING_METRICS)
                                    for metric_name in metrics:
                                        metric_key = (*run_key[:4], task.name, metric_name, run_key[4])
                                        metric_rows[metric_key] = {
                                            **run_rows[run_key],
                                            "task": task.name,
                                            "metric_name": metric_name,
                                            "expected_nonfinite": bool(
                                                not has_labels
                                                and task.name == "integration_eval"
                                                and metric_name in UNLABELLED_INTEGRATION_METRICS
                                            ),
                                        }

    return pd.DataFrame(run_rows.values()), pd.DataFrame(metric_rows.values())


def selected_feature_audit(results_root: Path) -> tuple[pd.DataFrame, dict[tuple, set[str]]]:
    rows = []
    hashes: dict[tuple, set[str]] = {}
    for path in sorted(results_root.rglob("selected_features.json")):
        payload = read_json(path)
        metadata = payload.get("metadata", {})
        feature_names = [str(value) for value in payload.get("feature_names", [])]
        computed_hash = hashlib.md5("\n".join(feature_names).encode("utf-8")).hexdigest()
        stored_hash = str(metadata.get("feature_name_hash", ""))
        parts = path.parts
        feature_index = parts.index("feature_selection")
        dataset_key = parts[feature_index - 1]
        method = parts[feature_index + 1]
        n_features = int(parts[feature_index + 2].removeprefix("n"))
        seed = int(parts[feature_index + 3].removeprefix("seed"))
        key = (dataset_key, method, n_features)
        hashes.setdefault(key, set()).add(stored_hash)
        rows.append(
            {
                "path": str(path),
                "dataset_key": dataset_key,
                "fs_method": method,
                "n_features": n_features,
                "selection_seed": seed,
                "requested_n_features": metadata.get("requested_n_features", n_features),
                "effective_n_features": metadata.get("effective_n_features", len(feature_names)),
                "stored_feature_hash": stored_hash,
                "computed_feature_hash": computed_hash,
                "hash_matches": stored_hash == computed_hash,
                "feature_count_matches": int(metadata.get("effective_n_features", len(feature_names)))
                == len(feature_names),
            }
        )
    return pd.DataFrame(rows), hashes


def embedding_audit(
    expected_runs: pd.DataFrame,
    selection_hashes: dict[tuple, set[str]],
) -> pd.DataFrame:
    rows = []
    for run in expected_runs.itertuples(index=False):
        combination_dir = (
            ROOT_DIR
            / run.dataset_dir
            / run.integration_method
            / run.fs_method
            / f"n{run.n_features}"
            / f"seed{run.random_seed}"
        )
        meta_path = combination_dir / "embedding.meta.json"
        data_path = combination_dir / "embedding.npz"
        meta = read_json(meta_path) if meta_path.exists() else {}
        selection_hash = str(meta.get("selection_feature_hash", ""))
        dataset_key = run.dataset.lower()
        known_hashes = selection_hashes.get((dataset_key, run.fs_method, run.n_features), set())
        rows.append(
            {
                "dataset": run.dataset,
                "fs_method": run.fs_method,
                "n_features": run.n_features,
                "integration_method": run.integration_method,
                "random_seed": run.random_seed,
                "embedding_exists": data_path.exists(),
                "metadata_exists": meta_path.exists(),
                "integrator_version": meta.get("integrator_implementation_version", ""),
                "expected_integrator_version": INTEGRATOR_VERSIONS.get(run.integration_method, ""),
                "version_matches": meta.get("integrator_implementation_version")
                == INTEGRATOR_VERSIONS.get(run.integration_method),
                "selection_feature_hash": selection_hash,
                "selection_hash_found": bool(selection_hash) and selection_hash in known_hashes,
                "cache_signature": meta.get("cache_signature", ""),
                "path": str(data_path),
            }
        )
    return pd.DataFrame(rows)


def task_metadata_audit(expected_metrics: pd.DataFrame, embedding: pd.DataFrame) -> pd.DataFrame:
    expected_tasks = expected_metrics[
        ["dataset", "fs_method", "n_features", "integration_method", "random_seed", "task", "dataset_dir"]
    ].drop_duplicates()
    embedding_lookup = embedding.set_index(
        ["dataset", "fs_method", "n_features", "integration_method", "random_seed"]
    )
    rows = []
    for task in expected_tasks.itertuples(index=False):
        combination_dir = (
            ROOT_DIR
            / task.dataset_dir
            / task.integration_method
            / task.fs_method
            / f"n{task.n_features}"
            / f"seed{task.random_seed}"
        )
        records_path = combination_dir / f"{task.task}_records.json"
        meta_path = combination_dir / f"{task.task}_records.meta.json"
        meta = read_json(meta_path) if meta_path.exists() else {}
        embedding_key = (
            task.dataset,
            task.fs_method,
            task.n_features,
            task.integration_method,
            task.random_seed,
        )
        embedding_row = embedding_lookup.loc[embedding_key]
        rows.append(
            {
                "dataset": task.dataset,
                "fs_method": task.fs_method,
                "n_features": task.n_features,
                "integration_method": task.integration_method,
                "random_seed": task.random_seed,
                "task": task.task,
                "records_exist": records_path.exists(),
                "metadata_exists": meta_path.exists(),
                "selection_hash_matches_embedding": meta.get("selection_feature_hash")
                == embedding_row["selection_feature_hash"],
                "integration_signature_matches": meta.get("integration_cache_signature")
                == embedding_row["cache_signature"],
                "task_version_present": bool(meta.get("task_implementation_version")),
                "records_path": str(records_path),
            }
        )
    return pd.DataFrame(rows)


def status_row(check: str, expected: int | str, observed: int | str, passed: bool, details: str = "") -> dict:
    return {
        "check": check,
        "expected": expected,
        "observed": observed,
        "status": "PASS" if passed else "FAIL",
        "details": details,
    }


def warning_row(check: str, expected: int | str, observed: int | str, details: str) -> dict:
    return {
        "check": check,
        "expected": expected,
        "observed": observed,
        "status": "WARN",
        "details": details,
    }


def main() -> None:
    args = build_parser().parse_args()
    results_root = (ROOT_DIR / args.results_root).resolve()
    config_root_args = args.config_root or [
        "configs/rebuild_v1/downstream",
        "configs/rebuild_v1/non_spatial_downstream",
    ]
    config_roots = [(ROOT_DIR / value).resolve() for value in config_root_args]
    output_dir = Path(args.output_dir).resolve() if args.output_dir else results_root / "audit"
    output_dir.mkdir(parents=True, exist_ok=True)

    merged_path = results_root / "merged_results.csv"
    results = pd.read_csv(merged_path)
    results["n_features"] = results["n_features"].replace("all", 1)
    results["n_features"] = pd.to_numeric(results["n_features"], errors="raise").astype(int)
    results["random_seed"] = pd.to_numeric(results["random_seed"], errors="raise").astype(int)
    expected_runs, expected_metrics = expected_tables(config_roots)

    observed_keys = results[KEY_COLUMNS].copy()
    expected_keys = expected_metrics[KEY_COLUMNS].copy()
    coverage = expected_keys.merge(
        results[KEY_COLUMNS + ["metric_value"]],
        how="outer",
        on=KEY_COLUMNS,
        indicator=True,
    ).merge(
        expected_metrics[KEY_COLUMNS + ["expected_nonfinite"]],
        how="left",
        on=KEY_COLUMNS,
    )
    coverage["coverage_status"] = coverage["_merge"].map(
        {"left_only": "missing", "right_only": "unexpected", "both": "present"}
    )
    coverage = coverage.drop(columns="_merge")

    duplicates = results.loc[results.duplicated(KEY_COLUMNS, keep=False)].sort_values(KEY_COLUMNS)
    metric_values = pd.to_numeric(results["metric_value"], errors="coerce")
    results_with_flags = results.copy()
    results_with_flags["is_finite"] = np.isfinite(metric_values)
    results_with_flags["is_zero"] = metric_values.eq(0) & results_with_flags["is_finite"]
    results_with_flags = results_with_flags.merge(
        expected_metrics[KEY_COLUMNS + ["expected_nonfinite"]],
        how="left",
        on=KEY_COLUMNS,
    )
    results_with_flags["expected_nonfinite"] = results_with_flags["expected_nonfinite"].fillna(False)
    unexpected_nonfinite = results_with_flags.loc[
        ~results_with_flags["is_finite"] & ~results_with_flags["expected_nonfinite"]
    ]
    expected_nonfinite_but_finite = results_with_flags.loc[
        results_with_flags["is_finite"] & results_with_flags["expected_nonfinite"]
    ]

    selection, selection_hashes = selected_feature_audit(results_root)
    embedding = embedding_audit(expected_runs, selection_hashes)
    task_metadata = task_metadata_audit(expected_metrics, embedding)

    metric_summary = (
        results_with_flags.groupby(["dataset", "task", "metric_name"], dropna=False)
        .agg(
            n=("metric_value", "size"),
            n_finite=("is_finite", "sum"),
            n_nonfinite=("is_finite", lambda values: int((~values).sum())),
            n_zero=("is_zero", "sum"),
            minimum=("metric_value", "min"),
            maximum=("metric_value", "max"),
            mean=("metric_value", "mean"),
        )
        .reset_index()
    )
    coverage_summary = (
        coverage.groupby(["dataset", "task", "coverage_status"], observed=False)
        .size()
        .rename("n")
        .reset_index()
    )

    missing = coverage.loc[coverage["coverage_status"] == "missing"]
    unexpected = coverage.loc[coverage["coverage_status"] == "unexpected"]
    embedding_failures = embedding.loc[
        ~embedding[
            ["embedding_exists", "metadata_exists", "version_matches", "selection_hash_found"]
        ].all(axis=1)
    ]
    task_link_columns = [
        "records_exist",
        "metadata_exists",
        "selection_hash_matches_embedding",
        "integration_signature_matches",
    ]
    task_failures = task_metadata.loc[
        ~task_metadata[
            task_link_columns
        ].all(axis=1)
    ]
    selection_failures = selection.loc[~selection[["hash_matches", "feature_count_matches"]].all(axis=1)]

    alignment_datasets = sorted(results.loc[results["task"] == "alignment_eval", "dataset"].unique())
    expected_alignment_datasets = sorted(
        expected_metrics.loc[expected_metrics["task"] == "alignment_eval", "dataset"].unique()
    )
    expected_task_columns = [
        "dataset",
        "fs_method",
        "n_features",
        "integration_method",
        "task",
        "random_seed",
    ]
    expected_task_count = len(expected_metrics[expected_task_columns].drop_duplicates())
    task_versions_present = int(task_metadata["task_version_present"].sum())
    missing_task_versions = (
        task_metadata.loc[~task_metadata["task_version_present"]]
        .groupby("task")
        .size()
        .sort_index()
    )
    missing_task_version_details = ", ".join(
        f"{task}={count}" for task, count in missing_task_versions.items()
    )
    summary = pd.DataFrame(
        [
            status_row("unique embeddings", len(expected_runs), len(embedding), len(embedding) == len(expected_runs)),
            status_row(
                "task record payloads",
                expected_task_count,
                len(task_metadata),
                len(task_metadata) == expected_task_count,
            ),
            status_row("merged metric rows", len(expected_metrics), len(results), len(results) == len(expected_metrics)),
            status_row("missing metric keys", 0, len(missing), len(missing) == 0),
            status_row("unexpected metric keys", 0, len(unexpected), len(unexpected) == 0),
            status_row("duplicate metric keys", 0, len(duplicates), len(duplicates) == 0),
            status_row("unexpected non-finite values", 0, len(unexpected_nonfinite), len(unexpected_nonfinite) == 0),
            status_row(
                "unlabelled metrics stored as NA",
                int(expected_metrics["expected_nonfinite"].sum()),
                int((~results_with_flags["is_finite"] & results_with_flags["expected_nonfinite"]).sum()),
                len(expected_nonfinite_but_finite) == 0,
                "E8p5Embryo/E9p5Embryo dASW, dLISI, ILL and GC",
            ),
            status_row("selected-feature payload integrity", 0, len(selection_failures), len(selection_failures) == 0),
            status_row("embedding metadata integrity", 0, len(embedding_failures), len(embedding_failures) == 0),
            status_row("task metadata linkage integrity", 0, len(task_failures), len(task_failures) == 0),
            warning_row(
                "task implementation version recorded",
                len(task_metadata),
                task_versions_present,
                f"Missing explicit implementation-version strings: {missing_task_version_details}",
            )
            if task_versions_present < len(task_metadata)
            else status_row(
                "task implementation version recorded",
                len(task_metadata),
                task_versions_present,
                True,
            ),
            status_row(
                "alignment dataset scope",
                ",".join(expected_alignment_datasets),
                ",".join(alignment_datasets),
                alignment_datasets == expected_alignment_datasets,
            ),
        ]
    )

    tables = {
        "audit_summary.tsv": summary,
        "metric_coverage.tsv": coverage,
        "coverage_summary.tsv": coverage_summary,
        "metric_quality_summary.tsv": metric_summary,
        "duplicate_metric_rows.tsv": duplicates,
        "unexpected_nonfinite_metrics.tsv": unexpected_nonfinite,
        "expected_nonfinite_but_finite.tsv": expected_nonfinite_but_finite,
        "selected_feature_audit.tsv": selection,
        "selected_feature_failures.tsv": selection_failures,
        "embedding_audit.tsv": embedding,
        "embedding_failures.tsv": embedding_failures,
        "task_metadata_audit.tsv": task_metadata,
        "task_metadata_failures.tsv": task_failures,
    }
    for filename, frame in tables.items():
        frame.to_csv(output_dir / filename, sep="\t", index=False)

    print(summary.to_string(index=False))
    print(f"\nAudit tables: {output_dir}")
    if (summary["status"] == "FAIL").any():
        raise SystemExit(1)


if __name__ == "__main__":
    main()

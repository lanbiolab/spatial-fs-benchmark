from __future__ import annotations

import argparse
from pathlib import Path

import anndata as ad
import pandas as pd
import yaml


TASK_TYPE_MAP = {
    "integration_eval": "Integration",
    "clustering_eval": "Clustering",
    "alignment_eval": "Alignment",
    "slice_representation_eval": "SliceRepresentation",
}

METRIC_METADATA = {
    "dASW": {"name": "dASW", "included": True, "higher_better": True, "type": "Integration"},
    "dLISI": {"name": "dLISI", "included": True, "higher_better": True, "type": "Integration"},
    "ILL": {"name": "ILL", "included": True, "higher_better": True, "type": "Integration"},
    "bASW": {"name": "bASW", "included": True, "higher_better": True, "type": "Integration"},
    "iLISI": {"name": "iLISI", "included": True, "higher_better": True, "type": "Integration"},
    "GC": {"name": "GC", "included": True, "higher_better": True, "type": "Integration"},
    "ari": {"name": "ARI", "included": True, "higher_better": True, "type": "Clustering"},
    "nmi": {"name": "NMI", "included": True, "higher_better": True, "type": "Clustering"},
    "CHAOS": {"name": "CHAOS", "included": True, "higher_better": False, "type": "Clustering"},
    "PAS": {"name": "PAS", "included": True, "higher_better": False, "type": "Clustering"},
    "Accuracy": {"name": "Accuracy", "included": True, "higher_better": True, "type": "Alignment"},
    "Ratio": {"name": "Ratio", "included": True, "higher_better": False, "type": "Alignment"},
    "slice_repr_ARI": {"name": "ARI", "included": True, "higher_better": True, "type": "SliceRepresentation"},
    "slice_repr_NMI": {"name": "NMI", "included": True, "higher_better": True, "type": "SliceRepresentation"},
    "slice_repr_distance_mean": {
        "name": "DistanceMean",
        "included": False,
        "higher_better": False,
        "type": "SliceRepresentation",
    },
    "silhouette": {"name": "Silhouette", "included": False, "higher_better": True, "type": "Clustering"},
}

METHOD_METADATA = {
    "all": {"name": "All features", "is_baseline": True},
    "all_features": {"name": "All features", "is_baseline": True},
    "hvg": {"name": "HVG", "is_baseline": False},
    "svg": {"name": "SVG", "is_baseline": False},
    "highly_expressed": {"name": "Highly Expressed", "is_baseline": True},
    "random": {"name": "Random", "is_baseline": True},
    "scanpy_seurat": {"name": "scanpy-Seurat", "is_baseline": False},
    "scanpy_seurat_batch": {"name": "scanpy-Seurat (batch)", "is_baseline": False},
    "scanpy_cell_ranger": {"name": "scanpy-CellRanger", "is_baseline": False},
    "scanpy_cell_ranger_batch": {"name": "scanpy-CellRanger (batch)", "is_baseline": False},
    "scanpy_pearson": {"name": "scanpy-Pearson", "is_baseline": False},
    "scanpy_pearson_batch": {"name": "scanpy-Pearson (batch)", "is_baseline": False},
    "scanpy_seurat_v3": {"name": "scanpy-SeuratV3", "is_baseline": False},
    "scanpy_seurat_v3_batch": {"name": "scanpy-SeuratV3 (batch)", "is_baseline": False},
    "seurat_vst": {"name": "Seurat-VST", "is_baseline": False},
    "seurat_mvp": {"name": "Seurat-MVP", "is_baseline": False},
    "seurat_disp": {"name": "Seurat-Dispersion", "is_baseline": False},
    "seurat_sct": {"name": "Seurat-scTransform", "is_baseline": False},
    "scsegindex": {"name": "scSEGIndex", "is_baseline": False},
    "dubstepr": {"name": "DUBStepR", "is_baseline": False},
    "nbumi": {"name": "NBumi", "is_baseline": False},
    "osca": {"name": "OSCA", "is_baseline": False},
    "scry": {"name": "scry", "is_baseline": False},
    "singleCellHaystack": {"name": "singleCellHaystack", "is_baseline": False},
    "Brennecke": {"name": "Brennecke", "is_baseline": False},
    "scPNMF": {"name": "scPNMF", "is_baseline": False},
    "TFs": {"name": "Transcription factors", "is_baseline": False},
    "statistic_mean": {"name": "Statistic mean", "is_baseline": False},
    "statistic_variance": {"name": "Statistic variance", "is_baseline": False},
    "triku": {"name": "triku", "is_baseline": False},
    "hotspot": {"name": "Hotspot", "is_baseline": False},
    "anticor": {"name": "Anticor", "is_baseline": False},
    "wilcoxon": {"name": "Wilcoxon", "is_baseline": False},
    "conserved": {"name": "Cross-slice conserved", "is_baseline": False},
    "cross_slice_conserved": {"name": "Cross-slice conserved", "is_baseline": False},
}

INTEGRATION_METADATA = {
    "pca": "PCA",
    "combat": "ComBat",
    "scanorama": "Scanorama",
    "scvi": "scVI",
    "symphony": "Symphony",
    "cellcharter": "CellCharter",
    "gpsa": "GPSA",
    "staligner": "STAligner",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Adapt spatial benchmark results to atlas-style plotting inputs.")
    parser.add_argument("--results", default="results/baristaseq_mvp/results.csv", help="Path to benchmark results.csv")
    parser.add_argument("--configs-dir", default="configs/datasets", help="Directory containing dataset YAML configs")
    parser.add_argument("--output-dir", default="results/atlas_style", help="Output directory for adapted TSV files")
    parser.add_argument(
        "--slice-repr-datasets",
        default="",
        help="Comma-separated dataset names to retain for SliceRepresentation metrics. Others are excluded only from adapted plotting tables.",
    )
    return parser.parse_args()


def load_dataset_metadata(configs_dir: Path, datasets: list[str]) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for config_path in sorted(configs_dir.glob("*.yaml")):
        with config_path.open("r", encoding="utf-8") as handle:
            config = yaml.safe_load(handle)
        if config["name"] not in datasets:
            continue
        adata = ad.read_h5ad(config["path"])
        slice_values = adata.obs[config["slice_key"]].astype(str)
        n_slices = int(slice_values.nunique())
        n_spots = int(adata.n_obs)
        label_key = config.get("label_key")
        if label_key and label_key in adata.obs:
            label_values = adata.obs[label_key]
            n_labels = int(label_values.dropna().astype(str).nunique())
        else:
            n_labels = None
        rows.append(
            {
                "Dataset": config["name"],
                "Name": config["name"],
                "Platform": config.get("platform"),
                "NSlices": n_slices,
                "Spots": n_spots,
                "Features": int(adata.n_vars),
                "Labels": n_labels,
                "SpotsPerSlice": float(n_spots / n_slices) if n_slices > 0 else None,
            }
        )
    return pd.DataFrame(rows)


def build_metrics_metadata() -> pd.DataFrame:
    rows = []
    for metric, meta in METRIC_METADATA.items():
        rows.append(
            {
                "Metric": metric,
                "Name": meta["name"],
                "Included": meta["included"],
                "HigherBetter": meta["higher_better"],
                "Type": meta["type"],
            }
        )
    return pd.DataFrame(rows)


def sort_feature_levels(values: list[object]) -> list[object]:
    def _key(value: object) -> tuple[int, int | str]:
        text = str(value)
        if text == "all":
            return (1, text)
        return (0, int(text))

    return sorted(values, key=_key)


def build_methods_metadata(methods: list[str], n_features: list[object]) -> pd.DataFrame:
    rows = []
    for method in methods:
        method_meta = METHOD_METADATA.get(method, {"name": method, "is_baseline": False})
        rows.append(
            {
                "Method": method,
                "Name": method_meta["name"],
                "IsBaseline": method_meta["is_baseline"],
                "Kind": "selector-family",
            }
        )
        for n_feature in sort_feature_levels(n_features):
            rows.append(
                {
                    "Method": f"{method}-N{n_feature}",
                    "Name": f"{method_meta['name']} (N={n_feature})",
                    "IsBaseline": method_meta["is_baseline"],
                    "Kind": "selector-setting",
                }
            )
    return pd.DataFrame(rows)


def build_integrations_metadata(integrations: list[str]) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "Integration": integrations,
            "Name": [INTEGRATION_METADATA.get(item, item) for item in integrations],
        }
    )


def adapt_results(results: pd.DataFrame) -> pd.DataFrame:
    metrics_meta = build_metrics_metadata().set_index("Metric")
    adapted = results.copy()
    adapted["Task"] = adapted["task"]
    adapted["Type"] = adapted["task"].map(TASK_TYPE_MAP)
    adapted["Metric"] = adapted["metric_name"]
    adapted["MethodBase"] = adapted["fs_method"]
    adapted["Method"] = adapted["fs_method"] + "-N" + adapted["n_features"].astype(str)
    adapted["SelFeatures"] = adapted["n_features"]
    adapted["Integration"] = adapted["integration_method"]
    adapted["ValueRaw"] = adapted["metric_value"]
    adapted["HigherBetter"] = adapted["Metric"].map(metrics_meta["HigherBetter"])
    adapted["MetricName"] = adapted["Metric"].map(metrics_meta["Name"])
    adapted["Value"] = adapted["ValueRaw"].where(adapted["HigherBetter"], -adapted["ValueRaw"])
    adapted["Seed"] = adapted["random_seed"]
    adapted["Notes"] = adapted["notes"].fillna("")
    adapted["Platform"] = adapted["platform"]
    adapted["NSlices"] = adapted["n_slices"]
    adapted["Runtime"] = adapted["runtime"]
    adapted["Dataset"] = adapted["dataset"]
    adapted["TaskLabel"] = adapted["Type"]
    adapted["IntegrationLabel"] = adapted["Integration"].map(lambda x: INTEGRATION_METADATA.get(x, x))
    adapted["EffectiveFeatures"] = adapted.get("effective_n_features", adapted["SelFeatures"])
    return adapted[
        [
            "Dataset",
            "Platform",
            "NSlices",
            "Method",
            "MethodBase",
            "SelFeatures",
            "Integration",
            "IntegrationLabel",
            "Task",
            "TaskLabel",
            "Type",
            "Metric",
            "MetricName",
            "Value",
            "ValueRaw",
            "HigherBetter",
            "Seed",
            "Runtime",
            "EffectiveFeatures",
            "Notes",
        ]
    ]


def compute_ranges(metrics: pd.DataFrame) -> pd.DataFrame:
    ranges = (
        metrics.groupby(["Dataset", "Metric", "Type"], dropna=False)["Value"]
        .agg(Lower="min", Upper="max")
        .reset_index()
    )
    ranges["Range"] = ranges["Upper"] - ranges["Lower"]
    ranges["Range"] = ranges["Range"].where(ranges["Range"] != 0, 1.0)
    return ranges


def main() -> None:
    args = parse_args()
    results_path = Path(args.results)
    output_dir = Path(args.output_dir)
    data_dir = output_dir / "data"
    output_subdir = output_dir / "output"
    data_dir.mkdir(parents=True, exist_ok=True)
    output_subdir.mkdir(parents=True, exist_ok=True)

    results = pd.read_csv(results_path)
    adapted = adapt_results(results)
    if args.slice_repr_datasets.strip():
        allowed = {item.strip() for item in args.slice_repr_datasets.split(",") if item.strip()}
        adapted = adapted[
            ~(
                (adapted["Type"] == "SliceRepresentation")
                & (~adapted["Dataset"].isin(allowed))
            )
        ].copy()

    adapted.to_csv(data_dir / "benchmark.tsv", sep="\t", index=False)
    adapted.to_csv(data_dir / "num-features.tsv", sep="\t", index=False)
    adapted.to_csv(data_dir / "spatial-benchmark-long.tsv", sep="\t", index=False)

    metric_ranges = compute_ranges(adapted)
    metric_ranges.to_csv(output_subdir / "baseline-ranges.tsv", sep="\t", index=False)

    datasets_meta = load_dataset_metadata(Path(args.configs_dir), sorted(adapted["Dataset"].unique().tolist()))
    datasets_meta.to_csv(data_dir / "datasets-metadata.tsv", sep="\t", index=False)

    build_metrics_metadata().to_csv(data_dir / "metrics-metadata.tsv", sep="\t", index=False)
    build_methods_metadata(
        sorted(adapted["MethodBase"].unique().tolist()),
        sort_feature_levels(adapted["SelFeatures"].dropna().astype(str).unique().tolist()),
    ).to_csv(data_dir / "methods-metadata.tsv", sep="\t", index=False)
    build_integrations_metadata(sorted(adapted["Integration"].unique().tolist())).to_csv(
        data_dir / "integrations-metadata.tsv", sep="\t", index=False
    )


if __name__ == "__main__":
    main()

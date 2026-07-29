from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from scipy.stats import rankdata, spearmanr


ROOT = Path(__file__).resolve().parents[1]
COMPONENT_METRICS = {
    "BatchMixing": ["bASW", "iLISI"],
    "BiologicalConservation": ["dASW", "dLISI", "ILL", "GC"],
    "LabelIndependentClustering": ["silhouette", "CHAOS", "PAS"],
    "LabelAgreement": ["ari", "nmi"],
    "Alignment": ["Accuracy", "Ratio"],
}
COMPETITIVE_GROUPS = {"expression_driven", "spatially_informed"}
LABELLED_DATASETS = {
    "DLPFC",
    "MouseBrainSerialSections",
    "STOmics0212",
    "STOmics0218",
    "STOmics0224",
}


METHOD_DETAILS = {
    "all_features": ("Control/reference", "All genes", "Project implementation", "NA"),
    "random": ("Control/reference", "Counts-derived gene universe", "Project implementation", "NA"),
    "TFs": ("Control/reference", "Predefined transcription-factor list", "Project resource", "NA"),
    "scsegindex": ("Control/reference", "Normalized expression", "scMerge 1.26.0", "lin2019scsegindex"),
    "scanpy_seurat": ("Expression-driven", "Log-normalized expression", "Scanpy 1.10.4", "wolf2018scanpy"),
    "scanpy_seurat_batch": ("Expression-driven", "Log-normalized expression + slice labels", "Scanpy 1.10.4", "wolf2018scanpy"),
    "scanpy_seurat_v3": ("Expression-driven", "Integer counts", "Scanpy 1.10.4", "stuart2019seurat"),
    "scanpy_seurat_v3_batch": ("Expression-driven", "Integer counts + slice labels", "Scanpy 1.10.4", "stuart2019seurat"),
    "scanpy_cell_ranger": ("Expression-driven", "Log-normalized expression", "Scanpy 1.10.4", "wolf2018scanpy"),
    "scanpy_cell_ranger_batch": ("Expression-driven", "Log-normalized expression + slice labels", "Scanpy 1.10.4", "wolf2018scanpy"),
    "scanpy_pearson": ("Expression-driven", "Integer counts", "Scanpy 1.10.4", "lause2021pearson"),
    "scanpy_pearson_batch": ("Expression-driven", "Integer counts + slice labels", "Scanpy 1.10.4", "lause2021pearson"),
    "seurat_vst": ("Expression-driven", "Counts", "Seurat 5.4.0", "stuart2019seurat"),
    "seurat_mvp": ("Expression-driven", "Counts", "Seurat 5.4.0", "stuart2019seurat"),
    "seurat_disp": ("Expression-driven", "Counts", "Seurat 5.4.0", "stuart2019seurat"),
    "seurat_sct": ("Expression-driven", "Counts", "Seurat 5.4.0; sctransform 0.4.3", "hafemeister2019sctransform"),
    "triku": ("Expression-driven", "Counts; expression-neighbor graph", "triku 2.1.4", "ascension2022triku"),
    "hotspot": ("Expression-driven", "Counts; PCA expression-neighbor graph (30-NN)", "Hotspot 1.0.0", "detomaso2021hotspot"),
    "nbumi": ("Expression-driven", "Integer counts", "M3Drop 1.36.0", "andrews2019m3drop"),
    "scry": ("Expression-driven", "Counts", "scry 1.22.0", "townes2019scry"),
    "osca": ("Expression-driven", "Counts + slice labels", "scran 1.38.0", "amezquita2020osca"),
    "dubstepr": ("Expression-driven", "Normalized expression", "DUBStepR 1.2.0", "ranjan2021dubstepr"),
    "anticor": ("Expression-driven", "Counts/expression", "anticor_features 0.1.8", "tyler2024anticor"),
    "Brennecke": ("Expression-driven", "Counts", "scran 1.38.0", "brennecke2013technical"),
    "scPNMF": ("Expression-driven", "Log-normalized expression", "scPNMF 1.0.1", "song2021scpnmf"),
    "singleCellHaystack": ("Expression-driven", "Detection matrix + PCA expression coordinates", "singleCellHaystack 0.3.4", "vandenbon2020singlecellhaystack"),
    "statistic_mean": ("Expression-driven", "Expression", "Project implementation", "NA"),
    "statistic_variance": ("Expression-driven", "Expression", "Project implementation", "NA"),
    "morans_i": ("Spatially informed", "Expression + spatial coordinates (8-NN)", "Project implementation", "moran1950notes"),
    "sparkx": ("Spatially informed", "Counts + spatial coordinates", "SPARK 1.1.1", "zhu2021sparkx"),
    "nnsvg": ("Spatially informed", "Counts + spatial coordinates", "nnSVG 1.14.0", "weber2023nnsvg"),
    "spatialde": ("Spatially informed", "Counts + spatial coordinates", "SpatialDE 1.1.3", "svensson2018spatialde"),
    "somde": ("Spatially informed", "Counts + spatial coordinates", "SOMDE 0.1.8", "hao2021somde"),
    "wilcoxon": ("Label-informed oracle", "Expression + benchmark labels", "Scanpy 1.10.4", "wilcoxon1945individual"),
}

CITATION_LABELS = {
    "NA": "Not applicable",
    "lin2019scsegindex": "Lin et al. (2019)",
    "wolf2018scanpy": "Wolf et al. (2018)",
    "stuart2019seurat": "Stuart et al. (2019)",
    "lause2021pearson": "Lause et al. (2021)",
    "hafemeister2019sctransform": "Hafemeister and Satija (2019)",
    "ascension2022triku": "Ascensi\u00f3n et al. (2022)",
    "detomaso2021hotspot": "DeTomaso and Yosef (2021)",
    "andrews2019m3drop": "Andrews and Hemberg (2019)",
    "townes2019scry": "Townes et al. (2019)",
    "amezquita2020osca": "Amezquita et al. (2020)",
    "ranjan2021dubstepr": "Ranjan et al. (2021)",
    "tyler2024anticor": "Tyler et al. (2024)",
    "brennecke2013technical": "Brennecke et al. (2013)",
    "song2021scpnmf": "Song et al. (2021)",
    "vandenbon2020singlecellhaystack": "Vandenbon and Diez (2020)",
    "moran1950notes": "Moran (1950)",
    "zhu2021sparkx": "Zhu et al. (2021)",
    "weber2023nnsvg": "Weber et al. (2023)",
    "svensson2018spatialde": "Svensson et al. (2018)",
    "hao2021somde": "Hao et al. (2021)",
    "wilcoxon1945individual": "Wilcoxon (1945)",
}

METHOD_LABELS = {
    "all_features": "All features",
    "random": "Random",
    "TFs": "Transcription factors",
    "scsegindex": "scSEGIndex",
    "scanpy_seurat": "Scanpy Seurat",
    "scanpy_seurat_batch": "Scanpy Seurat, batch",
    "scanpy_seurat_v3": "Scanpy Seurat v3",
    "scanpy_seurat_v3_batch": "Scanpy Seurat v3, batch",
    "scanpy_cell_ranger": "Scanpy Cell Ranger",
    "scanpy_cell_ranger_batch": "Scanpy Cell Ranger, batch",
    "scanpy_pearson": "Scanpy Pearson",
    "scanpy_pearson_batch": "Scanpy Pearson, batch",
    "seurat_vst": "Seurat VST",
    "seurat_mvp": "Seurat MVP",
    "seurat_disp": "Seurat Dispersion",
    "seurat_sct": "Seurat sctransform",
    "triku": "triku",
    "hotspot": "Hotspot",
    "nbumi": "NBumi",
    "scry": "scry",
    "osca": "OSCA",
    "dubstepr": "DUBStepR",
    "anticor": "Anticor",
    "Brennecke": "Brennecke",
    "scPNMF": "scPNMF",
    "singleCellHaystack": "singleCellHaystack",
    "statistic_mean": "Mean expression",
    "statistic_variance": "Expression variance",
    "morans_i": "Moran's I",
    "sparkx": "SPARK-X",
    "nnsvg": "nnSVG",
    "spatialde": "SpatialDE",
    "somde": "SOMDE",
    "wilcoxon": "Wilcoxon",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--frozen-dir",
        type=Path,
        default=ROOT / "results/spatial_svg_rebuild_v1/frozen_scores",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=ROOT / "results/submission_robustness_v1",
    )
    parser.add_argument("--bootstrap", type=int, default=5000)
    parser.add_argument("--weight-draws", type=int, default=2000)
    parser.add_argument("--seed", type=int, default=20260729)
    return parser.parse_args()


def normalized_n(value: object) -> str:
    text = str(value)
    if text.lower() == "all":
        return "all"
    return str(int(float(text)))


def strict_mean(frame: pd.DataFrame, columns: list[str]) -> pd.Series:
    values = frame.reindex(columns=columns)
    return values.mean(axis=1).where(values.notna().all(axis=1))


def weighted_available(a: pd.Series, b: pd.Series, weight_a: float) -> pd.Series:
    values = pd.concat([a, b], axis=1)
    weights = np.array([weight_a, 1.0 - weight_a])
    present = values.notna().to_numpy()
    numerator = np.nansum(values.to_numpy() * weights, axis=1)
    denominator = (present * weights).sum(axis=1)
    return pd.Series(np.divide(numerator, denominator, where=denominator > 0), index=values.index).where(
        denominator > 0
    )


def task_scores(
    metric_frame: pd.DataFrame,
    value: str,
    component_weights: tuple[float, float, float] = (0.5, 0.5, 0.5),
) -> pd.DataFrame:
    keys = ["dataset", "fs_method", "n_features", "integration_method"]
    wide = metric_frame.pivot(index=keys, columns="metric_name", values=value).reset_index()
    wide.columns.name = None
    for component, metrics in COMPONENT_METRICS.items():
        wide[component] = strict_mean(wide, metrics)
    integration_weight, clustering_weight, core_weight = component_weights
    wide["Integration"] = weighted_available(
        wide["BatchMixing"], wide["BiologicalConservation"], integration_weight
    )
    wide["Clustering"] = weighted_available(
        wide["LabelIndependentClustering"], wide["LabelAgreement"], clustering_weight
    )
    wide["CoreOverall"] = weighted_available(wide["Integration"], wide["Clustering"], core_weight)
    return wide


def representative_metrics(frozen_dir: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    metrics = pd.read_csv(
        frozen_dir / "scaled_setting_metrics.tsv",
        sep="\t",
        dtype={"n_features": str},
        low_memory=False,
    )
    reps = pd.read_csv(frozen_dir / "representative_settings.tsv", sep="\t")
    metrics["n_features"] = metrics["n_features"].map(normalized_n)
    reps["n_features"] = reps["n_features"].map(normalized_n)
    keys = ["dataset", "fs_method", "n_features"]
    selected = metrics.merge(reps[keys + ["MethodGroup"]], on=keys, how="inner", validate="many_to_one")
    return metrics, selected


def add_alternative_scaling(all_metrics: pd.DataFrame) -> pd.DataFrame:
    grouping = ["dataset", "task", "metric_name"]
    frame = all_metrics.copy()
    frame["RankScaledValue"] = frame.groupby(grouping, observed=True)["OrientedMean"].transform(
        lambda x: (x.rank(method="average") - 1) / max(x.notna().sum() - 1, 1)
    )
    frame["ZScaledValue"] = frame.groupby(grouping, observed=True)["OrientedMean"].transform(
        lambda x: (x - x.mean()) / x.std(ddof=0) if x.std(ddof=0) > 0 else np.nan
    )
    return frame


def global_ranks(scores: pd.DataFrame, scenario: str, labelled_only: bool = False) -> pd.DataFrame:
    frame = scores.copy()
    if labelled_only:
        frame = frame[frame["dataset"].isin(LABELLED_DATASETS)]
    frame["DatasetRank"] = frame.groupby(["dataset", "integration_method"], observed=True)[
        "CoreOverall"
    ].rank(method="average", ascending=False)
    result = (
        frame.groupby(["fs_method", "integration_method"], observed=True)
        .agg(MeanDatasetRank=("DatasetRank", "mean"), MeanScore=("CoreOverall", "mean"), NDatasets=("dataset", "nunique"))
        .reset_index()
    )
    result["GlobalRank"] = result.groupby("integration_method", observed=True)["MeanDatasetRank"].rank(
        method="average"
    )
    result["Scenario"] = scenario
    return result


def kendalls_w(matrix: np.ndarray) -> float:
    matrix = np.asarray(matrix, dtype=float)
    m, n = matrix.shape
    sums = matrix.sum(axis=0)
    s_value = np.square(sums - sums.mean()).sum()
    tie_total = 0.0
    for row in matrix:
        _, counts = np.unique(row, return_counts=True)
        tie_total += np.sum(counts**3 - counts)
    denominator = m * m * (n**3 - n) - m * tie_total
    return float(12 * s_value / denominator) if denominator > 0 else float("nan")


def concordance_table(ranks: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for integrator, part in ranks.groupby("integration_method", observed=True):
        wide = part.pivot(index="Scenario", columns="fs_method", values="GlobalRank").dropna(axis=1)
        rows.append(
            {
                "integration_method": integrator,
                "NScenarios": wide.shape[0],
                "NMethods": wide.shape[1],
                "KendallsW": kendalls_w(wide.to_numpy()),
            }
        )
        for first_idx, first in enumerate(wide.index):
            for second in wide.index[first_idx + 1 :]:
                rho = spearmanr(wide.loc[first], wide.loc[second]).statistic
                rows.append(
                    {
                        "integration_method": integrator,
                        "NScenarios": 2,
                        "NMethods": wide.shape[1],
                        "KendallsW": np.nan,
                        "Scenario1": first,
                        "Scenario2": second,
                        "SpearmanRho": rho,
                    }
                )
    return pd.DataFrame(rows)


def bootstrap_ranks(
    base_scores: pd.DataFrame, n_bootstrap: int, rng: np.random.Generator
) -> pd.DataFrame:
    rows = []
    for integrator, frame in base_scores.groupby("integration_method", observed=True):
        frame = frame.copy()
        frame["DatasetRank"] = frame.groupby("dataset", observed=True)["CoreOverall"].rank(
            method="average", ascending=False
        )
        datasets = sorted(frame["dataset"].unique())
        pivot = frame.pivot(index="dataset", columns="fs_method", values="DatasetRank").loc[datasets]
        draws = np.empty((n_bootstrap, pivot.shape[1]), dtype=float)
        values = pivot.to_numpy()
        for index in range(n_bootstrap):
            sampled = rng.integers(0, len(datasets), size=len(datasets))
            mean_dataset_rank = np.nanmean(values[sampled], axis=0)
            draws[index] = rankdata(mean_dataset_rank, method="average")
        observed_mean = pivot.mean(axis=0)
        observed_global = pd.Series(rankdata(observed_mean, method="average"), index=pivot.columns)
        for method_idx, method in enumerate(pivot.columns):
            rows.append(
                {
                    "integration_method": integrator,
                    "fs_method": method,
                    "NDataSets": len(datasets),
                    "ObservedMeanDatasetRank": observed_mean[method],
                    "ObservedGlobalRank": observed_global[method],
                    "BootstrapGlobalRankMedian": np.median(draws[:, method_idx]),
                    "BootstrapGlobalRankLowerQuartile": np.quantile(draws[:, method_idx], 0.25),
                    "BootstrapGlobalRankUpperQuartile": np.quantile(draws[:, method_idx], 0.75),
                    "BootstrapGlobalRankLower95": np.quantile(draws[:, method_idx], 0.025),
                    "BootstrapGlobalRankUpper95": np.quantile(draws[:, method_idx], 0.975),
                    "Top5Probability": np.mean(draws[:, method_idx] <= 5),
                }
            )
    return pd.DataFrame(rows)


def weight_perturbation(
    metrics: pd.DataFrame, n_draws: int, rng: np.random.Generator
) -> tuple[pd.DataFrame, pd.DataFrame]:
    draws = []
    for draw in range(n_draws):
        weights = tuple(rng.uniform(0.4, 0.6, size=3))
        score = task_scores(metrics, "ScaledValue", weights)
        rank = global_ranks(score, f"weight_draw_{draw}")
        rank["Draw"] = draw
        rank["BatchMixingWeight"] = weights[0]
        rank["LabelIndependentWeight"] = weights[1]
        rank["IntegrationWeight"] = weights[2]
        draws.append(rank)
    all_draws = pd.concat(draws, ignore_index=True)
    summary = (
        all_draws.groupby(["fs_method", "integration_method"], observed=True)
        .agg(
            MedianRank=("GlobalRank", "median"),
            Lower95Rank=("GlobalRank", lambda x: x.quantile(0.025)),
            Upper95Rank=("GlobalRank", lambda x: x.quantile(0.975)),
            Top5Probability=("GlobalRank", lambda x: (x <= 5).mean()),
        )
        .reset_index()
    )
    return all_draws, summary


def bootstrap_spearman(
    x: np.ndarray,
    y: np.ndarray,
    clusters: np.ndarray | None,
    n_bootstrap: int,
    rng: np.random.Generator,
) -> tuple[int, float, float, float]:
    valid = np.isfinite(x) & np.isfinite(y)
    x, y = x[valid], y[valid]
    clusters = clusters[valid] if clusters is not None else None
    observed = float(spearmanr(x, y).statistic)
    estimates = []
    if clusters is None:
        for _ in range(n_bootstrap):
            idx = rng.integers(0, len(x), size=len(x))
            if np.unique(x[idx]).size > 1 and np.unique(y[idx]).size > 1:
                estimates.append(spearmanr(x[idx], y[idx]).statistic)
    else:
        unique = np.unique(clusters)
        for _ in range(n_bootstrap):
            sampled = rng.choice(unique, size=len(unique), replace=True)
            idx = np.concatenate([np.flatnonzero(clusters == value) for value in sampled])
            if np.unique(x[idx]).size > 1 and np.unique(y[idx]).size > 1:
                estimates.append(spearmanr(x[idx], y[idx]).statistic)
    return len(x), observed, float(np.quantile(estimates, 0.025)), float(np.quantile(estimates, 0.975))


def correlation_summary(n_bootstrap: int, rng: np.random.Generator) -> pd.DataFrame:
    rows = []
    held_path = ROOT / "manuscript_genome_research_spatial_omics_v2/source_data/Supplemental_Fig_S2d_association.tsv"
    held = pd.read_csv(held_path, sep="\t")
    fold = (held["dataset"].astype(str) + "::" + held["heldout_slice"].astype(str)).to_numpy()
    result = bootstrap_spearman(
        held["JaccardPercentile"].to_numpy(float),
        held["ARIPercentile"].to_numpy(float),
        fold,
        n_bootstrap,
        rng,
    )
    rows.append(dict(zip(["Analysis", "N", "SpearmanRho", "Lower95", "Upper95"], ["Held-out Jaccard percentile vs ARI percentile", *result], strict=True)))

    ranks_path = ROOT / "manuscript_genome_research_spatial_omics_v2/source_data/Supplemental_Fig_S5c_global_rank_agreement.tsv"
    ranks = pd.read_csv(ranks_path, sep="\t").dropna(subset=["scvi", "cellcharter"])
    result = bootstrap_spearman(
        ranks["scvi"].to_numpy(float),
        ranks["cellcharter"].to_numpy(float),
        None,
        n_bootstrap,
        rng,
    )
    rows.append(dict(zip(["Analysis", "N", "SpearmanRho", "Lower95", "Upper95"], ["scVI vs CellCharter global mean ranks", *result], strict=True)))
    return pd.DataFrame(rows)


def runtime_summary(frozen_dir: Path) -> pd.DataFrame:
    metrics = pd.read_csv(frozen_dir / "scaled_seed_metrics.tsv", sep="\t")
    keys = ["dataset", "fs_method", "n_features", "integration_method", "task", "random_seed"]
    records = metrics.drop_duplicates(keys).copy()
    records["MemoryGiB"] = records["memory_usage"] / 1024
    summary = (
        records.groupby(["fs_method", "integration_method", "task"], observed=True)
        .agg(
            NRuns=("runtime", "size"),
            RuntimeMedianSeconds=("runtime", "median"),
            RuntimeIQRSeconds=("runtime", lambda x: x.quantile(0.75) - x.quantile(0.25)),
            RuntimeMaximumSeconds=("runtime", "max"),
            PeakMemoryMedianGiB=("MemoryGiB", "median"),
            PeakMemoryMaximumGiB=("MemoryGiB", "max"),
        )
        .reset_index()
    )
    return summary


def methods_inventory(reps: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for method in sorted(reps["fs_method"].unique()):
        category, inputs, version, citation = METHOD_DETAILS[method]
        settings = reps.loc[reps["fs_method"].eq(method), "n_features"].map(normalized_n).unique()
        rows.append(
            {
                "Method": method,
                "DisplayMethod": METHOD_LABELS[method],
                "Category": category,
                "UsesSpatialCoordinates": category == "Spatially informed",
                "Inputs": inputs,
                "RepresentativeFeatureCount": ",".join(sorted(settings)),
                "SoftwareVersion": version,
                "CitationKey": citation,
                "PrimaryReference": CITATION_LABELS[citation],
            }
        )
    return pd.DataFrame(rows)


def latex_escape(value: object) -> str:
    text = str(value)
    replacements = {
        "\\": r"\textbackslash{}",
        "&": r"\&",
        "%": r"\%",
        "$": r"\$",
        "#": r"\#",
        "_": r"\_",
        "{": r"\{",
        "}": r"\}",
    }
    for old, new in replacements.items():
        text = text.replace(old, new)
    return text


def write_latex_tables(output_dir: Path, methods: pd.DataFrame, datasets: pd.DataFrame) -> None:
    table_dir = ROOT / "manuscript_genome_research_spatial_omics_v2/supplemental_tables"
    table_dir.mkdir(parents=True, exist_ok=True)
    method_lines = [
        r"\begin{landscape}",
        r"\scriptsize",
        r"\setlength{\tabcolsep}{3pt}",
        r"\begin{longtable}{p{0.16\linewidth}p{0.15\linewidth}p{0.29\linewidth}p{0.18\linewidth}p{0.16\linewidth}}",
        r"\caption{\textbf{Complete representative feature-selection panel.} Coordinate use refers to the ranking step itself; downstream integrators may separately use spatial coordinates.}\label{tab:methods}\\",
        r"\toprule Method & Classification & Ranking input & Software/version & Primary reference \\",
        r"\midrule\endfirsthead",
        r"\toprule Method & Classification & Ranking input & Software/version & Primary reference \\",
        r"\midrule\endhead",
    ]
    for row in methods.sort_values(["Category", "Method"]).itertuples(index=False):
        method_lines.append(
            " & ".join(
                latex_escape(value)
                for value in [row.DisplayMethod, row.Category, row.Inputs, row.SoftwareVersion, row.PrimaryReference]
            )
            + r" \\" 
        )
    method_lines.extend([r"\bottomrule", r"\end{longtable}", r"\end{landscape}"])
    (table_dir / "Supplemental_Table_S1.tex").write_text("\n".join(method_lines) + "\n")

    dataset_lines = [
        r"\begin{landscape}",
        r"\scriptsize",
        r"\setlength{\tabcolsep}{2.5pt}",
        r"\begin{longtable}{p{0.16\linewidth}p{0.07\linewidth}p{0.08\linewidth}p{0.05\linewidth}p{0.07\linewidth}p{0.07\linewidth}p{0.23\linewidth}p{0.17\linewidth}}",
        r"\caption{\textbf{Datasets and benchmark-level quality control.} Spot and gene counts are those entering the rebuilt count-dependent workflow.}\label{tab:datasets}\\",
        r"\toprule Dataset & Species & Platform & Slices & Spots & Genes & Spot-level QC & Gene-level QC \\",
        r"\midrule\endfirsthead",
        r"\toprule Dataset & Species & Platform & Slices & Spots & Genes & Spot-level QC & Gene-level QC \\",
        r"\midrule\endhead",
    ]
    for row in datasets.itertuples(index=False):
        dataset_lines.append(
            " & ".join(
                latex_escape(value)
                for value in [
                    row.Dataset,
                    row.Species,
                    row.Platform,
                    row.Slices,
                    row.Spots,
                    row.GenesAfterBenchmarkFiltering,
                    row.SpotQC,
                    row.GeneQC,
                ]
            )
            + r" \\" 
        )
    dataset_lines.extend(
        [
            r"\bottomrule",
            r"\end{longtable}",
            r"\noindent No Stereo-seq dataset entered the rebuilt benchmark; a Stereo-seq bin size is therefore not applicable. Mitochondrial genes were not removed globally and were filtered only within method-specific official workflows where required.",
            r"\end{landscape}",
        ]
    )
    (table_dir / "Supplemental_Table_S2.tex").write_text("\n".join(dataset_lines) + "\n")


def dataset_inventory() -> pd.DataFrame:
    config_paths = [
        "configs/datasets/dlpfc.yaml",
        "configs/datasets/mouse_brain_serial_sections.yaml",
        "configs/datasets/stomics_0212_wilcoxon.yaml",
        "configs/datasets/stomics_0218_wilcoxon.yaml",
        "configs/datasets/stomics_0224_wilcoxon.yaml",
        "configs/datasets/e8p5_embryo.yaml",
        "configs/datasets/e9p5_embryo.yaml",
    ]
    rows = []
    for relative in config_paths:
        config = yaml.safe_load((ROOT / relative).read_text())
        dataset = config["name"]
        metadata_candidates = list(
            (ROOT / "results/spatial_svg_rebuild_v1").glob(
                "**/feature_selection/**/selected_features.meta.json"
            )
        )
        found = None
        for candidate in metadata_candidates:
            payload = json.loads(candidate.read_text())
            if payload.get("dataset", {}).get("dataset_name") == dataset:
                found = payload["dataset"]
                break
        if found is None:
            raise RuntimeError(f"No run metadata found for {dataset}")
        preprocessing = config.get("preprocess", {})
        min_cells = int(preprocessing.get("min_cells_per_gene", 0) or 0)
        rows.append(
            {
                "Dataset": dataset,
                "Species": config.get("species"),
                "Platform": config.get("platform"),
                "Slices": np.nan,
                "Spots": found["n_obs"],
                "GenesAfterBenchmarkFiltering": found["n_vars"],
                "LabelsAvailable": bool(config.get("label_key")),
                "AlignmentPairs": len(config.get("alignment_pairs", [])),
                "SpotQC": "source-provided tissue spots; none additional",
                "GeneQC": f"detected in >= {min_cells} spots" if min_cells else "none additional",
                "MitochondrialGeneFilter": "method-specific only (SPARK-X/nnSVG)",
                "StereoSeqBinSize": "NA (no Stereo-seq dataset in rebuilt benchmark)",
                "Config": relative,
            }
        )
    slice_counts = {
        "DLPFC": 4,
        "MouseBrainSerialSections": 2,
        "STOmics0212": 10,
        "STOmics0218": 7,
        "STOmics0224": 8,
        "E8p5Embryo": 6,
        "E9p5Embryo": 6,
    }
    frame = pd.DataFrame(rows)
    frame["Slices"] = frame["Dataset"].map(slice_counts)
    return frame


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(args.seed)

    all_metrics, representative = representative_metrics(args.frozen_dir)
    all_metrics = add_alternative_scaling(all_metrics)
    keys = ["dataset", "fs_method", "n_features"]
    alternative_representative = all_metrics.merge(
        representative[keys + ["MethodGroup"]].drop_duplicates(),
        on=keys,
        how="inner",
        validate="many_to_one",
    )
    competitive = alternative_representative[
        alternative_representative["MethodGroup"].isin(COMPETITIVE_GROUPS)
    ].copy()

    scenarios = []
    observed_scores = task_scores(competitive, "ScaledValue")
    scenarios.append(global_ranks(observed_scores, "observed_range"))
    scenarios.append(global_ranks(task_scores(competitive, "RankScaledValue"), "rank_percentile"))
    scenarios.append(global_ranks(task_scores(competitive, "ZScaledValue"), "z_score"))
    scenarios.append(global_ranks(task_scores(competitive, "ScaledValue", (0.5, 0.5, 0.6)), "integration_weight_0.6"))
    scenarios.append(global_ranks(task_scores(competitive, "ScaledValue", (0.5, 0.5, 0.4)), "integration_weight_0.4"))
    scenarios.append(global_ranks(observed_scores, "labelled_datasets_only", labelled_only=True))
    scenario_ranks = pd.concat(scenarios, ignore_index=True)

    weight_draws, weight_summary = weight_perturbation(competitive, args.weight_draws, rng)
    bootstrap = bootstrap_ranks(observed_scores, args.bootstrap, rng)
    correlations = correlation_summary(args.bootstrap, rng)

    method_table = methods_inventory(representative)
    dataset_table = dataset_inventory()
    outputs = {
        "methods_inventory.tsv": method_table,
        "datasets_inventory.tsv": dataset_table,
        "scaling_weight_sensitivity_ranks.tsv": scenario_ranks,
        "scaling_weight_concordance.tsv": concordance_table(scenario_ranks),
        "weight_perturbation_summary.tsv": weight_summary,
        "dataset_bootstrap_rank_intervals.tsv": bootstrap,
        "spearman_correlations_with_ci.tsv": correlations,
        "runtime_memory_summary.tsv": runtime_summary(args.frozen_dir),
    }
    for name, frame in outputs.items():
        frame.to_csv(args.output_dir / name, sep="\t", index=False)
    write_latex_tables(args.output_dir, method_table, dataset_table)
    weight_draws.to_csv(args.output_dir / "weight_perturbation_draws.tsv.gz", sep="\t", index=False)
    manifest = {
        "input": str(args.frozen_dir),
        "bootstrap_draws": args.bootstrap,
        "weight_draws": args.weight_draws,
        "random_seed": args.seed,
        "labelled_datasets": sorted(LABELLED_DATASETS),
        "method_count": int(representative["fs_method"].nunique()),
        "competitive_method_count": int(competitive["fs_method"].nunique()),
        "outputs": {name: len(frame) for name, frame in outputs.items()},
    }
    (args.output_dir / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")


if __name__ == "__main__":
    main()

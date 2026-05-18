from __future__ import annotations

import argparse
from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml
from matplotlib.colors import ListedColormap
from sklearn.cluster import KMeans
from umap import UMAP


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Generate a spatial case study figure.")
    parser.add_argument("--results", default="results/dlpfc_mvp/results.csv")
    parser.add_argument("--results-dir", default="results/dlpfc_mvp")
    parser.add_argument("--dataset", default="DLPFC")
    parser.add_argument("--dataset-config", default="configs/datasets/dlpfc.yaml")
    parser.add_argument("--output", default="results/atlas_style/figures/figure5_case_study.png")
    return parser.parse_args()


def choose_methods(results: pd.DataFrame, dataset_name: str) -> tuple[str, int, list[str]]:
    metric_types = {
        "dASW", "dLISI", "ILL", "bASW", "iLISI", "GC",
        "ari", "nmi", "CHAOS", "PAS",
        "Accuracy", "Ratio",
    }
    direction = {"CHAOS": -1, "PAS": -1, "Ratio": -1}
    subset = results[results["dataset"] == dataset_name].copy()
    subset = subset[subset["metric_name"].isin(metric_types)]
    subset["score"] = subset["metric_value"] * subset["metric_name"].map(lambda x: direction.get(x, 1))
    aggregated = (
        subset.groupby(["integration_method", "n_features", "fs_method"], as_index=False)["score"]
        .mean()
    )
    best_combo = (
        aggregated.groupby(["integration_method", "n_features"], as_index=False)["score"]
        .mean()
        .sort_values("score", ascending=False)
        .iloc[0]
    )
    pair = (
        aggregated[
            (aggregated["integration_method"] == best_combo["integration_method"])
            & (aggregated["n_features"] == best_combo["n_features"])
        ]
        .sort_values("score", ascending=False)
    )
    methods = [pair.iloc[0]["fs_method"], pair.iloc[-1]["fs_method"]]
    return str(best_combo["integration_method"]), int(best_combo["n_features"]), methods


def load_embedding(base_dir: Path, dataset_name: str, integration: str, fs_method: str, n_features: int) -> np.ndarray:
    path = base_dir / dataset_name.lower() / integration / fs_method / f"n{n_features}" / "seed0" / "embedding.npz"
    data = np.load(path)
    return data["embedding"]


def format_method_name(name: str) -> str:
    mapping = {
        "hvg": "HVG",
        "svg": "SVG",
        "highly_expressed": "Highly expressed",
        "random": "Random",
        "TFs": "Transcription factors",
        "scanpy_seurat": "scanpy-Seurat",
        "scanpy_seurat_batch": "scanpy-Seurat (batch)",
        "scanpy_seurat_v3": "scanpy-SeuratV3",
        "scanpy_seurat_v3_batch": "scanpy-SeuratV3 (batch)",
        "scanpy_cell_ranger": "scanpy-CellRanger",
        "scanpy_cell_ranger_batch": "scanpy-CellRanger (batch)",
        "scanpy_pearson": "scanpy-Pearson",
        "scanpy_pearson_batch": "scanpy-Pearson (batch)",
        "seurat_vst": "Seurat-VST",
        "seurat_mvp": "Seurat-MVP",
        "seurat_disp": "Seurat-Dispersion",
        "seurat_sct": "Seurat-scTransform",
        "scsegindex": "scSEGIndex",
        "dubstepr": "DUBStepR",
        "nbumi": "NBumi",
        "osca": "OSCA",
        "scry": "scry",
        "singleCellHaystack": "singleCellHaystack",
        "Brennecke": "Brennecke",
        "scPNMF": "scPNMF",
        "statistic_mean": "Statistic mean",
        "statistic_variance": "Statistic variance",
        "triku": "triku",
        "hotspot": "Hotspot",
        "anticor": "Anticor",
        "wilcoxon": "Wilcoxon",
        "all": "All features",
    }
    return mapping.get(name, name)


def make_spatial_panel(axs, coords, labels, slices) -> None:
    unique_slices = np.unique(slices)
    label_codes = pd.Categorical(labels).codes
    cmap = ListedColormap(plt.cm.tab20(np.linspace(0, 1, max(20, len(np.unique(label_codes))))))
    for idx, slice_id in enumerate(unique_slices):
        ax = axs[idx]
        mask = slices == slice_id
        ax.scatter(coords[mask, 0], coords[mask, 1], c=label_codes[mask], cmap=cmap, s=4)
        ax.set_title(str(slice_id), fontsize=8, pad=3)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.invert_yaxis()
        for spine in ax.spines.values():
            spine.set_visible(False)


def main() -> None:
    args = parse_args()
    results = pd.read_csv(args.results)
    integration, n_features, fs_methods = choose_methods(results, args.dataset)
    with Path(args.dataset_config).open("r", encoding="utf-8") as handle:
        dataset_config = yaml.safe_load(handle)
    adata = ad.read_h5ad(dataset_config["path"])
    coords = np.asarray(adata.obsm["spatial"])
    slices = adata.obs[dataset_config["slice_key"]].astype(str).to_numpy()
    if dataset_config.get("label_key") is None:
        raise ValueError("Case study requires a dataset config with label_key.")
    truth = adata.obs[dataset_config["label_key"]].astype(str).to_numpy()
    n_clusters = len(np.unique(truth))
    unique_slices = np.unique(slices)

    fig = plt.figure(figsize=(12.5, 9.5), constrained_layout=True)
    outer = fig.add_gridspec(3, 2, height_ratios=[1.0, 1.0, 1.55], hspace=0.08, wspace=0.08)

    for col, fs_method in enumerate(fs_methods):
        embedding = load_embedding(Path(args.results_dir), args.dataset, integration, fs_method, n_features)
        embedding_2d = UMAP(random_state=0).fit_transform(embedding)
        pred = KMeans(n_clusters=n_clusters, n_init=20, random_state=0).fit_predict(embedding)

        ax1 = fig.add_subplot(outer[0, col])
        truth_codes = pd.Categorical(truth).codes
        ax1.scatter(embedding_2d[:, 0], embedding_2d[:, 1], c=truth_codes, cmap="tab20", s=4)
        ax1.set_title(
            f"{format_method_name(fs_method)}\n{integration.upper()} | N={n_features}",
            fontsize=10,
            pad=6,
            fontweight="bold",
        )
        ax1.set_xticks([])
        ax1.set_yticks([])
        for spine in ax1.spines.values():
            spine.set_visible(False)
        if col == 0:
            ax1.set_ylabel("Embedding\nby label", fontsize=10)

        ax2 = fig.add_subplot(outer[1, col])
        slice_codes = pd.Categorical(slices).codes
        ax2.scatter(embedding_2d[:, 0], embedding_2d[:, 1], c=slice_codes, cmap="Set2", s=4)
        ax2.set_xticks([])
        ax2.set_yticks([])
        for spine in ax2.spines.values():
            spine.set_visible(False)
        if col == 0:
            ax2.set_ylabel("Embedding\nby slice", fontsize=10)

        subgrid = outer[2, col].subgridspec(1, len(unique_slices), wspace=0.03)
        sub_axes = [fig.add_subplot(subgrid[0, idx]) for idx in range(len(unique_slices))]
        make_spatial_panel(sub_axes, coords, pred.astype(str), slices)
        if col == 0:
            sub_axes[0].set_ylabel("Spatial\ndomains", fontsize=10)

    fig.suptitle(
        f"Figure 5. Spatial case study on {args.dataset}: best vs. worst feature selection",
        fontsize=12,
        fontweight="bold",
    )
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=300, bbox_inches="tight")
    fig.savefig(output.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()

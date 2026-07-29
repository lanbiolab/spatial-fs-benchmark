from __future__ import annotations

from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / "results/resource_profile_v1/feature_selection_resource_profile.tsv"
OUTPUT = ROOT / "manuscript_genome_research_spatial_omics_v2/supplemental_tables/Supplemental_Table_S3.tex"

LABELS = {
    "all_features": "All features",
    "random": "Random",
    "TFs": "Transcription factors",
    "scsegindex": "scSEGIndex",
    "scanpy_cell_ranger": "Scanpy Cell Ranger",
    "scanpy_cell_ranger_batch": "Scanpy Cell Ranger, batch",
    "scanpy_pearson": "Scanpy Pearson",
    "scanpy_pearson_batch": "Scanpy Pearson, batch",
    "scanpy_seurat": "Scanpy Seurat",
    "scanpy_seurat_batch": "Scanpy Seurat, batch",
    "scanpy_seurat_v3": "Scanpy Seurat v3",
    "scanpy_seurat_v3_batch": "Scanpy Seurat v3, batch",
    "seurat_disp": "Seurat Dispersion",
    "seurat_mvp": "Seurat MVP",
    "seurat_sct": "Seurat sctransform",
    "seurat_vst": "Seurat VST",
    "singleCellHaystack": "singleCellHaystack",
    "statistic_mean": "Statistic mean",
    "statistic_variance": "Statistic variance",
    "morans_i": "Moran's I",
    "sparkx": "SPARK-X",
    "spatialde": "SpatialDE",
    "nnsvg": "nnSVG",
    "somde": "SOMDE",
    "scPNMF": "scPNMF",
    "dubstepr": "DUBStepR",
    "nbumi": "NBumi",
    "anticor": "Anticor",
    "hotspot": "Hotspot",
    "osca": "OSCA",
}


def escape(value: object) -> str:
    return str(value).replace("_", r"\_").replace("&", r"\&").replace("%", r"\%")


def main() -> None:
    frame = pd.read_csv(INPUT, sep="\t")
    frame["Display"] = frame["method"].map(LABELS).fillna(frame["method"])
    frame = frame.sort_values("Display", key=lambda values: values.str.lower())
    lines = [
        r"\clearpage",
        r"\section*{Supplemental Table S3}",
        r"\small",
        r"\begin{longtable}{p{0.27\textwidth}rrrrr}",
        r"\caption{\textbf{Standardized feature-selection computational profile.} Each method was run in an isolated wrapper process on 2,000 Mouse Brain spots (1,000 per slice) and the 5,000 most frequently detected genes. Process wall time includes data loading and preprocessing; selector time measures only the selection call. Peak resident memory is the maximum RSS reported by GNU time for the isolated wrapper command. These values characterize one fixed profile and are not full-pipeline runtimes.}\\",
        r"\toprule",
        r"Method & Requested & Returned & Selector (s) & Process (s) & Peak RSS (MiB) \\",
        r"\midrule",
        r"\endfirsthead",
        r"\toprule",
        r"Method & Requested & Returned & Selector (s) & Process (s) & Peak RSS (MiB) \\",
        r"\midrule",
        r"\endhead",
    ]
    for row in frame.itertuples(index=False):
        lines.append(
            f"{escape(row.Display)} & {int(row.requested_n_features):,} & "
            f"{int(row.effective_n_features):,} & {row.wall_seconds_internal:.2f} & "
            f"{row.wall_seconds_process:.2f} & {row.peak_rss_mib:.0f} \\\\"
        )
    lines.extend([r"\bottomrule", r"\end{longtable}", r"\normalsize", ""])
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text("\n".join(lines), encoding="utf-8")


if __name__ == "__main__":
    main()

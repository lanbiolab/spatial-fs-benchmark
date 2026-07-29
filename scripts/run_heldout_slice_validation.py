#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import pandas as pd

from spatial_fs_benchmark.config import load_dataset_config
from spatial_fs_benchmark.data.io import load_dataset
from spatial_fs_benchmark.data.preprocess import preprocess_dataset
from spatial_fs_benchmark.feature_selection import build_feature_selector
from spatial_fs_benchmark.feature_selection.base import FeatureSelectionResult
from spatial_fs_benchmark.validation import (
    HeldoutFold,
    heldout_metrics_frame,
    heldout_moran_reference,
    score_heldout_clustering,
    score_spatial_reproducibility,
    split_heldout_slice,
)


METHODS = {
    "all_features": {},
    "random": {},
    "scanpy_seurat_v3_batch": {},
    "scanpy_cell_ranger_batch": {},
    "seurat_vst": {},
    "triku": {"max_cells": 5000},
    "singleCellHaystack": {"max_cells": 5000},
    "scsegindex": {},
    "morans_i": {"n_neighbors": 8},
    "sparkx": {"n_threads": 4},
    "nnsvg": {"n_threads": 4},
    "spatialde": {"max_cells_per_slice": 800},
    "somde": {"spots_per_node": 20},
    "wilcoxon": {},
}


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser()
    result.add_argument("--dataset-config", required=True)
    result.add_argument("--heldout-slice", required=True)
    result.add_argument("--output-root", default="results/validation_v1/heldout_slice")
    result.add_argument("--n-features", default=2000, type=int)
    result.add_argument("--seed", default=0, type=int)
    result.add_argument("--methods", default=",".join(METHODS))
    result.add_argument(
        "--reference-method",
        choices=("morans_i", "nnsvg"),
        default="morans_i",
        help="Independent held-out gene-ranking statistic used as the spatial reference.",
    )
    return result


def _atomic_write_text(path: Path, value: str) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(value, encoding="utf-8")
    temporary.replace(path)


def _write_selection(path: Path, method: str, selector, selection: FeatureSelectionResult) -> None:
    _atomic_write_text(
        path,
        json.dumps(
            {
                "method": method,
                "implementation_version": getattr(selector, "implementation_version", None),
                "feature_names": selection.feature_names,
                "scores": selection.scores,
                "metadata": selection.metadata,
            },
            indent=2,
            default=str,
        ),
    )


def _read_selection(path: Path, gene_names: list[str]) -> FeatureSelectionResult:
    payload = json.loads(path.read_text(encoding="utf-8"))
    index_by_gene = {gene: index for index, gene in enumerate(gene_names)}
    names = [str(gene) for gene in payload["feature_names"] if str(gene) in index_by_gene]
    return FeatureSelectionResult(
        method_name=str(payload["method"]),
        feature_names=names,
        feature_indices=[index_by_gene[gene] for gene in names],
        scores=[float(value) for value in payload["scores"][: len(names)]],
        metadata=dict(payload.get("metadata", {})),
    )


def main() -> None:
    args = parser().parse_args()
    requested = [value.strip() for value in args.methods.split(",") if value.strip()]
    unknown = sorted(set(requested) - set(METHODS))
    if unknown:
        raise ValueError(f"Unknown methods: {unknown}")

    config = load_dataset_config(args.dataset_config)
    dataset = preprocess_dataset(load_dataset(config), config.preprocess)
    fold = HeldoutFold(dataset=dataset.name, heldout_slice=args.heldout_slice)
    training, heldout = split_heldout_slice(dataset, args.heldout_slice)

    fold_dir = Path(args.output_root) / fold.name
    ranking_dir = fold_dir / "rankings"
    metric_dir = fold_dir / "metrics"
    status_dir = fold_dir / "status"
    for path in (ranking_dir, metric_dir, status_dir):
        path.mkdir(parents=True, exist_ok=True)

    reference_name = "moran" if args.reference_method == "morans_i" else "nnsvg"
    reference_path = fold_dir / f"heldout_{reference_name}_reference.json"
    if reference_path.exists():
        reference = _read_selection(reference_path, heldout.gene_names)
    else:
        if args.reference_method == "morans_i":
            reference_selector = build_feature_selector("morans_i", n_neighbors=8)
            reference = heldout_moran_reference(heldout, random_seed=args.seed)
        else:
            reference_selector = build_feature_selector("nnsvg", n_threads=4)
            reference = reference_selector.select(
                heldout,
                n_features=heldout.n_vars,
                random_seed=args.seed,
            )
        _write_selection(reference_path, args.reference_method, reference_selector, reference)

    manifest = {
        "dataset": dataset.name,
        "dataset_config": str(Path(args.dataset_config).resolve()),
        "heldout_slice": args.heldout_slice,
        "training_slices": sorted(set(training.slice_ids.tolist())),
        "n_training_spots": training.n_obs,
        "n_heldout_spots": heldout.n_obs,
        "n_genes": dataset.n_vars,
        "n_features": args.n_features,
        "seed": args.seed,
        "reference_method": args.reference_method,
    }
    _atomic_write_text(fold_dir / "fold.json", json.dumps(manifest, indent=2))

    for method in requested:
        ranking_path = ranking_dir / f"{method}.json"
        metric_path = metric_dir / f"{method}.tsv"
        status_path = status_dir / f"{method}.json"
        if metric_path.exists() and status_path.exists():
            status = json.loads(status_path.read_text(encoding="utf-8"))
            if status.get("status") == "complete":
                print(f"fold={fold.name} method={method} status=cached", flush=True)
                continue

        started = time.perf_counter()
        try:
            if ranking_path.exists():
                selection = _read_selection(ranking_path, training.gene_names)
                cache_state = "ranking_cache"
                runtime = 0.0
            else:
                selector = build_feature_selector(method, **METHODS[method])
                requested_n = training.n_vars if method == "all_features" else args.n_features
                selection = selector.select(training, n_features=requested_n, random_seed=args.seed)
                runtime = time.perf_counter() - started
                cache_state = "computed"
                _write_selection(ranking_path, method, selector, selection)

            spatial_scores = score_spatial_reproducibility(
                heldout,
                selection,
                reference,
                n_features=args.n_features,
            )
            if args.reference_method != "morans_i":
                spatial_scores = {
                    key.replace("heldout_moran", f"heldout_{reference_name}").replace(
                        "heldout_topn_jaccard", f"heldout_{reference_name}_topn_jaccard"
                    ).replace(
                        "heldout_rank_spearman", f"heldout_{reference_name}_rank_spearman"
                    ): value
                    for key, value in spatial_scores.items()
                }
            clustering = score_heldout_clustering(
                heldout,
                selection,
                n_features=args.n_features,
                random_seed=args.seed,
            )
            metrics = heldout_metrics_frame(
                fold,
                method,
                spatial_scores,
                clustering,
                n_features=args.n_features,
                seed=args.seed,
                runtime_seconds=runtime,
            )
            temporary = metric_path.with_suffix(".tsv.tmp")
            metrics.to_csv(temporary, sep="\t", index=False)
            temporary.replace(metric_path)
            status = {
                "method": method,
                "status": "complete",
                "runtime_seconds": runtime,
                "n_selected": len(selection.feature_names),
                "cache_state": cache_state,
                "error": "",
            }
            _atomic_write_text(status_path, json.dumps(status, indent=2))
            print(
                f"fold={fold.name} method={method} status=complete runtime={runtime:.1f}s",
                flush=True,
            )
        except Exception as exc:
            runtime = time.perf_counter() - started
            status = {
                "method": method,
                "status": "failed",
                "runtime_seconds": runtime,
                "n_selected": 0,
                "cache_state": "failed",
                "error": repr(exc),
            }
            _atomic_write_text(status_path, json.dumps(status, indent=2))
            print(f"fold={fold.name} method={method} status=failed error={exc!r}", flush=True)

    metric_frames = [
        pd.read_csv(metric_dir / f"{method}.tsv", sep="\t")
        for method in requested
        if (metric_dir / f"{method}.tsv").exists()
    ]
    status_rows = [
        json.loads((status_dir / f"{method}.json").read_text(encoding="utf-8"))
        for method in requested
        if (status_dir / f"{method}.json").exists()
    ]
    merged = pd.concat(metric_frames, ignore_index=True) if metric_frames else pd.DataFrame()
    merged.to_csv(fold_dir / "metrics.tsv", sep="\t", index=False)
    pd.DataFrame(status_rows).to_csv(fold_dir / "status.tsv", sep="\t", index=False)
    failures = sum(row["status"] != "complete" for row in status_rows)
    print(f"fold={fold.name} methods={len(status_rows)} failures={failures} metrics={len(merged)}")


if __name__ == "__main__":
    main()

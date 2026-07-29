#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import pandas as pd

from spatial_fs_benchmark.feature_selection import build_feature_selector
from spatial_fs_benchmark.feature_selection.base import FeatureSelectionResult
from spatial_fs_benchmark.validation import (
    SemiSyntheticSpec,
    build_semi_synthetic_dataset,
    score_feature_ranking,
)


METHODS = {
    "all_features": {},
    "random": {},
    "scanpy_seurat": {},
    "scanpy_seurat_batch": {},
    "scanpy_seurat_v3": {},
    "scanpy_seurat_v3_batch": {},
    "scanpy_cell_ranger": {},
    "scanpy_cell_ranger_batch": {},
    "scanpy_pearson": {},
    "scanpy_pearson_batch": {},
    "seurat_vst": {},
    "seurat_mvp": {},
    "seurat_disp": {},
    "seurat_sct": {"max_cells": 3000},
    "scsegindex": {},
    "dubstepr": {"max_cells": 3000},
    "nbumi": {"max_cells": 3000},
    "osca": {},
    "scry": {},
    "triku": {"max_cells": 3000},
    "hotspot": {"max_cells": 3000},
    "singleCellHaystack": {"max_cells": 3000},
    "Brennecke": {},
    "scPNMF": {"max_cells": 3000},
    "anticor": {"max_cells": 3000},
    "statistic_mean": {},
    "statistic_variance": {},
    "morans_i": {"n_neighbors": 8},
    "sparkx": {"n_threads": 4},
    "nnsvg": {"n_threads": 4},
    "spatialde": {"max_cells_per_slice": 200},
    "somde": {"spots_per_node": 20},
    "wilcoxon": {},
}


def _atomic_write_text(path: Path, value: str) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(value, encoding="utf-8")
    temporary.replace(path)


def _selection_from_cache(path: Path, gene_names: list[str]) -> FeatureSelectionResult:
    payload = json.loads(path.read_text(encoding="utf-8"))
    index_by_gene = {gene: index for index, gene in enumerate(gene_names)}
    selected_names = [str(gene) for gene in payload["feature_names"] if str(gene) in index_by_gene]
    return FeatureSelectionResult(
        method_name=str(payload["method"]),
        feature_names=selected_names,
        feature_indices=[index_by_gene[gene] for gene in selected_names],
        scores=[float(value) for value in payload["scores"][: len(selected_names)]],
        metadata=dict(payload.get("metadata", {})),
    )


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser()
    result.add_argument("--prevalence", required=True, type=float)
    result.add_argument("--seed", required=True, type=int)
    result.add_argument("--output-root", default="results/validation_v1/semi_synthetic")
    result.add_argument("--n-genes", default=3000, type=int)
    result.add_argument("--grid-size", default=15, type=int)
    result.add_argument("--n-slices", default=3, type=int)
    result.add_argument("--methods", default=",".join(METHODS))
    result.add_argument("--save-dataset", action="store_true")
    return result


def main() -> None:
    args = parser().parse_args()
    requested_methods = [value.strip() for value in args.methods.split(",") if value.strip()]
    unknown = sorted(set(requested_methods) - set(METHODS))
    if unknown:
        raise ValueError(f"Unknown methods: {unknown}")

    spec = SemiSyntheticSpec(
        prevalence=args.prevalence,
        seed=args.seed,
        n_genes=args.n_genes,
        grid_size=args.grid_size,
        n_slices=args.n_slices,
    )
    dataset = build_semi_synthetic_dataset(spec)
    scenario_dir = Path(args.output_root) / dataset.name
    ranking_dir = scenario_dir / "rankings"
    metric_dir = scenario_dir / "metrics"
    status_dir = scenario_dir / "status"
    ranking_dir.mkdir(parents=True, exist_ok=True)
    metric_dir.mkdir(parents=True, exist_ok=True)
    status_dir.mkdir(parents=True, exist_ok=True)
    if args.save_dataset:
        dataset.save_h5ad(scenario_dir / "dataset.h5ad")
    dataset.adata.var.reset_index(names="gene").to_csv(
        scenario_dir / "truth.tsv", sep="\t", index=False
    )
    (scenario_dir / "scenario.json").write_text(
        json.dumps(dataset.adata.uns["semi_synthetic_spec"], indent=2, sort_keys=True),
        encoding="utf-8",
    )

    for method in requested_methods:
        cache_path = ranking_dir / f"{method}.json"
        metric_path = metric_dir / f"{method}.tsv"
        status_path = status_dir / f"{method}.json"
        if metric_path.exists() and status_path.exists():
            status = json.loads(status_path.read_text(encoding="utf-8"))
            if status.get("status") == "complete":
                print(f"scenario={dataset.name} method={method} status=cached", flush=True)
                continue
        started = time.perf_counter()
        try:
            selector = None
            if cache_path.exists():
                selection = _selection_from_cache(cache_path, dataset.gene_names)
                runtime = 0.0
                cache_state = "ranking_cache"
            else:
                selector = build_feature_selector(method, **METHODS[method])
                selection = selector.select(dataset, n_features=dataset.n_vars, random_seed=args.seed)
                runtime = time.perf_counter() - started
                cache_state = "computed"
                _atomic_write_text(
                    cache_path,
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
            metrics = score_feature_ranking(dataset, selection)
            metrics.insert(0, "method", method)
            metrics.insert(0, "seed", args.seed)
            metrics.insert(0, "prevalence", args.prevalence)
            metrics["runtime_seconds"] = runtime
            metrics.to_csv(metric_path.with_suffix(".tsv.tmp"), sep="\t", index=False)
            metric_path.with_suffix(".tsv.tmp").replace(metric_path)
            status = {
                "method": method,
                "status": "complete",
                "runtime_seconds": runtime,
                "n_ranked": len(selection.feature_names),
                "cache_state": cache_state,
                "error": "",
            }
            _atomic_write_text(status_path, json.dumps(status, indent=2))
            print(
                f"scenario={dataset.name} method={method} status=complete "
                f"runtime={runtime:.1f}s source={cache_state}",
                flush=True,
            )
        except Exception as exc:
            runtime = time.perf_counter() - started
            status = {
                "method": method,
                "status": "failed",
                "runtime_seconds": runtime,
                "n_ranked": 0,
                "cache_state": "failed",
                "error": repr(exc),
            }
            _atomic_write_text(status_path, json.dumps(status, indent=2))
            print(
                f"scenario={dataset.name} method={method} status=failed error={exc!r}",
                flush=True,
            )

    metric_frames = [
        pd.read_csv(metric_dir / f"{method}.tsv", sep="\t")
        for method in requested_methods
        if (metric_dir / f"{method}.tsv").exists()
    ]
    status_rows = [
        json.loads((status_dir / f"{method}.json").read_text(encoding="utf-8"))
        for method in requested_methods
        if (status_dir / f"{method}.json").exists()
    ]
    metrics = pd.concat(metric_frames, ignore_index=True) if metric_frames else pd.DataFrame()
    metrics.to_csv(scenario_dir / "metrics.tsv", sep="\t", index=False)
    pd.DataFrame(status_rows).to_csv(scenario_dir / "status.tsv", sep="\t", index=False)
    failures = sum(row["status"] != "complete" for row in status_rows)
    print(
        f"scenario={dataset.name} methods={len(status_rows)} failures={failures} "
        f"metrics={len(metrics)}"
    )


if __name__ == "__main__":
    main()

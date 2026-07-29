from __future__ import annotations

import argparse
import multiprocessing as mp
from pathlib import Path

import numpy as np
import pandas as pd
import scipy
from scipy import io as spio
from scipy import sparse
from scipy.stats import rankdata


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--max-cells-per-slice", type=int, default=0)
    parser.add_argument("--max-genes", type=int, default=2000)
    parser.add_argument("--n-threads", type=int, default=8)
    parser.add_argument("--seed", type=int, default=0)
    return parser.parse_args()


def patch_legacy_dependencies() -> None:
    scipy.inf = np.inf
    np.infty = np.inf
    np.int = int
    np.float = float
    np.bool = bool


def row_variances(matrix) -> np.ndarray:
    means = np.asarray(matrix.mean(axis=1)).ravel()
    means_sq = np.asarray(matrix.multiply(matrix).mean(axis=1)).ravel()
    return np.maximum(means_sq - means**2, 0)


def main() -> None:
    args = parse_args()
    patch_legacy_dependencies()
    import scGCO

    input_dir = Path(args.input_dir)
    counts = spio.mmread(input_dir / "counts.mtx").tocsr()
    genes = pd.read_csv(input_dir / "genes.tsv", sep="\t", header=None)[0].astype(str).to_numpy()
    observations = pd.read_csv(input_dir / "observations.tsv", sep="\t")
    score_rows = []
    for slice_number, slice_id in enumerate(sorted(observations["slice"].astype(str).unique())):
        indices = np.flatnonzero(observations["slice"].astype(str).to_numpy() == slice_id)
        if args.max_cells_per_slice > 0 and len(indices) > args.max_cells_per_slice:
            rng = np.random.default_rng(args.seed + slice_number * 1009)
            indices = np.sort(rng.choice(indices, args.max_cells_per_slice, replace=False))
        if len(indices) < 10:
            continue
        matrix = counts[:, indices]
        variances = row_variances(matrix)
        candidates = np.flatnonzero(np.isfinite(variances) & (variances > 0))
        if len(candidates) > args.max_genes:
            candidates = candidates[np.argsort(variances[candidates])[::-1][: args.max_genes]]
        coords = observations.iloc[indices][["x", "y"]].to_numpy(dtype=float)
        dense = matrix[candidates].T.toarray()
        expression = pd.DataFrame(dense, columns=genes[candidates])
        library_sizes = expression.to_numpy().sum(axis=1)
        positive_sizes = library_sizes[library_sizes > 0]
        median_size = np.median(positive_sizes) if len(positive_sizes) else 1.0
        factors = np.maximum(library_sizes / max(median_size, 1.0), 1e-12)
        normalized = pd.DataFrame(
            np.log1p(expression.to_numpy() / factors[:, None]),
            columns=expression.columns,
            index=expression.index,
        )
        graph = scGCO.create_graph_with_weight(coords, normalized.sum(axis=1).to_numpy())
        mixtures = scGCO.multiGMM(normalized, random_state=args.seed, ncores=args.n_threads)
        chunks = [chunk for chunk in np.array_split(normalized, args.n_threads, axis=1) if chunk.shape[1] > 0]
        worker_args = [
            (coords, chunk, graph, mixtures, 30, 100, 10, "expansion", "gmm")
            for chunk in chunks
        ]
        if len(worker_args) == 1:
            chunk_results = [scGCO.compute_single_fixed_sf(*worker_args[0])]
        else:
            with mp.get_context("fork").Pool(processes=len(worker_args)) as pool:
                chunk_results = pool.starmap(scGCO.compute_single_fixed_sf, worker_args)
        result_genes = []
        pvalues = []
        for chunk_result in chunk_results:
            result_genes.extend(chunk_result[2])
            pvalues.extend(min(values) if len(values) else 1.0 for values in chunk_result[1])
        pvalues = np.asarray(pvalues, dtype=float)
        statistic = -np.log10(np.maximum(pvalues, np.finfo(float).tiny))
        percentiles = rankdata(statistic, method="average") / len(statistic)
        row = np.zeros(len(genes), dtype=float)
        score_map = dict(zip(map(str, result_genes), percentiles, strict=True))
        for candidate in candidates:
            row[candidate] = score_map.get(genes[candidate], 0.0)
        score_rows.append(row)
    if not score_rows:
        raise RuntimeError("No slice produced a valid scGCO ranking")
    scores = np.mean(np.vstack(score_rows), axis=0)
    output = pd.DataFrame({"Feature": genes, "Score": scores}).sort_values("Score", ascending=False)
    output.to_csv(args.output, sep="\t", index=False)


if __name__ == "__main__":
    main()

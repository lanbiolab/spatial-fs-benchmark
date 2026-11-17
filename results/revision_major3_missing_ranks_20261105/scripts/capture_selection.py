from __future__ import annotations

import argparse
import json
import os
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse

from spatial_fs_benchmark.config import load_benchmark_config
from spatial_fs_benchmark.data.io import load_dataset
from spatial_fs_benchmark.data.preprocess import preprocess_dataset
from spatial_fs_benchmark.feature_selection import build_feature_selector
from spatial_fs_benchmark.benchmark.experiment import write_feature_selection

ROOT = Path('/home/xd/y/Benchmark')
OUT = ROOT / 'results/revision_major3_missing_ranks_20261105'


def prepare_r_script():
    text = (ROOT / 'scripts/run_spatial_r_feature_selector.R').read_text()
    old = '    slice_id <- slice_ids[[slice_number]]\n    indices <- which(as.character(observations$slice) == slice_id)'
    new = '''    slice_id <- slice_ids[[slice_number]]
    slice_cache <- file.path(dirname(output_path), sprintf("slice_%04d.rds", slice_number))
    if (file.exists(slice_cache)) {
        slice_scores[[slice_id]] <- readRDS(slice_cache)
        next
    }
    indices <- which(as.character(observations$slice) == slice_id)'''
    assert text.count(old) == 1
    text = text.replace(old, new)
    assert text.count('    scores <- rep(0, length(genes))') == 1
    text = text.replace('    scores <- rep(0, length(genes))', '    scores <- rep(NA_real_, length(genes))')
    old = '    slice_scores[[slice_id]] <- scores\n}'
    new = '    saveRDS(scores, slice_cache)\n    slice_scores[[slice_id]] <- scores\n}'
    assert text.count(old) == 1
    text = text.replace(old, new)
    old = 'score_matrix <- do.call(cbind, slice_scores)\naggregate_score <- rowMeans(score_matrix, na.rm = TRUE)'
    new = '''score_matrix <- do.call(cbind, slice_scores)
colnames(score_matrix) <- names(slice_scores)
slice_table <- data.frame(Feature = genes, score_matrix, check.names = FALSE)
con <- gzfile(paste0(output_path, ".slices.tsv.gz"), "wt")
write.table(slice_table, con, sep = "\\t", quote = FALSE, row.names = FALSE, na = "NA")
close(con)
aggregate_score <- rowMeans(score_matrix, na.rm = TRUE)'''
    assert text.count(old) == 1
    target = OUT / 'scripts/run_spatial_r_capture.R'
    target.parent.mkdir(parents=True, exist_ok=True)
    modified = text.replace(old, new)
    if not target.exists() or target.read_text() != modified:
        temporary = target.with_name(f'{target.name}.{os.getpid()}.tmp')
        temporary.write_text(modified)
        os.replace(temporary, target)
    return target


def python_scores(dataset, selector, seed, job):
    genes = dataset.adata.var_names.to_numpy().astype(str)
    slices = np.unique(dataset.slice_ids)
    rows, used = [], []
    for slice_number, slice_id in enumerate(slices):
        indices = np.flatnonzero(dataset.slice_ids == slice_id)
        if selector.name == 'morans_i':
            indices = selector._select_slice_indices(indices, seed, slice_number, len(slices))
            if len(indices) < 3:
                continue
        else:
            if selector.max_cells_per_slice is not None and len(indices) > selector.max_cells_per_slice:
                rng = np.random.default_rng(seed + slice_number * 1009)
                indices = np.sort(rng.choice(indices, selector.max_cells_per_slice, replace=False))
            if len(indices) < selector._minimum_cells_per_slice():
                continue
        cache = job / f'slice_{slice_number:04d}.npy'
        if cache.exists():
            row = np.load(cache)
        elif selector.name == 'morans_i':
            row = selector._percentile_ranks(selector._moran_scores(dataset.adata.X[indices], dataset.coords[indices]))
            row = np.where(row > 0, row, np.nan)
            np.save(cache, row)
        else:
            matrix = selector._matrix(dataset)[indices].T
            variances = selector._row_variances(matrix)
            candidates = np.flatnonzero(np.isfinite(variances) & (variances > 0))
            if selector.max_genes is not None and len(candidates) > selector.max_genes:
                candidates = candidates[np.argsort(variances[candidates])[::-1][:selector.max_genes]]
            score_map = selector._score_slice(matrix[candidates], dataset.coords[indices], genes[candidates])
            if not score_map:
                continue
            row = np.full(len(genes), np.nan)
            for candidate in candidates:
                value = score_map.get(genes[candidate], 0.0)
                if np.isfinite(value) and value > 0:
                    row[candidate] = value
            np.save(cache, row)
        assert row.shape == (len(genes),)
        rows.append(row)
        used.append(str(slice_id))
        print(f'{selector.name} seed={seed} slice={slice_number+1}/{len(slices)} ranked={np.isfinite(row).sum()}', flush=True)
    if not rows:
        raise RuntimeError('No eligible slices')
    return np.stack(rows), used


def r_scores(dataset, selector, seed, job):
    genes = dataset.adata.var_names.to_numpy().astype(str)
    runner = prepare_r_script()
    matrix_path = job / 'matrix.mtx'
    selector._write_matrix(matrix_path, dataset.adata.layers['counts'])
    pd.Series(genes).to_csv(job / 'genes.tsv', sep='\t', index=False, header=False)
    pd.DataFrame({'slice': dataset.slice_ids, 'x': dataset.coords[:, 0], 'y': dataset.coords[:, 1]}).to_csv(job / 'observations.tsv', sep='\t', index=False)
    output = job / 'ranked_features.tsv'
    cmd = [str(ROOT / '.conda-env/bin/Rscript'), str(runner), str(job), str(output), selector.name,
           str(int(selector.max_cells_per_slice or 0)), str(int(selector.max_genes or 0)), str(selector.n_threads), str(seed)]
    env = os.environ.copy(); env.update(OMP_NUM_THREADS=str(selector.n_threads), OPENBLAS_NUM_THREADS=str(selector.n_threads))
    subprocess.run(cmd, cwd=ROOT, env=env, check=True)
    table = pd.read_csv(str(output) + '.slices.tsv.gz', sep='\t')
    assert np.array_equal(table['Feature'].astype(str).to_numpy(), genes)
    matrix_path.unlink()
    return table.iloc[:, 1:].to_numpy(dtype=float).T, table.columns[1:].astype(str).tolist()


def main():
    p = argparse.ArgumentParser(); p.add_argument('dataset'); p.add_argument('method'); p.add_argument('seed', type=int); args=p.parse_args()
    cfg = load_benchmark_config(ROOT / f'configs/rebuild_v1/{args.dataset}_svg_features.yaml')
    dc = cfg.datasets[0]
    mc = next(x for x in cfg.feature_selection_methods if x.name == args.method)
    selector = build_feature_selector(mc.name, **mc.params)
    assert args.seed == 0 or selector.stochastic_selection
    job = OUT / 'selection' / args.dataset / args.method / f'seed{args.seed}'
    job.mkdir(parents=True, exist_ok=True)
    dataset = preprocess_dataset(load_dataset(dc), dc.preprocess)
    print('Loaded', dataset.name, dataset.n_obs, dataset.n_vars, args.method, args.seed, flush=True)
    matrix, slices = r_scores(dataset, selector, args.seed, job) if args.method in {'nnsvg','sparkx'} else python_scores(dataset, selector, args.seed, job)
    genes = dataset.adata.var_names.to_numpy().astype(str)
    assert matrix.shape == (len(slices), len(genes)) and np.all((matrix[np.isfinite(matrix)] > 0) & (matrix[np.isfinite(matrix)] <= 1))
    np.savez_compressed(job/'slice_percentiles.npz', scores=matrix, genes=genes, slices=np.asarray(slices))
    valid = np.isfinite(matrix); count=valid.sum(axis=0); zero_mean=np.nan_to_num(matrix, nan=0.0).mean(axis=0)
    missing_mean=np.divide(np.nansum(matrix,axis=0), count, out=np.full(len(genes), np.nan), where=count>0)
    for version, scores in [('original_zero', zero_mean), ('missing_excluded', missing_mean)]:
        if args.method == 'morans_i':
            order=np.argsort(np.nan_to_num(scores,nan=-np.inf))[::-1]
            if version == 'missing_excluded': order=order[np.isfinite(scores[order])]
            result=selector._build_result(args.method, genes, np.nan_to_num(scores,nan=-np.inf), 2000) if len(order)>=2000 else selector._build_ranked_result(args.method,genes,genes[order].tolist(),scores[order],2000)
        else:
            eligible=np.flatnonzero(np.isfinite(scores) & (scores>0)); order=eligible[np.argsort(scores[eligible])[::-1]]
            result=selector._build_ranked_result(args.method,genes,genes[order].tolist(),scores[order],2000)
        write_feature_selection(job / f'{version}_n2000.json', result)
    original_path=ROOT / f'results/spatial_svg_rebuild_v1/{args.dataset}/{dataset.name.lower()}/feature_selection/{args.method}/n2000/seed{args.seed}/selected_features.json'
    original=json.loads(original_path.read_text()); rebuilt=json.loads((job/'original_zero_n2000.json').read_text()); alternative=json.loads((job/'missing_excluded_n2000.json').read_text())
    old, ref, new=original['feature_names'],rebuilt['feature_names'],alternative['feature_names']
    summary={'dataset':dataset.name,'dataset_config':args.dataset,'method':args.method,'seed':args.seed,'total_genes':len(genes),'slices_used':len(slices),'union_ranked_genes':int((count>0).sum()),'intersection_ranked_genes':int((count==len(slices)).sum()),'total_gene_slice_cells':int(valid.size),'ranked_gene_slice_cells':int(valid.sum()),'per_slice_ranked_counts':dict(zip(slices,valid.sum(axis=1).astype(int).tolist())),'old_reproduction_exact_order':old==ref,'old_reproduction_top2000_overlap':len(set(old)&set(ref)),'new_vs_original_top2000_overlap':len(set(old)&set(new)),'original_selected':len(old),'new_selected':len(new)}
    (job/'summary.json').write_text(json.dumps(summary,indent=2)); (job/'DONE').write_text('ok\n')
    print(json.dumps({k:v for k,v in summary.items() if k!='per_slice_ranked_counts'}),flush=True)

if __name__=='__main__':
    import fcntl
    import sys
    if len(sys.argv) != 4:
        raise SystemExit('Expected dataset method seed')
    lock_dir=OUT/'selection'/sys.argv[1]/sys.argv[2]/f'seed{sys.argv[3]}'
    lock_dir.mkdir(parents=True,exist_ok=True)
    with (lock_dir/'job.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX)
        if (lock_dir/'DONE').exists():
            print('Job already complete; no recomputation',flush=True)
        else:
            main()

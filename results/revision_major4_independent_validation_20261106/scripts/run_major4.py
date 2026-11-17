"""Supplementary independent-reference validation; original files are read only."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse
from sklearn.decomposition import TruncatedSVD
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score

from spatial_fs_benchmark.config import load_dataset_config
from spatial_fs_benchmark.data.io import load_dataset
from spatial_fs_benchmark.data.preprocess import preprocess_dataset
from spatial_fs_benchmark.metrics.clustering_metrics import cluster_embedding
from spatial_fs_benchmark.validation.heldout_slice import split_heldout_slice, _normalized_log_matrix

ROOT = Path('/home/xd/y/Benchmark')
OUT = ROOT / 'results/revision_major4_independent_validation_20261106'
ORIGINAL = ROOT / 'results/validation_v1/heldout_slice'
METHODS = ['all_features', 'random', 'scanpy_seurat_v3_batch', 'scanpy_cell_ranger_batch',
           'seurat_vst', 'triku', 'singleCellHaystack', 'scsegindex', 'morans_i',
           'sparkx', 'nnsvg', 'spatialde', 'somde', 'wilcoxon']
REFERENCES = {'all_features', 'random', 'wilcoxon'}

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def group(method):
    return 'reference' if method in REFERENCES else 'competitive'

def atomic_json(path, value):
    temp = path.with_suffix('.tmp')
    temp.write_text(json.dumps(value, indent=2, default=str))
    temp.replace(path)

def freeze():
    OUT.mkdir(exist_ok=True)
    spec = {
        'reviewer': 'R1 Major 4', 'purpose': 'Independent annotation ARI and externally defined cortical-layer marker recovery',
        'original_fold_root': str(ORIGINAL), 'methods': METHODS,
        'folds': sorted(p.name for p in ORIGINAL.iterdir() if (p / 'fold.json').exists()),
        'feature_selection': 'Original seed-0 training-only selections, unchanged; no feature reselection',
        'n_features': 2000, 'all_features_reference': 'All input genes, not budget-matched',
        'evaluation_seeds': [0, 1, 2], 'primary_ari': 'Existing six-fold seed-0 ARI is reproduced and tabulated directly',
        'valid_annotation_sensitivity': 'Exclude genuinely missing/unknown reference labels from scoring and from defining K; cluster all spots',
        'clustering': 'Original selected-count library normalization to 10000, log1p, 30-component TruncatedSVD, original KMeans/MiniBatchKMeans implementation',
        'marker_primary': 'All Zeng_Markers.csv rows with Species H or MH, fixed before overlap evaluation',
        'marker_mapping': 'Unique case-insensitive exact symbol match to DLPFC gene_name; no guesses for aliases and no blind use of source Ensmbl column',
        'marker_background': 'Original preprocessed DLPFC gene universe',
        'marker_endpoints': ['all-panel recall', 'per-layer recall', 'uniform-random expected overlap and enrichment', 'training-expression/detection-matched expected overlap and null 95% range'],
        'matched_null': {'expression_quantiles': 10, 'detection_quantiles': 10, 'draws': 10000, 'seed': 20261106,
                         'matching_data': 'Training slices only, mean raw counts and fraction detected; conditional hypergeometric draws in joint bins'},
        'inference': 'All folds and methods retained; no method-by-fold independence assumption or claim of new biological ground truth',
        'source_manifest_sha256': sha(OUT / 'sources/source_manifest.json'),
    }
    path = OUT / 'analysis_protocol.json'
    if path.exists():
        assert json.loads(path.read_text()) == spec, 'Protocol changed'
    else:
        atomic_json(path, spec)
    snapshots = []
    for fold in spec['folds']:
        for method in METHODS:
            for section, suffix in [('rankings', '.json'), ('metrics', '.tsv'), ('status', '.json')]:
                p = ORIGINAL / fold / section / (method + suffix)
                assert p.exists(), p
                if section == 'status':
                    assert json.loads(p.read_text())['status'] == 'complete'
                snapshots.append({'path': str(p), 'sha256': sha(p)})
    frame = pd.DataFrame(snapshots)
    p = OUT / 'original_input_hashes.tsv'
    if p.exists():
        assert frame.equals(pd.read_csv(p, sep='\t'))
    else:
        frame.to_csv(p, sep='\t', index=False)
    print('Frozen six folds, 14 methods including 3 references, 3 evaluation seeds; original inputs hashed.', flush=True)

def marker_analysis(dataset, training, fold, selections, job):
    genes = np.asarray(dataset.gene_names)
    symbols = dataset.adata.var['gene_name'].astype(str).str.upper()
    lookup = {}
    for g, symbol in zip(genes, symbols):
        lookup.setdefault(symbol, []).append(str(g))
    source = pd.read_csv(OUT / 'sources/Zeng_Markers.csv')
    assert len(source) == 81 and source.Paper.eq('Zeng').all()
    mapping = []
    for r in source.to_dict('records'):
        matches = lookup.get(str(r['Gene']).upper(), [])
        include = str(r['Species']) in {'H', 'MH'}
        mapping.append({**r, 'human_supported': include, 'symbol_matches': len(matches),
                        'dataset_gene_id': matches[0] if len(matches) == 1 else '',
                        'mapping_status': 'unique_exact_symbol' if len(matches) == 1 else 'not_uniquely_matched',
                        'source_id_agrees': len(matches) == 1 and str(r['Ensmbl']) == matches[0]})
    mapped = pd.DataFrame(mapping)
    mapped.to_csv(job / 'external_marker_mapping.tsv', sep='\t', index=False)
    panel = mapped[mapped.human_supported & mapped.dataset_gene_id.ne('')].copy()
    assert len(panel) > 50 and panel.dataset_gene_id.is_unique
    marker_sets = {'All human-supported layers': set(panel.dataset_gene_id)}
    for layer in ['1', '2', '3', '4', '5', '6', '6b', 'WM']:
        marker_sets['L' + layer if layer != 'WM' else 'WM'] = set(panel.loc[panel[layer].eq(1), 'dataset_gene_id'])
    marker_details = []
    background = set(genes)
    counts = sparse.csr_matrix(training.adata.layers['counts'])
    mean = np.asarray(counts.mean(axis=0)).ravel()
    detection = np.asarray((counts > 0).mean(axis=0)).ravel()
    mean_bins = pd.qcut(mean, q=10, labels=False, duplicates='drop')
    detection_bins = pd.qcut(detection, q=10, labels=False, duplicates='drop')
    bins = pd.Series(list(zip(mean_bins, detection_bins))).astype(str).to_numpy()
    bin_names = sorted(set(bins))
    bgroups = {b: set(genes[bins == b]) for b in bin_names}
    pd.DataFrame({'gene': genes, 'mean_training_counts': mean, 'training_detection_fraction': detection,
                  'joint_bin': bins}).to_csv(job / 'training_background_bins.tsv', sep='\t', index=False)
    seed = 20261106 + int(hashlib.sha256(fold.encode()).hexdigest()[:8], 16)
    for method, selected in selections.items():
        chosen = set(selected)
        for layer, markers in marker_sets.items():
            hits = chosen & markers
            uniform = len(chosen) * len(markers) / len(background)
            rng = np.random.default_rng(seed)
            null = np.zeros(10000, dtype=int)
            expected = 0.
            for b, genes_b in bgroups.items():
                n = len(chosen & genes_b)
                if n == 0:
                    continue
                m = len(markers & genes_b)
                expected += n * m / len(genes_b)
                null += rng.hypergeometric(m, len(genes_b) - m, n, size=len(null))
            marker_details.append({'fold': fold, 'dataset': dataset.name, 'method': method, 'method_role': group(method),
                                   'marker_set': layer, 'background_genes': len(background), 'selected_genes': len(chosen),
                                   'source_human_supported_genes': int(mapped.human_supported.sum()),
                                   'available_markers': len(markers), 'overlap': len(hits),
                                   'marker_recall': len(hits) / len(markers) if markers else np.nan,
                                   'uniform_expected_overlap': uniform,
                                   'uniform_fold_enrichment': len(hits) / uniform if uniform else np.nan,
                                   'matched_expected_overlap': expected,
                                   'matched_fold_enrichment': len(hits) / expected if expected else np.nan,
                                   'matched_null_lower': np.quantile(null, .025),
                                   'matched_null_upper': np.quantile(null, .975),
                                   'recovered_marker_ids': ';'.join(sorted(hits))})
    pd.DataFrame(marker_details).to_csv(job / 'external_marker_recovery.tsv', sep='\t', index=False)
    print('Marker recovery finished', fold, 'available markers', len(panel), flush=True)

def run_fold(fold, only=None):
    manifest = json.loads((ORIGINAL / fold / 'fold.json').read_text())
    config = load_dataset_config(manifest['dataset_config'])
    dataset = preprocess_dataset(load_dataset(config), config.preprocess)
    training, heldout = split_heldout_slice(dataset, str(manifest['heldout_slice']))
    assert sorted(np.unique(training.slice_ids).tolist()) == sorted(manifest['training_slices'])
    assert str(manifest['heldout_slice']) not in set(training.slice_ids)
    job = OUT / 'folds' / fold
    job.mkdir(parents=True, exist_ok=True)
    label_column = heldout.adata.obs[heldout.label_key]
    truth = np.asarray([str(v) for v in heldout.labels])
    valid = label_column.notna().to_numpy() & ~pd.Series(truth).str.strip().str.lower().isin(['', 'nan', 'none', 'unknown', 'unassigned', 'na', 'n/a']).to_numpy()
    legacy_k = len(np.unique(truth))
    valid_k = len(np.unique(truth[valid]))
    assert valid.sum() > valid_k > 1
    pd.DataFrame({'spot': heldout.adata.obs_names, 'benchmark_label': truth, 'valid_annotation': valid}).to_csv(job / 'spot_annotations.tsv', sep='\t', index=False)
    atomic_json(job / 'fold_audit.json', {**manifest, 'legacy_k': legacy_k, 'valid_annotation_k': valid_k,
                                         'annotated_spots': int(valid.sum()), 'excluded_unannotated_spots': int((~valid).sum()),
                                         'no_heldout_slice_in_training': True})
    print('Loaded', fold, 'spots', heldout.n_obs, 'missing annotation', int((~valid).sum()), 'K', legacy_k, valid_k, flush=True)
    selections = {}
    for method in METHODS:
        if only and method not in only:
            continue
        p = ORIGINAL / fold / 'rankings' / f'{method}.json'
        saved = json.loads(p.read_text())
        selected = [str(g) for g in saved['feature_names'] if str(g) in set(heldout.gene_names)]
        if method != 'all_features':
            selected = selected[:2000]
            assert len(selected) == 2000
        assert len(set(selected)) == len(selected)
        selections[method] = selected
        matrix = None
        for seed in [0, 1, 2]:
            target = job / f'{method}_seed{seed}.json'
            predictions = job / f'{method}_seed{seed}.npz'
            if target.exists() and predictions.exists():
                print('Cached', fold, method, seed, flush=True)
                continue
            if matrix is None:
                matrix = _normalized_log_matrix(heldout, selected)
            components = min(30, matrix.shape[0] - 1, matrix.shape[1] - 1)
            embedding = TruncatedSVD(n_components=components, random_state=seed).fit_transform(matrix)
            original_pred = cluster_embedding(embedding, n_clusters=legacy_k, random_seed=seed)
            annotation_pred = original_pred if valid_k == legacy_k else cluster_embedding(embedding, n_clusters=valid_k, random_seed=seed)
            legacy_ari = adjusted_rand_score(truth, original_pred)
            valid_ari = adjusted_rand_score(truth[valid], annotation_pred[valid])
            old = pd.read_csv(ORIGINAL / fold / 'metrics' / f'{method}.tsv', sep='\t')
            old_ari = float(old.loc[old.metric.eq('heldout_ari'), 'value'].iloc[0])
            row = {'fold': fold, 'dataset': dataset.name, 'heldout_slice': str(manifest['heldout_slice']),
                   'method': method, 'method_role': group(method), 'evaluation_seed': seed, 'selection_seed': 0,
                   'actual_features': len(selected), 'legacy_ari': float(legacy_ari),
                   'legacy_nmi': float(normalized_mutual_info_score(truth, original_pred)),
                   'valid_annotation_ari': float(valid_ari),
                   'valid_annotation_nmi': float(normalized_mutual_info_score(truth[valid], annotation_pred[valid])),
                   'frozen_original_seed0_ari': old_ari,
                   'seed0_absolute_ari_difference': abs(float(legacy_ari) - old_ari) if seed == 0 else None,
                   'selection_file_sha256': sha(p), 'feature_order_sha256': hashlib.sha256('\n'.join(selected).encode()).hexdigest()}
            np.savez_compressed(predictions, spot_ids=heldout.adata.obs_names.to_numpy().astype(str),
                                true_labels=truth, valid_annotation=valid,
                                original_predictions=original_pred, valid_annotation_predictions=annotation_pred)
            atomic_json(target, row)
            print('Finished', fold, method, seed, 'legacy ARI', round(legacy_ari, 4), 'valid ARI', round(valid_ari, 4),
                  'seed0 reproduction', row['seed0_absolute_ari_difference'], flush=True)
        del matrix
    if dataset.name == 'DLPFC' and only is None:
        marker_analysis(dataset, training, fold, selections, job)
    if only is None:
        assert len(list(job.glob('*_seed*.json'))) == 42
        (job / 'DONE').write_text('complete\n')

if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--freeze', action='store_true')
    p.add_argument('--fold')
    p.add_argument('--methods')
    a = p.parse_args()
    if a.freeze:
        freeze()
    if a.fold:
        run_fold(a.fold, a.methods.split(',') if a.methods else None)

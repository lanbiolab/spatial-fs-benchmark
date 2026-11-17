"""Read-only audit of the completed M3 experiments; write audit results only."""
import json
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path('/home/xd/y/Benchmark')
OUT = ROOT / 'results/revision_major3_missing_ranks_20261105'
rows = []
for path in sorted((OUT / 'selection').glob('*/*/seed*/summary.json')):
    job = path.parent
    info = json.loads(path.read_text())
    assert (job / 'DONE').exists()
    with np.load(job / 'slice_percentiles.npz') as z:
        x, genes, slices = z['scores'], z['genes'], z['slices']
    valid = np.isfinite(x)
    counts = valid.sum(axis=0)
    assert x.shape == (info['slices_used'], info['total_genes'])
    assert (counts > 0).sum() == info['union_ranked_genes']
    assert (counts == len(slices)).sum() == info['intersection_ranked_genes']
    assert valid.sum() == info['ranked_gene_slice_cells']
    assert all(info['per_slice_ranked_counts'][str(s)] == n for s, n in zip(slices, valid.sum(axis=1)))
    means = {
        'original_zero': np.nan_to_num(x, nan=0).mean(axis=0),
        'missing_excluded': np.divide(np.nansum(x, axis=0), counts,
                                     out=np.full(len(genes), np.nan), where=counts > 0),
    }
    saved = {}
    for version, expected in means.items():
        sel = json.loads((job / f'{version}_n2000.json').read_text())
        idx = np.array(sel['feature_indices'])
        assert len(idx) == len(set(idx)) == 2000
        assert np.array_equal(genes[idx], sel['feature_names'])
        assert np.isfinite(expected[idx]).all()
        assert np.allclose(expected[idx], sel['scores'], rtol=1e-12, atol=1e-12)
        assert np.all(np.diff(expected[idx]) <= 1e-12)
        excluded = np.ones(len(genes), dtype=bool)
        excluded[idx] = False
        assert np.nanmax(expected[excluded]) <= expected[idx[-1]] + 1e-12
        saved[version] = sel['feature_names']
    old_path = ROOT / 'results/spatial_svg_rebuild_v1' / info['dataset_config'] / info['dataset'].lower() / 'feature_selection' / info['method'] / 'n2000' / f"seed{info['seed']}" / 'selected_features.json'
    old = json.loads(old_path.read_text())['feature_names']
    assert set(old) == set(saved['original_zero'])
    assert (old == saved['original_zero']) == info['old_reproduction_exact_order']
    assert len(set(old) & set(saved['missing_excluded'])) == info['new_vs_original_top2000_overlap']
    rows.append({'dataset': info['dataset'], 'method': info['method'], 'seed': info['seed'],
                 'formula_and_top2000_verified': True, 'original_gene_set_verified': True,
                 'original_order_exact': info['old_reproduction_exact_order']})
assert len(rows) == 49
pd.DataFrame(rows).to_csv(OUT / 'analysis/final_selection_audit.tsv', sep='\t', index=False)

jobs = []
for path in sorted((OUT / 'downstream').glob('*/*/*/seed*/*/manifest.json')):
    job = path.parent
    assert (job / 'DONE').exists() and not (job / 'FAILED').exists()
    manifest = json.loads(path.read_text())
    result = Path(manifest.get('result_dir', manifest.get('original_dir')))
    nrecords = 0
    for task in manifest['tasks']:
        records = json.loads((result / f'{task}_records.json').read_text())
        assert records
        nrecords += len(records)
        for r in records:
            assert str(r['n_features']) == '2000'
            assert r['fs_method'] == job.parents[1].name
            assert r['integration_method'] == job.parents[2].name
            assert int(r['random_seed']) == int(job.parent.name.removeprefix('seed'))
    jobs.append({'job': str(job.relative_to(OUT / 'downstream')),
                 'version': job.name, 'mode': manifest['mode'], 'records_verified': nrecords})
assert sum(r['version'] == 'missing_excluded' for r in jobs) == 225
assert not list((OUT / 'downstream').rglob('FAILED'))
pd.DataFrame(jobs).to_csv(OUT / 'analysis/final_downstream_audit.tsv', sep='\t', index=False)
frozen = pd.read_csv(OUT / 'analysis/original_frozen_score_reproduction_audit.tsv', sep='\t')
assert len(frozen) == 75 and frozen.abs_difference.max() < 1e-7
result = {'selection_jobs_verified': len(rows), 'selection_original_order_exact': sum(r['original_order_exact'] for r in rows),
          'downstream_jobs_verified': len(jobs), 'downstream_versions': dict(Counter(r['version'] for r in jobs)),
          'downstream_modes': dict(Counter(r['mode'] for r in jobs)),
          'frozen_scores_verified': len(frozen), 'max_frozen_score_difference': float(frozen.abs_difference.max()),
          'failed_markers': 0, 'all_checks_passed': True}
(OUT / 'analysis/final_audit_summary.json').write_text(json.dumps(result, indent=2))
print(json.dumps(result, indent=2))

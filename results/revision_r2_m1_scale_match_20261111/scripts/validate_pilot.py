from pathlib import Path
import json,math
import numpy as np,anndata as ad
OUT=Path('/home/xd/y/Benchmark/results/revision_r2_m1_scale_match_20261111')
done=json.loads((OUT/'pilot_queue_finished.json').read_text());assert done['all_completed']
records=[]
for job in done['jobs']:
    for flow in ['scvi','cellcharter']:
        base=Path(job['output_dir'])/job['dataset_name'].lower()/flow/job['method']/'n2000/seed0'
        rows=[]
        for task in ['integration_eval','clustering_eval']:rows+=json.loads((base/f'{task}_records.json').read_text())
        vals={r['metric_name']:r['metric_value'] for r in rows}
        for m in ['bASW','iLISI','silhouette','CHAOS','PAS']:assert m in vals and math.isfinite(vals[m]),(job,flow,m,vals)
        arr=np.load(base/'embedding.npz',allow_pickle=True);assert all(np.isfinite(arr[k]).all() for k in arr.files if np.issubdtype(arr[k].dtype,np.number))
        records.extend(rows)
a=ad.read_h5ad(OUT/'data/fine_footprint.h5ad');b=ad.read_h5ad(OUT/'data/pseudo_visium.h5ad')
assert a.var_names.equals(b.var_names)
assert np.array_equal(np.asarray(a.X.sum(axis=0)),np.asarray(b.X.sum(axis=0)))
assert a.obs.sample_id.nunique()==b.obs.sample_id.nunique()==6
result={'pilot_jobs':len(done['jobs']),'workflow_runs':8,'required_finite_metric_cells':40,'same_gene_universe':True,'exact_per_gene_count_conservation':True,'all_six_slices_retained':True,'all_checks_passed':True}
(OUT/'pilot_validation.json').write_text(json.dumps(result,indent=2))
print(json.dumps(result,indent=2))

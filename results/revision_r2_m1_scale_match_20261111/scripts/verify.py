from pathlib import Path
import json,hashlib
import numpy as np,pandas as pd,anndata as ad
from scipy import sparse

o=Path('/home/xd/y/Benchmark/results/revision_r2_m1_scale_match_20261111')
p=json.loads((o/'protocol.json').read_text())
fine=ad.read_h5ad(o/'data/fine_footprint.h5ad');coarse=ad.read_h5ad(o/'data/pseudo_visium.h5ad')
mapping=pd.read_csv(o/'analysis/bead_to_pseudo_spot.tsv',sep='\t')
assert fine.var_names.equals(coarse.var_names)
assert np.array_equal(fine.obs_names.to_numpy(),mapping.bead_id.to_numpy())
matrix=sparse.csr_matrix((np.ones(len(mapping),dtype=np.int64),(mapping.pseudo_index,np.arange(len(mapping)))),shape=(coarse.n_obs,fine.n_obs))
assert (matrix@fine.layers['counts']!=coarse.layers['counts']).nnz==0
assert fine.layers['counts'].sum()==coarse.layers['counts'].sum()==p['total_counts_per_arm']
assert len(set(fine.obs.sample_id))==len(set(coarse.obs.sample_id))==6
assert (mapping.distance_to_center_um_assumed<=27.5+1e-9).all()
source=ad.read_h5ad(p['source']);source_layer=next((k for k in ['counts','raw_count'] if k in source.layers),None)
raw=sparse.csr_matrix(source.layers[source_layer] if source_layer else source.X)
varidx=source.var_names.get_indexer(fine.var_names);assert (varidx>=0).all()
assert (raw[mapping.source_index.to_numpy()][:,varidx]!=fine.layers['counts']).nnz==0
expected=sparse.csr_matrix(raw[mapping.source_index.to_numpy()]);mask=np.asarray((expected>0).sum(axis=0)).ravel()>=3
assert np.array_equal(np.flatnonzero(mask),varidx)
coverage=pd.read_csv(o/'analysis/run_coverage.tsv',sep='\t')
embeddings=[]
for r in coverage.itertuples():
    if not r.complete:continue
    folder=o/'runs'/r.arm/r.method/f'e95scale_{r.arm}'/r.workflow/r.method/'n2000'/f'seed{r.seed}'
    v=np.load(folder/'embedding.npz',allow_pickle=True)
    shapes={}
    for key in v.files:
        ar=v[key]
        if np.issubdtype(ar.dtype,np.number):assert np.isfinite(ar).all(),(folder,key)
        shapes[key]=list(ar.shape)
    expected_n=fine.n_obs if r.arm=='fine_footprint' else coarse.n_obs
    assert any(len(shape)==2 and shape[0]==expected_n for shape in shapes.values()),(folder,shapes)
    embeddings.append({'arm':r.arm,'method':r.method,'workflow':r.workflow,'seed':int(r.seed),'shapes':shapes})
assert len(embeddings)==int(coverage.complete.sum())
metrics=pd.read_csv(o/'analysis/seed_level_metrics.tsv',sep='\t')
assert len(metrics)==len(embeddings)*5
assert np.isfinite(metrics.metric_value).all()
for (arm,method,seed),part in metrics.groupby(['arm','fs_method','random_seed']):
    sp=o/'runs'/arm/method/f'e95scale_{arm}'/'feature_selection'/method/'n2000'/f'seed{seed}'/'selected_features.json'
    if not sp.exists():sp=sp.parent.parent/'seed0'/'selected_features.json'
    names=json.loads(sp.read_text())['feature_names']
    assert len(names)==len(set(names))
    assert set(names)<=set(fine.var_names)
    assert part.effective_n_features.eq(len(names)).all(),(arm,method,seed)
result={'source_to_fine_counts_identical':True,'fine_to_pooled_counts_identical':True,'gene_universe_fixed_before_pooling':True,'all_six_pucks_retained_in_both_input_arms':True,'all_beads_within_disks':True,'finite_embeddings_verified':len(embeddings),'finite_metric_rows_verified':len(metrics),'embedding_shapes':embeddings}
(o/'analysis/independent_verification.json').write_text(json.dumps(result,indent=2))
print({k:v for k,v in result.items() if k!='embedding_shapes'})

from pathlib import Path
import json,yaml,copy,hashlib,time
import numpy as np,anndata as ad
from scipy import sparse
R=Path('/home/xd/y/Benchmark');O=R/'results/revision_r3_comment5_20261107';P=O/'prepared';D=O/'candidate_data'
status=json.loads((D/'openst_adult_mouse_by_cell.h5ad.tar.gz.status.json').read_text())
original=json.loads((O/'stereo_benchmark_protocol.json').read_text())['source_code_sha256']
assert all(hashlib.sha256((R/p).read_bytes()).hexdigest()==h for p,h in original.items()),'Original benchmark code differs from frozen protocol'
assert status['state']=='downloaded' and status['md5']=='c837f3184d25fb4842a0e27bde1b3941'
a=ad.read_h5ad(D/'openst_by_cell.h5ad');shape=list(a.shape)
assert a.var_names.is_unique and a.obs_names.is_unique
a.X=sparse.csr_matrix(a.X);v=a.X.data
assert np.isfinite(v).all() and np.min(v)>=0 and np.equal(v,np.rint(v)).all()
total=np.asarray(a.X.sum(axis=1)).ravel();mt=a.var_names.str.lower().str.startswith('mt-');assert mt.sum()>0
mt_counts=np.asarray(a.X[:,mt].sum(axis=1)).ravel();pct=100*mt_counts/np.maximum(total,1)
positive=total>0
assert np.allclose(total,a.obs.total_counts)
assert np.allclose(pct[positive],a.obs.pct_counts_mt.to_numpy()[positive],atol=1e-4)
assert np.all(np.isfinite(a.obs.pct_counts_mt.to_numpy()) | ~positive)
print('QC checked: source values match on nonzero-count cells; zero-count rows:',int((~positive).sum()),flush=True)
rules=[('nonbackground',a.obs.cell_ID_mask.to_numpy()!=0),('min_counts_75',total>=75),('max_counts_10000',total<=10000),('pct_mt_below_30',pct<30)]
keep=np.ones(a.n_obs,dtype=bool);steps=[]
for name,mask in rules:
 before=int(keep.sum());keep &= mask;steps.append({'rule':name,'removed_at_step':before-int(keep.sum()),'remaining':int(keep.sum())})
a=a[keep].copy();a.obs['sample_id']='OpenST_adult_hippocampus';a.uns.clear();a.layers.clear()
assert np.isfinite(a.obsm['spatial']).all() and (np.asarray(a.X.sum(axis=1)).ravel()>0).all()
eligible=int((np.asarray((a.X>0).sum(axis=0)).ravel()>=3).sum());assert eligible>=2000
path=P/'openst_adult_brain_counts.h5ad';assert not path.exists(),'Prepared file already exists; inspect before overwriting'
a.write_h5ad(path,compression='gzip')
ds={'name':'OpenSTAdultBrain','path':str(path),'platform':'Open-ST','species':'Mouse','slice_key':'sample_id','label_key':None,'coord_key':'spatial','alignment_pairs':[],'preprocess':{'normalize_total':True,'target_sum':10000.0,'log1p':True,'scale':False,'min_cells_per_gene':3,'count_layer':'X','require_integer_counts':True}}
(P/'openst.yaml').write_text(yaml.safe_dump(ds,sort_keys=False))
manifest={'prepared_at':time.strftime('%Y-%m-%dT%H:%M:%S%z'),'raw_shape':shape,'prepared_shape':list(a.shape),'eligible_genes_min3':eligible,'counts_full_integer_audit_passed':True,'qc_steps':steps,'cell_qc_source':'https://rajewsky-lab.github.io/openst/latest/static/examples/notebooks/hypo_adult_eda/','gene_qc':'Original benchmark min_cells_per_gene=3, consistently with ST/Stereo-seq; not tutorial min_cells=10.','label_key':None,'reference_labels_available':False,'cluster_number':12,'cluster_number_basis':'Fixed before results; same K=12 as original label-free E8.5/E9.5 canonical configurations, not a known biological class count.','metrics':['silhouette','CHAOS','PAS'],'ineligible':['ARI','NMI','Wilcoxon supervised reference','batch mixing','alignment'],'source_archive_md5':status['md5'],'source_archive_sha256':status['sha256']}
(O/'openst_preparation_audit.json').write_text(json.dumps(manifest,indent=2))
jobs=[];C=O/'configs/openst';C.mkdir(parents=True,exist_ok=True)
for j in json.loads((O/'stereo_jobs.json').read_text()):
 if j['method']=='wilcoxon':continue
 cfg=yaml.safe_load(Path(j['config']).read_text());cfg.update(name='r3c5_openst_'+j['method'],datasets=[str(P/'openst.yaml')],output_dir=str(O/'benchmark/openst'/j['method']));cfg['tasks'][0]['params']['n_clusters']=12
 f=C/(j['method']+'.yaml');f.write_text(yaml.safe_dump(cfg,sort_keys=False));jobs.append({**j,'config':str(f)})
assert len(jobs)==33 and sum(j['expected_records'] for j in jobs)==285
(O/'openst_jobs.json').write_text(json.dumps(jobs,indent=2))
protocol=json.loads((O/'stereo_benchmark_protocol.json').read_text());protocol.update(dataset='Open-ST adult mouse hippocampus',preparation_manifest='openst_preparation_audit.json',feature_settings=33,references=['TFs','random','all_features','scsegindex'],expected_combinations=285,task_parameters=cfg['tasks'],scope='Separate label-free platform extension; no ARI/NMI or biological accuracy claim. K fixed at 12 from original label-free protocol; rankings are conditional on this resolution.',excluded_settings={'wilcoxon':'Requires reference labels absent from this object.'},consumers='Original scVI and CellCharter for all 33 eligible settings; original GraphST competitive 29-setting scope.')
(O/'openst_benchmark_protocol.json').write_text(json.dumps(protocol,indent=2))
pilot=yaml.safe_load((C/'morans_i.yaml').read_text());pilot.update(name='r3c5_openst_pilot',output_dir=str(O/'openst_pilot'),integration_methods=[x for x in pilot['integration_methods'] if x['name']=='scvi'],seeds=[0]);(P/'openst_pilot.yaml').write_text(yaml.safe_dump(pilot,sort_keys=False))
print(json.dumps(manifest,indent=2),flush=True);print('Frozen 33 jobs, 285 combinations',flush=True)

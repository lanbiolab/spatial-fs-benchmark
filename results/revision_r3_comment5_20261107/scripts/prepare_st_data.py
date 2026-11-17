from pathlib import Path
import pandas as pd,numpy as np,anndata as ad,json,yaml,hashlib
from scipy import sparse
ROOT=Path('/home/xd/y/Benchmark');O=ROOT/'results/revision_r3_comment5_20261107';D=O/'candidate_data';P=O/'prepared';P.mkdir(exist_ok=True)
SECTIONS=['20A','20B','21A','21B','22A','22B']
metadata=pd.read_csv(D/'st_metadata.tsv.gz',sep='\t',index_col=0)
keep=metadata['section_index'].isin(SECTIONS)&metadata['passed_QC'].eq(True)&~metadata['ABA_parent'].isin(['Undefined areas','ventricular systems'])
m=metadata.loc[keep].copy();assert set(m.section_index)==set(SECTIONS);ids=set(m.index);pieces=[];obs=[];genes=None
for chunk in pd.read_csv(D/'st_counts.tsv.gz',sep='\t',index_col=0,chunksize=500):
 if genes is None:genes=chunk.columns.astype(str).tolist();assert len(set(genes))==len(genes)
 else:assert chunk.columns.astype(str).tolist()==genes
 hit=chunk.loc[chunk.index.isin(ids)]
 if len(hit):
  v=hit.to_numpy(dtype=np.float32);assert np.all(np.isfinite(v)) and np.min(v)>=0 and np.all(v==np.rint(v));pieces.append(sparse.csr_matrix(v));obs.extend(hit.index.tolist())
assert len(obs)==len(ids) and set(obs)==ids
a=ad.AnnData(X=sparse.vstack(pieces,format='csr'),obs=m.loc[obs].copy(),var=pd.DataFrame(index=genes));a.obs['sample_id']=a.obs['section_index'].astype(str);a.obs['annotation']=a.obs['ABA_parent'].astype(str);a.obsm['spatial']=a.obs[['HE_X','HE_Y']].to_numpy(dtype=float)
assert np.isfinite(a.obsm['spatial']).all() and np.asarray(a.X.sum(axis=1)).min()>0
eligible=int((np.asarray((a.X>0).sum(axis=0)).ravel()>=3).sum());assert eligible>=2000
path=P/'st_adult_brain_counts.h5ad';a.write_h5ad(path,compression='gzip')
labels=a.obs['annotation'].value_counts().to_dict();audit={'source_accession':'GSE147747','source_count_file':'GSE147747_expr_raw_counts_table.tsv.gz','sections':SECTIONS,'section_selection':'Six adjacent central coronal sections at levels 20-22, chosen from metadata before outcome inspection; anatomical overlap with cortical/hippocampal structures. No selection by selector results.','inclusion':'Author passed_QC=True; exclude undefined atlas areas and ventricular systems, retaining other source anatomical annotations.','prepared_shape':list(a.shape),'eligible_genes_min3':eligible,'per_section_counts':a.obs.sample_id.value_counts().to_dict(),'atlas_region_counts':labels,'AP_by_section':a.obs.groupby('sample_id')['stereo_AP'].first().to_dict(),'annotation_source':'Source ABA_parent from histological registration to Allen Brain Atlas; do not use source expression cluster_id.','animal_mapping':'Publication reports three mice; mapping of these selected sections to individuals not yet verified. Do not call six sections six animals.','task_scope':'Original clustering endpoints under three consumers; cross-platform comparisons remain separate from batch mixing/alignment.','new_reference_comparison':True}
(O/'st_preparation_audit.json').write_text(json.dumps(audit,indent=2));print(json.dumps(audit),flush=True)
ds={'name':'STAdultBrain','path':str(path),'platform':'ST','species':'Mouse','slice_key':'sample_id','label_key':'annotation','coord_key':'spatial','alignment_pairs':[],'preprocess':{'normalize_total':True,'target_sum':10000.0,'log1p':True,'scale':False,'min_cells_per_gene':3,'count_layer':'X','require_integer_counts':True}}
(P/'st.yaml').write_text(yaml.safe_dump(ds,sort_keys=False))
jobs=json.loads((O/'stereo_jobs.json').read_text());outjobs=[];C=O/'configs/st';C.mkdir(parents=True,exist_ok=True)
for j in jobs:
 cfg=yaml.safe_load(Path(j['config']).read_text());cfg.update(name='r3c5_st_'+j['method'],datasets=[str(P/'st.yaml')],output_dir=str(O/'benchmark/st'/j['method']));cfg['tasks'][0]['params']['n_clusters']=len(labels);f=C/(j['method']+'.yaml');f.write_text(yaml.safe_dump(cfg,sort_keys=False));outjobs.append({**j,'config':str(f)})
(O/'st_jobs.json').write_text(json.dumps(outjobs,indent=2))
protocol=json.loads((O/'stereo_benchmark_protocol.json').read_text());protocol.update(dataset='ST adult brain',preparation_manifest='st_preparation_audit.json',task_parameters=cfg['tasks'],unsupported_endpoints=['Batch mixing/alignment not added to this cross-platform clustering comparison'],scope='Anatomical reference agreement using histology-to-atlas registration; six selected sections are not six independent biological replicates.');(O/'st_benchmark_protocol.json').write_text(json.dumps(protocol,indent=2))

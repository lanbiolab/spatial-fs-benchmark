from pathlib import Path
import anndata as ad,numpy as np,json,yaml,hashlib,copy
from scipy import sparse
ROOT=Path('/home/xd/y/Benchmark');OUT=ROOT/'results/revision_r3_comment5_20261107';D=OUT/'candidate_data';P=OUT/'prepared';P.mkdir(exist_ok=True)
a=ad.read_h5ad(D/'stereo_mouse_brain.h5ad');assert 'count' in a.layers and 'annotation' in a.obs and 'spatial' in a.obsm
counts=a.layers['count'];v=counts.data if sparse.issparse(counts) else np.asarray(counts).ravel()
assert np.all(np.isfinite(v)) and np.min(v)>=0 and np.all(np.abs(v-np.rint(v))<=1e-6),'Source count layer is not nonnegative integer counts'
assert a.var_names.is_unique and a.obs_names.is_unique
labels=a.obs['annotation'].astype(str);keep=labels.ne('Cavity').to_numpy();rawshape=list(a.shape)
a.X=counts.copy();a=a[keep].copy();a.layers.clear();a.obsm['spatial']=np.asarray(a.obsm['spatial'],dtype=float);a.obs['sample_id']='Mouse_brain_Adult';a.uns.clear()
assert np.isfinite(a.obsm['spatial']).all()
detected=np.asarray((a.X>0).sum(axis=0)).ravel();eligible=int((detected>=3).sum());assert eligible>=2000
f=P/'stereo_seq_adult_brain_counts.h5ad'
if not f.exists():a.write_h5ad(f,compression='gzip')
label_counts=a.obs['annotation'].astype(str).value_counts().to_dict()
manifest={'source':'https://ftp.cngb.org/pub/SciRAID/stomics/STDS0000058/stomics/Mouse_brain.h5ad','source_sha256':json.loads((D/'stereo_mouse_brain.h5ad.status.json').read_text())['sha256'],'source_shape':rawshape,'source_count_layer':'count','X_was_log_transformed':True,'counts_full_integer_audit_passed':True,'excluded_observations':int((~keep).sum()),'exclusion_rule':'Remove only source-annotated Cavity (non-tissue); retain all other source annotations; no outcome-dependent filtering.','prepared_shape':list(a.shape),'source_annotation_counts':label_counts,'eligible_genes_min3':eligible,'true_slices':1,'batch_mixing_eligible':False,'alignment_eligible':False,'reference_annotation_provenance_pending':'Source-provided anatomical names; exact annotation construction must be checked before interpreting ARI/NMI as independent validation.','pilot_purpose':'Original implementation compatibility and resource measurement, not final platform comparison','feature_budget':2000,'selection_seed':0,'evaluation_seed':0,'task_n_clusters':len(label_counts),'main_manuscript_unchanged':True}
(OUT/'stereo_preparation_audit.json').write_text(json.dumps(manifest,indent=2));print(json.dumps(manifest),flush=True)
cfg=yaml.safe_load((ROOT/'configs/rebuild_v1/downstream/canonical/mouse_brain_serial_sections.yaml').read_text())
ds={'name':'StereoSeqAdultBrain','path':str(f),'platform':'Stereo-seq','species':'Mouse','slice_key':'sample_id','label_key':'annotation','coord_key':'spatial','alignment_pairs':[],'preprocess':{'normalize_total':True,'target_sum':10000.0,'log1p':True,'scale':False,'min_cells_per_gene':3,'count_layer':'X','require_integer_counts':True}}
(P/'stereo_seq.yaml').write_text(yaml.safe_dump(ds,sort_keys=False))
cfg.update(name='r3c5_stereo_compatibility_pilot',output_dir=str(OUT/'pilot'),datasets=[str(P/'stereo_seq.yaml')],feature_selection_methods=[{'name':'morans_i','params':{'n_neighbors':8},'n_features':[2000]}],integration_methods=[x for x in cfg['integration_methods'] if x['name']=='scvi'],tasks=[copy.deepcopy(x) for x in cfg['tasks'] if x['name']=='clustering_eval'],n_features=[2000],seeds=[0]);cfg['tasks'][0]['params']['n_clusters']=len(label_counts)
(P/'stereo_pilot.yaml').write_text(yaml.safe_dump(cfg,sort_keys=False))

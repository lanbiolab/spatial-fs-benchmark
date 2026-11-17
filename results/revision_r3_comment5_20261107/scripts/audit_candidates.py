from pathlib import Path
import json,tarfile,shutil
import anndata as ad,numpy as np,pandas as pd
from scipy import sparse
ROOT=Path('/home/xd/y/Benchmark');OUT=ROOT/'results/revision_r3_comment5_20261107';D=OUT/'candidate_data';rows=[]
p=D/'openst_adult_mouse_by_cell.h5ad.tar.gz'
if p.exists() and not (D/'openst_by_cell.h5ad').exists():
 with tarfile.open(p) as tar:
  members=[x for x in tar.getmembers() if x.isfile() and x.name.endswith('.h5ad')];assert len(members)==1
  assert members[0].size<5_000_000_000
  with tar.extractfile(members[0]) as src,(D/'openst_by_cell.h5ad').open('wb') as dst:shutil.copyfileobj(src,dst)
for name,p in [('original_visium',ROOT/'data/raw/mouse_brain_serial_sections/mouse_brain_serial_sections.h5ad'),('stereo_seq',D/'stereo_mouse_brain.h5ad'),('open_st',D/'openst_by_cell.h5ad')]:
 if not p.exists():continue
 a=ad.read_h5ad(p,backed='r');row={'dataset':name,'path':str(p),'shape':list(a.shape),'obs_columns':list(a.obs),'var_columns':list(a.var),'layers':list(a.layers),'obsm':list(a.obsm),'uns':list(a.uns),'obs_examples':a.obs.head(3).astype(str).to_dict(),'var_examples':a.var.head(3).astype(str).to_dict(),'gene_examples':a.var_names[:10].tolist()}
 matrix=a.X.to_memory() if hasattr(a.X,'to_memory') else a.X
 x=matrix[:min(1000,a.n_obs),:];v=x.data if sparse.issparse(x) else np.asarray(x).ravel();row['sample_X']={'nonnegative':bool(np.min(v)>=0),'integer_fraction':float(np.mean(np.isclose(v,np.rint(v)))),'max':float(np.max(v))};del matrix,x
 row['candidate_categories']={k:a.obs[k].astype(str).value_counts().head(45).to_dict() for k in a.obs if (a.obs[k].nunique()<100 or any(w in k.lower() for w in ['anno','type','region','cluster','label','sample','batch','slice']))}
 for k in a.obsm:
  v=np.asarray(a.obsm[k]);
  if any(w in k.lower() for w in ['spatial','coord']):row.setdefault('coordinate_bounds',{})[k]={'shape':list(v.shape),'min':np.nanmin(v,axis=0).tolist(),'max':np.nanmax(v,axis=0).tolist()}
 a.file.close();rows.append(row);print(json.dumps(row,default=str),flush=True)
if (D/'st_metadata.tsv.gz').exists():
 m=pd.read_csv(D/'st_metadata.tsv.gz',sep='\t',index_col=0);row={'dataset':'st','shape':list(m.shape),'columns':list(m),'examples':m.head(3).astype(str).to_dict(),'categories':{k:m[k].astype(str).value_counts().head(80).to_dict() for k in m if m[k].nunique()<100 or any(w in k.lower() for w in ['region','label','type','sample','section'])}};rows.append(row);print(json.dumps(row,default=str),flush=True)
(OUT/'candidate_audit.json').write_text(json.dumps(rows,indent=2,default=str))

from pathlib import Path
import anndata as ad,numpy as np,pandas as pd
from scipy.spatial import cKDTree
o=Path('/home/xd/y/Benchmark/results/revision_r2_m1_scale_match_20261111')
a=ad.read_h5ad(o/'data/fine_footprint.h5ad',backed='r');b=ad.read_h5ad(o/'data/pseudo_visium.h5ad',backed='r')
sid=sorted(b.obs.sample_id.astype(str).unique())[0]
x=np.asarray(b.obsm['spatial'])[b.obs.sample_id.astype(str).to_numpy()==sid];center=x[np.argmin(np.linalg.norm(x-np.median(x,axis=0),axis=1))]
for arm,d in [('fine_footprint',a),('pseudo_visium',b)]:
    xy=np.asarray(d.obsm['spatial'])[d.obs.sample_id.astype(str).to_numpy()==sid]-center
    keep=np.max(np.abs(xy),axis=1)<=220
    pd.DataFrame({'x':xy[keep,0],'y':xy[keep,1],'arm':arm,'sample_id':sid}).to_csv(o/'analysis'/f'geometry_{arm}.tsv',sep='\t',index=False)
a.file.close();b.file.close()

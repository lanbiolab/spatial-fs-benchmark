from pathlib import Path
import json
import anndata as ad
import numpy as np
import pandas as pd
import yaml
from scipy.spatial import cKDTree

ROOT=Path('/home/xd/y/Benchmark');OUT=ROOT/'results/revision_r2_major1_20261106/analysis'
rows=[]
for cf in ['dlpfc','mouse_brain_serial_sections','stomics_0212_wilcoxon','stomics_0218_wilcoxon','stomics_0224_wilcoxon','e9p5_embryo']:
    cfg=yaml.safe_load((ROOT/'configs/datasets'/f'{cf}.yaml').read_text());a=ad.read_h5ad(ROOT/cfg['path'],backed='r')
    is_bead=cf=='e9p5_embryo'
    for sid,part in a.obs.groupby(cfg['slice_key'],observed=True):
        inds=np.flatnonzero(a.obs[cfg['slice_key']].astype(str).to_numpy()==str(sid))
        if is_bead:
            xy=np.asarray(a.obsm[cfg['coord_key']])[inds]*0.65
            basis='conditional author convention 0.65 um/native-unit; not an independently verified per-puck calibration'
        else:
            xy=np.column_stack([part.array_col.to_numpy()*50.,part.array_row.to_numpy()*50*np.sqrt(3)])
            basis='standard Visium array-grid pitch 100 um'
        tree=cKDTree(xy);dist=tree.query(xy,k=9)[0]
        for radius in [100,150,200]:
            neighbors=tree.query_ball_point(xy,radius+1e-8,return_length=True)-1
            rows.append({'dataset':cfg['name'],'sample_id':sid,'radius_um_assumed':radius,'n_retained_observations':len(xy),
                         'degree_median':np.median(neighbors),'degree_q25':np.quantile(neighbors,.25),'degree_q75':np.quantile(neighbors,.75),
                         'knn8_support_median_um_assumed':np.median(dist[:,8]),'calibration_basis':basis,
                         'conditional_author_convention':is_bead,'new_model_fit':False})
    a.file.close()
df=pd.DataFrame(rows);df.to_csv(OUT/'conditional_physical_radius_neighborhood_audit.tsv',sep='\t',index=False)
print(df[df.radius_um_assumed.eq(150)].groupby('dataset')[['degree_median','knn8_support_median_um_assumed']].agg(['min','max']).to_string())
audit=json.loads((OUT/'final_audit.json').read_text());audit['conditional_radius_diagnostic']=True
audit['conditional_radius_diagnostic_scope']='5 Visium datasets + E9.5 under author 0.65 convention; no E8.5 physical-unit claim; no refit or support-equivalence claim'
(OUT/'final_audit.json').write_text(json.dumps(audit,indent=2))

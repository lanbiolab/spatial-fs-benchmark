from pathlib import Path
import anndata as ad,pandas as pd,json
o=Path('/home/xd/y/Benchmark/results/revision_r2_m1_scale_match_20261111')
a=ad.read_h5ad(o/'data/pseudo_visium.h5ad',backed='r')
d=pd.DataFrame(a.obsm['spatial'],columns=['x','y']);d['sample_id']=a.obs.sample_id.astype(str).to_numpy()
d.to_csv(o/'analysis/pooled_coordinates.tsv',sep='\t',index=False);a.file.close()
records=[{'arm':'pseudo_visium','method':'sparkx','stage':'feature selection','seeds':[0,1,2],'failure':'SPARK-X mixture cosine kernel 1: exactly singular crossprod of centred transformed coordinates (Lapack dgesv U[2,2]=0)','policy':'Original implementation and coordinates preserved; absent results are unavailable, never zero or a worst rank. Exclude SPARK-X from BOTH arms of paired rank analyses; preserve its fine-arm results separately.','log':'logs/accelerate__pseudo_visium__sparkx.log'}]
(o/'analysis/documented_method_failures.json').write_text(json.dumps(records,indent=2))

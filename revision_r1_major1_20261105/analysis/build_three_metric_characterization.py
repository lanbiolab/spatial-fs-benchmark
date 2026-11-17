"""Figure-1 metric diagnostics on all competitive method-budget settings.

Correlations are computed within each dataset x integrator across 29 x 7
setting rows, then averaged; these are descriptive, not independent replicates.
"""
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

root=Path('/home/xd/y/Benchmark/revision_r1_major1_20261105');out=root/'figures_generated/primary_metric_corrected';out.mkdir(parents=True,exist_ok=True)
scores=pd.read_csv(root/'frozen_three/scaled_setting_metrics.tsv',sep='\t',dtype={'n_features':str});rep=pd.read_csv(root/'frozen_three/representative_settings.tsv',sep='\t',dtype={'n_features':str});groups=rep[rep.MethodGroup.isin(['expression_driven','spatially_informed'])][['dataset','fs_method']].drop_duplicates();assert len(groups)==87
v=scores.merge(groups,on=['dataset','fs_method'],validate='many_to_one');v=v[v.task.isin(['integration_eval','clustering_eval'])].copy();v['n_features']=pd.to_numeric(v.n_features,errors='raise').astype(int)
levels=['bASW','dASW','iLISI','dLISI','ILL','GC','ari','nmi','CHAOS','PAS','silhouette']
assert set(v.metric_name)==set(levels) and not v.duplicated(['dataset','integration_method','fs_method','n_features','metric_name']).any()
assert v.groupby(['dataset','integration_method','metric_name']).size().eq(203).all()
v.to_csv(out/'metric_characterization_values.tsv',sep='\t',index=False,na_rep='NA')
features=[];pairs=[]
for (dataset,integ),frame in v.groupby(['dataset','integration_method']):
 wide=frame.pivot(index=['fs_method','n_features'],columns='metric_name',values='ScaledValue').reindex(columns=levels);assert wide.shape==(203,11) and wide.notna().all().all()
 n=wide.index.get_level_values('n_features').to_numpy()
 for metric in levels:
  features.append({'dataset':dataset,'integrator':integ,'metric':metric,'rho':spearmanr(n,wide[metric]).statistic})
 corr=wide.corr(method='spearman')
 for a in levels:
  for b in levels:pairs.append({'dataset':dataset,'integrator':integ,'metric_a':a,'metric_b':b,'rho':corr.loc[a,b]})
feature=pd.DataFrame(features);pair=pd.DataFrame(pairs)
feature.to_csv(out/'feature_budget_correlations.tsv',sep='\t',index=False)
pair.to_csv(out/'metric_pair_correlations.tsv',sep='\t',index=False)
summary=pair.groupby(['metric_a','metric_b'],as_index=False).agg(mean_rho=('rho','mean'),sd_rho=('rho','std'),n=('rho','count'));assert summary.n.eq(6).all();summary.to_csv(out/'metric_pair_correlations_summary.tsv',sep='\t',index=False)
print('rows',len(v),'feature correlations',len(feature),'metric correlations',len(pair),'all 6 strata x 203 competitive settings')

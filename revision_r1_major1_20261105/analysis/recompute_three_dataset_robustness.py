"""R1 M1: recalculate primary sensitivity with three named-label datasets.

The original Dataset x Metric ranges are retained; three datasets are selected
before ranking and bootstrap. This imports the original scoring implementation
and writes only into the isolated revision directory.
"""
import importlib.util
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import rankdata,spearmanr

ROOT=Path('/home/xd/y/Benchmark');REV=ROOT/'revision_r1_major1_20261105';OUT=REV/'analysis';OUT.mkdir(exist_ok=True)
spec=importlib.util.spec_from_file_location('original_robustness',ROOT/'scripts/build_submission_robustness_tables.py');m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
datasets=['DLPFC','STOmics0212','STOmics0218'];seed=20261105;rng=np.random.default_rng(seed)
allmetrics,rep=m.representative_metrics(ROOT/'results/spatial_svg_rebuild_v1/frozen_scores')
allmetrics=allmetrics[allmetrics.dataset.isin(datasets)].copy();rep=rep[rep.dataset.isin(datasets)].copy()
alt=m.add_alternative_scaling(allmetrics);keys=['dataset','fs_method','n_features']
alt=alt.merge(rep[keys+['MethodGroup']].drop_duplicates(),on=keys,how='inner',validate='many_to_one')
comp=alt[alt.MethodGroup.isin(m.COMPETITIVE_GROUPS)].copy()
base=m.task_scores(comp,'ScaledValue')
scenarios=[m.global_ranks(base,'observed_range'),m.global_ranks(m.task_scores(comp,'RankScaledValue'),'rank_percentile'),m.global_ranks(m.task_scores(comp,'ZScaledValue'),'z_score'),m.global_ranks(m.task_scores(comp,'ScaledValue',(0.5,0.5,0.6)),'integration_weight_0.6'),m.global_ranks(m.task_scores(comp,'ScaledValue',(0.5,0.5,0.4)),'integration_weight_0.4')]
ranks=pd.concat(scenarios,ignore_index=True);ranks.to_csv(OUT/'three_dataset_sensitivity_ranks.tsv',sep='\t',index=False)
concordance=m.concordance_table(ranks);concordance.to_csv(OUT/'three_dataset_sensitivity_concordance.tsv',sep='\t',index=False,na_rep='NA')
bootstrap=m.bootstrap_ranks(base,5000,rng);bootstrap.to_csv(OUT/'three_dataset_bootstrap_ranks.tsv',sep='\t',index=False)

# Same three independent weights are used across datasets in each draw.
w=rng.uniform(.4,.6,size=(2000,3));results=[]
for integrator,p in base.groupby('integration_method'):
 methods=sorted(p.fs_method.unique());assert len(methods)==29
 def arr(col):
  matrix=p.pivot(index='dataset',columns='fs_method',values=col).loc[datasets,methods].to_numpy();assert np.isfinite(matrix).all();return matrix
 bm,bc,li,la=(arr(k) for k in ['BatchMixing','BiologicalConservation','LabelIndependentClustering','LabelAgreement'])
 scores=w[:,2,None,None]*(w[:,0,None,None]*bm+(1-w[:,0,None,None])*bc)+(1-w[:,2,None,None])*(w[:,1,None,None]*li+(1-w[:,1,None,None])*la)
 dataset_ranks=rankdata(-scores,method='average',axis=2);global_ranks=rankdata(dataset_ranks.mean(axis=1),method='average',axis=1)
 for j,method in enumerate(methods):results.append(dict(integration_method=integrator,fs_method=method,draws=len(w),seed=seed,median_rank=np.median(global_ranks[:,j]),lower95=np.quantile(global_ranks[:,j],.025),upper95=np.quantile(global_ranks[:,j],.975),top1_frequency=np.mean(global_ranks[:,j]==1),top5_frequency=np.mean(global_ranks[:,j]<=5)))
weight=pd.DataFrame(results);weight.to_csv(OUT/'three_dataset_weight_perturbation.tsv',sep='\t',index=False)

# Budget scores by dataset then equal-weight averages, preserving requested N.
frame=pd.read_csv(REV/'frozen_three/setting_task_scores.tsv',sep='\t',dtype={'n_features':str});frame=frame[frame.MethodGroup.isin(m.COMPETITIVE_GROUPS)].copy();frame['N']=pd.to_numeric(frame.n_features,errors='coerce');assert set(frame.N.dropna().unique())=={100,200,500,1000,2000,5000,10000}
dataset=frame.groupby(['dataset','integration_method','N'],as_index=False).agg(competitive_methods=('fs_method','nunique'),core_mean=('CoreOverallMean','mean'),integration_mean=('IntegrationMean','mean'),clustering_mean=('ClusteringMean','mean'));assert (dataset.competitive_methods==29).all();dataset.to_csv(OUT/'three_dataset_budget_per_dataset.tsv',sep='\t',index=False)
summary=dataset.groupby(['integration_method','N'],as_index=False).agg(ndatasets=('dataset','nunique'),core_mean=('core_mean','mean'),core_sd_datasets=('core_mean','std'),integration_mean=('integration_mean','mean'),clustering_mean=('clustering_mean','mean'));summary.to_csv(OUT/'three_dataset_budget_summary.tsv',sep='\t',index=False)

# Canonical method and dataset pairs, matched across all three integrators.
r=pd.read_csv(REV/'frozen_three/combined_integrator_scores.tsv',sep='\t');r=r[r.MethodGroup.isin(m.COMPETITIVE_GROUPS)]
pair=r.pivot(index=['dataset','fs_method'],columns='integration_method',values=['CoreOverallMean','IntegrationMean','ClusteringMean','AlignmentMean'])
rows=[]
for integ in ['cellcharter','graphst']:
 for metric in ['CoreOverallMean','IntegrationMean','ClusteringMean','AlignmentMean']:
  delta=(pair[(metric,integ)]-pair[(metric,'scvi')]).dropna();rows.append(dict(integration=integ,metric=metric,n_combinations=len(delta),n_datasets=len(delta.index.get_level_values('dataset').unique()),mean_difference=delta.mean(),median_difference=delta.median()))
pd.DataFrame(rows).to_csv(OUT/'three_dataset_integrator_differences.tsv',sep='\t',index=False)
print('Concordance:');print(concordance[concordance.NScenarios==5].to_string(index=False))
print('Top bootstrap:');print(bootstrap[bootstrap.ObservedGlobalRank<=2].to_string(index=False))
print('Weight perturbation leaders:');print(weight[weight.top1_frequency>.3].to_string(index=False))
print('Budget:');print(summary.to_string(index=False))
print('Integrator differences:');print(pd.DataFrame(rows).to_string(index=False))

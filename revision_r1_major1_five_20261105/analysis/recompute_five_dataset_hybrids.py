"""Rerank hybrid controls on the three-dataset primary scope."""
from pathlib import Path
import pandas as pd
root=Path('/home/xd/y/Benchmark');out=root/'revision_r1_major1_five_20261105/analysis/hybrids';out.mkdir(parents=True,exist_ok=True)
src=root/'results/submission_validations_v1';keep=['DLPFC','MouseBrainSerialSections','STOmics0212','STOmics0218','STOmics0224']
for file in ['hybrid_effective_features.tsv','hybrid_differences_from_best_parent.tsv','hybrid_comparison_task_scores.tsv']:
 d=pd.read_csv(src/file,sep='\t');d=d[d.dataset.isin(keep)].copy();d.to_csv(out/file,sep='\t',index=False,na_rep='NA')
hybrid=pd.read_csv(out/'hybrid_comparison_task_scores.tsv',sep='\t');hybrid=hybrid[hybrid.MethodGroup=='hybrid']
base=pd.read_csv(root/'revision_r1_major1_five_20261105/frozen_five/representative_task_scores.tsv',sep='\t');base=base[base.MethodGroup.isin(['expression_driven','spatially_informed'])]
scores=pd.concat([base,hybrid],ignore_index=True);assert (scores.groupby(['dataset','integration_method']).fs_method.nunique()==31).all()
scores['DatasetRank']=scores.groupby(['dataset','integration_method']).CoreOverallMean.rank(method='average',ascending=False)
rank=scores.groupby(['fs_method','integration_method','MethodGroup'],as_index=False).agg(MeanCoreOverall=('CoreOverallMean','mean'),SDCoreOverall=('CoreOverallMean','std'),MeanDatasetRank=('DatasetRank','mean'),NDatasets=('dataset','nunique'))
rank['GlobalRankWithHybrids']=rank.groupby('integration_method').MeanDatasetRank.rank(method='average');eff=pd.read_csv(out/'hybrid_effective_features.tsv',sep='\t').groupby('fs_method',as_index=False).agg(EffectiveFeaturesMin=('effective_n_features','min'),EffectiveFeaturesMedian=('effective_n_features','median'),EffectiveFeaturesMax=('effective_n_features','max'));rank=rank.merge(eff,on='fs_method',how='left');rank.to_csv(out/'hybrid_global_rank_summary.tsv',sep='\t',index=False,na_rep='NA')
d=pd.read_csv(out/'hybrid_differences_from_best_parent.tsv',sep='\t');print('hybrid delta:');print(d.groupby(['fs_method','integration_method']).DifferenceFromBestParent.agg(['mean','min','max']).to_string());print('hybrid rank:');print(rank[rank.MethodGroup=='hybrid'].to_string(index=False))

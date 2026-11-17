from pathlib import Path
import importlib.util,json
import numpy as np
import pandas as pd

root=Path('/home/xd/y/Benchmark');out=root/'results/revision_major3_missing_ranks_20261105';base=root/'revision_r1_major1_five_20261105/frozen_five'
spec=importlib.util.spec_from_file_location('s',root/'scripts/build_frozen_scores.py');s=importlib.util.module_from_spec(spec);spec.loader.exec_module(s)
rows=[]
for scenario in ['original_frozen','original_zero','missing_excluded']:
 for integ in ['scvi','cellcharter','graphst']:
  for seed in range(3):
   if scenario=='original_frozen':
    d=root/f'results/spatial_svg_rebuild_v1/dlpfc/dlpfc/{integ}/nnsvg/n2000/seed{seed}'
   else:
    job=out/f'downstream/dlpfc/{integ}/nnsvg/seed{seed}/{scenario}';assert (job/'DONE').exists();manifest=json.loads((job/'manifest.json').read_text());d=Path(manifest['original_dir'] if manifest['mode']=='reuse_unchanged_original' else manifest['result_dir'])
   for task in ['integration_eval','clustering_eval','alignment_eval']:
    r=json.loads((d/f'{task}_records.json').read_text())
    for x in r:x['scenario']=scenario
    rows.extend(r)
allrows=pd.DataFrame(rows);allrows.n_features=allrows.n_features.astype(str);allrows['OrientedValue']=allrows.apply(lambda x:float(x.metric_value)*(1 if s.METRIC_METADATA[x.metric_name][2] else -1) if pd.notna(x.metric_value) else np.nan,axis=1)
ranges=pd.read_csv(base/'frozen_metric_ranges.tsv',sep='\t');res=[]
for name,frame in allrows.groupby('scenario'):
 scaled=s.scale_metrics(frame,ranges,'OrientedValue');scores=s.aggregate_setting_scores(s.build_seed_task_scores(scaled));scores['scenario']=name;res.append(scores)
scores=pd.concat(res,ignore_index=True);output=out/'analysis_interim_four_methods';output.mkdir(exist_ok=True);scores.to_csv(output/'dlpfc_nnsvg_only_task_scores.tsv',sep='\t',index=False)
old=pd.read_csv(base/'combined_integrator_scores.tsv',sep='\t');old=old[(old.dataset=='DLPFC')&(old.fs_method=='nnsvg')]
check=scores[scores.scenario=='original_frozen'].merge(old,on=['dataset','fs_method','integration_method'],suffixes=('_computed','_frozen'))
assert max(abs(check.CoreOverallMean_computed-check.CoreOverallMean_frozen))<1e-7
print(scores[['integration_method','scenario','CoreOverallMean','IntegrationMean','ClusteringMean']].to_string(index=False))

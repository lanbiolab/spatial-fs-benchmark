import importlib.util,json
from pathlib import Path
import pandas as pd
import numpy as np
r=Path('/home/xd/y/Benchmark');f=r/'revision_r1_major1_five_20261105/frozen_five'
spec=importlib.util.spec_from_file_location('frozen',r/'scripts/build_frozen_scores.py');m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
records=[]
for integ in ['scvi','cellcharter','graphst']:
 for seed in range(3):
  p=r/f'results/spatial_svg_rebuild_v1/mouse_brain_serial_sections/mousebrainserialsections/{integ}/morans_i/n2000/seed{seed}'
  for task in ['integration_eval','clustering_eval','alignment_eval']:records.extend(json.loads((p/f'{task}_records.json').read_text()))
d=pd.DataFrame(records);d['n_features']=d.n_features.astype(str);d['OrientedValue']=d.apply(lambda x:float(x.metric_value)*(1 if m.METRIC_METADATA[x.metric_name][2] else -1) if pd.notna(x.metric_value) else np.nan,axis=1)
ranges=pd.read_csv(f/'frozen_metric_ranges.tsv',sep='\t');scaled=m.scale_metrics(d,ranges,'OrientedValue');score=m.aggregate_setting_scores(m.build_seed_task_scores(scaled))
orig=pd.read_csv(f/'combined_integrator_scores.tsv',sep='\t')
fields=['BatchMixingMean','BiologicalConservationMean','IntegrationMean','ClusteringMean','AlignmentMean','CoreOverallMean']
for integ in ['scvi','cellcharter','graphst']:
 ref=orig[(orig.dataset=='MouseBrainSerialSections')&(orig.fs_method=='morans_i')&(orig.integration_method==integ)].iloc[0];row=score[score.integration_method==integ].iloc[0]
 diff={x:abs(float(row[x])-float(ref[x])) for x in fields};print(integ,'frozen-score differences:',diff);assert max(diff.values())<1e-8

import importlib.util
import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT=Path('/home/xd/y/Benchmark');OUT=ROOT/'results/revision_major3_missing_ranks_20261105';FIVE=ROOT/'revision_r1_major1_five_20261105/frozen_five'
DATASETS=['mouse_brain_serial_sections','dlpfc','stomics_0212','stomics_0218','stomics_0224']
ALL_DATASETS=DATASETS
METHODS=['morans_i','somde','spatialde','sparkx']
INTEGRATORS=['scvi','cellcharter','graphst']
ANALYSIS=OUT/'analysis_interim_four_methods';ANALYSIS.mkdir(parents=True,exist_ok=True)
spec=importlib.util.spec_from_file_location('scoring',ROOT/'scripts/build_frozen_scores.py')
scoring=importlib.util.module_from_spec(spec);spec.loader.exec_module(scoring)

summaries=[];slice_rows=[]
for dataset in ALL_DATASETS:
    for method in METHODS:
        for selection_seed in (range(3) if method=='spatialde' else [0]):
            job=OUT/'selection'/dataset/method/f'seed{selection_seed}'
            assert (job/'DONE').exists(),f'Missing selection {job}'
            data=json.loads((job/'summary.json').read_text())
            counts=data.pop('per_slice_ranked_counts')
            summaries.append(data)
            for sid,n in counts.items():
                slice_rows.append({'dataset':data['dataset'],'method':method,'seed':selection_seed,'slice':sid,'total_genes':data['total_genes'],'ranked_genes':n,'unranked_genes':data['total_genes']-n})
pd.DataFrame(summaries).to_csv(ANALYSIS/'method_gene_universe_by_dataset.tsv',sep='\t',index=False)
pd.DataFrame(slice_rows).to_csv(ANALYSIS/'method_gene_universe_by_slice.tsv',sep='\t',index=False)

def records_for(dataset, method, seed, integrator, version):
    ss=seed if method=='spatialde' else 0
    selection=OUT/'selection'/dataset/method/f'seed{ss}'
    summary=json.loads((selection/'summary.json').read_text())
    job=OUT/'downstream'/dataset/integrator/method/f'seed{seed}'/version
    if version=='original_frozen' or (version=='original_zero' and summary['old_reproduction_exact_order']):
        job=OUT/'downstream'/dataset/integrator/method/f'seed{seed}'/'missing_excluded'
        assert (job/'DONE').exists()
        manifest=json.loads((job/'manifest.json').read_text())
        # For an exact original selection the frozen original is the paired baseline,
        # independent of whether the alternative was reused or newly evaluated.
        orig_dataset=summary['dataset'].lower()
        result=ROOT/f'results/spatial_svg_rebuild_v1/{dataset}/{orig_dataset}/{integrator}/{method}/n2000/seed{seed}'
        tasks=manifest['tasks']
    else:
        assert (job/'DONE').exists(),f'Missing downstream {job}'
        manifest=json.loads((job/'manifest.json').read_text())
        result=Path(manifest['original_dir'] if manifest['mode']=='reuse_unchanged_original' else manifest['result_dir'])
        tasks=manifest['tasks']
    rows=[]
    for task in tasks:
        path=result/f'{task}_records.json';assert path.exists(),path
        rows.extend(json.loads(path.read_text()))
    for row in rows:
        assert row['fs_method']==method and row['integration_method']==integrator and int(row['random_seed'])==seed
    return rows

ranges=pd.read_csv(FIVE/'frozen_metric_ranges.tsv',sep='\t')
base=pd.read_csv(FIVE/'combined_integrator_scores.tsv',sep='\t',dtype={'n_features':str})
base=base[base.MethodGroup.isin(['expression_driven','spatially_informed']) & base.integration_method.isin(INTEGRATORS)].copy()
assert len(base)==29*5*3 and not base.duplicated(['dataset','fs_method','integration_method']).any()
outputs=[]
for version in ['original_frozen','original_zero','missing_excluded']:
    data=[]
    for dataset in DATASETS:
        for method in METHODS:
            for seed in range(3):
                for integ in INTEGRATORS:data.extend(records_for(dataset,method,seed,integ,version))
    frame=pd.DataFrame(data);frame['n_features']=frame.n_features.astype(str)
    assert not frame.duplicated(scoring.METRIC_KEY).any()
    frame['OrientedValue']=frame.apply(lambda r:float(r.metric_value)*(1 if scoring.METRIC_METADATA[r.metric_name][2] else -1) if pd.notna(r.metric_value) else np.nan,axis=1)
    scaled=scoring.scale_metrics(frame,ranges,'OrientedValue')
    task=scoring.aggregate_setting_scores(scoring.build_seed_task_scores(scaled))
    assert len(task)==5*4*3
    keys=['dataset','fs_method','integration_method']
    if version in ['original_frozen','original_zero']:
        exact_by_method=pd.DataFrame(summaries).groupby(['dataset','method']).old_reproduction_exact_order.all()
        check=task.merge(base[keys+['CoreOverallMean']],on=keys,suffixes=('_recomputed','_frozen'),validate='one_to_one')
        check['exact_selection']=[version=='original_frozen' or bool(exact_by_method.loc[(r.dataset,r.fs_method)]) for r in check.itertuples()]
        check['abs_difference']=(check.CoreOverallMean_recomputed-check.CoreOverallMean_frozen).abs()
        exact=check[check.exact_selection]
        assert exact.abs_difference.notna().all() and exact.abs_difference.max()<1e-7, exact.sort_values('abs_difference').tail()
        check.to_csv(ANALYSIS/f'{version}_score_reproduction_audit.tsv',sep='\t',index=False)
    cols=[c for c in task.columns if c not in keys and c!='n_features']
    combined=base.set_index(keys).copy();recomputed=task.set_index(keys)
    assert recomputed.index.is_unique and recomputed.index.isin(combined.index).all()
    for col in cols:
        if col in combined.columns:combined.loc[recomputed.index,col]=recomputed[col]
    combined=combined.reset_index()
    combined['DatasetRank']=combined.groupby(['dataset','integration_method']).CoreOverallMean.rank(ascending=False,method='average')
    rank=combined.groupby(['integration_method','fs_method'],as_index=False).agg(MeanDatasetRank=('DatasetRank','mean'),MeanCoreOverall=('CoreOverallMean','mean'),NDatasets=('dataset','nunique'))
    assert rank.NDatasets.eq(5).all()
    rank['GlobalRank']=rank.groupby('integration_method').MeanDatasetRank.rank(method='average')
    rank['version']=version;combined['version']=version
    frame.to_csv(ANALYSIS/f'{version}_seed_metrics.tsv',sep='\t',index=False,na_rep='NA')
    combined.to_csv(ANALYSIS/f'{version}_dataset_ranks.tsv',sep='\t',index=False,na_rep='NA')
    rank.to_csv(ANALYSIS/f'{version}_global_ranks.tsv',sep='\t',index=False,na_rep='NA')
    outputs.append(rank)
before,paired,after=outputs
changes=before.merge(after,on=['integration_method','fs_method'],suffixes=('_zero','_missing'),validate='one_to_one')
changes['GlobalRankChange']=changes.GlobalRank_missing-changes.GlobalRank_zero
changes['MeanDatasetRankChange']=changes.MeanDatasetRank_missing-changes.MeanDatasetRank_zero
changes[changes.fs_method.isin(METHODS)].sort_values(['integration_method','GlobalRank_missing']).to_csv(ANALYSIS/'four_spatial_methods_rank_changes.tsv',sep='\t',index=False)
paired_changes=paired.merge(after,on=['integration_method','fs_method'],suffixes=('_zero_refit','_missing'),validate='one_to_one')
paired_changes['GlobalRankChange']=paired_changes.GlobalRank_missing-paired_changes.GlobalRank_zero_refit
paired_changes[paired_changes.fs_method.isin(METHODS)].sort_values(['integration_method','GlobalRank_missing']).to_csv(ANALYSIS/'four_spatial_methods_paired_refit_rank_changes.tsv',sep='\t',index=False)
print('INTERIM ONLY: four spatial methods recalculated; nnSVG remains at original scores.',flush=True)

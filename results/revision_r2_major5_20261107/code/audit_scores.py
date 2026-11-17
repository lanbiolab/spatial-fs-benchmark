from pathlib import Path
import importlib.util,json,hashlib,itertools,shutil
import numpy as np
import pandas as pd
from scipy.stats import rankdata
ROOT=Path('/home/xd/y/Benchmark');OUT=ROOT/'results/revision_r2_major5_20261107';OUT.mkdir(parents=True,exist_ok=True)
FROZEN=ROOT/'results/spatial_svg_rebuild_v1/frozen_scores'
spec=importlib.util.spec_from_file_location('robust',ROOT/'scripts/build_submission_robustness_tables.py');m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
allmetrics,selected=m.representative_metrics(FROZEN)
selected=selected[selected.MethodGroup.isin(m.COMPETITIVE_GROUPS)].copy()
ds=['DLPFC','MouseBrainSerialSections','STOmics0212','STOmics0218','STOmics0224','E8p5Embryo','E9p5Embryo']
primary=ds[:5];align=ds[:2]
task=m.task_scores(selected,'ScaledValue')
assert len(task)==406 and set(task.integration_method)=={'scvi','cellcharter'}
original=pd.read_csv(FROZEN/'representative_task_scores.tsv',sep='\t',dtype={'n_features':str})
check=task.merge(original,on=['dataset','fs_method','n_features','integration_method'],validate='one_to_one')
cols=list(m.COMPONENT_METRICS)+['Integration','Clustering','CoreOverall']
for c in cols:np.testing.assert_allclose(check[c],check[c+'Mean'],rtol=1e-10,atol=1e-12,equal_nan=True)
availability=[];composition=[];contributions=[]
for dataset in ds:
    annotated=dataset in primary
    for comp,metrics in m.COMPONENT_METRICS.items():
        expected=comp in ['BatchMixing','LabelIndependentClustering'] or (comp in ['BiologicalConservation','LabelAgreement'] and annotated) or(comp=='Alignment' and dataset in align)
        for metric in metrics:
            part=selected[selected.dataset.eq(dataset)&selected.metric_name.eq(metric)]
            finite=int(np.isfinite(part.ScaledValue).sum());assert finite==(58 if expected else 0),(dataset,metric,finite)
            allp=allmetrics[allmetrics.dataset.eq(dataset)&allmetrics.metric_name.eq(metric)]
            allp=allp[allp.fs_method.isin(task.fs_method.unique())]
            assert np.isfinite(allp.ScaledValue).all() if expected else not np.isfinite(allp.ScaledValue).any()
            availability.append({'dataset':dataset,'component':comp,'metric':metric,'eligible':expected,'canonical_metric_rows':len(part),'canonical_finite_rows':finite,'all_budget_metric_rows':len(allp),'all_budget_finite_rows':int(np.isfinite(allp.ScaledValue).sum())})
            coreweight=0 if not expected or comp=='Alignment' else (.25 if annotated else .5)/len(metrics)
            contributions.append({'dataset':dataset,'component':comp,'metric':metric,'eligible':expected,'within_component_weight':1/len(metrics) if expected else 0,'effective_weight_in_CoreOverall':coreweight,'interpretation':'not included in revised primary composite' if not annotated else 'label-available Visium composite; dataset-specific annotation meaning'})
    role='primary_label_available_Visium' if annotated else 'separate_label_free_Slide_seqV2'
    for t,c,formula,retained in [('Integration','BatchMixing;BiologicalConservation' if annotated else 'BatchMixing','0.5 BM + 0.5 BC' if annotated else 'BM',True),('Clustering','LabelIndependentClustering;LabelAgreement' if annotated else 'LabelIndependentClustering','0.5 LI + 0.5 LA' if annotated else 'LI',True),('Alignment','Alignment' if dataset in align else '', 'mean(Accuracy, oriented Ratio)' if dataset in align else 'NA',dataset in align),('CoreOverall','BatchMixing;BiologicalConservation;LabelIndependentClustering;LabelAgreement' if annotated else 'BatchMixing;LabelIndependentClustering','0.25 BM + 0.25 BC + 0.25 LI + 0.25 LA' if annotated else 'historical only: 0.5 BM + 0.5 LI',annotated)]:
        composition.append({'dataset':dataset,'analysis_role':role,'task':t,'components':c,'formula':formula,'retained_in_revised_presentation':retained,'alignment_enters_CoreOverall':False})
pd.DataFrame(availability).to_csv(OUT/'metric_availability_audit.tsv',sep='\t',index=False)
pd.DataFrame(contributions).to_csv(OUT/'metric_contribution_weights.tsv',sep='\t',index=False)
pd.DataFrame(composition).to_csv(OUT/'dataset_task_composition.tsv',sep='\t',index=False)
task['analysis_role']=np.where(task.dataset.isin(primary),'primary_label_available_Visium','separate_label_free_Slide_seqV2')
task.to_csv(OUT/'audited_canonical_scores.tsv',sep='\t',index=False)
long=task.melt(id_vars=['dataset','fs_method','n_features','integration_method','analysis_role'],value_vars=cols,var_name='endpoint',value_name='score')
long=long[long.score.notna()].copy()
# Embryo historical aggregates are retained in the audit, not promoted as revised biological endpoints.
long=long[~((~long.dataset.isin(primary))&long.endpoint.isin(['CoreOverall','Integration','Clustering']))]
long['within_dataset_rank']=long.groupby(['dataset','integration_method','endpoint']).score.rank(ascending=False,method='average')
long['method_label']=long.fs_method.map(m.METHOD_LABELS)
long.to_csv(OUT/'task_dataset_scores_and_ranks.tsv',sep='\t',index=False)
leaders=long[long.dataset.isin(primary)&long.endpoint.isin(['Integration','Clustering','Alignment'])&long.within_dataset_rank.eq(1)].copy()
leaders.to_csv(OUT/'task_dataset_observed_leaders.tsv',sep='\t',index=False)
base=task[task.dataset.isin(primary)]
ranks=m.global_ranks(base,'five_label_available_Visium_secondary_summary')
ranks.to_csv(OUT/'five_dataset_secondary_global_ranks.tsv',sep='\t',index=False)
prior=pd.read_csv(ROOT/'results/revision_minor1_3_20261106/analysis/study_bootstrap_and_leave_one_study_out.tsv',sep='\t')
exact=np.array(list(itertools.product(range(5),repeat=5)))
uncert=[]
for integ,p in base.groupby('integration_method'):
    mat=p.pivot(index='dataset',columns='fs_method',values='CoreOverall').loc[primary]
    v=rankdata(-mat.to_numpy(),method='average',axis=1)
    observed=rankdata(v.mean(axis=0),method='average')
    boot=rankdata(v[exact].mean(axis=1),method='average',axis=1)
    loo=np.stack([rankdata(v[np.arange(5)!=i].mean(axis=0),method='average') for i in range(5)])
    for j,method in enumerate(mat.columns):
        row={'integration_method':integ,'fs_method':method,'method_label':m.METHOD_LABELS[method],'ObservedMeanDatasetRank':v[:,j].mean(),'ObservedGlobalRank':observed[j],'Lower95':np.quantile(boot[:,j],.025),'Upper95':np.quantile(boot[:,j],.975),'LOSO_MinRank':loo[:,j].min(),'LOSO_MaxRank':loo[:,j].max(),'NStudies':5,'NExactResamples':len(exact)}
        old=prior[prior.integration_method.eq(integ)&prior.fs_method.eq(method)].iloc[0]
        for k in ['ObservedMeanDatasetRank','ObservedGlobalRank','Lower95','Upper95','LOSO_MinRank','LOSO_MaxRank']:assert abs(row[k]-old[k])<1e-10,(method,k,row[k],old[k])
        uncert.append(row)
u=pd.DataFrame(uncert);u.to_csv(OUT/'primary_rank_uncertainty.tsv',sep='\t',index=False)
with (OUT/'component_definitions.json').open('w') as f:json.dump(m.COMPONENT_METRICS,f,indent=2)
sources=[ROOT/'scripts/build_frozen_scores.py',ROOT/'scripts/build_submission_robustness_tables.py',FROZEN/'scaled_setting_metrics.tsv',FROZEN/'representative_settings.tsv',FROZEN/'representative_task_scores.tsv',FROZEN/'frozen_metric_ranges.tsv',ROOT/'results/revision_minor1_3_20261106/analysis/study_bootstrap_and_leave_one_study_out.tsv']
(OUT/'source_hashes.json').write_text(json.dumps({str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources},indent=2))
audit={'canonical_settings_reproduced':len(check),'components_and_scores_per_setting':len(cols),'dataset_metric_combinations_audited':len(availability),'dataset_task_rows':len(composition),'primary_rank_intervals_reproduced':len(u),'exact_study_resamples':len(exact),'new_model_fits':0,'alignment_excluded_from_CoreOverall':True,'all_checks_passed':True}
(OUT/'analysis_audit.json').write_text(json.dumps(audit,indent=2))
print(json.dumps(audit));print(u[u.ObservedGlobalRank.le(5)].to_string(index=False));print(leaders[['dataset','integration_method','endpoint','method_label']].to_string(index=False))

from pathlib import Path
import json,hashlib,importlib.util,itertools,shutil
import numpy as np
import pandas as pd
from scipy.stats import rankdata,spearmanr
ROOT=Path('/home/xd/y/Benchmark');OUT=ROOT/'results/revision_r2_major6_20261107';OUT.mkdir(parents=True,exist_ok=True)
FROZEN=ROOT/'results/spatial_svg_rebuild_v1/frozen_scores';PRIOR=ROOT/'results/revision_minor1_3_20261106/analysis'
spec=importlib.util.spec_from_file_location('robust',ROOT/'scripts/build_submission_robustness_tables.py');m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
ds=['DLPFC','MouseBrainSerialSections','STOmics0212','STOmics0218','STOmics0224','E8p5Embryo','E9p5Embryo'];primary=ds[:5]
configs=dict(zip(ds,['dlpfc','mouse_brain_serial_sections','stomics_0212','stomics_0218','stomics_0224','e8p5_embryo','e9p5_embryo']))
meta=pd.read_csv(FROZEN/'setting_metric_seed_summary.tsv',sep='\t',dtype={'n_features':str})
groups=pd.read_csv(FROZEN/'representative_settings.tsv',sep='\t').set_index('fs_method').MethodGroup.to_dict()
counts=meta.groupby(['dataset','fs_method','n_features'],as_index=False).agg(actual_min=('EffectiveFeaturesMin','min'),actual_max=('EffectiveFeaturesMax','max'),evaluation_seed_count_max=('NSeeds','max'))
records=[]
for r in counts.itertuples():
    folder=ROOT/'results/spatial_svg_rebuild_v1'/configs[r.dataset]/r.dataset.lower()/'feature_selection'/r.fs_method/('n1' if r.n_features=='all' else 'n'+r.n_features)
    files=sorted(folder.glob('seed*/selected_features.json'));assert files,(r,folder)
    sizes=[]
    for p in files:
        names=json.loads(p.read_text())['feature_names'];n=len(names);assert n==len(set(names)),p
        sizes.append(n);records.append({'dataset':r.dataset,'fs_method':r.fs_method,'requested_features':r.n_features,'selection_seed':int(p.parent.name.removeprefix('seed')),'actual_features':n,'path':str(p.relative_to(ROOT)),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()})
    assert min(sizes)==r.actual_min and max(sizes)==r.actual_max,(r,sizes)
counts=counts.rename(columns={'n_features':'requested_features'})
counts['method_group']=counts.fs_method.map(groups)
counts['role']=np.where(counts.method_group.isin(m.COMPETITIVE_GROUPS),'competitive','descriptive_reference')
counts['numeric_request']=pd.to_numeric(counts.requested_features,errors='coerce')
counts['exact_request']=counts.actual_min.eq(counts.numeric_request)&counts.actual_max.eq(counts.numeric_request)
counts['scope']=np.where(counts.dataset.isin(primary),'primary_Visium','separate_embryo')
counts.to_csv(OUT/'requested_returned_counts_all_datasets.tsv',sep='\t',index=False)
pd.DataFrame(records).to_csv(OUT/'selected_gene_file_audit.tsv',sep='\t',index=False)
coverage=[]
for dataset in ds:
    for method in sorted(groups):
        for request in ['100','200','500','1000','2000','5000','10000','all']:
            q=counts[counts.dataset.eq(dataset)&counts.fs_method.eq(method)&counts.requested_features.eq(request)]
            if len(q):status='evaluated_original_grid'
            elif method=='wilcoxon' and dataset not in primary:status='not_applicable_no_labels'
            elif method=='all_features' or request=='all':status='not_applicable_budget_definition'
            else:status='not_run_original_grid'
            coverage.append({'dataset':dataset,'fs_method':method,'requested_features':request,'status':status,'actual_min':q.iloc[0].actual_min if len(q) else np.nan,'actual_max':q.iloc[0].actual_max if len(q) else np.nan})
pd.DataFrame(coverage).to_csv(OUT/'complete_grid_coverage.tsv',sep='\t',index=False)
comp=counts[counts.role.eq('competitive')&counts.dataset.isin(primary)]
assert len(comp)==1015
old=pd.read_csv(PRIOR/'requested_vs_effective_feature_counts.tsv',sep='\t')
check=comp.merge(old,left_on=['dataset','fs_method','numeric_request'],right_on=['dataset','fs_method','N'],validate='one_to_one')
assert len(check)==1015
np.testing.assert_array_equal(check.actual_min,check.EffectiveFeaturesMin);np.testing.assert_array_equal(check.actual_max,check.EffectiveFeaturesMax)
elig=comp.groupby('fs_method',as_index=False).agg(all_budget_exact=('exact_request','all'),nonexact_settings=('exact_request',lambda x:int((~x).sum())))
at2000=comp[comp.numeric_request.eq(2000)].groupby('fs_method').exact_request.all()
elig['all_primary_exact_at_2000']=elig.fs_method.map(at2000)
elig.to_csv(OUT/'competitive_panel_eligibility.tsv',sep='\t',index=False)
fixed=elig.loc[elig.all_budget_exact,'fs_method'].tolist();exact2000=elig.loc[elig.all_primary_exact_at_2000,'fs_method'].tolist()
assert len(fixed)==25 and len(exact2000)==27
comp[~comp.exact_request].to_csv(OUT/'primary_competitive_shortfalls.tsv',sep='\t',index=False)
canonical=comp[comp.numeric_request.eq(2000)].copy();canonical.to_csv(OUT/'primary_2000_counts.tsv',sep='\t',index=False)

# Verify the separately completed Random2000 reference using its original unchanged selector.
randomrows=[]
for dataset in primary:
    folder=ROOT/'results/revision_major2_random2000_20261105/runs'/configs[dataset]/dataset.lower()/'feature_selection/random/n2000'
    files=sorted(folder.glob('seed*/selected_features.json'));assert len(files)==3
    for p in files:
        names=json.loads(p.read_text())['feature_names'];assert len(names)==len(set(names))==2000
        randomrows.append({'dataset':dataset,'fs_method':'random','requested_features':2000,'selection_seed':p.parent.name,'actual_features':2000,'role':'budget_matched_reference_not_competitive','path':str(p.relative_to(ROOT)),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()})
pd.DataFrame(randomrows).to_csv(OUT/'random2000_count_audit.tsv',sep='\t',index=False)

rep=pd.read_csv(FROZEN/'representative_task_scores.tsv',sep='\t',dtype={'n_features':str})
rep=rep[rep.dataset.isin(primary)&rep.MethodGroup.isin(m.COMPETITIVE_GROUPS)]
rep=rep.rename(columns={c+'Mean':c for c in ['BatchMixing','BiologicalConservation','LabelIndependentClustering','LabelAgreement','Integration','Clustering','Alignment','CoreOverall']})
assert len(rep)==290 and rep.n_features.eq('2000').all()
allmethods=sorted(rep.fs_method.unique());panels={'all_29_nominal_2000':allmethods,'exact_27_at_2000':exact2000,'fixed_25_all_budgets':fixed}
rankrows=[];taskrows=[];exact=np.array(list(itertools.product(range(5),repeat=5)))
for panel,methods in panels.items():
    part=rep[rep.fs_method.isin(methods)].copy()
    long=part.melt(id_vars=['dataset','fs_method','integration_method'],value_vars=['Integration','Clustering','Alignment','CoreOverall'],var_name='endpoint',value_name='score').dropna(subset=['score'])
    long['panel']=panel;long['within_dataset_rank']=long.groupby(['dataset','integration_method','endpoint']).score.rank(ascending=False,method='average');taskrows.append(long)
    for integ,q in part.groupby('integration_method'):
        mat=q.pivot(index='dataset',columns='fs_method',values='CoreOverall').loc[primary,sorted(methods)]
        v=rankdata(-mat.to_numpy(),axis=1,method='average');mean=v.mean(axis=0);obs=rankdata(mean,method='average')
        boot=rankdata(v[exact].mean(axis=1),axis=1,method='average')
        for j,method in enumerate(mat.columns):rankrows.append({'panel':panel,'integration_method':integ,'fs_method':method,'method_label':m.METHOD_LABELS[method],'n_competitors':len(methods),'mean_dataset_rank':mean[j],'global_rank':obs[j],'lower95':np.quantile(boot[:,j],.025),'upper95':np.quantile(boot[:,j],.975),'n_studies':5})
ranks=pd.DataFrame(rankrows);ranks.to_csv(OUT/'canonical_rank_sensitivity.tsv',sep='\t',index=False);pd.concat(taskrows).to_csv(OUT/'task_rank_sensitivity.tsv',sep='\t',index=False)
prior=pd.read_csv(ROOT/'results/revision_r2_major5_20261107/primary_rank_uncertainty.tsv',sep='\t')
test=ranks[ranks.panel.eq('all_29_nominal_2000')].merge(prior,on=['integration_method','fs_method'])
np.testing.assert_allclose(test.global_rank,test.ObservedGlobalRank);np.testing.assert_allclose(test.lower95,test.Lower95);np.testing.assert_allclose(test.upper95,test.Upper95)
concord=[]
for panel in ['exact_27_at_2000','fixed_25_all_budgets']:
    for integ in ['scvi','cellcharter']:
        new=ranks[ranks.panel.eq(panel)&ranks.integration_method.eq(integ)]
        original=ranks[ranks.panel.eq('all_29_nominal_2000')&ranks.integration_method.eq(integ)]
        pair=new.merge(original,on='fs_method',suffixes=('_new','_old'))
        concord.append({'panel':panel,'integration_method':integ,'n_retained_methods':len(pair),'spearman_retained_order':spearmanr(pair.global_rank_new,pair.global_rank_old).statistic,'interpretation':'ordering among retained methods; rank numbers have different candidate-set denominators'})
pd.DataFrame(concord).to_csv(OUT/'rank_order_concordance.tsv',sep='\t',index=False)

scores=pd.read_csv(FROZEN/'setting_task_scores.tsv',sep='\t');scores=scores[scores.dataset.isin(primary)&scores.MethodGroup.isin(m.COMPETITIVE_GROUPS)].copy();scores['N']=pd.to_numeric(scores.n_features)
budgets=[100,200,500,1000,2000,5000,10000];resample=np.random.default_rng(20261106).integers(0,5,size=(5000,5))
curves=[];contrasts=[]
for panel,methods in [('all_29_requested_budget',allmethods),('fixed_exact_budget_panel',fixed)]:
    part=scores[scores.fs_method.isin(methods)]
    bystudy=part.groupby(['dataset','integration_method','N'],as_index=False).agg(score=('CoreOverallMean','mean'),methods=('fs_method','nunique'));assert bystudy.methods.eq(len(methods)).all()
    for integ,q in bystudy.groupby('integration_method'):
        mat=q.pivot(index='dataset',columns='N',values='score').loc[primary,budgets].to_numpy();boot=mat[resample].mean(axis=1)
        for j,n in enumerate(budgets):curves.append({'panel':panel,'integration_method':integ,'N':n,'NMethods':len(methods),'Mean':mat[:,j].mean(),'Lower95':np.quantile(boot[:,j],.025),'Upper95':np.quantile(boot[:,j],.975)})
        for a,b in [(500,1000),(500,2000),(1000,2000),(2000,10000),(500,10000)]:
            delta=boot[:,budgets.index(b)]-boot[:,budgets.index(a)]
            contrasts.append({'panel':panel,'integration_method':integ,'N_from':a,'N_to':b,'MeanDifference':mat[:,budgets.index(b)].mean()-mat[:,budgets.index(a)].mean(),'Lower95':np.quantile(delta,.025),'Upper95':np.quantile(delta,.975),'NMethods':len(methods)})
curves=pd.DataFrame(curves);contrasts=pd.DataFrame(contrasts)
for frame,file,keys,vals in [(curves,'budget_curves_with_study_uncertainty.tsv',['panel','integration_method','N'],['Mean','Lower95','Upper95']),(contrasts,'paired_budget_contrasts.tsv',['panel','integration_method','N_from','N_to'],['MeanDifference','Lower95','Upper95'])]:
    old=pd.read_csv(PRIOR/file,sep='\t');pair=frame.merge(old,on=keys,suffixes=('_new','_old'));assert len(pair)==len(frame)
    for v in vals:np.testing.assert_allclose(pair[v+'_new'],pair[v+'_old'],rtol=1e-10,atol=1e-12)
    frame.to_csv(OUT/file,sep='\t',index=False)
(OUT/'figures').mkdir(exist_ok=True)
shutil.copy2(ROOT/'results/revision_minor1_3_20261106/figures/budget_exact_panel_sensitivity.pdf',OUT/'figures/budget_exact_panel_sensitivity.pdf')
profile=ROOT/'results/resource_profile_v1/feature_selection_resource_profile.tsv'
shutil.copy2(profile,OUT/'standardized_profile_counts_reference.tsv')
sourcefiles=[FROZEN/'setting_metric_seed_summary.tsv',FROZEN/'representative_settings.tsv',FROZEN/'representative_task_scores.tsv',FROZEN/'setting_task_scores.tsv',PRIOR/'requested_vs_effective_feature_counts.tsv',PRIOR/'budget_curves_with_study_uncertainty.tsv',PRIOR/'paired_budget_contrasts.tsv',ROOT/'scripts/build_submission_robustness_tables.py',profile]
(OUT/'source_hashes.json').write_text(json.dumps({str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in sourcefiles},indent=2))
audit={'all_original_dataset_method_budget_conditions':len(counts),'selected_gene_files_checked':len(records),'primary_competitive_conditions_reproduced':len(check),'primary_nonexact_conditions':int((~comp.exact_request).sum()),'exact_at2000_methods':len(exact2000),'fixed_all_budget_exact_methods':len(fixed),'random2000_selection_files_checked':len(randomrows),'original_Major5_rank_rows_reproduced':len(test),'budget_curve_rows_reproduced':len(curves),'paired_contrasts_reproduced':len(contrasts),'new_model_fits':0,'all_checks_passed':True}
(OUT/'analysis_audit.json').write_text(json.dumps(audit,indent=2))
(OUT/'protocol.json').write_text(json.dumps({'primary_datasets':primary,'budgets':budgets,'exact_2000_panel':exact2000,'fixed_all_budget_panel':fixed,'ranking':'Within-dataset ranks recomputed among each fixed eligible panel, then averaged equally across five datasets/studies. Exact 3125 study resamples.','budget':'Same fixed 25 methods at all budgets; 5000 paired study resamples, seed 20261106; original frozen scaling; 3 evaluation seeds at2000 and1 atotherbudgets.','references':'Original references and updated Random2000 separately audited; not promoted to competitive selectors.','missing':"Unrun combinations are not created; all-features has request 'all' and storage folder n1.",'no_selector_padding_or_model_changes':True},indent=2))
print(json.dumps(audit));print(elig[~elig.all_budget_exact].to_string(index=False));print(pd.DataFrame(concord).to_string(index=False));print(ranks[ranks.global_rank.le(3)].to_string(index=False));print(contrasts[contrasts.panel.eq('fixed_exact_budget_panel')].to_string(index=False))

from pathlib import Path
import json,hashlib,itertools
import numpy as np,pandas as pd
from scipy.stats import spearmanr

ROOT=Path('/home/xd/y/Benchmark');OUT=ROOT/'results/revision_r2_m1_scale_match_20261111';A=OUT/'analysis'
jobs=json.loads((OUT/'jobs.json').read_text());required=['bASW','iLISI','silhouette','CHAOS','PAS'];rows=[];coverage=[];sets={}
for j in jobs:
    base=Path(j['output_dir'])/j['dataset_name'].lower()
    for seed in [0,1,2]:
        sp=base/'feature_selection'/j['method']/'n2000'/f'seed{seed}'/'selected_features.json'
        if not sp.exists():sp=base/'feature_selection'/j['method']/'n2000/seed0/selected_features.json'
        if sp.exists():
            payload=json.loads(sp.read_text());sets[(j['arm'],j['method'],seed)]=payload
        for flow in ['scvi','cellcharter']:
            folder=base/flow/j['method']/'n2000'/f'seed{seed}';values=[]
            for task in ['integration_eval','clustering_eval']:
                p=folder/f'{task}_records.json'
                if p.exists():
                    for row in json.loads(p.read_text()):values.append({**row,'arm':j['arm'],'source_path':str(p.relative_to(ROOT))})
            present={r['metric_name'] for r in values if r['metric_name'] in required and r['metric_value'] is not None and np.isfinite(r['metric_value'])}
            coverage.append({'arm':j['arm'],'method':j['method'],'workflow':flow,'seed':seed,'complete':set(required)==present,'missing':';'.join(sorted(set(required)-present))})
            rows+=values
coverage=pd.DataFrame(coverage);coverage.to_csv(A/'run_coverage.tsv',sep='\t',index=False)
# Only an explicitly diagnosed terminal failure can be omitted. Missing/running
# jobs must still stop analysis; unavailable methods are never assigned a rank.
fp=A/'documented_method_failures.json'
failures=json.loads(fp.read_text()) if fp.exists() else []
allowed={(x['arm'],x['method']) for x in failures}
unexpected=coverage[~coverage.complete & ~coverage.apply(lambda r:(r.arm,r.method) in allowed,axis=1)]
if len(unexpected):
    print(coverage.groupby('arm').complete.agg(['sum','count']).to_string());print('Incomplete',unexpected.to_dict('records'));raise SystemExit(2)
complete_methods=coverage.groupby('method').complete.all();complete_methods=complete_methods.index[complete_methods].tolist()
d=pd.DataFrame(rows);d=d[d.metric_name.isin(required)].copy()
key=['arm','fs_method','integration_method','random_seed','metric_name'];assert len(d)==5*int(coverage.complete.sum()) and not d.duplicated(key).any()
assert d.metric_value.notna().all() and np.isfinite(d.metric_value).all()
d.to_csv(A/'seed_level_metrics.tsv',sep='\t',index=False)
refs=['random','all_features'];assert set(refs)<=set(complete_methods)
s=d.groupby(['arm','fs_method','integration_method','metric_name'],as_index=False).agg(mean=('metric_value','mean'),seed_sd=('metric_value','std'),n_seeds=('random_seed','nunique'),actual_min=('effective_n_features','min'),actual_max=('effective_n_features','max'))
assert (s.n_seeds==3).all();s['role']=np.where(s.fs_method.isin(refs),'reference','competitive');s['direction']=np.where(s.metric_name.isin(['CHAOS','PAS']),'lower','higher')
s.to_csv(A/'metric_summary.tsv',sep='\t',index=False)
c=s[s.role.eq('competitive') & s.fs_method.isin(complete_methods)].copy();group=['arm','integration_method','metric_name']
c['rank']=c.groupby(group)['mean'].transform(lambda x:x.rank(ascending=(c.loc[x.index,'direction'].iloc[0]=='lower'),method='average'))
c.to_csv(A/'endpoint_ranks.tsv',sep='\t',index=False)
counts=s.groupby(['arm','fs_method','role'],as_index=False).agg(actual_min=('actual_min','min'),actual_max=('actual_max','max'))
counts.to_csv(A/'actual_gene_counts.tsv',sep='\t',index=False)
valid=counts[counts.role.eq('competitive')].groupby('fs_method').agg(lo=('actual_min','min'),hi=('actual_max','max'),arms=('arm','nunique'))
exact=valid.index[(valid.lo==2000)&(valid.hi==2000)&(valid.arms==2)].tolist()
cor=[];pair=[]
panels=[('paired_complete',sorted(c.fs_method.unique())),('exact_2000_both',exact),('paired_excluding_nnsvg',sorted(set(c.fs_method)-{'nnsvg'}))]
for panel,methods in panels:
    cc=c[c.fs_method.isin(methods)].copy()
    cc['panel_rank']=cc.groupby(group)['mean'].transform(lambda x:x.rank(ascending=(cc.loc[x.index,'direction'].iloc[0]=='lower'),method='average'))
    cc['panel']=panel
    if panel=='exact_2000_both':cc.to_csv(A/'exact_budget_ranks.tsv',sep='\t',index=False)
    for (flow,metric),part in cc.groupby(['integration_method','metric_name']):
        w=part.pivot(index='fs_method',columns='arm',values='panel_rank');rho=float(spearmanr(w.fine_footprint,w.pseudo_visium).statistic)
        topfine=set(w.nsmallest(5,'fine_footprint',keep='all').index);topcoarse=set(w.nsmallest(5,'pseudo_visium',keep='all').index)
        cor.append({'panel':panel,'n_methods':len(methods),'workflow':flow,'metric':metric,'rank_rho':rho,'median_absolute_rank_shift':float((w.pseudo_visium-w.fine_footprint).abs().median()),'top5_tie_inclusive_overlap':len(topfine&topcoarse),'fine_top5_size':len(topfine),'pooled_top5_size':len(topcoarse)})
        for method,row in w.iterrows():pair.append({'panel':panel,'workflow':flow,'metric':metric,'method':method,'fine_rank':row.fine_footprint,'pseudo_rank':row.pseudo_visium,'rank_change':row.pseudo_visium-row.fine_footprint})
pd.DataFrame(cor).to_csv(A/'paired_rank_correlations.tsv',sep='\t',index=False);pd.DataFrame(pair).to_csv(A/'paired_rank_changes.tsv',sep='\t',index=False)
# Separate scale sensitivity from ordinary technical-seed rank variation.
seed_agreement=[]
for panel,methods in panels[:2]:
    sd=d[d.fs_method.isin(methods)].copy()
    seed_group=['arm','integration_method','metric_name','random_seed']
    sd['seed_rank']=sd.groupby(seed_group).metric_value.transform(lambda x:x.rank(ascending=sd.loc[x.index,'metric_name'].iloc[0] in ['CHAOS','PAS']))
    for (flow,metric),part in sd.groupby(['integration_method','metric_name']):
        w=part.pivot(index='fs_method',columns=['arm','random_seed'],values='seed_rank')
        for seed in [0,1,2]:
            seed_agreement.append({'panel':panel,'workflow':flow,'metric':metric,'comparison':'fine_vs_pooled_same_seed','arm':'paired','seed_a':seed,'seed_b':seed,'rho':float(spearmanr(w[('fine_footprint',seed)],w[('pseudo_visium',seed)]).statistic)})
        for arm in ['fine_footprint','pseudo_visium']:
            for sa,sb in itertools.combinations([0,1,2],2):
                seed_agreement.append({'panel':panel,'workflow':flow,'metric':metric,'comparison':'within_arm_technical_seeds','arm':arm,'seed_a':sa,'seed_b':sb,'rho':float(spearmanr(w[(arm,sa)],w[(arm,sb)]).statistic)})
pd.DataFrame(seed_agreement).to_csv(A/'seed_rank_agreement.tsv',sep='\t',index=False)
# Requested sets are compared seed by seed; deterministic selections reuse seed 0.
overlaps=[]
for method in sorted(complete_methods):
    for seed in [0,1,2]:
        def names(payload):
            if isinstance(payload,list):return set(payload)
            return set(payload['feature_names'])
        x=names(sets[('fine_footprint',method,seed)]);y=names(sets[('pseudo_visium',method,seed)])
        if method in refs:assert x==y,(method,seed,'reference gene sets must match across paired arms')
        overlaps.append({'method':method,'seed':seed,'fine_count':len(x),'pseudo_count':len(y),'intersection':len(x&y),'union':len(x|y),'jaccard':len(x&y)/len(x|y)})
pd.DataFrame(overlaps).to_csv(A/'selection_overlap.tsv',sep='\t',index=False)
# Same frozen scaling is a secondary component summary, never a primary embryo score.
frozen=ROOT/'results/spatial_svg_rebuild_v1/frozen_scores/frozen_metric_ranges.tsv'
fr=pd.read_csv(frozen,sep='\t');(A/'frozen_range_schema.json').write_text(json.dumps({'columns':fr.columns.tolist(),'path':str(frozen),'sha256':hashlib.sha256(frozen.read_bytes()).hexdigest()},indent=2))
# Descriptive reference: within-dataset endpoint ranks from the existing five Visium studies.
orig=pd.read_csv(ROOT/'results/revision_r2_major1_20261106/analysis/common_label_free_metric_scores.tsv',sep='\t')
orig=orig[orig.metric_name.isin(required)].copy()
orig['rank']=orig.groupby(['dataset','integration_method','metric_name']).ScaledValue.rank(ascending=False,method='average')
context=[]
for panel,methods in panels:
    vis0=orig[(~orig.dataset.str.contains('Embryo'))&orig.fs_method.isin(methods)].copy()
    vis0['rank']=vis0.groupby(['dataset','integration_method','metric_name']).ScaledValue.rank(ascending=False,method='average')
    vis=vis0.groupby(['fs_method','integration_method','metric_name'],as_index=False)['rank'].mean().rename(columns={'rank':'mean_dataset_rank'})
    for (arm,flow,metric),part in c[c.fs_method.isin(methods)].groupby(group):
        v=vis[(vis.integration_method==flow)&(vis.metric_name==metric)&vis.fs_method.isin(methods)]
        w=part.merge(v,on=['fs_method','integration_method','metric_name'],validate='one_to_one')
        assert len(w)==len(methods)
        context.append({'panel':panel,'arm':arm,'workflow':flow,'metric':metric,'n_methods':len(w),'rho_with_five_visium_mean_ranks':float(spearmanr(w['rank'],w.mean_dataset_rank).statistic),'interpretation':'descriptive reference; tissue, development and annotation differ'})
pd.DataFrame(context).to_csv(A/'visium_reference_rank_correlations.tsv',sep='\t',index=False)
protocol=json.loads((OUT/'protocol.json').read_text());top=c[c['rank']==1].to_dict('records')
hashes=pd.read_csv(A/'pipeline_source_hashes.tsv',sep='\t')
pipeline_unchanged=all(hashlib.sha256((ROOT/r.path).read_bytes()).hexdigest()==r.sha256 for r in hashes.itertuples())
assert pipeline_unchanged
summary={'planned_workflow_seed_runs':len(coverage),'complete_workflow_seed_runs':int(coverage.complete.sum()),'finite_required_metric_rows':len(d),'planned_competitive_methods':29,'paired_complete_competitive_methods':len(c.fs_method.unique()),'documented_failures':failures,'reference_methods':refs,'exact_count_methods':exact,'n_exact_count_methods':len(exact),'count_shortfalls':counts[(counts.role=='competitive')&((counts.actual_min!=2000)|(counts.actual_max!=2000))].to_dict('records'),'rank_correlations':cor,'endpoint_winners':top,'counts_conserved':protocol['counts_conserved_per_gene'],'original_input_unchanged':hashlib.sha256(Path(protocol['source']).read_bytes()).hexdigest()==protocol['source_sha256'],'pipeline_source_unchanged':pipeline_unchanged,'three_seeds_are_technical_repeats':True,'coordinate_convention_conditional':True,'primary_ranking_unchanged':True,'nnsvg_support_note':'Original wrapper skips pucks with fewer than 65 observations: pooled GSM5915060_201112_03 has 60, versus 741 fine beads. Additional paired_excluding_nnsvg sensitivity is reported.'}
assert summary['original_input_unchanged']
(A/'evidence_summary.json').write_text(json.dumps(summary,indent=2))
print(json.dumps({'runs':summary['complete_workflow_seed_runs'],'exact_methods':len(exact),'correlations':cor,'shortfalls':summary['count_shortfalls']},indent=2))

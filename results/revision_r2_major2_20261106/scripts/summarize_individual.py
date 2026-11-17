"""Audit every expected result, independently rescore predictions, aggregate by individual."""
from pathlib import Path
import collections,hashlib,json
import numpy as np
import pandas as pd
from sklearn.metrics import adjusted_rand_score,normalized_mutual_info_score

ROOT=Path('/home/xd/y/Benchmark');OUT=ROOT/'results/revision_r2_major2_20261106'
AN=OUT/'analysis';AN.mkdir(exist_ok=True)
spec=json.loads((OUT/'analysis_protocol.json').read_text())
metrics=[];statuses=[];missing=[];manifests=[];prediction_checks=[]
for fold in spec['folds']:
    job=OUT/'folds'/(fold['dataset']+'__individual_'+fold['individual'])
    if (job/'fold.json').exists():
        d=json.loads((job/'fold.json').read_text());assert not(set(d['training_slices'])&set(d['heldout_slices']));assert d['heldout_individual'] not in d['training_individuals'];manifests.append(d)
    for method in spec['methods']:
        p=job/'status'/(method+'.json');m=job/(method+'_section_metrics.tsv')
        if not p.exists():missing.append({'dataset':fold['dataset'],'individual':fold['individual'],'method':method,'status':'not_started'});continue
        status=json.loads(p.read_text());statuses.append(status)
        if status['status']=='not_applicable':
            assert fold['dataset']=='E9p5Embryo' and method=='wilcoxon'
            continue
        if status['status']!='complete' or not m.exists():missing.append(status);continue
        frame=pd.read_csv(m,sep='\t');assert set(frame.heldout_slice.astype(str))==set(fold['heldout_slices'])
        assert frame.heldout_individual.astype(str).eq(fold['individual']).all()
        if method=='random':assert frame.actual_features.eq(2000).all()
        if method!='all_features':assert frame.actual_features.le(2000).all()
        metrics.append(frame)
        if fold['dataset'] in ['DLPFC','STOmics0212']:
            assert len(frame)==3*len(fold['heldout_slices'])
            for r in frame.to_dict('records'):
                pred=job/'predictions'/f"{method}__{r['heldout_slice']}__seed{r['evaluation_seed']}.npz"
                x=np.load(pred);valid=x['valid'];ari=adjusted_rand_score(x['truth'][valid],x['pred'][valid]);nmi=normalized_mutual_info_score(x['truth'][valid],x['pred'][valid])
                assert abs(ari-r['heldout_ari'])<1e-12 and abs(nmi-r['heldout_nmi'])<1e-12
                prediction_checks.append({'dataset':fold['dataset'],'individual':fold['individual'],'section':r['heldout_slice'],'method':method,'seed':r['evaluation_seed'],'ari_absolute_error':abs(ari-r['heldout_ari']),'nmi_absolute_error':abs(nmi-r['heldout_nmi'])})
if metrics:
    data=pd.concat(metrics,ignore_index=True);data.to_csv(AN/'all_section_seed_results.tsv',sep='\t',index=False)
    group=['dataset','heldout_individual','heldout_slice','method','role']
    numeric=[c for c in ['actual_features','test_spots','heldout_ari','heldout_nmi','heldout_topn_jaccard','heldout_moran_percentile','heldout_rank_spearman'] if c in data]
    section=data.groupby(group,dropna=False,observed=True)[numeric].mean().reset_index();section.to_csv(AN/'section_summary.tsv',sep='\t',index=False)
    individual=section.groupby(['dataset','heldout_individual','method','role'],observed=True)[numeric].mean().reset_index();individual.to_csv(AN/'individual_summary.tsv',sep='\t',index=False)
    datasets=individual.groupby(['dataset','method','role'],observed=True)[numeric].mean().reset_index();datasets.to_csv(AN/'dataset_equal_individual_summary.tsv',sep='\t',index=False)
    # Display scores only for complete method-by-dataset blocks, with all donors present.
    complete=[]
    for (ds,method),g in individual.groupby(['dataset','method']):
        expected=sum(f['dataset']==ds for f in spec['folds'])
        if len(g)==expected:complete.append({'dataset':ds,'method':method,'individuals':expected})
    pd.DataFrame(complete).to_csv(AN/'complete_dataset_method_blocks.tsv',sep='\t',index=False)
pd.DataFrame(statuses).to_csv(AN/'method_status.tsv',sep='\t',index=False)
pd.DataFrame(prediction_checks).to_csv(AN/'independent_prediction_rescoring.tsv',sep='\t',index=False)
original=json.loads((OUT/'original_source_hashes.json').read_text());changed=[p for p,h in original.items() if hashlib.sha256(Path(p).read_bytes()).hexdigest()!=h]
assert not changed,changed
audit={'expected_individual_folds':len(spec['folds']),'expected_method_fold_status_records':len(spec['folds'])*len(spec['methods']),'expected_applicable_method_fold_jobs':len(spec['folds'])*len(spec['methods'])-sum(f['dataset']=='E9p5Embryo' for f in spec['folds']),'not_applicable_policy':'E9.5 Wilcoxon reference requires biological labels absent from original config; no fabricated labels','status_counts':dict(collections.Counter(r['status'] for r in statuses)),'incomplete':missing,'completed':not missing,'individual_fold_manifests':len(manifests),'independently_rescored_predictions':len(prediction_checks),'original_code_config_manuscript_changed':changed,'mapping_sha256':hashlib.sha256((OUT/'section_individual_mapping.tsv').read_bytes()).hexdigest()}
(AN/'completion_audit.json').write_text(json.dumps(audit,indent=2,default=str));print(json.dumps({k:v for k,v in audit.items() if k!='incomplete'},indent=2))
if missing:print('Incomplete jobs:',len(missing))

"""Additive, fixed-budget external-target retention case study; no original writes."""
from pathlib import Path
import hashlib,json,time,shutil
import numpy as np
import pandas as pd
from sklearn.decomposition import TruncatedSVD
from sklearn.metrics import adjusted_rand_score,normalized_mutual_info_score
from spatial_fs_benchmark.config import load_dataset_config
from spatial_fs_benchmark.data.io import load_dataset
from spatial_fs_benchmark.data.preprocess import preprocess_dataset
from spatial_fs_benchmark.validation.heldout_slice import subset_slices,_normalized_log_matrix,score_spatial_reproducibility
from spatial_fs_benchmark.feature_selection.base import FeatureSelectionResult
from spatial_fs_benchmark.metrics.clustering_metrics import cluster_embedding

ROOT=Path('/home/xd/y/Benchmark')
OUT=ROOT/'results/revision_r2_major4_20261107'
PRIOR=ROOT/'results/revision_r2_major2_20261106'
SOURCE=ROOT/'results/revision_major4_independent_validation_20261106/sources'
BUDGET=2000
PARENTS={'HVG':'scanpy_seurat_v3_batch','SVG':'morans_i'}

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def dump(p,v):
    p.parent.mkdir(parents=True,exist_ok=True)
    p.write_text(json.dumps(v,indent=2,default=str))

def fill(h,s,n):
    """Original balanced-union policy, applied after removal of mandatory targets."""
    out=[];seen=set()
    def add(g):
        if g not in seen and len(out)<n:out.append(g);seen.add(g)
    for g in h[:(n+1)//2]:add(g)
    for g in s[:n//2]:add(g)
    for a,b in zip(h,s):add(a);add(b)
    assert len(out)==n
    return out

def select(strategy,h,s,mandatory):
    m=set(mandatory)
    if len(m)>BUDGET:raise ValueError('Mandatory targets exceed the budget; explicitly revise the budget or target list.')
    if len(m)==BUDGET:return sorted(m)
    n=BUDGET-len(m)
    hh=[g for g in h if g not in m];ss=[g for g in s if g not in m]
    chosen=hh[:n] if strategy=='HVG' else ss[:n] if strategy=='SVG' else fill(hh,ss,n)
    chosen+=sorted(m)
    assert len(chosen)==len(set(chosen))==BUDGET and m.issubset(chosen)
    return chosen

def main():
    OUT.mkdir(parents=True,exist_ok=True)
    spec={'reviewer':'R2 Major 4','case_study':'DLPFC external cortical-marker retention under fixed computational budget','donors':['Br5292','Br8100'],'budget':BUDGET,'strategies':['HVG','SVG','Hybrid'],'conditions':['Unconstrained','Retain targets'],'parents':PARENTS,'hybrid':'Original balanced union: half from each eligible parent ranking; alternate down rankings to replace duplicates. Apply same rule to residual budget after removing mandatory targets. Odd residual: HVG receives the extra initial slot.','mandatory':'All H/MH Zeng markers, fixed external list already archived for R1 Major 4; exact case-insensitive unique gene_name mapping within each training-only gene universe. No outcome-dependent gene selection. Unavailable/ambiguous targets reported, never silently padded.','split':'Reuse R2 Major 2 leave-one-donor-out folds, complete donor excluded including both sections. Frozen parent rankings fitted solely to training donor.','evaluation':'Unchanged R2 Major 2 selected-count normalization, log1p, SVD30, original clustering, seeds 0/1/2, valid source layer annotations. Held-out Moran percentile separately; no artificial rank interpretation for mandatory genes.','aggregation':'Seeds within section, sections within donor, equal donor means. Both donor results shown, no significance test with two donors.','interpretation':'Target retention is satisfied by construction, not independent evidence of biological recovery. ARI/NMI and held-out spatial percentiles quantify downstream consequences in this case study; no universal improvement claim.','source_sha256':sha(SOURCE/'Zeng_Markers.csv')}
    if (OUT/'analysis_protocol.json').exists():assert json.loads((OUT/'analysis_protocol.json').read_text())==spec
    else:dump(OUT/'analysis_protocol.json',spec)
    (OUT/'sources').mkdir(exist_ok=True)
    for f in ['Zeng_Markers.csv','source_manifest.json']:shutil.copy2(SOURCE/f,OUT/'sources'/f)
    source=pd.read_csv(SOURCE/'Zeng_Markers.csv')
    assert len(source)==81 and source.Paper.eq('Zeng').all()
    source=source[source.Species.isin(['H','MH'])].copy()
    cfg=load_dataset_config(ROOT/'configs/datasets/dlpfc.yaml');raw=load_dataset(cfg)
    rows=[];selection_rows=[];mapping_rows=[];hashes={}
    for donor in spec['donors']:
        old=PRIOR/'folds'/('DLPFC__individual_'+donor)
        manifest=json.loads((old/'fold.json').read_text());job=OUT/'folds'/donor;job.mkdir(parents=True,exist_ok=True)
        train=preprocess_dataset(subset_slices(raw,manifest['training_slices'],'DLPFC_train_'+donor),cfg.preprocess)
        universe=train.gene_names
        assert hashlib.sha256('\n'.join(universe).encode()).hexdigest()==manifest['training_gene_hash']
        lookup={}
        for g,sym in zip(universe,train.adata.var.gene_name.astype(str).str.upper()):lookup.setdefault(sym,[]).append(g)
        mandatory=[]
        for r in source.to_dict('records'):
            matches=lookup.get(str(r['Gene']).upper(),[])
            mapped=matches[0] if len(matches)==1 else ''
            if mapped:mandatory.append(mapped)
            mapping_rows.append({'heldout_donor':donor,**r,'mapped_gene':mapped,'mapping_status':'unique_training_gene' if mapped else 'unavailable_or_ambiguous'})
        assert len(mandatory)==len(set(mandatory))
        h=json.loads((old/'rankings'/f'{PARENTS["HVG"]}.json').read_text())['feature_names']
        s=json.loads((old/'rankings'/f'{PARENTS["SVG"]}.json').read_text())['feature_names']
        assert len(h)==len(s)==BUDGET
        for p in [old/'fold.json',old/'training_gene_universe.txt',old/'rankings'/f'{PARENTS["HVG"]}.json',old/'rankings'/f'{PARENTS["SVG"]}.json']:
            hashes[str(p)]=sha(p)
        selections={}
        for strategy in spec['strategies']:
            for retain in [False,True]:
                condition='Retain targets' if retain else 'Unconstrained'
                genes=select(strategy,h,s,mandatory if retain else [])
                selections[(strategy,condition)]=genes
                overlap=len(set(genes)&set(mandatory))
                selection_rows.append({'heldout_donor':donor,'strategy':strategy,'condition':condition,'requested_features':BUDGET,'actual_features':len(genes),'source_targets':len(source),'available_targets':len(mandatory),'retained_targets':overlap,'target_retention':overlap/len(mandatory),'training_genes':len(universe)})
                dump(job/(strategy+('_retain' if retain else '_base')+'_genes.json'),genes)
        for sid in manifest['heldout_slices']:
            test=subset_slices(raw,[sid],'DLPFC_test_'+sid);test.adata=test.adata[:,universe].copy()
            opts=dict(cfg.preprocess);opts['min_cells_per_gene']=0;test=preprocess_dataset(test,opts)
            lab=test.adata.obs[test.label_key];truth=np.array([str(v) for v in test.labels])
            valid=lab.notna().to_numpy() & ~pd.Series(truth).str.strip().str.lower().isin(['','nan','none','unknown','unassigned','na','n/a']).to_numpy()
            k=len(np.unique(truth[valid]));assert valid.sum()>k>1
            refpath=old/'spatial_references'/(sid+'.json');hashes[str(refpath)]=sha(refpath)
            ref=FeatureSelectionResult(**json.loads(refpath.read_text()))
            for (strategy,condition),genes in selections.items():
                sel=FeatureSelectionResult(method_name=strategy,feature_names=genes,feature_indices=[universe.index(g) for g in genes],scores=[0.0]*len(genes))
                spatial=score_spatial_reproducibility(test,sel,ref,BUDGET)
                matrix=_normalized_log_matrix(test,genes)
                for seed in [0,1,2]:
                    emb=TruncatedSVD(n_components=min(30,matrix.shape[0]-1,matrix.shape[1]-1),random_state=seed).fit_transform(matrix)
                    pred=cluster_embedding(emb,n_clusters=k,random_seed=seed)
                    ari=adjusted_rand_score(truth[valid],pred[valid]);nmi=normalized_mutual_info_score(truth[valid],pred[valid])
                    if condition=='Unconstrained' and strategy in PARENTS:
                        oldmetrics=pd.read_csv(old/(PARENTS[strategy]+'_section_metrics.tsv'),sep='\t',dtype={'heldout_slice':str})
                        row=oldmetrics[oldmetrics.heldout_slice.eq(sid)&oldmetrics.evaluation_seed.eq(seed)].iloc[0]
                        assert abs(ari-row.heldout_ari)<1e-10 and abs(nmi-row.heldout_nmi)<1e-10,(donor,sid,strategy,seed,ari,row.heldout_ari)
                    tag=strategy+('_retain' if condition=='Retain targets' else '_base')
                    pp=job/'predictions'/f'{tag}__{sid}__seed{seed}.npz';pp.parent.mkdir(exist_ok=True)
                    np.savez_compressed(pp,truth=truth,pred=pred,valid=valid,spot=np.asarray(test.adata.obs_names).astype(str))
                    rows.append({'heldout_donor':donor,'heldout_slice':sid,'strategy':strategy,'condition':condition,'seed':seed,'actual_features':len(genes),'heldout_ari':ari,'heldout_nmi':nmi,'heldout_moran_percentile':spatial['heldout_moran_percentile'],'heldout_topn_jaccard':spatial['heldout_topn_jaccard']})
            print('Complete',donor,sid,'6 conditions x 3 seeds',flush=True)
    frame=pd.DataFrame(rows);sf=pd.DataFrame(selection_rows)
    frame.to_csv(OUT/'section_seed_metrics.tsv',sep='\t',index=False)
    sf.to_csv(OUT/'selection_audit.tsv',sep='\t',index=False)
    pd.DataFrame(mapping_rows).to_csv(OUT/'target_mapping.tsv',sep='\t',index=False)
    metrics=['heldout_ari','heldout_nmi','heldout_moran_percentile','heldout_topn_jaccard']
    sections=frame.groupby(['heldout_donor','heldout_slice','strategy','condition'],as_index=False)[metrics].mean()
    donors=sections.groupby(['heldout_donor','strategy','condition'],as_index=False)[metrics].mean().merge(sf,on=['heldout_donor','strategy','condition'],validate='one_to_one')
    donors.to_csv(OUT/'donor_summary.tsv',sep='\t',index=False)
    summary=donors.groupby(['strategy','condition'],as_index=False)[metrics+['target_retention','retained_targets','available_targets']].mean()
    summary.to_csv(OUT/'mean_summary.tsv',sep='\t',index=False)
    delta=donors[donors.condition.eq('Retain targets')].merge(donors[donors.condition.eq('Unconstrained')],on=['heldout_donor','strategy'],suffixes=('_retained','_base'))
    for m in metrics+['target_retention','retained_targets']:delta['delta_'+m]=delta[m+'_retained']-delta[m+'_base']
    delta.to_csv(OUT/'paired_donor_differences.tsv',sep='\t',index=False)
    assert all(sha(Path(p))==h for p,h in hashes.items())
    dump(OUT/'input_hashes.json',hashes)
    audit={'status':'complete','section_seed_runs':len(frame),'donors':2,'sections':4,'conditions_per_section':6,'all_actual_counts':sorted(sf.actual_features.unique().tolist()),'parent_baseline_seed_results_reproduced':24,'original_inputs_unchanged':True,'unavailable_target_rows':sum(not r['mapped_gene'] for r in mapping_rows)}
    dump(OUT/'completion.json',audit)
    print(summary.to_string(index=False),flush=True)
    print(json.dumps(audit),flush=True)

if __name__=='__main__':main()

"""Additive donor-blocked feature-transfer validation; frozen originals are read only."""
import argparse,hashlib,json,time,traceback,os
import fcntl
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.decomposition import TruncatedSVD
from sklearn.metrics import adjusted_rand_score,normalized_mutual_info_score
from spatial_fs_benchmark.config import load_dataset_config
from spatial_fs_benchmark.data.io import load_dataset
from spatial_fs_benchmark.data.preprocess import preprocess_dataset
from spatial_fs_benchmark.feature_selection import build_feature_selector
from spatial_fs_benchmark.feature_selection.base import FeatureSelectionResult
from spatial_fs_benchmark.validation.heldout_slice import subset_slices,heldout_moran_reference,score_spatial_reproducibility,_normalized_log_matrix
from spatial_fs_benchmark.metrics.clustering_metrics import cluster_embedding

ROOT=Path('/home/xd/y/Benchmark')
OUT=ROOT/'results/revision_r2_major2_20261106'
CONFIGS={'DLPFC':'dlpfc','STOmics0212':'stomics_0212_wilcoxon','STOmics0218':'stomics_0218_wilcoxon','STOmics0224':'stomics_0224_wilcoxon','E9p5Embryo':'e9p5_embryo'}
# Identical selector implementations and settings to the original held-out panel.
METHODS={'all_features':{},'random':{},'scanpy_seurat_v3_batch':{},'scanpy_cell_ranger_batch':{},'seurat_vst':{},'triku':{'max_cells':5000},'singleCellHaystack':{'max_cells':5000},'scsegindex':{},'morans_i':{'n_neighbors':8},'sparkx':{'n_threads':4},'wilcoxon':{},'nnsvg':{'n_threads':4},'spatialde':{'max_cells_per_slice':800},'somde':{'spots_per_node':20}}
REFERENCES={'all_features','random','scsegindex','wilcoxon'}
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def atomic(p,x):
    p.parent.mkdir(parents=True,exist_ok=True);tmp=p.with_suffix(p.suffix+f'.{os.getpid()}.tmp');tmp.write_text(json.dumps(x,indent=2,default=str));tmp.replace(p)
def table(p,rows):
    tmp=p.with_suffix('.tsv.tmp');pd.DataFrame(rows).to_csv(tmp,sep='\t',index=False);tmp.replace(p)
def freeze():
    mapping=pd.read_csv(OUT/'section_individual_mapping.tsv',sep='\t',dtype=str)
    folds=[]
    for ds,cfg in CONFIGS.items():
        rows=mapping[mapping.dataset.eq(ds)]
        assert rows.individual_id.notna().all()
        ids=sorted(rows.individual_id.unique())
        for donor in ids:
            folds.append({'dataset':ds,'individual':donor,'config':cfg,'heldout_slices':rows.loc[rows.individual_id.eq(donor),'sample_id'].tolist()})
    spec={'reviewer':'R2 Major 2','split':'leave one biological individual out; all retained sections of that individual excluded together','folds':folds,'n_features':2000,'selection_seed':0,'evaluation_seeds':[0,1,2],'methods':METHODS,'reference_methods':sorted(REFERENCES),'gene_filter':'Split raw input first; fit configured minimum-detected-cell filter on training individuals only. Retain the same frozen training gene universe in all test sections. No test expression or labels used for feature selection. Source-provided processed X is retained where specified by the original input contract.','test_processing':'No integration retraining. Original per-section selected-count normalization to 1e4, log1p, TruncatedSVD(30), KMeans/MiniBatchKMeans. Embedding/clustering fit on test section, so this is feature-set transfer, not a trained classifier prediction.','annotation_endpoints':'DLPFC and STOmics0212 only; valid source annotations; K determined from nonmissing annotation classes, as original protocol. Three evaluation seeds for fixed seed-0 features. Source cluster IDs in other datasets are not treated as independent biological truth.','spatial_endpoints':'Original held-out 8NN Moran rank agreement, top-2000 Jaccard and mean percentile; selection methods may themselves target spatial structure, so spatial transfer does not constitute independent biological validation.','aggregation':'Average seeds within section, then sections equally within held-out individual; report all individuals. Average individuals equally within dataset, and dataset summaries separately; do not pool spots/sections as biological replicates or treat overlapping training folds as independent experiments.','scope':'DLPFC 2, skin 7, treatment mouse brain 7, cortex 2, E9.5 3 folds. Serial MouseBrain and retained E8.5 each one individual: cross-section evidence only. E9.5 remains a separate label-free analysis. Condition/sex/time and subject can be confounded; no population-level donor-effect inference.','mapping_sha256':sha(OUT/'section_individual_mapping.tsv')}
    path=OUT/'analysis_protocol.json'
    if path.exists():assert json.loads(path.read_text())==spec,'Protocol has changed'
    else:atomic(path,spec)
    codefiles=list((ROOT/'src/spatial_fs_benchmark/feature_selection').glob('*.py'))+[ROOT/'src/spatial_fs_benchmark/data/preprocess.py',ROOT/'scripts/run_spatial_r_feature_selector.R']
    original=list((ROOT/'manuscript_genome_research_spatial_omics_v2').glob('*.tex'))+list((ROOT/'manuscript_genome_research_spatial_omics_v2/figures').glob('*.pdf'))
    hashes={str(p):sha(p) for p in codefiles+original+[ROOT/'configs/datasets'/f'{c}.yaml' for c in CONFIGS.values()]}
    hp=OUT/'original_source_hashes.json'
    if hp.exists():assert json.loads(hp.read_text())==hashes
    else:atomic(hp,hashes)
    print('Frozen',len(folds),'individual folds and',len(METHODS),'methods',flush=True)

def run(ds,donor,methods):
    mapping=pd.read_csv(OUT/'section_individual_mapping.tsv',sep='\t',dtype=str)
    mm=mapping[mapping.dataset.eq(ds)]
    test_ids=mm.loc[mm.individual_id.eq(donor),'sample_id'].tolist();train_ids=mm.loc[~mm.individual_id.eq(donor),'sample_id'].tolist()
    assert test_ids and train_ids and not(set(test_ids)&set(train_ids))
    cfg=load_dataset_config(ROOT/'configs/datasets'/f'{CONFIGS[ds]}.yaml')
    raw=load_dataset(cfg)
    training=preprocess_dataset(subset_slices(raw,train_ids,ds+'_train_without_'+donor),cfg.preprocess)
    job=OUT/'folds'/(ds+'__individual_'+donor);job.mkdir(parents=True,exist_ok=True)
    genes=training.gene_names;raw_nvars=raw.n_vars
    testopt=dict(cfg.preprocess);testopt['min_cells_per_gene']=0
    tests={}
    for sid in test_ids:
        test=subset_slices(raw,[sid],ds+'_test_'+sid)
        test.adata=test.adata[:,genes].copy()
        tests[sid]=preprocess_dataset(test,testopt)
    del raw
    manifest={'dataset':ds,'heldout_individual':donor,'heldout_slices':test_ids,'training_slices':train_ids,'training_individuals':sorted(mm.loc[~mm.individual_id.eq(donor),'individual_id'].unique()),'no_heldout_individual_in_training':True,'raw_input_genes':raw_nvars,'training_filtered_genes':len(genes),'training_spots':training.n_obs,'heldout_spots':{s:t.n_obs for s,t in tests.items()},'training_gene_hash':hashlib.sha256('\n'.join(genes).encode()).hexdigest()}
    atomic(job/'fold.json',manifest)
    (job/'training_gene_universe.txt').write_text('\n'.join(genes))
    print('Loaded',ds,donor,'train',training.n_obs,'test',manifest['heldout_spots'],'genes',len(genes),flush=True)
    references={}
    for sid,test in tests.items():
        p=job/'spatial_references'/(sid+'.json')
        if p.exists():r=json.loads(p.read_text());ref=FeatureSelectionResult(**r)
        else:ref=heldout_moran_reference(test,random_seed=0);atomic(p,{k:getattr(ref,k) for k in ['method_name','feature_names','feature_indices','scores','metadata']})
        references[sid]=ref
    for method in methods:
        (job/'locks').mkdir(exist_ok=True)
        with (job/'locks'/(method+'.lock')).open('a') as method_lock:
            fcntl.flock(method_lock,fcntl.LOCK_EX)
            statuspath=job/'status'/(method+'.json');metricpath=job/(method+'_section_metrics.tsv')
            if method=='wilcoxon' and training.label_key is None:
                atomic(statuspath,{'status':'not_applicable','dataset':ds,'individual':donor,'method':method,'reason':'Original dataset config provides no biological label_key; label-informed reference is not defined. No labels are fabricated.'})
                print('Not applicable',ds,donor,method,'no source labels',flush=True);continue
            if statuspath.exists() and json.loads(statuspath.read_text()).get('status')=='complete' and metricpath.exists():
                print('Cached',ds,donor,method,flush=True);continue
            start=time.perf_counter()
            atomic(statuspath,{'status':'running','dataset':ds,'individual':donor,'method':method,'started_at':time.time()})
            try:
                p=job/'rankings'/(method+'.json')
                if p.exists():selection=FeatureSelectionResult(**json.loads(p.read_text()))
                else:
                    selector=build_feature_selector(method,**METHODS[method]);selection=selector.select(training,n_features=training.n_vars if method=='all_features' else 2000,random_seed=0)
                    atomic(p,{k:getattr(selection,k) for k in ['method_name','feature_names','feature_indices','scores','metadata']})
                assert len(set(selection.feature_names))==len(selection.feature_names)
                if method!='all_features' and len(selection.feature_names)>2000:raise AssertionError('Budget exceeded')
                rows=[]
                for sid,test in tests.items():
                    spatial=score_spatial_reproducibility(test,selection,references[sid],n_features=2000)
                    base={'dataset':ds,'heldout_individual':donor,'heldout_slice':sid,'method':method,'role':'reference' if method in REFERENCES else 'competitive','actual_features':len(selection.feature_names),'test_spots':test.n_obs,**spatial}
                    if ds not in {'DLPFC','STOmics0212'}:
                        rows.append({**base,'evaluation_seed':-1});continue
                    label=test.adata.obs[test.label_key]
                    truth=np.array([str(v) for v in test.labels]);valid=label.notna().to_numpy() & ~pd.Series(truth).str.strip().str.lower().isin(['','nan','none','unknown','unassigned','na','n/a']).to_numpy()
                    k=len(np.unique(truth[valid]));assert valid.sum()>k>1
                    matrix=_normalized_log_matrix(test,selection.feature_names)
                    for seed in [0,1,2]:
                        embedding=TruncatedSVD(n_components=min(30,matrix.shape[0]-1,matrix.shape[1]-1),random_state=seed).fit_transform(matrix)
                        pred=cluster_embedding(embedding,n_clusters=k,random_seed=seed)
                        pp=job/'predictions'/f'{method}__{sid}__seed{seed}.npz';pp.parent.mkdir(exist_ok=True)
                        np.savez_compressed(pp,spot=np.asarray(test.adata.obs_names).astype(str),truth=truth,pred=pred,valid=valid)
                        rows.append({**base,'evaluation_seed':seed,'annotated_spots':int(valid.sum()),'reference_k':k,'heldout_ari':float(adjusted_rand_score(truth[valid],pred[valid])),'heldout_nmi':float(normalized_mutual_info_score(truth[valid],pred[valid]))})
                table(metricpath,rows)
                atomic(statuspath,{'status':'complete','dataset':ds,'individual':donor,'method':method,'runtime_seconds':time.perf_counter()-start,'actual_features':len(selection.feature_names),'sections_evaluated':len(test_ids)})
                print('Complete',ds,donor,method,round(time.perf_counter()-start,1),'sec',flush=True)
            except Exception as e:
                atomic(statuspath,{'status':'failed','dataset':ds,'individual':donor,'method':method,'error':repr(e),'traceback':traceback.format_exc(),'runtime_seconds':time.perf_counter()-start})
                print('FAILED',ds,donor,method,repr(e),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--freeze',action='store_true');p.add_argument('--dataset');p.add_argument('--individual');p.add_argument('--methods',default=','.join(METHODS));a=p.parse_args()
    if a.freeze:freeze()
    else:run(a.dataset,a.individual,a.methods.split(','))

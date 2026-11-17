from pathlib import Path
import hashlib,json,importlib.util,itertools
import anndata as ad
import numpy as np
import pandas as pd
import yaml
from scipy.spatial import cKDTree,Delaunay
from scipy import sparse
from scipy.stats import spearmanr

ROOT=Path('/home/xd/y/Benchmark')
OUT=ROOT/'results/revision_r2_major1_20261106'
AN=OUT/'analysis';AN.mkdir(exist_ok=True)
FROZEN=ROOT/'results/spatial_svg_rebuild_v1/frozen_scores'
cfgs=['dlpfc','mouse_brain_serial_sections','stomics_0212_wilcoxon','stomics_0218_wilcoxon','stomics_0224_wilcoxon','e8p5_embryo','e9p5_embryo']
paths=[ROOT/'configs/datasets'/f'{x}.yaml' for x in cfgs]
paths += [FROZEN/x for x in ['scaled_setting_metrics.tsv','representative_settings.tsv','setting_task_scores.tsv','frozen_metric_ranges.tsv']]
paths += [ROOT/x for x in ['src/spatial_fs_benchmark/feature_selection/svg.py','src/spatial_fs_benchmark/integration/graphst.py','scripts/run_cellcharter_native.py','src/spatial_fs_benchmark/metrics/spatial_clustering.py','manuscript_genome_research_spatial_omics_v2/main.tex','manuscript_genome_research_spatial_omics_v2/supplement.tex']]
hashes={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
coordinates=[];edges=[];cal=[];counts=[]
for cf in cfgs:
    cfg=yaml.safe_load((ROOT/'configs/datasets'/f'{cf}.yaml').read_text())
    a=ad.read_h5ad(ROOT/cfg['path'],backed='r')
    allxy=np.asarray(a.obsm[cfg['coord_key']],dtype=np.float64)
    platform='Slide-seqV2' if cf.startswith('e') else 'Visium standard/non-HD'
    counts.append({'dataset':cfg['name'],'platform':platform,'n_observations':a.n_obs,'n_slices':a.obs[cfg['slice_key']].nunique(),
                   'unit':'bead' if cf.startswith('e') else 'spot','n_genes_raw_input':a.n_vars})
    for sid,part in a.obs.groupby(cfg['slice_key'],observed=True):
        idx=np.flatnonzero(a.obs[cfg['slice_key']].astype(str).to_numpy()==str(sid));xy=allxy[idx]
        tree=cKDTree(xy);ds,ns=tree.query(xy,k=9)
        metric_xy=None
        if not cf.startswith('e'):
            # Standard Visium's staggered array uses 50 um per array_col step,
            # sqrt(3)*50 um per array_row step (100 um nearest-spot pitch).
            metric_xy=np.column_stack([part.array_col.to_numpy()*50.,part.array_row.to_numpy()*50*np.sqrt(3)])
            design=np.column_stack([xy,np.ones(len(xy))]);coef=np.linalg.lstsq(design,metric_xy,rcond=None)[0]
            fitted=design@coef;residual=np.linalg.norm(fitted-metric_xy,axis=1)
            sv=np.linalg.svd(coef[:2],compute_uv=False)
            cal.append({'dataset':cfg['name'],'sample_id':sid,'basis':'supplied array_row/array_col; standard Visium 100 um pitch',
                        'um_per_native_unit_mean':sv.mean(),'affine_anisotropy_ratio':sv.max()/sv.min(),
                        'median_fit_residual_um':np.median(residual),'max_fit_residual_um':residual.max(),
                        'n_observations':len(xy)})
        row={'dataset':cfg['name'],'sample_id':str(sid),'platform':platform,'n_observations':len(xy),
             'native_nn1_median':np.median(ds[:,1]),'native_nn3_median':np.median(ds[:,3]),'native_nn8_median':np.median(ds[:,8]),
             'metric_calibration_available':metric_xy is not None}
        if metric_xy is not None:
            md=cKDTree(metric_xy).query(metric_xy,k=9)[0]
            row.update({'grid_um_nn1_median':np.median(md[:,1]),'grid_um_nn3_median':np.median(md[:,3]),'grid_um_nn8_median':np.median(md[:,8])})
            for radius in [150,200]:
                neighbor_counts=cKDTree(metric_xy).query_ball_point(metric_xy,radius,return_length=True)-1
                row[f'radius_{radius}um_degree_median']=np.median(neighbor_counts)
        coordinates.append(row)
        # Connectivity is reconstructed for diagnostics only; no training or
        # selection outputs are changed. Distances refer to retained observations.
        for label,k in [('GraphST_3NN',3),('Moran_8NN',8)]:
            ii=np.repeat(np.arange(len(xy)),k);jj=ns[:,1:k+1].ravel()
            pairs=np.unique(np.sort(np.column_stack([ii,jj]),axis=1),axis=0)
            nd=np.linalg.norm(xy[pairs[:,0]]-xy[pairs[:,1]],axis=1)
            er={'dataset':cfg['name'],'sample_id':str(sid),'graph':label,'edges':len(pairs),
                'native_distance_median':np.median(nd),'native_distance_p95':np.quantile(nd,.95),'native_distance_max':nd.max()}
            if metric_xy is not None:
                mm=np.linalg.norm(metric_xy[pairs[:,0]]-metric_xy[pairs[:,1]],axis=1)
                er.update({'grid_um_distance_median':np.median(mm),'grid_um_distance_p95':np.quantile(mm,.95),'grid_um_distance_max':mm.max()})
            edges.append(er)
        tri=Delaunay(xy).simplices
        pairs=np.unique(np.sort(np.concatenate([tri[:,[0,1]],tri[:,[0,2]],tri[:,[1,2]]]),axis=1),axis=0)
        nd=np.linalg.norm(xy[pairs[:,0]]-xy[pairs[:,1]],axis=1)
        er={'dataset':cfg['name'],'sample_id':str(sid),'graph':'CellCharter_Delaunay_diagnostic','edges':len(pairs),
            'native_distance_median':np.median(nd),'native_distance_p95':np.quantile(nd,.95),'native_distance_max':nd.max()}
        if metric_xy is not None:
            mm=np.linalg.norm(metric_xy[pairs[:,0]]-metric_xy[pairs[:,1]],axis=1)
            er.update({'grid_um_distance_median':np.median(mm),'grid_um_distance_p95':np.quantile(mm,.95),'grid_um_distance_max':mm.max()})
        edges.append(er)
    a.file.close()
pd.DataFrame(counts).to_csv(AN/'dataset_observation_units.tsv',sep='\t',index=False)
pd.DataFrame(coordinates).to_csv(AN/'slice_neighborhood_support.tsv',sep='\t',index=False)
pd.DataFrame(edges).to_csv(AN/'graph_edge_support.tsv',sep='\t',index=False)
pd.DataFrame(cal).to_csv(AN/'visium_lattice_calibration_audit.tsv',sep='\t',index=False)

spec=importlib.util.spec_from_file_location('robust',ROOT/'scripts/build_submission_robustness_tables.py')
robust=importlib.util.module_from_spec(spec);spec.loader.exec_module(robust)
_,selected=robust.representative_metrics(FROZEN)
selected=selected[selected.MethodGroup.isin(['expression_driven','spatially_informed'])].copy()
assert selected.fs_method.nunique()==29
datasets=pd.DataFrame(counts).set_index('dataset')
selected['platform']=selected.dataset.map(datasets.platform)
selected.to_csv(AN/'platform_stratified_all_metric_scores.tsv',sep='\t',index=False)
tasks=robust.task_scores(selected,'ScaledValue')
tasks['platform']=tasks.dataset.map(datasets.platform)
common=['bASW','iLISI','silhouette','CHAOS','PAS']
raw=selected[selected.metric_name.isin(common)].copy()
assert raw.ScaledValue.notna().all()
raw.to_csv(AN/'common_label_free_metric_scores.tsv',sep='\t',index=False)
long=tasks.melt(id_vars=['dataset','fs_method','integration_method','n_features','platform'],
                value_vars=['BatchMixing','LabelIndependentClustering'],var_name='Component',value_name='Score')
assert long.Score.notna().all()
long['DatasetRank']=long.groupby(['dataset','integration_method','Component']).Score.rank(ascending=False,method='average')
long.to_csv(AN/'platform_common_component_dataset_scores.tsv',sep='\t',index=False)
summary=long.groupby(['platform','integration_method','Component','fs_method'],as_index=False).agg(
    MeanScore=('Score','mean'),MeanDatasetRank=('DatasetRank','mean'),NDatasets=('dataset','nunique'))
summary['PlatformRank']=summary.groupby(['platform','integration_method','Component']).MeanDatasetRank.rank(method='average')
summary['OriginalStudyUnits']=np.where(summary.platform.str.startswith('Visium'),5,1)
summary['Inference']='descriptive within-platform summary; no cross-platform technology effect estimate'
summary.to_csv(AN/'platform_common_component_method_summary.tsv',sep='\t',index=False)
primary=tasks[tasks.platform.str.startswith('Visium')]
ranks=robust.global_ranks(primary,'Visium_only_primary')
expected=pd.read_csv(ROOT/'results/revision_minor1_3_20261106/analysis/study_bootstrap_and_leave_one_study_out.tsv',sep='\t')
check=ranks.merge(expected[expected.analysis_role.eq('primary')],on=['fs_method','integration_method'])
assert len(check)==58
np.testing.assert_allclose(check.GlobalRank,check.ObservedGlobalRank,atol=1e-12)
np.testing.assert_allclose(check.MeanDatasetRank,check.ObservedMeanDatasetRank,atol=1e-12)
ranks.to_csv(AN/'visium_only_primary_rank_reproduction.tsv',sep='\t',index=False)
pair=[]
for (integ,comp),part in summary.groupby(['integration_method','Component']):
    wide=part.pivot(index='fs_method',columns='platform',values='PlatformRank')
    rho=spearmanr(wide.iloc[:,0],wide.iloc[:,1]).statistic
    pair.append({'integration_method':integ,'Component':comp,'Spearman_descriptive':rho,'NMethods':29,
                 'cross_platform_causal_interpretation':False})
pd.DataFrame(pair).to_csv(AN/'between_platform_rank_description.tsv',sep='\t',index=False)
support=pd.DataFrame(coordinates).groupby('dataset').agg(
    NativeNN1Min=('native_nn1_median','min'),NativeNN1Max=('native_nn1_median','max'),
    GridUmNN8Min=('grid_um_nn8_median','min'),GridUmNN8Max=('grid_um_nn8_median','max'))
support.to_csv(AN/'dataset_neighborhood_summary.tsv',sep='\t')
after={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
assert hashes==after
pd.DataFrame([{'path':p,'sha256_before':v,'sha256_after':after[p]} for p,v in hashes.items()]).to_csv(OUT/'original_input_hashes.tsv',sep='\t',index=False)
audit={'datasets':7,'slices':len(coordinates),'visium_slices_with_grid_calibration':len(cal),
       'slide_seq_metric_unit_calibration_claimed':False,'primary_visium_rank_rows_reproduced':len(check),
       'competitive_selectors':29,'platform_component_summary_rows':len(summary),'same_common_metrics_across_platforms':common,
       'new_model_fits':0,'original_inputs_unchanged':True,'all_checks_passed':True}
(AN/'final_audit.json').write_text(json.dumps(audit,indent=2))
print(json.dumps(audit,indent=2))
print(support.to_string())
print(pd.DataFrame(pair).to_string(index=False))

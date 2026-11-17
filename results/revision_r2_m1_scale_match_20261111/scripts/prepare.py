from pathlib import Path
import json,hashlib,copy
import numpy as np,pandas as pd,anndata as ad,yaml
from scipy import sparse
from scipy.spatial import cKDTree

ROOT=Path('/home/xd/y/Benchmark')
OUT=ROOT/'results/revision_r2_m1_scale_match_20261111'
for name in ['data','configs','scripts','logs','analysis']: (OUT/name).mkdir(parents=True,exist_ok=True)
cfg=yaml.safe_load((ROOT/'configs/datasets/e9p5_embryo.yaml').read_text())
source=ROOT/cfg['path'];a=ad.read_h5ad(source)
count_layer=next((k for k in ['counts','raw_count'] if k in a.layers),None)
X=sparse.csr_matrix(a.layers[count_layer] if count_layer else a.X)
assert np.all(np.isfinite(X.data)) and np.all(X.data>=0) and np.allclose(X.data,np.rint(X.data))
X=X.astype(np.int64)
xy=np.asarray(a.obsm[cfg['coord_key']],dtype=float)*0.65
ids=a.obs[cfg['slice_key']].astype(str).to_numpy();summary=[];centers_all=[];groups_all=[];mappings=[];offset=0
radius=27.5;pitch=100.;dy=np.sqrt(3)*pitch/2
for sid in sorted(set(ids)):
    ii=np.flatnonzero(ids==sid);pos=xy[ii];lo=pos.min(axis=0);hi=pos.max(axis=0)
    centers=[];names=[]
    for row in range(int(np.ceil((lo[1]+radius)/dy)),int(np.floor((hi[1]-radius)/dy))+1):
        shift=(row%2)*pitch/2
        for col in range(int(np.ceil((lo[0]+radius-shift)/pitch)),int(np.floor((hi[0]-radius-shift)/pitch))+1):
            centers.append([col*pitch+shift,row*dy]);names.append(f'{sid}:r{row}:c{col}')
    centers=np.asarray(centers);assert len(centers)>12,(sid,len(centers))
    dist,assignment=cKDTree(centers).query(pos)
    use=(dist<=radius+1e-9)&(np.asarray(X[ii].sum(axis=1)).ravel()>0)
    inds=ii[use];assigned=assignment[use];occupied=np.unique(assigned)
    assert len(occupied)>12,(sid,len(occupied))
    local=np.searchsorted(occupied,assigned)
    for g,k in enumerate(occupied):
        centers_all.append(centers[k]);groups_all.append({'sample_id':sid,'pseudo_spot':names[k],'n_beads':int((assigned==k).sum())})
    for ix,g,d in zip(inds,local,dist[use]):mappings.append({'source_index':int(ix),'bead_id':str(a.obs_names[ix]),'sample_id':sid,'pseudo_index':int(offset+g),'distance_to_center_um_assumed':float(d)})
    summary.append({'sample_id':sid,'input_beads':len(ii),'retained_beads':len(inds),'retained_fraction':len(inds)/len(ii),'pseudo_spots':len(occupied),'beads_per_pseudo_spot_median':float(np.median(np.bincount(local))),'library_retained_fraction':float(X[inds].sum()/X[ii].sum())})
    offset+=len(occupied)
mapping=pd.DataFrame(mappings).sort_values('source_index').reset_index(drop=True)
indices=mapping.source_index.to_numpy();native=X[indices]
# Fix the same measured gene universe before summation in both arms.
gene_mask=np.asarray((native>0).sum(axis=0)).ravel()>=3
native=native[:,gene_mask].tocsr();var=a.var.iloc[np.flatnonzero(gene_mask)].copy()
membership=sparse.csr_matrix((np.ones(len(indices),dtype=np.int64),(mapping.pseudo_index.to_numpy(),np.arange(len(indices)))),shape=(len(groups_all),len(indices)))
coarse=(membership@native).tocsr()
assert np.array_equal(np.asarray(native.sum(axis=0)),np.asarray(coarse.sum(axis=0)))
assert int(native.sum())==int(coarse.sum())
fineobs=a.obs.iloc[indices].copy();fineobs['sample_id']=fineobs[cfg['slice_key']].astype(str).astype('category')
fineobs['pseudo_spot']=[groups_all[k]['pseudo_spot'] for k in mapping.pseudo_index]
coarseobs=pd.DataFrame(groups_all).set_index('pseudo_spot',drop=False);coarseobs['sample_id']=coarseobs.sample_id.astype('category')
objects={'fine_footprint':ad.AnnData(native,obs=fineobs,var=var.copy()),'pseudo_visium':ad.AnnData(coarse,obs=coarseobs,var=var.copy())}
objects['fine_footprint'].obsm['spatial']=xy[indices]
objects['pseudo_visium'].obsm['spatial']=np.asarray(centers_all)
for arm,d in objects.items():
    d.layers['counts']=d.X.copy()
    d.uns['scale_sensitivity']={'coordinate_multiplier':0.65,'coordinate_basis':'original-author E9.5 convention; conditional, not independently calibrated per puck','diameter_um_assumed':55.,'pitch_um_assumed':100.,'arm':arm,'pairing':'exactly identical underlying beads and per-gene summed counts','gene_universe_rule':'detected in >=3 retained fine beads; same universe in both arms','lattice_origin':'global source-coordinate origin, alternating rows; not selected using expression or outcomes'}
    d.write_h5ad(OUT/'data'/f'{arm}.h5ad',compression='gzip')
    dc={'name':'E95Scale_'+arm,'path':str(OUT/'data'/f'{arm}.h5ad'),'platform':'Slide-seqV2_'+arm,'species':'Mouse','slice_key':'sample_id','coord_key':'spatial','preprocess':{'normalize_total':True,'target_sum':10000.,'log1p':True,'scale':False,'min_cells_per_gene':0,'count_layer':'counts','require_integer_counts':True}}
    (OUT/'configs'/f'data_{arm}.yaml').write_text(yaml.safe_dump(dc,sort_keys=False))
    for sid in d.obs.sample_id.cat.categories:
        coords=d.obsm['spatial'][d.obs.sample_id.astype(str).to_numpy()==sid]
        distances=cKDTree(coords).query(coords,k=9)[0][:,8]
        for row in summary:
            if row['sample_id']==sid:row[f'{arm}_median_knn8_um_assumed']=float(np.median(distances))
spatial=yaml.safe_load((ROOT/'configs/rebuild_v1/downstream/canonical/e9p5_embryo.yaml').read_text())
expr=yaml.safe_load((ROOT/'configs/rebuild_v1/non_spatial_downstream/canonical/e9p5_embryo.yaml').read_text())
excluded={'TFs','random','all','all_features','scsegindex','wilcoxon'}
methods=[x for x in expr['feature_selection_methods'] if x['name'] not in excluded]+spatial['feature_selection_methods']
assert len(methods)==29,(len(methods),[x['name'] for x in methods])
methods += [{'name':'random'},{'name':'all_features'}]
jobs=[]
for arm in objects:
    for method in methods:
        m=copy.deepcopy(method);m['n_features']=[2000]
        bc=copy.deepcopy(spatial);bc.update(name=f'scale_match_{arm}_{m["name"]}',datasets=[str(OUT/'configs'/f'data_{arm}.yaml')],feature_selection_methods=[m],output_dir=str(OUT/'runs'/arm/m['name']))
        path=OUT/'configs'/f'{arm}__{m["name"]}.yaml';path.write_text(yaml.safe_dump(bc,sort_keys=False))
        jobs.append({'arm':arm,'method':m['name'],'config':str(path),'output_dir':bc['output_dir'],'dataset_name':'E95Scale_'+arm})
mapping.to_csv(OUT/'analysis/bead_to_pseudo_spot.tsv',sep='\t',index=False)
pd.DataFrame(summary).to_csv(OUT/'analysis/aggregation_audit.tsv',sep='\t',index=False)
(OUT/'jobs.json').write_text(json.dumps(jobs,indent=2))
protocol={'source':str(source),'source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'source_shape':a.shape,'coordinate_multiplier':0.65,'coordinate_evidence':'Trunk_analysis_Fig2.Rmd original author commit 3d401f2e1cce973cbf89b219f279f5b3b36b03f0, line 710','coordinate_limit':'conditional author convention; not independently verified per puck; not applied to E8.5','scope':'paired scale-sensitivity on E9.5; five-dataset primary ranking unchanged','diameter_um_assumed':55,'pitch_um_assumed':100,'fine_shape':objects['fine_footprint'].shape,'coarse_shape':objects['pseudo_visium'].shape,'total_counts_per_arm':int(native.sum()),'counts_conserved_per_gene':True,'primary_contrast':'same retained beads before versus after summation on non-overlapping pseudo-Visium spots','genes_filtered_once_before_aggregation':True,'selectors':29,'references':['random_2000','all_features'],'workflows':['scvi','cellcharter'],'seeds':[0,1,2],'n_clusters':12,'requested_budget':2000,'expected_jobs':len(jobs),'expected_embedding_runs':len(jobs)*6,'interpretation':'label-free metric/rank sensitivity, not biological conservation or a technology causal effect','metric_policy':'bASW/iLISI and silhouette/CHAOS/PAS separately; no embryo CoreOverall; use pre-existing frozen E9.5 ranges only for descriptive components; raw metrics and within-arm ranks primary','planned_checks':['integer counts','per-gene exact count conservation','same feature universe and bead origin','all six pucks retained','minimum 13 spots per puck','actual output counts and matched eligible panel','complete metrics per workflow/seed','selection overlap and scale-specific ranks']}
(OUT/'protocol.json').write_text(json.dumps(protocol,indent=2))
print(json.dumps(protocol,indent=2));print(pd.DataFrame(summary).to_string(index=False))

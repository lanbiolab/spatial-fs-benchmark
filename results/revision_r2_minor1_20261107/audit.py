from pathlib import Path
import json,glob,tarfile,io,hashlib
import anndata as ad,numpy as np,pandas as pd
out=Path('results/revision_r2_minor1_20261107');out.mkdir(exist_ok=True)
records=json.loads((out/'input_inspection.json').read_text());rows=[]
for rec in records:
 ds=rec['dataset']; a=ad.read_h5ad(rec['config']['path'],backed='r');coords=np.asarray(a.obsm['spatial']);key=rec['config']['slice_key']
 for sid in a.obs[key].astype(str).unique():
  mask=a.obs[key].astype(str).eq(sid).to_numpy(); obs=a.obs.loc[mask];xy=coords[mask]; src=''; sf=''; unit=''; exact=False
  if ds=='DLPFC':
   src='data/raw/dlpfc/spatial.csv';p=pd.read_csv(src,index_col=0);o=pd.read_csv('data/raw/dlpfc/obs.csv',index_col=0);assert len(p)==len(o);assert [str(x).replace('.','-',1) for x in p.index]==list(o.index);p.index=o['key'];expected=p.loc[obs.index,['pxl_col_in_fullres','pxl_row_in_fullres']].to_numpy();unit='full-resolution image pixels';exact=np.array_equal(xy,expected)
  elif ds=='MouseBrainSerialSections':
   tarpaths=glob.glob('data/raw/Mouse Brain */*/*spatial.tar.gz');src=next(p for p in tarpaths if ('Section_2' in p)==(sid=='section2'))
   with tarfile.open(src) as t:
    m=next(m for m in t.getmembers() if 'tissue_positions' in m.name);data=t.extractfile(m).read();has_header=data.startswith(b'barcode,');p=pd.read_csv(io.BytesIO(data),header=0 if has_header else None)
    p.columns=['barcode','in_tissue','array_row','array_col','pxl_row_in_fullres','pxl_col_in_fullres'];p=p.set_index('barcode');bar=['-'.join(str(x).split(':')[-1].split('-')[:2]) for x in obs.index];expected=p.loc[bar,['pxl_col_in_fullres','pxl_row_in_fullres']].to_numpy();exact=np.array_equal(xy,expected);unit='full-resolution image pixels'
  elif ds.startswith('STOmics'):
   src=next(str(p) for p in Path(rec['config']['path']).parent.glob(sid+'*processed.h5ad'));b=ad.read_h5ad(src,backed='r');bar=[str(x).split(':')[-1] for x in obs.index];idx=b.obs_names.get_indexer(bar);assert (idx>=0).all();exact=np.array_equal(xy,np.asarray(b.obsm['spatial'])[idx]);sp=b.uns['spatial'];sf=json.dumps({str(k):{str(s):float(v) for s,v in val['scalefactors'].items()} for k,val in sp.items()});unit='source Visium image pixels (spatial/scalefactors convention)';b.file.close()
  else:
   src=rec['config']['path'];exact=np.array_equal(xy,obs[['x','y']].to_numpy());unit='native bead-localisation units; physical calibration not supplied'
  assert exact,(ds,sid)
  rows.append(dict(dataset=ds,section=sid,n_observations=len(obs),coordinate_unit=unit,source=src,source_coordinate_match=bool(exact),coordinate_sha256=hashlib.sha256(np.ascontiguousarray(xy).tobytes()).hexdigest(),scalefactors=sf))
 a.file.close()
pd.DataFrame(rows).to_csv(out/'coordinate_provenance_audit.tsv',sep='\t',index=False)
original=json.loads(Path('results/revision_r2_major2_20261106/manuscript_original_hashes.json').read_text());changed=[p for p,h in original.items() if hashlib.sha256(Path(p).read_bytes()).hexdigest()!=h];assert not changed
(out/'audit.json').write_text(json.dumps(dict(datasets=len(records),sections=len(rows),all_coordinates_match=True,original_manuscript_and_figures_checked=len(original),changed=changed),indent=2))
print((out/'audit.json').read_text())



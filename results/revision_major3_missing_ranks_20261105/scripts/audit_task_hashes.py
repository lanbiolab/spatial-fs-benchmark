import json,pathlib
root=pathlib.Path('/home/xd/y/Benchmark/results/spatial_svg_rebuild_v1');ds={'dlpfc':'dlpfc','mouse_brain_serial_sections':'mousebrainserialsections','stomics_0212':'stomics0212','stomics_0218':'stomics0218','stomics_0224':'stomics0224'};methods=['morans_i','somde','spatialde','sparkx','nnsvg'];integs=['scvi','cellcharter','graphst'];good=0;bad=[];missing=[]
for d,name in ds.items():
 for m in methods:
  for seed in range(3):
   for i in integs:
    path=root/d/name/i/m/'n2000'/('seed'+str(seed)); meta=json.load(open(path/'embedding.meta.json'))
    for task in (['integration_eval','clustering_eval','alignment_eval'] if d in ['dlpfc','mouse_brain_serial_sections'] else ['integration_eval','clustering_eval']):
     p=path/(task+'_records.meta.json')
     if not p.exists():missing.append(str(p));continue
     obj=json.load(open(p))
     if obj.get('selection_feature_hash')!=meta.get('selection_feature_hash') or obj.get('integration_cache_signature')!=meta.get('cache_signature'):bad.append((str(p),obj.get('integration_cache_signature'),meta.get('cache_signature')))
     else:good+=1
print('original task provenance good:',good,'missing:',len(missing),'mismatch:',len(bad));print('examples:',bad[:4]);print('missing examples',missing[:4])

import json, pathlib
root=pathlib.Path('/home/xd/y/Benchmark/results/spatial_svg_rebuild_v1');ds={'dlpfc':'dlpfc','mouse_brain_serial_sections':'mousebrainserialsections','stomics_0212':'stomics0212','stomics_0218':'stomics0218','stomics_0224':'stomics0224'};methods=['morans_i','somde','spatialde','sparkx','nnsvg'];integrators=['scvi','cellcharter','graphst'];mismatch=[];missing=[];ok=0
for d,name in ds.items():
 for m in methods:
  for seed in range(3):
   sel=root/d/name/'feature_selection'/m/'n2000'/('seed'+str(seed if m=='spatialde' else 0))/'selected_features.json'
   if not sel.exists():missing.append(str(sel));continue
   h=json.load(open(sel))['metadata']['feature_name_hash']
   for i in integrators:
    meta=root/d/name/i/m/'n2000'/('seed'+str(seed))/'embedding.meta.json'
    if not meta.exists():missing.append(str(meta));continue
    j=json.load(open(meta));h2=j.get('selection_feature_hash')
    if h!=h2:mismatch.append((d,m,seed,i,h,h2))
    else:ok+=1
print('original embedding provenance: matches',ok,'missing',len(missing),'mismatch',len(mismatch));print('missing examples',missing[:3]);print('mismatch examples',mismatch[:5])

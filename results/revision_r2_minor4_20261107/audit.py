from pathlib import Path
import json,re,hashlib,csv,shutil,subprocess,platform
ROOT=Path('/home/xd/y/Benchmark');OUT=ROOT/'results/revision_r2_minor4_20261107';src=ROOT/'results/resource_profile_v1';e=OUT/'evidence';e.mkdir(exist_ok=True)
rows=[]
for p in sorted((src/'json').glob('*.json')):
 d=json.loads(p.read_text());method=d['method'];t=(src/'time'/f'{method}.txt').read_text();lines={k.strip():v.strip() for line in t.splitlines() if ':' in line for k,v in [line.rsplit(': ',1)] if ': ' in line}
 # GNU time labels contain colons within the unit description, so parse by final colon-space.
 user=float(lines['User time (seconds)']);system=float(lines['System time (seconds)']);pct=int(lines['Percent of CPU this job got'].rstrip('%'));clock=lines['Elapsed (wall clock) time (h:mm:ss or m:ss)'];parts=clock.split(':');wall=sum(float(x)*60**i for i,x in enumerate(reversed(parts)));rss=int(lines['Maximum resident set size (kbytes)'])/1024
 if method in ['sparkx','nnsvg']:control='4 requested method threads; R child OMP_NUM_THREADS=4 and OPENBLAS_NUM_THREADS=4'
 elif method=='seurat_sct':control='future::sequential in R; BLAS/OpenMP effective thread count not recorded'
 elif method=='scPNMF':control='mc.cores=1 for feature-selection call; native library thread count not recorded'
 else:control='no explicit uniform per-method thread cap in profiling launcher; effective library threads not recorded'
 rows.append(dict(method=method,n_spots=d['n_spots'],n_genes=d['n_genes'],requested_n_features=d['requested_n_features'],returned_n_features=d['effective_n_features'],selector_seconds=d['wall_seconds_internal'],process_wall_seconds=wall,peak_rss_mib=rss,cpu_user_seconds=user,cpu_system_seconds=system,whole_process_cpu_percent=pct,configured_parallelism=control,selector_params=json.dumps(d['selector_params'],sort_keys=True),selector_implementation_version=d['selector_implementation_version'],feature_hash=d['feature_hash']))
 for q in [p,src/'time'/f'{method}.txt']:
  target=e/q.relative_to(src);target.parent.mkdir(exist_ok=True);shutil.copy2(q,target)
assert len(rows)==34 and all(r['n_spots']==2000 and r['n_genes']==5000 for r in rows)
with (OUT/'resource_profile_with_execution_metadata.tsv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
original=list(csv.DictReader((src/'feature_selection_resource_profile.tsv').open(),delimiter='\t'));by={r['method']:r for r in original}
for r in rows:
 b=by[r['method']];assert abs(float(b['wall_seconds_internal'])-r['selector_seconds'])<1e-9;assert abs(float(b['wall_seconds_process'])-r['process_wall_seconds'])<.011;assert abs(float(b['peak_rss_mib'])-r['peak_rss_mib'])<1e-8
scriptpaths=['scripts/run_feature_selection_resource_profile.sh','scripts/profile_feature_selector.py','scripts/run_spatial_r_feature_selector.R','scripts/run_r_feature_selector.R','src/spatial_fs_benchmark/feature_selection/spatial_r_method.py','src/spatial_fs_benchmark/feature_selection/r_method.py','src/spatial_fs_benchmark/feature_selection/svg.py','configs/rebuild_v1/non_spatial_downstream/canonical/mouse_brain_serial_sections.yaml','configs/rebuild_v1/downstream/canonical/mouse_brain_serial_sections.yaml']
for rel in scriptpaths:
 target=e/'code'/rel;target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(ROOT/rel,target)
p=Path('/home/xd/.codex/sessions/2026/03/14/rollout-2026-03-14T23-07-50-019cece3-af0d-7fa1-a671-7afa3736a3d6.jsonl');historical=[]
def extract_text(v):
 if isinstance(v,str):return v
 if isinstance(v,list):return '\n'.join(extract_text(x) for x in v)
 if isinstance(v,dict):return str(v.get('text',''))
 return ''
for i,line in enumerate(p.open(),1):
 if i!=37540:continue
 d=json.loads(line);s=extract_text(d['payload']['output'])
 for line in s.splitlines():
  if re.search(r'^CPU 112|^CPU\(s\):|^Model name:|^Thread\(s\)|^Core\(s\)|^Socket\(s\)|^Mem:',line):historical.append(line)
assert any('6330' in s for s in historical),historical
(e/'historical_hardware_snapshot.txt').write_text('Archived server snapshot: 2026-07-22T06:04:07Z, before 2026-07-29 profiling launch.\n'+'\n'.join(historical))
current={'checked_date':'2026-11-07','cpu_info':subprocess.check_output(['lscpu'],text=True,env={'LC_ALL':'C','PATH':'/usr/bin:/bin'}),'memory_total_bytes':int(re.search(r'MemTotal:\s+(\d+)',Path('/proc/meminfo').read_text()).group(1))*1024,'os_release':Path('/etc/os-release').read_text(),'python_current':platform.python_version()};(e/'current_host_snapshot.json').write_text(json.dumps(current,indent=2))
original_hashes=json.loads((ROOT/'results/revision_r2_major2_20261106/manuscript_original_hashes.json').read_text());changed=[p for p,h in original_hashes.items() if hashlib.sha256(Path(p).read_bytes()).hexdigest()!=h];assert not changed
(OUT/'audit.json').write_text(json.dumps(dict(profiled_methods=34,all_profile_values_reproduced_from_original_logs=True,uniform_input_spots=2000,uniform_input_genes=5000,launcher_concurrent_process_limit=2,matched_total_thread_budget=False,global_effective_thread_counts_recorded=False,new_timing_experiments=0,original_manuscript_and_figures_checked=len(original_hashes),changed=changed),indent=2))
(OUT/'source_hashes.json').write_text(json.dumps({str(p.relative_to(e)):hashlib.sha256(p.read_bytes()).hexdigest() for p in e.rglob('*') if p.is_file()},indent=2))
import zipfile
with zipfile.ZipFile(OUT/'audit_outputs.zip','w',zipfile.ZIP_DEFLATED) as z:
 for p in [OUT/'resource_profile_with_execution_metadata.tsv',OUT/'audit.json',OUT/'source_hashes.json',*e.rglob('*')]:
  if p.is_file():z.write(p,p.relative_to(OUT))
print((OUT/'audit.json').read_text());print((e/'historical_hardware_snapshot.txt').read_text())

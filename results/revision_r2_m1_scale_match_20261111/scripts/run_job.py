from pathlib import Path
import sys,json,time,math,fcntl,yaml
raw=yaml.safe_load(Path(sys.argv[1]).read_text());jobout=Path(raw['output_dir']);jobout.mkdir(parents=True,exist_ok=True)
lock_handle=(jobout/'execution.lock').open('w');fcntl.flock(lock_handle,fcntl.LOCK_EX)
if '--pilot' not in sys.argv and (jobout/'job_completion.json').exists():
    prior=json.loads((jobout/'job_completion.json').read_text())
    if prior.get('state')=='completed' and all(Path(p).exists() for p in prior.get('records',[])+prior.get('embeddings',[])) and len(prior.get('records',[]))==12:
        print('Verified existing completed job',raw['name'],flush=True);raise SystemExit(0)
from spatial_fs_benchmark.config import load_benchmark_config
from spatial_fs_benchmark.benchmark.runner import BenchmarkRunner
cfg=load_benchmark_config(sys.argv[1]);pilot='--pilot' in sys.argv
if pilot:cfg.seeds=[0]
out=Path(cfg.output_dir);out.mkdir(parents=True,exist_ok=True);started=time.monotonic()
BenchmarkRunner(cfg).run()
missing=[];records=[];embeddings=[]
for seed in cfg.seeds:
    for workflow in cfg.integration_methods:
        folder=out/cfg.datasets[0].name.lower()/workflow.name/cfg.feature_selection_methods[0].name/'n2000'/f'seed{seed}'
        for task in cfg.tasks:
            p=folder/f'{task.name}_records.json'
            if not p.exists():missing.append(str(p));continue
            data=json.loads(p.read_text())
            if not data:missing.append(str(p)+' EMPTY')
            else:records.append(str(p))
        if not (folder/'embedding.npz').exists():missing.append(str(folder/'embedding.npz'))
        else:embeddings.append(str(folder/'embedding.npz'))
result={'state':'completed' if not missing else 'incomplete','pilot':pilot,'elapsed_seconds':time.monotonic()-started,'records':records,'embeddings':embeddings,'missing':missing}
(out/('pilot_completion.json' if pilot else 'job_completion.json')).write_text(json.dumps(result,indent=2))
print(json.dumps(result),flush=True)
sys.exit(0 if not missing else 2)

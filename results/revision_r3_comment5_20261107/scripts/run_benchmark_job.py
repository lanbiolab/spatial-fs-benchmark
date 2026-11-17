from pathlib import Path
import sys,json,time
from spatial_fs_benchmark.config import load_benchmark_config
from spatial_fs_benchmark.benchmark.runner import BenchmarkRunner
cfg=load_benchmark_config(sys.argv[1]);out=Path(cfg.output_dir);out.mkdir(parents=True,exist_ok=True);t=time.monotonic();BenchmarkRunner(cfg).run()
records=[];missing=[]
for s in cfg.seeds:
 for method in cfg.integration_methods:
  p=out/cfg.datasets[0].name.lower()/method.name/cfg.feature_selection_methods[0].name/'n2000'/f'seed{s}'/'clustering_eval_records.json'
  if p.exists() and json.loads(p.read_text()):records.append(str(p))
  else:missing.append(str(p))
(out/'job_completion.json').write_text(json.dumps({'state':'completed' if not missing else 'incomplete','elapsed_seconds':time.monotonic()-t,'valid_records':records,'missing_records':missing},indent=2));print('JOB',cfg.feature_selection_methods[0].name,len(records),'valid,',len(missing),'missing',flush=True)
sys.exit(0 if not missing else 2)

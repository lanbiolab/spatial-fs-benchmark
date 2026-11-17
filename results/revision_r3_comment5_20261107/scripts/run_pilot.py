from pathlib import Path
import json,time
from spatial_fs_benchmark.config import load_benchmark_config
from spatial_fs_benchmark.benchmark.runner import BenchmarkRunner
O=Path('/home/xd/y/Benchmark/results/revision_r3_comment5_20261107');t=time.monotonic()
BenchmarkRunner(load_benchmark_config(O/'prepared/stereo_pilot.yaml')).run()
p=O/'pilot/stereoseqadultbrain/scvi/morans_i/n2000/seed0/clustering_eval_records.json'
assert p.exists() and json.loads(p.read_text()),'Pilot task did not produce records; inspect benchmark log'
(O/'pilot_verified.json').write_text(json.dumps({'state':'pilot_completed','elapsed_seconds':time.monotonic()-t,'records':str(p),'not_full_platform_experiment':True},indent=2))
print('PILOT verified',time.monotonic()-t,flush=True)

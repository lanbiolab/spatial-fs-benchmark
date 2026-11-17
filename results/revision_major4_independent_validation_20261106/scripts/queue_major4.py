import concurrent.futures
import json
import os
from pathlib import Path
import subprocess

ROOT = Path('/home/xd/y/Benchmark')
OUT = ROOT / 'results/revision_major4_independent_validation_20261106'
env = os.environ.copy()
env.update(PYTHONPATH='src', OPENBLAS_NUM_THREADS='4', OMP_NUM_THREADS='4', MKL_NUM_THREADS='4', NUMBA_NUM_THREADS='4')
folds = json.loads((OUT / 'analysis_protocol.json').read_text())['folds']

def run(fold):
    with (OUT / f'{fold}.log').open('a') as f:
        code = subprocess.run([str(ROOT / '.conda-env/bin/python'), '-u', str(OUT / 'scripts/run_major4.py'), '--fold', fold],
                              cwd=ROOT, env=env, stdout=f, stderr=subprocess.STDOUT).returncode
    print(fold, 'exit', code, flush=True)
    return code

with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
    results = list(pool.map(run, folds))
if any(results):
    raise SystemExit('A fold failed; inspect individual logs')
with (OUT / 'summary.log').open('a') as f:
    code = subprocess.run([str(ROOT / '.conda-env/bin/python'), '-u', str(OUT / 'scripts/summarize_major4.py')],
                          cwd=ROOT, env=env, stdout=f, stderr=subprocess.STDOUT).returncode
print('Summary exit', code, flush=True)
raise SystemExit(code)

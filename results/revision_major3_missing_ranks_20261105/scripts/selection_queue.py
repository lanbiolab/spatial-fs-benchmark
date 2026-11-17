import concurrent.futures
import os
import subprocess
from pathlib import Path

root = Path('/home/xd/y/Benchmark')
out = root/'results/revision_major3_missing_ranks_20261105'
datasets = ['mouse_brain_serial_sections','dlpfc','stomics_0212','stomics_0218','stomics_0224','e8p5_embryo','e9p5_embryo']
methods = ['morans_i','somde','spatialde','sparkx','nnsvg']
jobs = [(d,m,s) for d in datasets for m in methods for s in (range(3) if m=='spatialde' else [0])]
env = os.environ.copy()
env.update(PYTHONPATH='src', OPENBLAS_NUM_THREADS='2', OMP_NUM_THREADS='2', MKL_NUM_THREADS='2', NUMBA_NUM_THREADS='2')

def run(job):
    d,m,s=job
    folder=out/'selection'/d/m/f'seed{s}'
    folder.mkdir(parents=True,exist_ok=True)
    if (folder/'DONE').exists():
        return job,0,'already complete'
    command=[str(root/'.conda-env/bin/python'),'-u',str(out/'scripts/capture_selection.py'),d,m,str(s)]
    with (folder/'run.log').open('a') as log:
        result=subprocess.run(command,cwd=root,env=env,stdout=log,stderr=subprocess.STDOUT)
    return job,result.returncode,'completed' if result.returncode==0 else 'FAILED'

print('Selection queue:',len(jobs),'jobs; max 4 concurrent',flush=True)
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
    fs={pool.submit(run,j):j for j in jobs}
    for f in concurrent.futures.as_completed(fs):
        print(f.result(),flush=True)
print('Queue finished',flush=True)

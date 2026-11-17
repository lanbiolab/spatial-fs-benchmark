import concurrent.futures
import os
import subprocess
from pathlib import Path
root=Path('/home/xd/y/Benchmark');out=root/'results/revision_major3_missing_ranks_20261105'
datasets=['stomics_0224','stomics_0218','stomics_0212']
jobs=[(d,m,s) for d in datasets for m in ['somde','spatialde','sparkx','nnsvg'] for s in (range(3) if m=='spatialde' else [0])]
env=os.environ.copy();env.update(PYTHONPATH='src',OPENBLAS_NUM_THREADS='2',OMP_NUM_THREADS='2',MKL_NUM_THREADS='2',NUMBA_NUM_THREADS='2')
def run(job):
    d,m,s=job;folder=out/'selection'/d/m/f'seed{s}';folder.mkdir(parents=True,exist_ok=True)
    if (folder/'DONE').exists():return job,0,'already complete'
    command=[str(root/'.conda-env/bin/python'),'-u',str(out/'scripts/capture_selection.py'),d,m,str(s)]
    with (folder/'run.log').open('a') as log:
        ret=subprocess.run(command,cwd=root,env=env,stdout=log,stderr=subprocess.STDOUT).returncode
    return job,ret
print('Additional independent selection jobs:',len(jobs),flush=True)
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
    for x in concurrent.futures.as_completed([pool.submit(run,j) for j in jobs]):print(x.result(),flush=True)
print('Additional selection queue finished',flush=True)

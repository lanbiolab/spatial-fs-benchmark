import concurrent.futures
import os
import subprocess
from pathlib import Path

root=Path('/home/xd/y/Benchmark');out=root/'results/revision_major3_missing_ranks_20261105'
datasets=['dlpfc','stomics_0212','stomics_0218','stomics_0224','e8p5_embryo','e9p5_embryo']
jobs=[(d,'morans_i',0) for d in datasets]+[(d,'nnsvg',0) for d in ['e8p5_embryo','e9p5_embryo']]
env=os.environ.copy();env.update(PYTHONPATH='src',OPENBLAS_NUM_THREADS='2',OMP_NUM_THREADS='2',MKL_NUM_THREADS='2',NUMBA_NUM_THREADS='2')
def run(job):
    d,m,s=job;folder=out/'selection'/d/m/f'seed{s}';folder.mkdir(parents=True,exist_ok=True)
    if (folder/'DONE').exists():return job,0,'already complete'
    command=[str(root/'.conda-env/bin/python'),'-u',str(out/'scripts/capture_selection.py'),d,m,str(s)]
    with (folder/'run.log').open('a') as log:
        ret=subprocess.run(command,cwd=root,env=env,stdout=log,stderr=subprocess.STDOUT).returncode
    return job,ret
print('Independent priority jobs:',len(jobs),flush=True)
with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:
    for x in concurrent.futures.as_completed([pool.submit(run,j) for j in jobs]):print(x.result(),flush=True)
print('Priority queue finished',flush=True)

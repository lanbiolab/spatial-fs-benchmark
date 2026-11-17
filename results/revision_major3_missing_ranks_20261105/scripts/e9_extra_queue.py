import concurrent.futures,os,subprocess
from pathlib import Path
root=Path('/home/xd/y/Benchmark');out=root/'results/revision_major3_missing_ranks_20261105';ds='e9p5_embryo'
jobs=[(m,s) for m in ['somde','spatialde','sparkx'] for s in (range(3) if m=='spatialde' else [0])]
env=os.environ.copy();env.update(PYTHONPATH='src',OPENBLAS_NUM_THREADS='2',OMP_NUM_THREADS='2',MKL_NUM_THREADS='2',NUMBA_NUM_THREADS='2')
os.nice(10)
def run(job):
 m,s=job;folder=out/'selection'/ds/m/f'seed{s}';folder.mkdir(parents=True,exist_ok=True)
 if (folder/'DONE').exists():return job,0,'already complete'
 cmd=[str(root/'.conda-env/bin/python'),'-u',str(out/'scripts/capture_selection.py'),ds,m,str(s)]
 with (folder/'run.log').open('a') as log:ret=subprocess.run(cmd,cwd=root,env=env,stdout=log,stderr=subprocess.STDOUT).returncode
 return job,ret
print('E9.5 remaining gene-universe jobs:',len(jobs),flush=True)
with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:
 for f in concurrent.futures.as_completed([pool.submit(run,j) for j in jobs]):print(f.result(),flush=True)
print('E9.5 queue finished',flush=True)

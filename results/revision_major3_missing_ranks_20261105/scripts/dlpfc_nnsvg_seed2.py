import os,subprocess
from pathlib import Path
root=Path('/home/xd/y/Benchmark');out=root/'results/revision_major3_missing_ranks_20261105';ds='dlpfc';method='nnsvg'
env=os.environ.copy();env.update(PYTHONPATH='src',CUDA_VISIBLE_DEVICES='3',OPENBLAS_NUM_THREADS='2',OMP_NUM_THREADS='2',MKL_NUM_THREADS='2',NUMBA_NUM_THREADS='2')
print('Additional DLPFC nnSVG paired downstream seeds 1 and 2',flush=True)
for seed in [2]:
 for integ in ['scvi','cellcharter','graphst']:
  for version in ['missing_excluded','original_zero']:
   job=out/'downstream'/ds/integ/method/f'seed{seed}'/version
   if (job/'DONE').exists():continue
   job.mkdir(parents=True,exist_ok=True);cmd=[str(root/'.conda-env/bin/python'),'-u',str(out/'scripts/run_downstream.py'),ds,method,str(seed),integ,version]
   print('START',seed,integ,version,flush=True)
   with (job/'run.log').open('a') as log:ret=subprocess.run(cmd,cwd=root,env=env,stdout=log,stderr=subprocess.STDOUT).returncode
   print('FINISH',seed,integ,version,'exit',ret,flush=True)
   if ret:(job/'FAILED').write_text(f'exit {ret}\n')
print('DLPFC nnSVG extra worker finished',flush=True)

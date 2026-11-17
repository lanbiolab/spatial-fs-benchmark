import concurrent.futures
import os
import shutil
import subprocess
import tempfile
from pathlib import Path

root=Path('/home/xd/y/Benchmark');out=root/'results/revision_major3_missing_ranks_20261105';runner=out/'scripts/run_spatial_r_one_slice.R'
jobs=[('mouse_brain_serial_sections',2),('stomics_0212',10),('stomics_0218',7),('stomics_0224',8)]
env=os.environ.copy();env.update(OMP_NUM_THREADS='4',OPENBLAS_NUM_THREADS='4')
def run(job):
 ds,n=job;shared=out/'selection'/ds/'nnsvg'/'seed0';target=shared/f'slice_{n:04d}.rds'
 if target.exists() or (shared/'DONE').exists():return job,0,'already finished'
 helper=out/'parallel_slices'/ds/f'slice{n:04d}';helper.mkdir(parents=True,exist_ok=True)
 required=['matrix.mtx','genes.tsv','observations.tsv'];assert all((shared/x).exists() for x in required),f'Missing shared input {shared}'
 e=env.copy();e['M3_SLICE_NUMBER']=str(n)
 cmd=[str(root/'.conda-env/bin/Rscript'),str(runner),str(shared),str(helper/'ranks.tsv'),'nnsvg','0','0','4','0']
 with (helper/'run.log').open('a') as f:ret=subprocess.run(cmd,cwd=root,env=e,stdout=f,stderr=subprocess.STDOUT).returncode
 if ret:return job,ret,'failed'
 saved=helper/f'slice_{n:04d}.rds';assert saved.exists(),saved
 check='x<-readRDS(commandArgs(TRUE)[1]);stopifnot(is.numeric(x),length(x)==length(readLines(commandArgs(TRUE)[2])),any(is.finite(x)),all(x[is.finite(x)]>0),all(x[is.finite(x)]<=1))'
 subprocess.run([str(root/'.conda-env/bin/Rscript'),'-e',check,str(saved),str(shared/'genes.tsv')],cwd=root,check=True,stdout=subprocess.DEVNULL)
 temp=shared/f'.parallel_slice_{n:04d}_{os.getpid()}.tmp';shutil.copyfile(saved,temp)
 try:
  os.link(temp,target);message='published atomically'
 except FileExistsError:message='original process already completed slice'
 finally:temp.unlink()
 return job,0,message
print('Parallel exact slice jobs:',jobs,flush=True)
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
 for future in concurrent.futures.as_completed([pool.submit(run,j) for j in jobs]):print(future.result(),flush=True)
print('Parallel slice queue finished',flush=True)

import os,subprocess
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmwrite

root=Path('/home/xd/y/Benchmark');out=root/'results/revision_major3_missing_ranks_20261105';runner=out/'scripts/run_spatial_r_capture.R';helper=out/'scripts/run_spatial_r_one_slice.R'
text=runner.read_text();old='    slice_id <- slice_ids[[slice_number]]\n    slice_cache <-';new='    slice_id <- slice_ids[[slice_number]]\n    if (slice_number != as.integer(Sys.getenv("M3_SLICE_NUMBER"))) next\n    slice_cache <-';assert text.count(old)==1
helper.write_text(text.replace(old,new))
test=out/'parallel_validation';data=test/'input';base=test/'serial';extra=test/'single_slice';data.mkdir(parents=True,exist_ok=True);base.mkdir(exist_ok=True);extra.mkdir(exist_ok=True)
rng=np.random.default_rng(20261105);genes=np.array([f'gene{i:02}' for i in range(24)]);counts=rng.poisson(3,size=(24,160));mmwrite(data/'matrix.mtx',sparse.coo_matrix(counts));pd.Series(genes).to_csv(data/'genes.tsv',sep='\t',index=False,header=False)
pd.DataFrame({'slice':['A']*80+['B']*80,'x':rng.uniform(0,10,160),'y':rng.uniform(0,10,160)}).to_csv(data/'observations.tsv',sep='\t',index=False)
env=os.environ.copy();env.update(OMP_NUM_THREADS='4',OPENBLAS_NUM_THREADS='4')
def run(script,folder,number=None):
 e=env.copy()
 if number is not None:e['M3_SLICE_NUMBER']=str(number)
 cmd=[str(root/'.conda-env/bin/Rscript'),str(script),str(data),str(folder/'ranks.tsv'),'nnsvg','0','0','4','0']
 with (folder/'run.log').open('w') as f:return subprocess.run(cmd,cwd=root,env=e,stdout=f,stderr=subprocess.STDOUT).returncode
with ThreadPoolExecutor(max_workers=2) as ex:
 a=ex.submit(run,runner,base);b=ex.submit(run,helper,extra,2);ra,rb=a.result(),b.result()
print('Serial exit',ra,'independent slice exit',rb,flush=True)
assert ra==rb==0
code='a<-readRDS(commandArgs(TRUE)[1]);b<-readRDS(commandArgs(TRUE)[2]);cat("identical=",identical(a,b),";max_difference=",max(abs(a-b),na.rm=TRUE),";same_missing=",identical(is.na(a),is.na(b)),"\\n");stopifnot(identical(a,b))'
subprocess.run([str(root/'.conda-env/bin/Rscript'),'-e',code,str(base/'slice_0002.rds'),str(extra/'slice_0002.rds')],check=True,cwd=root)
print('VALIDATED exact per-gene slice percentiles and missing masks',flush=True)

from pathlib import Path
import subprocess,os,json,time,queue,threading,concurrent.futures,math
R=Path('/home/xd/y/Benchmark');O=R/'results/revision_r3_comment5_20261107';L=O/'logs_openst';L.mkdir(exist_ok=True)
PY=R/'.conda-env/bin/python';runner=O/'scripts/run_benchmark_job.py'
def save(p,v):
 t=p.with_suffix(p.suffix+'.tmp');t.write_text(json.dumps(v,indent=2));t.replace(p)
def env_for(gpu):
 e=os.environ.copy();e.update(CUDA_VISIBLE_DEVICES=str(gpu),OMP_NUM_THREADS='4',OPENBLAS_NUM_THREADS='4',MKL_NUM_THREADS='4',NUMBA_NUM_THREADS='4');return e
with (L/'preparation.log').open('w') as f:subprocess.run([str(PY),'-u',str(O/'scripts/prepare_openst.py')],cwd=R,stdout=f,stderr=subprocess.STDOUT,check=True)
print('Preparation complete; starting pilot',flush=True)
with (L/'pilot.log').open('w') as f:subprocess.run([str(PY),'-u',str(runner),str(O/'prepared/openst_pilot.yaml')],cwd=R,env=env_for(1),stdout=f,stderr=subprocess.STDOUT,check=True)
records=list((O/'openst_pilot').rglob('clustering_eval_records.json'));assert len(records)==1
rows=json.loads(records[0].read_text());assert len(rows)==3 and all(math.isfinite(x['metric_value']) for x in rows)
assert {x['metric_name'] for x in rows}=={'silhouette','CHAOS','PAS'}
save(O/'openst_pilot_verified.json',{'verified_at':time.strftime('%Y-%m-%dT%H:%M:%S%z'),'records':str(records[0])})
print('Pilot verified; starting full queue',flush=True)
jobs=json.loads((O/'openst_jobs.json').read_text());q=queue.Queue();[q.put(j) for j in jobs];lock=threading.Lock();summary=[]
def worker(gpu):
 while True:
  try:j=q.get_nowait()
  except queue.Empty:return
  e={'method':j['method'],'gpu':gpu,'state':'running','started_at':time.strftime('%Y-%m-%dT%H:%M:%S%z')};sf=L/(j['method']+'.status.json');save(sf,e);print('START',j['method'],gpu,flush=True)
  try:
   with (L/(j['method']+'.log')).open('w') as f:
    p=subprocess.Popen([str(PY),'-u',str(runner),j['config']],cwd=R,env=env_for(gpu),stdout=f,stderr=subprocess.STDOUT,start_new_session=True)
    e['pid']=p.pid;save(sf,e);rc=p.wait()
   e.update(state='completed' if rc==0 else 'failed_or_incomplete',exit_code=rc)
  except Exception as err:e.update(state='failed_or_incomplete',error=repr(err))
  e['ended_at']=time.strftime('%Y-%m-%dT%H:%M:%S%z');save(sf,e)
  with lock:summary.append(e);save(O/'openst_queue_status.json',summary)
  print('END',j['method'],e['state'],flush=True);q.task_done()
with concurrent.futures.ThreadPoolExecutor(3) as ex:list(ex.map(worker,[1,2,3]))
save(O/'openst_queue_finished.json',{'jobs':summary,'all_jobs_completed':len(summary)==33 and all(x['state']=='completed' for x in summary)})

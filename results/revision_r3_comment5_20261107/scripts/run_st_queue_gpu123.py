from pathlib import Path
import os,json,time,subprocess,queue,threading,concurrent.futures,signal
ROOT=Path('/home/xd/y/Benchmark');OUT=ROOT/'results/revision_r3_comment5_20261107';L=OUT/'logs_st';L.mkdir(exist_ok=True)
jobs=json.loads((OUT/'st_jobs.json').read_text());q=queue.Queue();[q.put(j) for j in jobs];lock=threading.Lock();summary=[]
start=time.monotonic()
while not (OUT/'pilot_verified.json').exists():
 if time.monotonic()-start>7200:raise RuntimeError('Pilot verification absent after 2h; full benchmark not started')
 time.sleep(15)
def worker(gpu):
 while True:
  try:j=q.get_nowait()
  except queue.Empty:return
  env=os.environ.copy();env.update(CUDA_VISIBLE_DEVICES=str(gpu),OMP_NUM_THREADS='4',OPENBLAS_NUM_THREADS='4',MKL_NUM_THREADS='4',NUMBA_NUM_THREADS='4')
  entry={'method':j['method'],'gpu':gpu,'state':'running','started_at':time.strftime('%Y-%m-%dT%H:%M:%S%z')};statefile=L/(j['method']+'.status.json');statefile.write_text(json.dumps(entry,indent=2));print('START',j['method'],gpu,flush=True)
  try:
   with (L/(j['method']+'.log')).open('w') as f:
    p=subprocess.Popen([str(ROOT/'.conda-env/bin/python'),'-u',str(OUT/'scripts/run_benchmark_job.py'),j['config']],cwd=ROOT,env=env,stdout=f,stderr=subprocess.STDOUT,start_new_session=True)
    try:rc=p.wait(timeout=12*3600)
    except subprocess.TimeoutExpired:os.killpg(p.pid,signal.SIGKILL);p.wait();raise
   entry.update(state='completed' if rc==0 else 'failed_or_incomplete',exit_code=rc)
  except Exception as e:entry.update(state='failed_or_incomplete',error=repr(e))
  entry['ended_at']=time.strftime('%Y-%m-%dT%H:%M:%S%z');statefile.write_text(json.dumps(entry,indent=2))
  with lock:summary.append(entry);(OUT/'st_queue_status.json').write_text(json.dumps(summary,indent=2))
  print('END',j['method'],entry['state'],flush=True);q.task_done()
with concurrent.futures.ThreadPoolExecutor(3) as ex:list(ex.map(worker,[1,2,3]))
(OUT/'st_queue_finished.json').write_text(json.dumps({'jobs':summary,'all_jobs_completed':all(x['state']=='completed' for x in summary)},indent=2))

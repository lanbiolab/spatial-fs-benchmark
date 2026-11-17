from pathlib import Path
import json,time,psutil,math
O=Path('/home/xd/y/Benchmark/results/revision_r3_comment5_20261107')
def save(p,obj):
 t=p.with_suffix(p.suffix+'.tmp');t.write_text(json.dumps(obj,indent=2));t.replace(p)
supervisors=[]
for p in psutil.process_iter(['cmdline']):
 c=p.info['cmdline'] or []
 if len(c)==3 and c[1]=='-u' and c[2]=='results/revision_r3_comment5_20261107/scripts/run_stereo_queue.py':supervisors.append(p)
assert len(supervisors)==1, 'Expected one live supervisor'
parent=supervisors[0];children=parent.children();assert len(children)==1
job=children[0];cmd=job.cmdline()
assert cmd[-1]==str(O/'configs/stereo_seq/nnsvg.yaml')
assert cmd[-2]==str(O/'scripts/run_benchmark_job.py')
summary=json.loads((O/'stereo_queue_status.json').read_text())
assert len(summary)==33 and all(x['state']=='completed' and x['method']!='nnsvg' for x in summary)
entry=json.loads((O/'logs/nnsvg.status.json').read_text());assert entry['state']=='running'
ident=job.create_time();desc=job.children(recursive=True)
before={p.pid:sum(p.cpu_times()[:2]) for p in desc}
# Only terminate the idle queue parent. The benchmark owns a separate process session.
parent.terminate();parent.wait(timeout=10)
assert job.is_running() and job.create_time()==ident
assert all(p.is_running() for p in desc)
entry.update(supervision='passive monitor without wall-time kill; existing benchmark preserved',monitor_started_at=time.strftime('%Y-%m-%dT%H:%M:%S%z'))
save(O/'logs/nnsvg.status.json',entry)
save(O/'nnsvg_supervisor_transition.json',{'old_supervisor_pid':parent.pid,'preserved_job_pid':job.pid,'preserved_job_create_time':ident,'descendant_cpu_before':before,'changed_experiment_parameters':False})
print('Preserved running job',job.pid,'and descendants',list(before),flush=True)
while job.is_running() and job.status()!=psutil.STATUS_ZOMBIE and job.create_time()==ident:
 save(O/'nnsvg_monitor_heartbeat.json',{'at':time.strftime('%Y-%m-%dT%H:%M:%S%z'),'pid':job.pid,'state':'running'})
 time.sleep(30)
completion=O/'benchmark/stereo_seq/nnsvg/job_completion.json'
result=json.loads(completion.read_text()) if completion.exists() else {'state':'missing_completion'}
records=list((O/'benchmark/stereo_seq/nnsvg').rglob('clustering_eval_records.json'))
valid=[]
for p in records:
 rows=json.loads(p.read_text())
 if rows and all(isinstance(r.get('metric_value'),(int,float)) and math.isfinite(r['metric_value']) for r in rows):valid.append(str(p))
ok=result.get('state')=='completed' and len(valid)==9
entry.update(state='completed' if ok else 'failed_or_incomplete',ended_at=time.strftime('%Y-%m-%dT%H:%M:%S%z'),validated_record_files=len(valid),completion_state=result.get('state'))
save(O/'logs/nnsvg.status.json',entry);summary.append(entry)
save(O/'stereo_queue_status.json',summary)
save(O/'stereo_queue_finished.json',{'jobs':summary,'all_jobs_completed':all(x['state']=='completed' for x in summary)})
print('END',entry['state'],len(valid),flush=True)

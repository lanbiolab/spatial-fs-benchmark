import psutil
roots=[]
for p in psutil.process_iter(['pid','username','cmdline']):
 try:
  a=p.info['cmdline'] or []
  if p.info['username']=='xd' and any('capture_selection.py' in x for x in a) and any(x in ['e8p5_embryo','e9p5_embryo'] for x in a):roots.append(p)
 except (psutil.NoSuchProcess,psutil.AccessDenied):pass
changed=set()
for p in roots:
 for q in [p]+p.children(recursive=True):
  try:q.nice(10);changed.add(q.pid)
  except (psutil.NoSuchProcess,psutil.AccessDenied):pass
print('Deprioritized embryo-specific capture processes:',len(changed),'roots:',[(p.pid,p.cmdline()[-3:]) for p in roots if p.is_running()])

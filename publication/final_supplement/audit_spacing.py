from pathlib import Path
from fontTools.ttLib import TTFont
import pymupdf as f
import json,math,collections
root=Path.cwd()
cache={}
def get_metrics(bold):
 if bold not in cache:
  font=TTFont('C:/Windows/Fonts/'+('arialbd.ttf' if bold else 'arial.ttf'))
  cmap=font.getBestCmap();units=font['head'].unitsPerEm
  kern={}
  if 'kern' in font:
   for table in font['kern'].kernTables:kern.update(table.kernTable)
  cache[bold]=(cmap,font['hmtx'].metrics,units,kern)
 return cache[bold]
def audit(path):
 d=f.open(path);counts=collections.Counter();bad=[];errors=[]
 for page in d:
  for block in page.get_text('rawdict')['blocks']:
   for line in block.get('lines',[]):
    dx,dy=line['dir']
    for s in line['spans']:
     if 'Arial' not in s['font']:continue
     chars=s['chars'];counts['text_runs']+=1;counts['characters']+=len(chars)
     counts['bold_runs']+=int('Bold' in s['font']);counts['rotated_runs']+=int(abs(dy)>.01)
     cmap,metrics,units,kern=get_metrics('Bold' in s['font'])
     for a,b in zip(chars,chars[1:]):
      if a['c'].isspace() or b['c'].isspace() or a['c']== '\ufffd' or b['c']=='\ufffd':continue
      ga,gb=cmap.get(ord(a['c'])),cmap.get(ord(b['c']))
      if not ga or not gb:continue
      advance=(b['origin'][0]-a['origin'][0])*dx+(b['origin'][1]-a['origin'][1])*dy
      expected=(metrics[ga][0]+kern.get((ga,gb),0))*s['size']/units
      error=abs(advance-expected);errors.append(error)
      if error>.04:bad.append({'text':''.join(c['c'] for c in chars),'pair':a['c']+b['c'],'advance':advance,'expected':expected,'error':error,'origin':a['origin']})
 return {'counts':dict(counts),'pairs_checked':len(errors),'max_error_pt':max(errors,default=0),'pairs_above_0.04pt':len(bad),'flagged_pairs':bad}
if __name__=='__main__':
 import sys
 dirs=[Path(p) for p in sys.argv[1:]] if len(sys.argv)>1 else [root/'outputs/GR_Main_Figures_Arial_Spacing_20261109/figures']
 report={str(p):audit(p) for directory in dirs for p in sorted(directory.glob('*.pdf'))}
 target=root/'work/all_typography/all_spacing_audit.json';target.write_text(json.dumps(report,indent=2,ensure_ascii=False),encoding='utf8')
 for n,r in report.items():print(Path(n).name,r['counts'],r['pairs_checked'],r['pairs_above_0.04pt'],round(r['max_error_pt'],4))

from pathlib import Path
from zipfile import ZipFile,ZIP_DEFLATED
import pymupdf as f,shutil,collections,json,re
R=Path.cwd();W=R/'work/font_spacing';O=R/'outputs/GR_Main_Figures_Arial_Spacing_20261109';P=O/'main_source';V=W/'plots';OLD=R/'outputs/GR_Main_Figures_Arial_20261109/main_source'
O.mkdir(exist_ok=True);shutil.copytree(OLD,P,dirs_exist_ok=True);F=O/'figures';F.mkdir(exist_ok=True);Q=O/'qa';Q.mkdir(exist_ok=True)
def compose(name,w,h,parts,tags=()):
 d=f.open();p=d.new_page(width=w,height=h)
 for fn,rect in parts:
  with f.open(V/fn) as src:p.show_pdf_page(f.Rect(rect),src,0)
 if tags:
  p.insert_font(fontname='ArialBold',fontfile='C:/Windows/Fonts/arialbd.ttf')
  for label,x,y in tags:p.insert_text((x,y),label,fontsize=9,fontname='ArialBold')
 d.subset_fonts();d.save(F/name,garbage=4,deflate=True)
compose('Figure1.pdf',540,585,[('Figure1a.pdf',(0,3,540,223)),('Figure1b.pdf',(0,225,540,585))],[('a',2,10),('b',2,234)])
for n in [2,3,5]:shutil.copy2(V/f'Figure{n}.pdf',F/f'Figure{n}.pdf')
compose('Figure4.pdf',540,658,[('Figure4a.pdf',(0,2,540,292)),('Figure4b.pdf',(0,294,252,524)),('Figure4c.pdf',(252,294,378,524)),('Figure4d.pdf',(378,294,540,524)),('Figure4e.pdf',(0,526,540,658))])
compose('Figure6.pdf',590,604,[('Figure6abcd.pdf',(0,0,590,438)),('Figure6e.pdf',(0,450,590,604))])
normalize=lambda t:re.sub(r'\s+','',t.translate(str.maketrans({'−':'-','–':'-','‑':'-','‐':'-','’':"'"})))
report={}
for p in F.glob('*.pdf'):
 shutil.copy2(p,P/'figures'/p.name)
 d=f.open(p);old=f.open(OLD/'figures'/p.name)
 fonts=collections.Counter(s['font'] for b in d[0].get_text('dict')['blocks'] for l in b.get('lines',[]) for s in l['spans'])
 assert all('Arial' in n for n in fonts),fonts
 for row in d[0].get_fonts(full=True):
  assert d.extract_font(row[0])[3],(p,row)
 a=collections.Counter(normalize(old[0].get_text()));b=collections.Counter(normalize(d[0].get_text()))
 report[p.name]={'fonts':dict(fonts),'all_fonts_embedded':True,'old_dimensions':list(old[0].rect),'new_dimensions':list(d[0].rect),'text_character_difference_old':dict(a-b),'text_character_difference_new':dict(b-a)}
 d[0].get_pixmap(matrix=f.Matrix(1.8,1.8),alpha=False).save(str(Q/(p.stem+'.png')))
(O/'Figure_Font_Audit.json').write_text(json.dumps(report,indent=2,ensure_ascii=False),encoding='utf8')
with ZipFile(W/'build.zip','w',ZIP_DEFLATED) as z:
 for n in ['main.tex','main_clean.tex','references.bib','authblk.sty','plainnat.bst']:z.write(P/n,n)
 for p in (P/'figures').glob('*.pdf'):z.write(p,'figures/'+p.name)
print(json.dumps(report,indent=2,ensure_ascii=False))

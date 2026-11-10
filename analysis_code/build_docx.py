
import re, docx
from docx.shared import Pt, Inches, RGBColor
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.oxml.ns import qn
from docx.oxml import OxmlElement

FONT="Times New Roman"
md=open("Manuscript_Frontiers_format.md").read()

# (a) strip ONLY a leading blockquote drafting note, never truncate at a marker
lines=md.split("\n"); k=0
while k<len(lines) and not lines[k].strip(): k+=1
if k<len(lines) and lines[k].startswith(">"):
    while k<len(lines) and (lines[k].startswith(">") or not lines[k].strip()): k+=1
    lines=lines[k:]

doc=docx.Document()
sec=doc.sections[0]
for a in ("top_margin","bottom_margin","left_margin","right_margin"): setattr(sec,a,Inches(1))

def inline(p, text):
    text=re.sub(r"\\(.)", r"\1", text)
    for tok in re.split(r"(\*\*[^*]+\*\*|\*[^*]+\*)", text):
        if not tok: continue
        if tok.startswith("**") and tok.endswith("**") and len(tok)>4:
            r=p.add_run(tok[2:-2]); r.bold=True
        elif tok.startswith("*") and tok.endswith("*") and len(tok)>2:
            r=p.add_run(tok[1:-1]); r.italic=True
        else: r=p.add_run(tok)
        r.font.name=FONT; r.font.size=Pt(11)

def body(txt, indent=0.0, hang=0.0, space=6):
    p=doc.add_paragraph(); p.paragraph_format.space_after=Pt(space)
    p.paragraph_format.line_spacing=1.5
    if indent: p.paragraph_format.left_indent=Inches(indent)
    if hang: p.paragraph_format.first_line_indent=Inches(-hang)
    inline(p, txt); return p

# group into blocks
blocks=[]; buf=[]
for ln in lines:
    s=ln.rstrip()
    if s.strip()=="---": continue
    if s.startswith("#") or s.startswith("|") or re.match(r"^\s*(\d+\.|[-*])\s+", s) or not s.strip():
        if buf: blocks.append(("p"," ".join(buf))); buf=[]
        if s.startswith("#"): blocks.append(("h",s))
        elif s.startswith("|"): blocks.append(("t",s))
        elif s.strip(): blocks.append(("l",s))
    else:
        # continuation of a list item keeps its item
        if blocks and blocks[-1][0]=="l" and ln.startswith("   ") and not buf:
            blocks[-1]=("l", blocks[-1][1].rstrip()+" "+s.strip())
        else: buf.append(s.strip())
if buf: blocks.append(("p"," ".join(buf)))

# merge consecutive table rows
merged=[]
for kind,txt in blocks:
    if kind=="t" and merged and merged[-1][0]=="t": merged[-1][1].append(txt)
    elif kind=="t": merged.append(["t",[txt]])
    else: merged.append([kind,txt])

FIGW=6.5
_title_done=[False]
for kind,txt in merged:
    if kind=="h":
        lvl=len(txt)-len(txt.lstrip("#")); t=txt.lstrip("# ").strip()
        is_title = (lvl==1 and not _title_done[0])
        if is_title: _title_done[0]=True
        if is_title:
            p=doc.add_paragraph(); p.paragraph_format.space_after=Pt(12)
            r=p.add_run(t); r.bold=True; r.font.size=Pt(16); r.font.name=FONT
        else:
            h=doc.add_heading(level=min(lvl,2)); h.paragraph_format.space_before=Pt(12)
            h.paragraph_format.space_after=Pt(6)
            r=h.add_run(t); r.font.name=FONT; r.font.color.rgb=RGBColor(0,0,0)
            r.font.size=Pt(14 if lvl==1 else 12); r.bold=True
    elif kind=="l":
        s=txt.strip()
        m=re.match(r"^(\d+\.)\s+(.*)$", s)
        if m: body(m.group(1)+"\t"+m.group(2), indent=0.35, hang=0.35)
        else: body("\u2022\t"+re.sub(r"^[-*]\s+","",s), indent=0.35, hang=0.35)
    elif kind=="t":
        rows=[r for r in txt if not re.match(r"^\|[\s:\-|]+\|$", r)]
        cells=[[c.strip() for c in r.strip().strip("|").split("|")] for r in rows]
        tb=doc.add_table(rows=len(cells), cols=len(cells[0])); tb.style="Table Grid"
        for i,row in enumerate(cells):
            for j,c in enumerate(row):
                cp=tb.cell(i,j).paragraphs[0]; cp.paragraph_format.space_after=Pt(2)
                inline(cp,c)
                for r in cp.runs:
                    r.font.size=Pt(9)
                    if i==0: r.bold=True
        doc.add_paragraph()
    else:
        p=body(txt)
        m=re.match(r"^\*\*Figure (\d+)\.\*\*", txt)
        if m:
            fp=f"submission/Figure{m.group(1)}.png"
            import os
            if os.path.exists(fp):
                ip=doc.add_paragraph(); ip.alignment=WD_ALIGN_PARAGRAPH.CENTER
                ip.paragraph_format.space_after=Pt(14)
                ip.add_run().add_picture(fp, width=Inches(FIGW))

# ---- font enforcement, all four places ----
def setfont(rPr):
    rf=rPr.find(qn('w:rFonts'))
    if rf is None:
        rf=OxmlElement('w:rFonts'); rPr.insert(0,rf)
    for a in ('w:ascii','w:hAnsi','w:cs','w:eastAsia'): rf.set(qn(a),FONT)
    for a in ('w:asciiTheme','w:hAnsiTheme','w:cstheme','w:eastAsiaTheme'):
        if rf.get(qn(a)) is not None: del rf.attrib[qn(a)]

st=doc.styles.element
dd=st.find(qn('w:docDefaults'))
rpd=dd.find(qn('w:rPrDefault'))
rPr=rpd.find(qn('w:rPr'))
if rPr is None: rPr=OxmlElement('w:rPr'); rpd.append(rPr)
setfont(rPr)
n_sty=0
for s in st.findall(qn('w:style')):
    sr=s.find(qn('w:rPr'))
    if sr is None: sr=OxmlElement('w:rPr'); s.append(sr)
    setfont(sr); n_sty+=1
n_tbl=0
for tsp in st.iter(qn('w:tblStylePr')):
    tr=tsp.find(qn('w:rPr'))
    if tr is None: tr=OxmlElement('w:rPr'); tsp.insert(0,tr)
    setfont(tr); n_tbl+=1
n_run=0
for r in doc.element.body.iter(qn('w:r')):
    rp=r.find(qn('w:rPr'))
    if rp is None: rp=OxmlElement('w:rPr'); r.insert(0,rp)
    setfont(rp); n_run+=1

doc.save("Manuscript_FrontiersImmunol_format.docx")
print("styles",n_sty,"tblStylePr",n_tbl,"runs",n_run,
      "imgs",len(docx.Document("Manuscript_FrontiersImmunol_format.docx").inline_shapes),
      "tables",len(docx.Document("Manuscript_FrontiersImmunol_format.docx").tables))

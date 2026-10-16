#!/usr/bin/env python3
"""
13_master_table_collapse.py
---------------------------
Rebuild pi_theta_MASTER_table under the "multiallelic = one site" (collapse)
convention, per complex + overall, for North / South / Between, with 95%
bootstrap CIs. Matches Table3_collapse; the grey-column figure can be redrawn
from this. Outputs pi_theta_MASTER_table_collapse.tsv/.xlsx
"""
import csv, random
from collections import defaultdict
random.seed(11); B=5000
def a_n(n): return sum(1.0/i for i in range(1,n))
AN={"Between":a_n(141),"North":a_n(77),"South":a_n(63)}
NOF={"Between":141,"North":77,"South":63}
def cls(vt): return "ns" if "missense" in vt else "syn"
COMPLEX={"I":["ND1","ND2","ND3","ND4","ND4L","ND5","ND6"],"III":["CYTB"],"IV":["COX1","COX2","COX3"],"V":["ATP6","ATP8"]}
G2C={g:c for c,gs in COMPLEX.items() for g in gs}

mem={}
for r in csv.DictReader(open("../vcf/pi_by_clade_persite.tsv.membership.tsv"),delimiter="\t"):
    mem[r["Sample"].replace("_MT","")]=r["Clade"]
MEMB={"North":set(s for s,c in mem.items() if c=="north"),
      "South":set(s for s,c in mem.items() if c=="south"),"Between":None}
rows=list(csv.DictReader(open("../vcf/haplotype_matrix.csv")))
idcols=[k for k in rows[0] if k.endswith("_MT")]
byPOS=defaultdict(list)
for r in rows: byPOS[r["POS"]].append(r)

def collapsed_sites(comp):
    """Return list of (gene, class, pi_site) — one entry per polymorphic position."""
    members=MEMB[comp]; ids=[k for k in idcols if (members is None or k[:-3] in members)]
    n=len(ids); out=[]
    for pos,rlist in byPOS.items():
        alts=[(sum(1 for k in ids if r[k]=="1"), cls(r["Variant_type"]), r["Gene"]) for r in rlist]
        tot=sum(a for a,_,_ in alts); ref=n-tot
        present=[a for a,_,_ in alts if a>0]
        if (1 if ref>0 else 0)+len(present) <= 1: continue      # monomorphic
        freqs=[ref]+present
        pis=(1.0-sum((f/n)**2 for f in freqs))*n/(n-1)
        maj=max([a for a in alts if a[0]>0],key=lambda x:x[0])
        out.append((maj[2], maj[1], pis))
    return out,AN[comp]

DATA={c:collapsed_sites(c) for c in ("North","South","Between")}
def sub(comp,group,cls_):
    sites,a=DATA[comp]
    return [pi for g,c,pi in sites if (group=="Overall" or G2C.get(g)==group.split("_")[-1]) and (cls_ is None or c==cls_)],a
def pt(pis,a):
    S=len(pis); return (a*sum(pis)/S) if S else None
def bootpt(pis,a):
    S=len(pis)
    if S<2: return (None,None)
    v=sorted(a*sum(pis[random.randrange(S)] for _ in range(S))/S for _ in range(B)); return (v[int(.025*B)],v[int(.975*B)])
def ratio_ci(comp,group):
    sites,a=DATA[comp]
    sub_=[(c,pi) for g,c,pi in sites if group=="Overall" or G2C.get(g)==group.split("_")[-1]]
    def r(s):
        sy=[pi for c,pi in s if c=="syn"]; ns=[pi for c,pi in s if c=="ns"]
        return (a*sum(ns)/len(ns))/(a*sum(sy)/len(sy)) if (sy and ns) else None
    pe=r(sub_); vals=[]
    for _ in range(B):
        bs=[sub_[random.randrange(len(sub_))] for _ in sub_]; v=r(bs)
        if v is not None: vals.append(v)
    vals.sort(); ci=(vals[int(.025*len(vals))],vals[int(.975*len(vals))]) if len(vals)>10 else (None,None)
    return pe,ci
def f(x,nd=3): return "" if x is None else "%.*f"%(nd,x)
groups=[("Complex_I","Complex I"),("Complex_III","Complex III"),("Complex_IV","Complex IV"),("Complex_V","Complex V"),("Overall","Overall")]
cols=["Comparison","Group","n","S_syn","piθ_Syn","Syn_CIlo","Syn_CIhi","S_NS","piθ_NS","NS_CIlo","NS_CIhi","piθ_NS/S","ratio_CIlo","ratio_CIhi"]
recs=[]
for gk,glab in groups:
    for comp in ("North","South","Between"):
        syn,a=sub(comp,gk,"syn"); ns,_=sub(comp,gk,"ns")
        peS=pt(syn,a); loS,hiS=bootpt(syn,a); peN=pt(ns,a); loN,hiN=bootpt(ns,a)
        rr,(rlo,rhi)=ratio_ci(comp,gk)
        recs.append([comp,glab,NOF[comp],len(syn),f(peS),f(loS),f(hiS),len(ns),f(peN),f(loN),f(hiN),f(rr),f(rlo),f(rhi)])

CAP=("Table. π/θ (mean, 95%% bootstrap CI) by complex and overall for North (n=77), South (n=63) and "
     "Between (whole panel, n=141), MULTIALLELIC POSITIONS COUNTED ONCE (collapse convention; π_site="
     "n/(n-1)(1-Σp²)). piθ_NS/S = nonsyn/syn ratio (1=neutral). Consistent with Table3_collapse. This is π/θ, not pN/pS.")
with open("pi_theta_MASTER_table_collapse.tsv","w") as o:
    o.write(CAP.replace("%%","%")+"\n"); o.write("\t".join(cols)+"\n")
    for r in recs: o.write("\t".join(str(x) for x in r)+"\n")

from openpyxl import Workbook
from openpyxl.styles import Font,PatternFill,Alignment,Border,Side
from openpyxl.utils import get_column_letter
wb=Workbook(); ws=wb.active; ws.title="pi_theta_master_collapse"
hf=Font(name="Arial",bold=True,color="FFFFFF",size=10)
fmeta=PatternFill("solid",fgColor="404040");fsyn=PatternFill("solid",fgColor="2E7D32");fns=PatternFill("solid",fgColor="C55A11");frat=PatternFill("solid",fgColor="595959")
base=Font(name="Arial",size=10);thin=Side(style="thin",color="D9D9D9");bd=Border(left=thin,right=thin,top=thin,bottom=thin);ctr=Alignment(horizontal="center")
cmpfill={"North":PatternFill("solid",fgColor="DDEBF7"),"South":PatternFill("solid",fgColor="FCE4E4"),"Between":PatternFill("solid",fgColor="E2EFDA")}
ws["A1"]="π/θ master table — multiallelic = one site (collapse)"; ws["A1"].font=Font(name="Arial",bold=True,size=12)
ws.cell(2,1,CAP.replace("%%","%")).font=Font(name="Arial",italic=True,size=9,color="444444"); ws.merge_cells("A2:N2")
grpfill=[fmeta,fmeta,fmeta,fsyn,fsyn,fsyn,fsyn,fns,fns,fns,fns,frat,frat,frat]
r0=4
for j,h in enumerate(cols,1):
    c=ws.cell(r0,j,h);c.font=hf;c.fill=grpfill[j-1];c.border=bd;c.alignment=ctr
for i,rec in enumerate(recs,r0+1):
    for j,x in enumerate(rec,1):
        try: xv=float(x) if ("." in str(x)) else (int(x) if str(x).lstrip("-").isdigit() else x)
        except: xv=x
        c=ws.cell(i,j,xv);c.border=bd;c.font=base
        if j>=3: c.alignment=ctr
        if isinstance(xv,float): c.number_format="0.000"
        if j<=2: c.fill=cmpfill[rec[0]]
    if rec[0]=="Between":
        for j in range(1,len(cols)+1): ws.cell(i,j).border=Border(left=thin,right=thin,top=thin,bottom=Side(style="medium",color="808080"))
ws.freeze_panes="C5"
for j,h in enumerate(cols,1): ws.column_dimensions[get_column_letter(j)].width=max(8,len(h)+1)
ws.column_dimensions["A"].width=11; ws.column_dimensions["B"].width=11
wb.save("pi_theta_MASTER_table_collapse.xlsx")
print("wrote pi_theta_MASTER_table_collapse.tsv/.xlsx")
print("\n=== collapse master (Overall rows) ===")
for r in recs:
    if r[1]=="Overall": print("  %-8s syn %s[%s-%s] NS %s[%s-%s] NS/S %s[%s-%s]"%(r[0],r[4],r[5],r[6],r[8],r[9],r[10],r[11],r[12],r[13]))

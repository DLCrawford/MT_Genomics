#!/usr/bin/env python3
"""
11_pitheta_master_table.py
--------------------------
THE single table behind figure_pitheta_greycolumn: mean pi/theta with 95% CIs.
One row per (Group x Comparison); everything the figure plots, nothing else.

For each of Complex I, III, IV, V, Overall  x  {North, South, Between}:
  S_syn, pi/theta synonymous (mean, 95% CI)
  S_NS,  pi/theta nonsynonymous (mean, 95% CI)
  pi/theta NonSyn/Syn ratio (mean, 95% CI)

Everything is pi/theta (NOT pN/pS). North n=77, South n=63 (within-clade,
segregating-within-clade sites); Between = whole panel n=141.
Sources: pi_theta_withinclade_by_complex.tsv, pi_theta_wholepanel_by_complex.tsv,
pi_theta_NS_over_S_ratio.tsv.
Outputs: pi_theta_MASTER_table.tsv + .xlsx
"""
import csv
from common import opath
from openpyxl import Workbook
from openpyxl.styles import Font, PatternFill, Alignment, Border, Side
from openpyxl.utils import get_column_letter
def R(p): return list(csv.DictReader(open(opath(p)), delimiter="\t"))
def num(x):
    try: return float(x)
    except: return None
WC=R("pi_theta_withinclade_by_complex.tsv"); WP=R("pi_theta_wholepanel_by_complex.tsv"); RT=R("pi_theta_NS_over_S_ratio.tsv")
def wc(cl,g,c): return next(r for r in WC if r["Clade"]==cl and r["Group"]==g and r["Class"]==c)
def wp(g,c):    return next(r for r in WP if r["Group"]==g and r["Class"]==c)
def rt(cmp,g):  return next(r for r in RT if r["Comparison"]==cmp and r["Group"]==g.replace("_"," "))
groups=[("Complex_I","Complex I"),("Complex_III","Complex III"),("Complex_IV","Complex IV"),("Complex_V","Complex V"),("Overall","Overall")]
COMPS=[("North","syn","nonsynonymous"),("South",None,None),("Between",None,None)]

cols=["Comparison","Group","n_clade","S_syn","piθ_Syn","Syn_CIlow","Syn_CIhi",
      "S_NS","piθ_NS","NS_CIlow","NS_CIhi","piθ_NS/S","ratio_CIlow","ratio_CIhi"]
N_OF={"North":77,"South":63,"Between":141}
def syn_row(cmp,gk):
    if cmp=="Between": return wp(gk,"synonymous")
    return wc(cmp,gk,"synonymous")
def ns_row(cmp,gk):
    if cmp=="Between": return wp(gk,"nonsynonymous")
    return wc(cmp,gk,"nonsynonymous")

rows=[]
for gk,glab in groups:
    for cmp in ("North","South","Between"):
        s=syn_row(cmp,gk); n=ns_row(cmp,gk); r=rt(cmp,gk)
        rows.append([cmp,glab,N_OF[cmp],
                     int(float(s["S_sites"])), num(s["pi_theta"]), num(s["CI_low"]), num(s["CI_high"]),
                     int(float(n["S_sites"])), num(n["pi_theta"]), num(n["CI_low"]), num(n["CI_high"]),
                     num(r["NS_over_S"]), num(r["CI_low"]), num(r["CI_high"])])

CAP=("Table. pi/theta (nucleotide diversity / Watterson theta) with 95%% bootstrap CIs (B=5000), "
     "for the 4 mtDNA-encoded OXPHOS complexes and overall. North (n=77) and South (n=63) use sites "
     "segregating within that clade; Between = whole panel (n=141). piθ_NS/S is the nonsynonymous-to-"
     "synonymous ratio of pi/theta (1 = neutral). This is pi/theta, not pN/pS. This table is the basis "
     "of the grey-column figure.")
with open(opath("pi_theta_MASTER_table.tsv"),"w") as o:
    o.write(CAP+"\n"); o.write("\t".join(cols)+"\n")
    for r in rows: o.write("\t".join("" if x is None else str(x) for x in r)+"\n")

# ---- xlsx ----
wb=Workbook(); ws=wb.active; ws.title="pi_theta_master"
hf=Font(name="Arial",bold=True,color="FFFFFF",size=10)
fmeta=PatternFill("solid",fgColor="404040"); fsyn=PatternFill("solid",fgColor="2E7D32"); fns=PatternFill("solid",fgColor="C55A11"); frat=PatternFill("solid",fgColor="595959")
base=Font(name="Arial",size=10); thin=Side(style="thin",color="D9D9D9"); bd=Border(left=thin,right=thin,top=thin,bottom=thin); ctr=Alignment(horizontal="center")
cmpfill={"North":PatternFill("solid",fgColor="DDEBF7"),"South":PatternFill("solid",fgColor="FCE4E4"),"Between":PatternFill("solid",fgColor="E2EFDA")}
ws["A1"]="π/θ master table (mean with 95% CI) — basis of the grey-column figure"; ws["A1"].font=Font(name="Arial",bold=True,size=12)
ws.cell(2,1,CAP.replace("%%","%")).font=Font(name="Arial",italic=True,size=9,color="444444"); ws.merge_cells("A2:N2")
r0=4
grpfill=[fmeta,fmeta,fmeta,fsyn,fsyn,fsyn,fsyn,fns,fns,fns,fns,frat,frat,frat]
for j,h in enumerate(cols,1):
    c=ws.cell(r0,j,h); c.font=hf; c.fill=grpfill[j-1]; c.border=bd; c.alignment=ctr
for i,row in enumerate(rows,r0+1):
    for j,x in enumerate(row,1):
        c=ws.cell(i,j,x); c.border=bd; c.font=base
        if j>=3: c.alignment=ctr
        if isinstance(x,float): c.number_format="0.000"
        if j<=2: c.fill=cmpfill[row[0]]
    if row[0]=="Between":  # thin separator under each group's Between row
        for j in range(1,len(cols)+1): ws.cell(i,j).border=Border(left=thin,right=thin,top=thin,bottom=Side(style="medium",color="808080"))
ws.freeze_panes="C5"
for j,h in enumerate(cols,1):
    ws.column_dimensions[get_column_letter(j)].width=max(8,len(h)+1.5)
ws.column_dimensions["A"].width=11; ws.column_dimensions["B"].width=11
wb.save(opath("pi_theta_MASTER_table.xlsx"))
print("wrote pi_theta_MASTER_table.tsv and .xlsx  (%d rows)"%len(rows))

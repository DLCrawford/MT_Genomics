#!/usr/bin/env python3
"""
06_tables_figures_Sept28.py
---------------------------
Deliverables requested Sept 28:
  (1) pN/pS Table 1 replacing the old dN/dS table (same layout: Codons, S_sites,
      N_sites, obs_S, obs_N) with pS, pN, pN/pS and 95% bootstrap CI + Overall row.
        -> pN_pS_Table1.txt  (tab, mirrors Old_dN_dS.txt)  + sheet in xlsx
  (2) Figure: A = pN/pS per gene (13) with CI; B = pN/pS per complex + overall with CI.
        -> figure_pnps_gene_complex.png/.pdf
  (3) Figure: A = within-clade pi/theta (North & South, syn/NS) per complex+overall;
              B = between-clade (whole-panel) pi/theta syn vs NS per complex+overall.
        -> figure_pitheta_within_between.png/.pdf

Reads the frozen TSVs produced by scripts 01-03 and ../vcf/dnds_per_gene.tsv.
"""
import csv
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from common import GENES, GENE2CX, vpath, opath
from openpyxl import Workbook
from openpyxl.styles import Font, PatternFill, Alignment, Border, Side
from openpyxl.utils import get_column_letter

def val(x):
    try: return float(x)
    except: return np.nan
def rd_dict(path):
    return list(csv.DictReader(open(path), delimiter="\t"))

# ---------- load site counts (codons, S_sites, N_sites, obs) ----------
DN = {r["Gene"]: r for r in rd_dict(vpath("dnds_per_gene.tsv")) if r["Gene"] and r["Gene"] != "Overall"}
DNall = next(r for r in rd_dict(vpath("dnds_per_gene.tsv")) if r["Gene"] == "Overall")
PN = {r["Gene"]: r for r in rd_dict(opath("pnps_bootstrap_by_gene.tsv"))}   # pN_pS, CI
PNC = {r["Group"]: r for r in rd_dict(opath("pnps_bootstrap_by_complex.tsv"))}

# =====================================================================
# (1) pN/pS Table 1
# =====================================================================
cols = ["Gene","Strand","Codons","S_sites","N_sites","obs_S","obs_N","pS","pN","pN_pS","CI_low","CI_high"]
def gene_row(g):
    d = DN[g]; p = PN[g]
    Ss=float(d["S_sites"]); Ns=float(d["N_sites"]); oS=int(d["obs_S"]); oN=int(d["obs_N"])
    return [g, d["Strand"], d["n_codons"], round(Ss,1), round(Ns,1), oS, oN,
            round(oS/Ss,4), round(oN/Ns,4), val(p["pN_pS"]), val(p["CI_low"]), val(p["CI_high"])]
def overall_row():
    Ss=float(DNall["S_sites"]); Ns=float(DNall["N_sites"]); oS=int(DNall["obs_S"]); oN=int(DNall["obs_N"])
    p=PNC["Overall"]
    return ["Overall","", DNall["n_codons"], round(Ss,1), round(Ns,1), oS, oN,
            round(oS/Ss,4), round(oN/Ns,4), val(p["pN_pS"]), val(p["CI_low"]), val(p["CI_high"])]

rows_tbl = [gene_row(g) for g in GENES] + [overall_row()]
CAPTION = ("Table 1. pN, pS and pN/pS for 13 mitochondrial proteins (141 individuals). "
           "Codons = codons per gene. obs_S / obs_N = observed synonymous / nonsynonymous "
           "segregating sites (polymorphisms). S_sites / N_sites = Nei-Gojobori synonymous / "
           "nonsynonymous site opportunities. pS = obs_S/S_sites, pN = obs_N/N_sites, pN/pS = pN/pS. "
           "CI = 95% percentile bootstrap (B=5000, sites resampled). pN/pS is a within-sample "
           "polymorphism ratio and does not depend on the reference sequence.")
with open(opath("pN_pS_Table1.txt"), "w") as o:
    o.write("pN_pS\n"+CAPTION+"\n")
    o.write("\t".join(cols)+"\n")
    for r in rows_tbl:
        o.write("\t".join(str(x) for x in r)+"\n")

# =====================================================================
# (2) Figure: pN/pS gene (A) + complex (B)
# =====================================================================
fig, ax = plt.subplots(1, 2, figsize=(14, 5), gridspec_kw={"width_ratios":[1.7,1]})
# A: 13 genes, ordered by complex
order = [g for cx in ("I","III","IV","V") for g in
         ["ND1","ND2","ND3","ND4","ND4L","ND5","ND6","CYTB","COX1","COX2","COX3","ATP6","ATP8"] if GENE2CX[g]==cx]
cxcol = {"I":"#2F5496","III":"#548235","IV":"#C55A11","V":"#7030A0"}
pe=[val(PN[g]["pN_pS"]) for g in order]; lo=[val(PN[g]["CI_low"]) for g in order]; hi=[val(PN[g]["CI_high"]) for g in order]
x=np.arange(len(order)); err=[np.array(pe)-np.array(lo),np.array(hi)-np.array(pe)]
bar_c=[cxcol[GENE2CX[g]] for g in order]
ax[0].bar(x,pe,yerr=err,capsize=3,color=bar_c,edgecolor="black",linewidth=0.5)
ax[0].set_xticks(x); ax[0].set_xticklabels(order,rotation=45,ha="right",fontsize=9)
ax[0].set_ylabel("pN/pS"); ax[0].set_title("A. pN/pS per gene (95% CI)",fontsize=11)
from matplotlib.patches import Patch
ax[0].legend(handles=[Patch(facecolor=cxcol[c],label="Complex "+c) for c in ("I","III","IV","V")],fontsize=8,title=None)
# B: complexes + overall
grp=["Complex_I","Complex_III","Complex_IV","Complex_V","Overall"]
pe=[val(PNC[g]["pN_pS"]) for g in grp]; lo=[val(PNC[g]["CI_low"]) for g in grp]; hi=[val(PNC[g]["CI_high"]) for g in grp]
x=np.arange(len(grp)); err=[np.array(pe)-np.array(lo),np.array(hi)-np.array(pe)]
cb=[cxcol.get(g.split("_")[1],"#c0392b") if g!="Overall" else "#c0392b" for g in grp]
ax[1].bar(x,pe,yerr=err,capsize=4,color=cb,edgecolor="black",linewidth=0.6)
ax[1].set_xticks(x); ax[1].set_xticklabels([g.replace("_"," ") for g in grp],rotation=25,ha="right",fontsize=9)
ax[1].set_title("B. pN/pS per complex (95% CI)",fontsize=11); ax[1].set_ylabel("pN/pS")
plt.tight_layout()
plt.savefig(opath("figure_pnps_gene_complex.png"),dpi=200,bbox_inches="tight")
plt.savefig(opath("figure_pnps_gene_complex.pdf"),bbox_inches="tight")
plt.savefig(opath("figure_pnps_gene_complex.svg"),bbox_inches="tight")
plt.close()

# =====================================================================
# (3) Figure: pi/theta within-clade (A) + between-clade (B)
# =====================================================================
WC = rd_dict(opath("pi_theta_withinclade_by_complex.tsv"))
WP = rd_dict(opath("pi_theta_wholepanel_by_complex.tsv"))
groups=["Complex_I","Complex_III","Complex_IV","Complex_V","Overall"]
def wc_pick(clade,g,cls): return next(r for r in WC if r["Clade"]==clade and r["Group"]==g and r["Class"]==cls)
def wp_pick(g,cls):       return next(r for r in WP if r["Group"]==g and r["Class"]==cls)

fig, ax = plt.subplots(1, 3, figsize=(19, 5))
def grouped(a, combos, title, ymax=None):
    xb=np.arange(len(groups)); w=0.8/len(combos)
    for k,(lab,color,pick) in enumerate(combos):
        pe=[val(pick(g)["pi_theta"]) for g in groups]; lo=[val(pick(g)["CI_low"]) for g in groups]; hi=[val(pick(g)["CI_high"]) for g in groups]
        err=[np.array(pe)-np.array(lo),np.array(hi)-np.array(pe)]
        off=(k-(len(combos)-1)/2)*w
        a.bar(xb+off,pe,w,yerr=err,capsize=3,label=lab,color=color,edgecolor="black",linewidth=0.4)
    a.set_xticks(xb); a.set_xticklabels([g.replace("_"," ") for g in groups],rotation=20,ha="right")
    a.set_ylabel("π / θ"); a.set_title(title,fontsize=11)
    if ymax: a.set_ylim(0,ymax)
# A: within-clade N & S, syn/NS
grouped(ax[0],
        [("North syn","#1f77b4",lambda g:wc_pick("North",g,"synonymous")),
         ("North NS","#aec7e8",lambda g:wc_pick("North",g,"nonsynonymous")),
         ("South syn","#d62728",lambda g:wc_pick("South",g,"synonymous")),
         ("South NS","#ff9896",lambda g:wc_pick("South",g,"nonsynonymous"))],
        "A. Within-clade π/θ  (North n=77, South n=63)")
ax[0].legend(fontsize=8,ncol=2)
# B: between-clade (whole panel) syn vs NS
grouped(ax[1],
        [("synonymous","#27ae60",lambda g:wp_pick(g,"synonymous")),
         ("nonsynonymous","#e67e22",lambda g:wp_pick(g,"nonsynonymous"))],
        "B. Between-clade (whole panel) π/θ")
ax[1].axhline(1.0,ls="--",color="grey",lw=0.8); ax[1].legend(fontsize=9)
plt.tight_layout()

# ---- Panel C: NS/S ratio of pi/theta, with bootstrap CI, per complex+overall ----
import random as _rnd
from common import A_N as _AN, GENE2CX as _G2C
_rnd.seed(7); _BR=5000
# raw per-site (gene,cls,pi): whole panel
_wp=[]
for _r in csv.DictReader(open(vpath("haplotype_matrix.csv"))):
    _gt=[_r[k] for k in _r if k.endswith("_MT")]; _n=len(_gt); _AC=sum(1 for x in _gt if x=="1")
    if _AC==0 or _AC==_n: continue
    _vt=_r["Variant_type"]; _c="syn" if ("synonymous" in _vt and "missense" not in _vt) else ("ns" if "missense" in _vt else None)
    if _c: _wp.append((_r["Gene"],_c,2.0*_AC*(_n-_AC)/(_n*(_n-1))))
_gm={(x["POS"],x["REF"],x["ALT"]):x["Gene"] for x in csv.DictReader(open(vpath("haplotype_matrix.csv")))}
_wc={"north":[],"south":[]}
_rr=csv.reader(open(vpath("pi_by_clade_persite.tsv"),encoding="latin-1"),delimiter="\t");_h=next(_rr);_cc=_h[:15]
for _row in _rr:
    if len(_row)<15: continue
    _d=dict(zip(_cc,_row[:15]))
    try:int(_d["north_n_called"])
    except:continue
    _g=_gm[(_d["POS"],_d["REF"],_d["ALT"])]; _c={"synonymous":"syn","nonsynonymous":"ns"}.get(_d["Class"])
    if not _c: continue
    for _cl,_f in (("north","north_pi_site"),("south","south_pi_site")):
        try:_pi=float(_d[_f])
        except:_pi=0.0
        if _pi>0: _wc[_cl].append((_g,_c,_pi))
def _sel(sites,cx): return [(c,pi) for g,c,pi in sites if cx=="Overall" or _G2C.get(g)==cx.split("_")[-1]]
def _R(sub,an):
    sy=[pi for c,pi in sub if c=="syn"]; ns=[pi for c,pi in sub if c=="ns"]
    return ((an*sum(ns)/len(ns))/(an*sum(sy)/len(sy))) if (sy and ns) else None
def _Rci(sub,an):
    v=[]
    for _ in range(_BR):
        bs=[sub[_rnd.randrange(len(sub))] for _ in sub]
        r=_R(bs,an)
        if r is not None: v.append(r)
    v.sort(); return (v[int(.025*len(v))],v[int(.975*len(v))]) if len(v)>10 else (None,None)
_series=[("Between","#e67e22",_wp,_AN["panel"]),
         ("North","#1f77b4",_wc["north"],_AN["north"]),
         ("South","#d62728",_wc["south"],_AN["south"])]
_ratio_rows=[]
xb=np.arange(len(groups)); w=0.8/len(_series)
for k,(lab,color,sites,an) in enumerate(_series):
    pe=[]; er_lo=[]; er_hi=[]
    for g in groups:
        sub=_sel(sites,g); r=_R(sub,an); lo,hi=_Rci(sub,an)
        pe.append(r if r is not None else np.nan)
        er_lo.append((r-lo) if (r is not None and lo is not None) else 0)
        er_hi.append((hi-r) if (r is not None and hi is not None) else 0)
        _ratio_rows.append([lab,g.replace("_"," "),
            "" if r is None else round(r,3),"" if lo is None else round(lo,3),"" if hi is None else round(hi,3)])
    off=(k-(len(_series)-1)/2)*w
    ax[2].bar(xb+off,pe,w,yerr=[er_lo,er_hi],capsize=3,label=lab,color=color,edgecolor="black",linewidth=0.4)
ax[2].axhline(1.0,ls="--",color="grey",lw=0.9)
ax[2].set_xticks(xb); ax[2].set_xticklabels([g.replace("_"," ") for g in groups],rotation=20,ha="right")
ax[2].set_ylabel("π/θ  NS / S"); ax[2].set_title("C. π/θ ratio  NS / S  (1 = neutral)",fontsize=11)
ax[2].legend(fontsize=8)
with open(opath("pi_theta_NS_over_S_ratio.tsv"),"w") as _o:
    _o.write("Comparison\tGroup\tNS_over_S\tCI_low\tCI_high\n")
    for row in _ratio_rows: _o.write("\t".join(str(x) for x in row)+"\n")

plt.savefig(opath("figure_pitheta_within_between.png"),dpi=200,bbox_inches="tight")
plt.savefig(opath("figure_pitheta_within_between.pdf"),bbox_inches="tight")
plt.savefig(opath("figure_pitheta_within_between.svg"),bbox_inches="tight")
plt.close()

# =====================================================================
# Workbook: pN/pS Table 1 as a formatted sheet
# =====================================================================
wb=Workbook(); ws=wb.active; ws.title="pN_pS_Table1"
hf=Font(name="Arial",bold=True,color="FFFFFF",size=10); fill=PatternFill("solid",fgColor="2F5496")
base=Font(name="Arial",size=10); thin=Side(style="thin",color="D9D9D9"); bd=Border(left=thin,right=thin,top=thin,bottom=thin); ctr=Alignment(horizontal="center")
ws["A1"]="Table 1. pN, pS, pN/pS for 13 mitochondrial proteins (95% bootstrap CI)"; ws["A1"].font=Font(name="Arial",bold=True,size=12)
ws.cell(2,1,CAPTION).font=Font(name="Arial",italic=True,size=9,color="444444"); ws.merge_cells("A2:L2")
r0=4
for j,h in enumerate(cols,1):
    c=ws.cell(r0,j,h); c.font=hf; c.fill=fill; c.border=bd; c.alignment=ctr
for i,row in enumerate(rows_tbl,r0+1):
    for j,v in enumerate(row,1):
        c=ws.cell(i,j,v); c.border=bd
        c.font=Font(name="Arial",bold=True,size=10) if row[0]=="Overall" else base
        if j>=3: c.alignment=ctr
        if isinstance(v,float): c.number_format="0.0000" if j>=8 else "0.0"
for j,h in enumerate(cols,1):
    ws.column_dimensions[get_column_letter(j)].width=max(8,len(h)+2)
wb.save(opath("pN_pS_Table1.xlsx"))

print("Wrote: pN_pS_Table1.txt, pN_pS_Table1.xlsx,")
print("       figure_pnps_gene_complex.png/.pdf, figure_pitheta_within_between.png/.pdf")

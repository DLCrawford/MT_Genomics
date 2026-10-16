#!/usr/bin/env python3
"""
10_pitheta_greycolumn.py
------------------------
One panel, 4 complexes + Overall. Per group, three comparisons, each shown as a
TRIPLET:  [ syn bar | grey NS/S ratio column | NS bar ].
  North   = blue pair    (dark/light)   ratio column grey
  South   = red pair      (dark/light)  ratio column grey
  Between = green / burnt-orange         ratio column grey
LEFT axis  : pi/theta (syn & NS bars), 95% CI whiskers.
RIGHT axis : pi/theta NonSyn/Syn ratio -- grey column from 0 to the MEAN ratio,
             with a 95% CI whisker on top; neutral dashed line at 1; capped at 2
             with a break through Complex V's Between-clade CI (to 4.39).
Outputs figure_pitheta_greycolumn.png/.pdf/.svg
"""
import csv
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from common import opath
def R(p): return list(csv.DictReader(open(opath(p)), delimiter="\t"))
def v(x):
    try: return float(x)
    except: return np.nan
WC=R("pi_theta_withinclade_by_complex.tsv"); WP=R("pi_theta_wholepanel_by_complex.tsv"); RT=R("pi_theta_NS_over_S_ratio.tsv")
groups=["Complex_I","Complex_III","Complex_IV","Complex_V","Overall"]; labels=["Complex I","Complex III","Complex IV","Complex V","Overall"]
def wc(cl,g,c): return next(r for r in WC if r["Clade"]==cl and r["Group"]==g and r["Class"]==c)
def wp(g,c):    return next(r for r in WP if r["Group"]==g and r["Class"]==c)
def rt(cmp,g):  return next(r for r in RT if r["Comparison"]==cmp and r["Group"]==g.replace("_"," "))

# comparison -> (syn source, NS source, syn color, NS color, ratio-Comparison-key, triplet-center-offset)
COMPS=[
 ("North", lambda g: wc("North",g,"synonymous"), lambda g: wc("North",g,"nonsynonymous"), "#1f4e9c","#9ecae1","North", -0.30),
 ("South", lambda g: wc("South",g,"synonymous"), lambda g: wc("South",g,"nonsynonymous"), "#b02418","#ff9d97","South",  0.00),
 ("Between",lambda g: wp(g,"synonymous"),        lambda g: wp(g,"nonsynonymous"),          "#27ae60","#d9701f","Between",0.30),
]
x=np.arange(len(groups))
fig,axL=plt.subplots(figsize=(15,6.2)); axR=axL.twinx()
bw=0.092; dsub=0.098; CAP=2.0; grey="#8c8c8c"

for name,synsrc,nssrc,csyn,cns,rkey,coff in COMPS:
    # syn bar (left)  at triplet-left
    for src,color,sub in ((synsrc,csyn,-dsub),(nssrc,cns,+dsub)):
        pe=[v(src(g)["pi_theta"]) for g in groups]; lo=[v(src(g)["CI_low"]) for g in groups]; hi=[v(src(g)["CI_high"]) for g in groups]
        err=[np.array(pe)-np.array(lo),np.array(hi)-np.array(pe)]
        axL.bar(x+coff+sub,pe,bw,yerr=err,capsize=2,color=color,edgecolor="black",linewidth=0.4,error_kw=dict(lw=0.7),zorder=2)
    # grey ratio column (right) in the middle of the triplet: 0 -> mean, whisker = CI
    for i,g in enumerate(groups):
        d=rt(rkey,g); r=v(d["NS_over_S"]); lo=v(d["CI_low"]); hi=v(d["CI_high"])
        xx=x[i]+coff
        axR.bar(xx,min(r,CAP),bw,color=grey,edgecolor="black",linewidth=0.4,zorder=1)
        if hi<=CAP:
            axR.errorbar(xx,r,yerr=[[r-lo],[hi-r]],fmt="none",ecolor="black",capsize=2.5,elinewidth=0.8,zorder=3)
        else:  # Complex V Between: whisker broken at cap
            axR.errorbar(xx,r,yerr=[[r-lo],[0]],fmt="none",ecolor="black",capsize=2.5,elinewidth=0.8,zorder=3)
            axR.plot([xx,xx],[r,1.74],color="black",lw=0.8,zorder=3)
            for dy in (1.74,1.82):
                axR.plot([xx-0.05,xx+0.05],[dy-0.03,dy+0.04],color="black",lw=1.0,zorder=4)
            axR.annotate("CI %.1f–%.2f"%(lo,hi),xy=(xx,1.88),ha="center",fontsize=7,zorder=4)

axL.set_ylim(0,1.0); axL.set_ylabel("π / θ   (syn & NS bars)",fontsize=11)
axR.set_ylim(0,CAP); axR.set_ylabel("π/θ  NonSyn / Syn   (grey columns; 1 = neutral)",fontsize=11,color=grey)
axR.tick_params(axis="y",colors=grey)
axR.axhline(1.0,ls="--",color=grey,lw=0.9,zorder=0)
axL.set_xticks(x); axL.set_xticklabels(labels,rotation=15,ha="right")

# legend
from matplotlib.patches import Patch
handles=[Patch(fc="#1f4e9c",ec="black",label="North syn"),Patch(fc="#9ecae1",ec="black",label="North NS"),
         Patch(fc="#b02418",ec="black",label="South syn"),Patch(fc="#ff9d97",ec="black",label="South NS"),
         Patch(fc="#27ae60",ec="black",label="Between syn"),Patch(fc="#d9701f",ec="black",label="Between NS"),
         Patch(fc=grey,ec="black",label="NS/S ratio (0→mean, whisker=95% CI)")]
axL.legend(handles=handles,loc="lower center",bbox_to_anchor=(0.5,1.01),ncol=4,fontsize=8,frameon=False,columnspacing=1.2)
fig.suptitle("π/θ by clade with NonSyn/Syn ratio columns — 4 complexes + Overall (95% CI)",fontsize=12,y=1.11)
axL.set_zorder(axR.get_zorder()+1); axL.patch.set_visible(False)
plt.tight_layout()
for ext in ("png","pdf","svg"):
    plt.savefig(opath("figure_pitheta_greycolumn."+ext),dpi=200 if ext=="png" else None,bbox_inches="tight")
print("wrote figure_pitheta_greycolumn.png/.pdf/.svg")

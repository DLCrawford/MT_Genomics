#!/usr/bin/env python3
"""
09_pitheta_merged.py
--------------------
All-in-one panel for the 4 complexes + Overall:
  LEFT axis  : pi/theta bars, 95% CI
      North  syn / NS  = dark blue / light blue
      South  syn / NS  = dark red  / light red
      Between syn / NS = green     / burnt orange
  RIGHT axis : pi/theta NonSyn/Syn ratio (grey), 95% CI, neutral line at 1
      North = circle, South = square, Between = diamond
      capped at 2 with a break through Complex V's Between-clade CI (to 4.39)
Reads pi_theta_withinclade_by_complex.tsv, pi_theta_wholepanel_by_complex.tsv,
pi_theta_NS_over_S_ratio.tsv.  Outputs figure_pitheta_merged.png/.pdf/.svg
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

x=np.arange(len(groups))
fig,axL=plt.subplots(figsize=(13,6)); axR=axL.twinx()

# ----- LEFT: 6 pi/theta bars per group -----
# (label, color, offset, source-callable returning row)
BARS=[
 ("North syn", "#1f4e9c", -0.375, lambda g: wc("North",g,"synonymous")),
 ("North NS",  "#9ecae1", -0.225, lambda g: wc("North",g,"nonsynonymous")),
 ("South syn", "#b02418", -0.075, lambda g: wc("South",g,"synonymous")),
 ("South NS",  "#ff9d97",  0.075, lambda g: wc("South",g,"nonsynonymous")),
 ("Between syn","#27ae60",  0.225, lambda g: wp(g,"synonymous")),
 ("Between NS", "#d9701f",  0.375, lambda g: wp(g,"nonsynonymous")),
]
w=0.14
for lab,color,off,src in BARS:
    pe=[v(src(g)["pi_theta"]) for g in groups]; lo=[v(src(g)["CI_low"]) for g in groups]; hi=[v(src(g)["CI_high"]) for g in groups]
    err=[np.array(pe)-np.array(lo),np.array(hi)-np.array(pe)]
    axL.bar(x+off,pe,w,yerr=err,capsize=2,color=color,edgecolor="black",linewidth=0.4,label=lab,error_kw=dict(lw=0.7),zorder=2)
axL.set_ylim(0,1.0); axL.set_ylabel("π / θ   (per site class)",fontsize=11)
axL.set_xticks(x); axL.set_xticklabels(labels,rotation=15,ha="right")

# ----- RIGHT: NS/S ratio, grey, 3 comparisons -----
CAP=2.0
RAT=[("North","o",-0.30),("South","s",0.0),("Between","D",0.30)]
grey="#5a5a5a"
for cmp,mk,off in RAT:
    for i,g in enumerate(groups):
        d=rt(cmp,g); r=v(d["NS_over_S"]); lo=v(d["CI_low"]); hi=v(d["CI_high"])
        xx=x[i]+off
        if hi<=CAP:
            axR.errorbar(xx,r,yerr=[[r-lo],[hi-r]],fmt=mk,ms=6,mfc=grey,mec="black",mew=0.5,
                         color=grey,ecolor=grey,capsize=3,elinewidth=1.0,zorder=5)
        else:  # Complex V Between: break at cap
            axR.errorbar(xx,r,yerr=[[r-lo],[0]],fmt=mk,ms=6,mfc=grey,mec="black",mew=0.5,
                         color=grey,ecolor=grey,capsize=3,elinewidth=1.0,zorder=5)
            axR.plot([xx,xx],[r,1.72],color=grey,lw=1.0,zorder=5)
            for dy in (1.72,1.80):
                axR.plot([xx-0.07,xx+0.07],[dy-0.03,dy+0.04],color=grey,lw=1.1,zorder=6)
            axR.annotate("CI %.1f–%.2f"%(lo,hi),xy=(xx,1.84),ha="center",fontsize=7,color=grey,zorder=6)
axR.axhline(1.0,ls="--",color=grey,lw=0.9,zorder=1)
axR.set_ylim(0,CAP); axR.set_ylabel("π/θ  NonSyn / Syn   (grey; 1 = neutral)",fontsize=11,color=grey)
axR.tick_params(axis="y",colors=grey)

# ----- legends -----
from matplotlib.lines import Line2D
h1,l1=axL.get_legend_handles_labels()
ratio_handles=[Line2D([],[],marker=mk,color=grey,mfc=grey,mec="black",ls="none",ms=6,label="NS/S "+cmp) for cmp,mk,off in RAT]
leg1=axL.legend(h1,l1,loc="lower center",bbox_to_anchor=(0.5,1.10),ncol=6,fontsize=8,frameon=False,columnspacing=1.0,handlelength=1.2)
axL.add_artist(leg1)
axR.legend(handles=ratio_handles,loc="lower center",bbox_to_anchor=(0.5,1.005),ncol=3,fontsize=8,frameon=False)
fig.suptitle("π/θ by clade and NonSyn/Syn ratio — 4 complexes + Overall (95% CI)",fontsize=12,y=1.16)
axL.set_zorder(axR.get_zorder()+1); axL.patch.set_visible(False)
plt.tight_layout()
for ext in ("png","pdf","svg"):
    plt.savefig(opath("figure_pitheta_merged."+ext),dpi=200 if ext=="png" else None,bbox_inches="tight")
print("wrote figure_pitheta_merged.png/.pdf/.svg")

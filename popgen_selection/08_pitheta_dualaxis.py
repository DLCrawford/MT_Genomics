#!/usr/bin/env python3
"""
08_pitheta_dualaxis.py
----------------------
Single-panel, dual-axis summary of the between-clade (whole-panel) result for
the 4 OXPHOS complexes + Overall:
  right y-axis : pi/theta for synonymous and nonsynonymous sites (bars, 95% CI)
  left  y-axis : pi/theta NonSyn/Syn ratio (diamonds, 95% CI), capped at 2 with a
                 break through Complex V's CI (true upper = 4.39), neutral line at 1.
Reads pi_theta_wholepanel_by_complex.tsv and pi_theta_NS_over_S_ratio.tsv.
Outputs figure_pitheta_dualaxis.png/.pdf/.svg
"""
import csv
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from common import opath

def rows(p): return list(csv.DictReader(open(opath(p)), delimiter="\t"))
def v(x):
    try: return float(x)
    except: return np.nan

WP = rows("pi_theta_wholepanel_by_complex.tsv")
RT = [r for r in rows("pi_theta_NS_over_S_ratio.tsv") if r["Comparison"] == "Between"]
groups = ["Complex_I","Complex_III","Complex_IV","Complex_V","Overall"]
labels = ["Complex I","Complex III","Complex IV","Complex V","Overall"]
def wp(g,cls): return next(r for r in WP if r["Group"]==g and r["Class"]==cls)
def rt(g):     return next(r for r in RT if r["Group"]==g.replace("_"," "))

x = np.arange(len(groups))
fig, axR = plt.subplots(figsize=(9,5.5))      # axR = right axis (pi/theta)
axL = axR.twinx()                              # axL = left axis (ratio) -- created 2nd, move to left
axL.spines["left"].set_position(("axes",0)); axL.yaxis.set_label_position("left"); axL.yaxis.set_ticks_position("left")
axR.spines["right"].set_position(("axes",1)); axR.yaxis.set_label_position("right"); axR.yaxis.set_ticks_position("right")

# ---- right axis: pi/theta syn & NS bars ----
w=0.26
for cls,color,off,lab in [("synonymous","#27ae60",-0.29,"π/θ synonymous"),
                          ("nonsynonymous","#e67e22",-0.02,"π/θ nonsynonymous")]:
    pe=[v(wp(g,cls)["pi_theta"]) for g in groups]
    lo=[v(wp(g,cls)["CI_low"]) for g in groups]; hi=[v(wp(g,cls)["CI_high"]) for g in groups]
    err=[np.array(pe)-np.array(lo),np.array(hi)-np.array(pe)]
    axR.bar(x+off,pe,w,yerr=err,capsize=3,color=color,edgecolor="black",linewidth=0.5,label=lab,zorder=2)
axR.set_ylim(0,1.0); axR.set_ylabel("π / θ  (synonymous, nonsynonymous)",fontsize=11)

# ---- left axis: NS/S ratio diamonds, capped at 2 with a break for Complex V ----
CAP=2.0
xr = x+0.30
for i,g in enumerate(groups):
    r=v(rt(g)["NS_over_S"]); lo=v(rt(g)["CI_low"]); hi=v(rt(g)["CI_high"])
    if hi<=CAP:
        axL.errorbar(xr[i],r,yerr=[[r-lo],[hi-r]],fmt="D",ms=7,color="black",
                     ecolor="black",capsize=4,elinewidth=1.3,zorder=5)
    else:  # Complex V: break the upper CI at the cap
        axL.errorbar(xr[i],r,yerr=[[r-lo],[0]],fmt="D",ms=7,color="black",
                     ecolor="black",capsize=4,elinewidth=1.3,zorder=5)
        axL.plot([xr[i],xr[i]],[r,1.72],color="black",lw=1.3,zorder=5)          # stem up to break
        for dy in (1.72,1.80):                                                   # break slashes (//)
            axL.plot([xr[i]-0.09,xr[i]+0.09],[dy-0.03,dy+0.04],color="black",lw=1.3,zorder=6)
        axL.annotate("CI %.1f–%.2f"%(lo,hi),xy=(xr[i],1.80),xytext=(xr[i]-0.02,1.90),
                     ha="center",fontsize=8,zorder=6)
axL.axhline(1.0,ls="--",color="grey",lw=1.0,zorder=1)                            # neutral null
axL.set_ylim(0,CAP); axL.set_ylabel("π/θ  NonSyn / Syn  (ratio; 1 = neutral)",fontsize=11,color="black")
# ratio proxy for legend
axL.scatter([],[],marker="D",color="black",label="π/θ NS/S ratio")

axR.set_xticks(x); axR.set_xticklabels(labels,rotation=15,ha="right")
fig.suptitle("Between-clade π/θ and NonSyn/Syn ratio (95% CI)",fontsize=12,y=1.09)
# combined legend
h1,l1=axR.get_legend_handles_labels(); h2,l2=axL.get_legend_handles_labels()
axR.legend(h1+h2,l1+l2,loc="lower center",bbox_to_anchor=(0.5,1.005),fontsize=9,ncol=3,frameon=False)
axR.set_zorder(axL.get_zorder()+1); axR.patch.set_visible(False)
plt.tight_layout()
for ext in ("png","pdf","svg"):
    plt.savefig(opath("figure_pitheta_dualaxis."+ext),dpi=200 if ext=="png" else None,bbox_inches="tight")
print("wrote figure_pitheta_dualaxis.png/.pdf/.svg")

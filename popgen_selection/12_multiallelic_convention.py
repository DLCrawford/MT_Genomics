#!/usr/bin/env python3
"""
12_multiallelic_convention.py
-----------------------------
Treat multiallelic positions consistently across S, theta_W, and pi, under three
conventions, for the reviewer's sensitivity analysis, and rebuild Table 3 on the
"multiallelic = one site" (collapse) convention.

Conventions (hard calls from haplotype_matrix; stop_retained counted synonymous):
  collapse : each POSITION is one segregating site; multiallelic pi_site uses the
             full multiallelic formula pi = n/(n-1) * (1 - sum_a p_a^2).  [REPORTED]
  split    : each ALT allele is its own biallelic site (previous Table 3).
  drop     : multiallelic positions excluded entirely.
Class of a position = its ALT's class; the one mixed position (2847) is assigned
to its MAJOR allele's class (synonymous).

Outputs: table3_multiallelic_sensitivity.tsv, Table3_collapse.tsv
"""
import csv
from collections import defaultdict
L_CDS=11417.0; L_S=2693.7; L_N=8694.3
def a_n(n): return sum(1.0/i for i in range(1,n))
def cls(vt): return "ns" if "missense" in vt else "syn"

mem={}
for r in csv.DictReader(open("../vcf/pi_by_clade_persite.tsv.membership.tsv"),delimiter="\t"):
    mem[r["Sample"].replace("_MT","")]=r["Clade"]
north=set(s for s,c in mem.items() if c=="north"); south=set(s for s,c in mem.items() if c=="south")

# group rows by POS; keep per-alt genotype vectors
rows=list(csv.DictReader(open("../vcf/haplotype_matrix.csv")))
idcols=[k for k in rows[0] if k.endswith("_MT")]
byPOS=defaultdict(list)
for r in rows:
    byPOS[r["POS"]].append(r)

def sample_stats(members):
    ids=[k for k in idcols if (members is None or k[:-3] in members)]
    n=len(ids); a=a_n(n)
    acc={c:{"S":0,"pi":0.0} for c in ("all","syn","ns")}
    accS={c:{"S":0,"pi":0.0} for c in ("all","syn","ns")}   # split
    accD={c:{"S":0,"pi":0.0} for c in ("all","syn","ns")}   # drop
    for pos,rlist in byPOS.items():
        alts=[]
        for r in rlist:
            AC=sum(1 for k in ids if r[k]=="1")
            alts.append((AC,cls(r["Variant_type"])))
        tot_alt=sum(AC for AC,_ in alts); ref=n-tot_alt
        present=[AC for AC,_ in alts if AC>0]
        # polymorphic in this sample?  (>1 distinct allele present)
        distinct = (1 if ref>0 else 0) + len(present)
        poly = distinct>1
        ismulti = len(present)>1
        # ---- collapse (position) ----
        if poly:
            freqs=[ref]+[AC for AC,_ in alts if AC>0]
            pis = 1.0 - sum((f/n)**2 for f in freqs); pis *= n/(n-1)
            # class: major alt's class
            maj = max([a for a in alts if a[0]>0], key=lambda x:x[0])
            c = maj[1]
            acc["all"]["S"]+=1; acc["all"]["pi"]+=pis
            acc[c]["S"]+=1;     acc[c]["pi"]+=pis
        # ---- split (each alt) ----
        for AC,c in alts:
            if 0<AC<n:
                pis=2.0*AC*(n-AC)/(n*(n-1))
                accS["all"]["S"]+=1; accS["all"]["pi"]+=pis
                accS[c]["S"]+=1;     accS[c]["pi"]+=pis
        # ---- drop multiallelic ----
        if not ismulti and poly:
            AC,c=[a for a in alts if a[0]>0][0]
            if 0<AC<n:
                pis=2.0*AC*(n-AC)/(n*(n-1))
                accD["all"]["S"]+=1; accD["all"]["pi"]+=pis
                accD[c]["S"]+=1;     accD[c]["pi"]+=pis
    def derive(acc):
        d={}
        S=acc["all"]["S"]; Ss=acc["syn"]["S"]; Sn=acc["ns"]["S"]
        d["S"]=S; d["Ssyn"]=Ss; d["Sns"]=Sn
        d["thT"]=S/a/L_CDS; d["thS"]=Ss/a/L_S; d["thN"]=Sn/a/L_N
        d["piT"]=acc["all"]["pi"]/L_CDS; d["piS"]=acc["syn"]["pi"]/L_S; d["piN"]=acc["ns"]["pi"]/L_N
        d["ptT"]=a*acc["all"]["pi"]/S if S else None
        d["ptS"]=a*acc["syn"]["pi"]/Ss if Ss else None
        d["ptN"]=a*acc["ns"]["pi"]/Sn if Sn else None
        d["pnps_th"]=d["thN"]/d["thS"]; d["pnps_pi"]=d["piN"]/d["piS"]
        d["nsS"]=d["ptN"]/d["ptS"]
        return d
    return n,derive(acc),derive(accS),derive(accD)

groups=[("Fhet_all",None),("North",north),("South",south)]
res={}
for name,mem_ in groups:
    n,col,spl,drp=sample_stats(mem_)
    res[name]=dict(n=n,collapse=col,split=spl,drop=drp)

# ---- sensitivity table (key totals) ----
with open("table3_multiallelic_sensitivity.tsv","w") as o:
    o.write("Group\tConvention\tS_total\tS_syn\tS_ns\ttheta_W_total\tpi_total\tpi_theta_overall\tpi_theta_syn\tpi_theta_ns\tpN_pS(theta)\tpN_pS(pi)\tNS/S(pi_theta)\n")
    for name,_ in groups:
        for conv in ("collapse","split","drop"):
            d=res[name][conv]
            o.write("%s\t%s\t%d\t%d\t%d\t%.4f\t%.4f\t%.3f\t%.3f\t%.3f\t%.4f\t%.4f\t%.3f\n"%(
                name,conv,d["S"],d["Ssyn"],d["Sns"],d["thT"],d["piT"],
                d["ptT"],d["ptS"],d["ptN"],d["pnps_th"],d["pnps_pi"],d["nsS"]))

# ---- console: sensitivity summary ----
print("=== SENSITIVITY (17 multiallelic positions; 1 mixed = POS2847) ===")
print("%-9s %-9s %5s %5s %5s | %-6s %-6s %-6s | %-6s %-6s"%("Group","Conv","S","Ssyn","Sns","π/θ_all","π/θsyn","π/θns","pNpSθ","NS/S"))
for name,_ in groups:
    for conv in ("collapse","split","drop"):
        d=res[name][conv]
        print("%-9s %-9s %5d %5d %5d | %6.3f %6.3f %6.3f | %6.4f %6.3f"%(
            name,conv,d["S"],d["Ssyn"],d["Sns"],d["ptT"],d["ptS"],d["ptN"],d["pnps_th"],d["nsS"]))
    print()

# ================= full Table 3 (collapse) + xlsx + sensitivity xlsx =================
def SN(b,c):
    return (c/b) if (b not in (0,None) and c not in (None,)) else None
def emit_table3():
    B=res["Fhet_all"]["collapse"]; No=res["North"]["collapse"]; So=res["South"]["collapse"]
    order=[("N (chromosomes)", lambda d:None, "n"),
           ("S total","S"),("S syn","Ssyn"),("S nonsyn","Sns"),
           ("θ_W total","thT"),("θ_W syn","thS"),("θ_W ns","thN"),("pN/pS (θ_W)","pnps_th"),
           ("π total","piT"),("π syn","piS"),("π ns","piN"),("pN/pS (π)","pnps_pi"),
           ("π/θ overall","ptT"),("π/θ syn","ptS"),("π/θ ns","ptN"),("NS/S (π/θ)","nsS")]
    Ns={"Fhet_all":141,"North":77,"South":63}
    lines=[]
    hdr=["Row","Fhet_all","North","South","South/North"]
    lines.append("\t".join(hdr))
    def val(d,k):
        if k=="n": return None
        return d[k]
    lines.append("\t".join(["N (chromosomes)","141","77","63","0.818"]))
    for label,key in [(x[0],x[1]) for x in order[1:]]:
        a=B[key]; b=No[key]; c=So[key]; r=SN(b,c)
        def f(x):
            if x is None: return ""
            if key in ("S","Ssyn","Sns"): return str(int(x))
            return "%.4f"%x if key.startswith("th") or key.startswith("pi") and key!="ptT" else "%.3f"%x
        # tidy formatting: counts int; theta/pi 4dp; ratios 3dp
        def ff(x):
            if x is None: return ""
            if key in ("S","Ssyn","Sns"): return str(int(round(x)))
            if key in ("thT","thS","thN","piT","piS","piN"): return "%.4f"%x
            return "%.3f"%x
        lines.append("\t".join([label,ff(a),ff(b),ff(c),("%.3f"%r if r is not None else "")]))
    return lines
lines=emit_table3()
with open("Table3_collapse.tsv","w") as o:
    o.write("Table 3. Segregating sites (S), Watterson's theta, nucleotide diversity (pi) and pi/theta "
            "among and within clades, multiallelic positions counted ONCE. North n=77, South n=63 "
            "(sites polymorphic within clade); Fhet_all = whole panel n=141. South/North = ratio of the "
            "two clade columns. L_CDS=11,417; L_Syn=2,694; L_NonSyn=8,694 bp. 17 multiallelic positions "
            "(1 mixed syn/NS, POS2847, assigned to its major allele=synonymous); stop_retained counted synonymous.\n")
    o.write("\n".join(lines)+"\n")

# xlsx for both
from openpyxl import Workbook
from openpyxl.styles import Font,PatternFill,Alignment,Border,Side
from openpyxl.utils import get_column_letter
wb=Workbook()
def style_sheet(ws,title,caption,rows2d,ncol):
    hf=Font(name="Arial",bold=True,color="FFFFFF",size=10); fill=PatternFill("solid",fgColor="404040")
    ws["A1"]=title; ws["A1"].font=Font(name="Arial",bold=True,size=12)
    ws.cell(2,1,caption).font=Font(name="Arial",italic=True,size=9,color="444444"); ws.merge_cells(start_row=2,start_column=1,end_row=2,end_column=ncol)
    thin=Side(style="thin",color="D9D9D9"); bd=Border(left=thin,right=thin,top=thin,bottom=thin); ctr=Alignment(horizontal="center")
    r0=4
    for i,row in enumerate(rows2d):
        for j,cellv in enumerate(row,1):
            c=ws.cell(r0+i,j,cellv); c.border=bd
            if i==0: c.font=hf; c.fill=fill; c.alignment=ctr
            else:
                c.font=Font(name="Arial",size=10)
                if j>=2: c.alignment=ctr
                try:
                    if "." in str(cellv): c.value=float(cellv)
                except: pass
    for j in range(1,ncol+1): ws.column_dimensions[get_column_letter(j)].width=13 if j==1 else 11
    ws.column_dimensions["A"].width=16
ws1=wb.active; ws1.title="Table3_collapse"
t3=[l.split("\t") for l in lines]
style_sheet(ws1,"Table 3 (multiallelic = one site)","S, θ_W, π, π/θ among & within clades; multiallelic positions counted once.",t3,5)
ws2=wb.create_sheet("multiallelic_sensitivity")
sens=[l.split("\t") for l in open("table3_multiallelic_sensitivity.tsv").read().strip().split("\n")]
style_sheet(ws2,"Multiallelic sensitivity (17 positions)","collapse = one site (reported); split = each ALT a biallelic site; drop = exclude multiallelic. Results stable across all three.",sens,13)
wb.save("Table3_collapse.xlsx")
print("\nwrote Table3_collapse.tsv/.xlsx and table3_multiallelic_sensitivity.tsv")
print("\n=== Table 3 (collapse) ===")
for l in lines: print(l)

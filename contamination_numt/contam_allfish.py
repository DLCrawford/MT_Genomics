import sys, math, statistics as st
from collections import Counter
ad_tsv, samples_txt, focals_txt, outdir = sys.argv[1:5]
groups_txt = sys.argv[5] if len(sys.argv)>5 else ""   # sample<TAB>hp|no_hp
samples=[l.strip() for l in open(samples_txt) if l.strip()]
focal=set(l.strip() for l in open(focals_txt) if l.strip())
grp={}
if groups_txt:
    for l in open(groups_txt):
        l=l.rstrip("\n")
        if not l or l.startswith("#"): continue
        p=l.split("\t") if "\t" in l else l.split()
        if len(p)>=2: grp[p[0]]=p[1]
MINDP,MAXDP=8,60; LO,HI=0.10,0.35
ns=[0]*len(samples); nm=[0]*len(samples); nt=[0]*len(samples)
def flush(rows):
    for si in range(len(samples)):
        refd=None; alts=[]; ok=True
        for r in rows:
            ad=r[si]
            if ad in (".",""): ok=False; break
            p=ad.split(",")
            try: rd=int(p[0])
            except: ok=False; break
            refd=rd if refd is None else refd
            if len(rows)>1: alts.append(int(p[1]) if len(p)>1 else 0)
            else: alts.extend(int(x) for x in p[1:])
        if not ok or refd is None: continue
        d=sorted([refd]+alts,reverse=True); dp=sum(d)
        if dp<MINDP or dp>MAXDP: continue
        ns[si]+=1
        if any(LO<=x/dp<HI for x in d[1:]): nm[si]+=1
        if any(LO<=d[r]/dp<HI for r in range(2,len(d))): nt[si]+=1
cur=None; buf=[]
for line in open(ad_tsv):
    c=line.rstrip("\n").split("\t"); pos=(c[0],c[1])
    if pos!=cur:
        if buf: flush(buf)
        buf=[]; cur=pos
    buf.append(c[4:])
if buf: flush(buf)
S=[]
for i,s in enumerate(samples):
    n=max(ns[i],1)
    S.append(dict(s=s,foc=s in focal,g=grp.get(s,"other"),n=ns[i],minor=nm[i]/n,third=nt[i]/n))
# baseline = no_hp group (clean reference). fall back to all non-focal if no groups.
base=[r for r in S if r["g"]=="no_hp"] or [r for r in S if not r["foc"]]
def zc(key):
    mu=st.mean(r[key] for r in base); sd=st.pstdev([r[key] for r in base]) or 1e-9
    return mu,sd
mm,sm=zc("minor"); mt,stt=zc("third")
for r in S: r["zminor"]=(r["minor"]-mm)/sm; r["zthird"]=(r["third"]-mt)/stt
def mwu(x,y):
    n1,n2=len(x),len(y)
    if n1==0 or n2==0: return float('nan')
    comb=sorted([(v,0) for v in x]+[(v,1) for v in y]); rk=[0]*len(comb); i=0
    while i<len(comb):
        j=i
        while j+1<len(comb) and comb[j+1][0]==comb[i][0]: j+=1
        rr=(i+j)/2+1
        for k in range(i,j+1): rk[k]=rr
        i=j+1
    R1=sum(rk[k] for k in range(len(comb)) if comb[k][1]==0)
    U1=R1-n1*(n1+1)/2; mu=n1*n2/2
    tie=Counter(v for v,_ in comb); T=sum(t**3-t for t in tie.values()); n=n1+n2
    sd=((n1*n2/12)*((n+1)-T/(n*(n-1))))**0.5 if n>1 else 0
    if sd==0: return 1.0
    z=(U1-mu)/sd; return 2*(1-0.5*(1+math.erf(abs(z)/2**0.5)))
S.sort(key=lambda r:-r["zthird"])
with open(f"{outdir}/per_sample_contam_allfish.tsv","w") as o:
    o.write("sample\tgroup\tfocal\tn_sites\trate_minor\tz_minor\trate_third\tz_third\n")
    for r in S:
        o.write(f'{r["s"]}\t{r["g"]}\t{int(r["foc"])}\t{r["n"]}\t{r["minor"]:.5f}\t{r["zminor"]:.1f}\t{r["third"]:.6f}\t{r["zthird"]:.1f}\n')
def grpvals(pred,key): return [r[key] for r in S if pred(r)]
hp=lambda r:r["g"]=="hp"; nohp=lambda r:r["g"]=="no_hp"; hpxf=lambda r:r["g"]=="hp" and not r["foc"]
with open(f"{outdir}/verdict.txt","w") as o:
    o.write("Nuclear contamination vs mitochondrial-heteroplasmy status (heteroplasmy-style metric, 0.10<=AD/DP<0.35)\n")
    o.write(f"metric per fish: rate_third (allele ranked>=3, unambiguous) and rate_minor (any non-major allele)\n")
    o.write(f"baseline for z = no_hp fish (n={len(base)})\n\n")
    o.write("GROUP                     n   median_rate_third  median_z_third   median_rate_minor median_z_minor\n")
    for name,pred in [("no_hp (no mt hp)",nohp),("hp (mt heteroplasmic)",hp),("hp excl. 4 focals",hpxf),("focal (4 high-load)",lambda r:r["foc"])]:
        v=[r for r in S if pred(r)]
        if not v: continue
        o.write(f"  {name:24s} {len(v):3d}   {st.median([r['third'] for r in v]):.6f}         {st.median([r['zthird'] for r in v]):5.1f}          {st.median([r['minor'] for r in v]):.5f}        {st.median([r['zminor'] for r in v]):5.1f}\n")
    o.write("\nRank-sum test (Mann-Whitney, two-sided) on rate_third:\n")
    o.write(f"  hp vs no_hp:            p = {mwu(grpvals(hp,'third'),grpvals(nohp,'third')):.3g}\n")
    o.write(f"  hp(excl focals) vs no_hp: p = {mwu(grpvals(hpxf,'third'),grpvals(nohp,'third')):.3g}\n")
    o.write("  (non-significant after excluding focals => mt heteroplasmy is NOT nuclear contamination)\n")
    # who is nuclear-contaminated
    cont=[r for r in S if r["zthird"]>=5]
    o.write(f"\nNuclear-contaminated fish (z_third >= 5): n={len(cont)}\n")
    for r in cont:
        o.write(f'  {r["s"]:8s} group={r["g"]:7s} focal={int(r["foc"])} z_third={r["zthird"]:.1f} z_minor={r["zminor"]:.1f}\n')
    inhp=sum(1 for r in cont if r["g"]=="hp"); innohp=sum(1 for r in cont if r["g"]=="no_hp")
    o.write(f"\n  of the contaminated: {inhp} are hp, {innohp} are no_hp; focals among them: {sum(1 for r in cont if r['foc'])}\n")
    o.write("  => if the ONLY contaminated hp fish are the known focals, the other ~49 hp fish are genuine heteroplasmy.\n")
print(open(f"{outdir}/verdict.txt").read())

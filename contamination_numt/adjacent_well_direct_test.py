#!/usr/bin/env python3
"""
Direct adjacent-well contamination test (contamination test #1).

Question: did the individual(s) in a physically ADJACENT well on the sequencing
plate contribute the minor (heteroplasmic) variants seen in an Hp individual?

For each Hp fish X:
  - adjacent wells = other panel members on the SAME plate within Chebyshev
    distance 1 (the <=8 surrounding wells). Donor pool is ALL such neighbours,
    Hp or not (scored on their MAJOR/consensus haplotype).
  - raw sharing      : fraction of X's minor variants whose base == neighbour major.
  - enrichment       : mean raw sharing for adjacent vs non-adjacent same-plate wells
                       (common clade variants inflate both; only an EXCESS in adjacent
                        wells would indicate well-to-well bleed).
  - rare-minor test  : X's minor variants that are carried as MAJOR by <=5 panel members
                       are the donor-identifying ("fingerprint") variants; a real
                       adjacent donor must supply these.
  - whole-haplotype  : is any adjacent neighbour's consensus carrying ALL of X's minors
                       with 0 contradicting sites (a near-clone in an adjacent well)?

Inputs : data_files_May/WGS_seq_plate.txt
         vcf/heteroplasmy_pileup_events_all.tsv
         vcf/haplotype_matrix.csv
Output : EE_revision_analyses/Hp_adjacent_well_direct_test.tsv
"""
import csv, re, statistics as st
from collections import Counter

PLATE="data_files_May/WGS_seq_plate.txt"
EVENTS="vcf/heteroplasmy_pileup_events_all.tsv"
MATRIX="vcf/haplotype_matrix.csv"
OUT="EE_revision_analyses/Hp_adjacent_well_direct_test.tsv"
RARE=5

coord={}
for r in csv.DictReader(open(PLATE),delimiter="\t"):
    r={k.strip():v for k,v in r.items()}
    wid=str(int(r["WGS_ID"])); lab=r["Plate_label"].strip()
    coord[wid]=(r["i5"].strip(),"ABCDEFGH".index(lab[0])+1,int(lab[1:]),lab)
def to_id(s):
    m=re.match(r'(\d+)',s.split("/")[-1]); return m.group(1) if m else None
M={}
for row in csv.DictReader(open(EVENTS),delimiter="\t"):
    M.setdefault(to_id(row["Individual"]),[]).append((row["POS"],row["Hp_allele"]))
rows=list(csv.DictReader(open(MATRIX)))
inds=[c[:-3] for c in rows[0].keys() if c.endswith("_MT")]
pos_ref={}; ind_alt={}; allpos=[]
for r in rows:
    p=r["POS"]
    if p not in pos_ref: pos_ref[p]=r["REF"]; allpos.append(p)
    for ind in inds:
        if r[ind+"_MT"]=="1": ind_alt[(ind,p)]=r["ALT"]
def major(ind,pos): return ind_alt.get((ind,pos),pos_ref.get(pos))
panel=set(inds)
freq=Counter((p,major(i,p)) for p in allpos for i in inds)

def adj(x):
    p,ro,co,_=coord[x]
    return [y for y in panel if y!=x and y in coord and coord[y][0]==p
            and max(abs(coord[y][1]-ro),abs(coord[y][2]-co))==1]
def nonadj(x):
    p=coord[x][0]; a=set(adj(x))
    return [y for y in panel if y!=x and y in coord and coord[y][0]==p and y not in a]
def extra(X,Y,hp):
    return sum(1 for p in allpos if p not in hp and major(X,p)!=major(Y,p))

hp=[x for x in M if x in coord and x in panel]
out=[["fish","well","plate","n_Hp","n_adj","n_rare_minor",
      "best_adj_raw_shared","frac_adj","frac_nonadj",
      "rare_minors_in_best_adj","adj_carry_all_n","min_extra_sites_adj","adjacent_contamination"]]
rare_tot=rare_adj=clean=0
for x in sorted(hp,key=lambda a:int(a)):
    mx=M[x]; hp_pos=set(p for p,_ in mx); nhp=len(hp_pos)
    A=adj(x); NA=nonadj(x)
    sh=lambda y: sum(1 for p,a in mx if major(y,p)==a)
    fa=st.mean([sh(y)/nhp for y in A]) if A else float('nan')
    fn=st.mean([sh(y)/nhp for y in NA]) if NA else float('nan')
    rare=[(p,a) for p,a in mx if freq[(p,a)]<=RARE]; rare_tot+=len(rare)
    rsh=lambda y: sum(1 for p,a in rare if major(y,p)==a)
    best_adj_rare=max((rsh(y) for y in A),default=0); rare_adj+=best_adj_rare
    full=[y for y in A if all(major(y,p)==a for p,a in mx)]
    mine=min((extra(x,y,hp_pos) for y in full),default=None)
    isclean = bool(full) and mine==0
    if isclean: clean+=1
    best=max(A,key=sh) if A else "-"
    out.append([x,coord[x][3],coord[x][0],nhp,len(A),len(rare),
                sh(best) if A else 0,round(fa,3),round(fn,3),
                best_adj_rare,len(full),mine if mine is not None else "NA",
                "YES" if isclean else "no"])
with open(OUT,"w",newline="") as f:
    csv.writer(f,delimiter="\t").writerows(out)
print(f"Hp fish: {len(hp)}  | rare minors total: {rare_tot}  carried by best adjacent well: {rare_adj}")
print(f"adjacent clean whole-haplotype matches: {clean}")
print(f"wrote {OUT}")

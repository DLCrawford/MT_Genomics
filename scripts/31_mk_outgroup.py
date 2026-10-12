#!/usr/bin/env python3
"""
31_mk_outgroup.py  —  Outgroup McDonald-Kreitman test for Fhet mt CDS.

Divergence (D) is measured against F. catenatus (NC_084079.1) — an external,
non-introgressed outgroup, so the test is reference-INDEPENDENT (unlike an
internal northern reference, whose D column flips with reference choice).

dN/dS of the divergence is corrected for multiple hits with a Kimura-2-parameter
(Ts/Tv-aware) model, because synonymous divergence to catenatus is high
(pS ~0.54, Ts/Tv ~3.2 — not saturated, but raw counts undercount Ds).  This is
the same logic as PAML yn00, implemented here with no external dependency.

Tests run:
  * Whole species  — D = Fhet reference (NC_012312.1) vs catenatus; P = all 141-panel variants.
  * North clade    — D = North-consensus vs catenatus;  P = within-North segregating sites.
  * South clade    — D = South-consensus vs catenatus;  P = within-South segregating sites.
  * Sensitivity    — whole species restricted to low-saturation genes
                     (COX1/COX2/COX3/ATP6/ATP8/CYTB).

Statistics per test: NI, alpha (raw NG and K2P-corrected), omega_a = alpha*dN/dS,
DoS, and a two-sided Fisher exact p (pure-Python).

Output: vcf/mk_catenatus.tsv  + stdout.
"""
import gzip, json, math, os, sys
from itertools import permutations

BASE = os.environ.get("MTG_BASE", os.path.expanduser("~/Projects/MT_Genomics_Cl_Ap2026"))
MTG  = f"{BASE}/MT_Genomics2"

OG_NAME  = os.environ.get("MTG_OG", "catenatus")   # set MTG_OG=grandis to switch
OG_FASTA = f"{MTG}/refs_outgroup/F_{OG_NAME}/ncbi_dataset/data/genome.fna"
OG_REPORT= f"{MTG}/refs_outgroup/F_{OG_NAME}/ncbi_dataset/data/organelle_gene_report.jsonl"
FHET_FA  = f"{MTG}/Missing_Files/SSM_MT_ref/Fhet_MT.fasta"
FHET_GFF = f"{MTG}/Missing_Files/SSM_MT_ref/Fhet_MT.gff"
VCF      = f"{MTG}/vcf/141_MT_variants.vcf.gz"

LOWSAT = {"COX1","COX2","COX3","ATP6","ATP8","CYTB"}

MT2 = {
 "TTT":"F","TTC":"F","TTA":"L","TTG":"L","TCT":"S","TCC":"S","TCA":"S","TCG":"S",
 "TAT":"Y","TAC":"Y","TAA":"*","TAG":"*","TGT":"C","TGC":"C","TGA":"W","TGG":"W",
 "CTT":"L","CTC":"L","CTA":"L","CTG":"L","CCT":"P","CCC":"P","CCA":"P","CCG":"P",
 "CAT":"H","CAC":"H","CAA":"Q","CAG":"Q","CGT":"R","CGC":"R","CGA":"R","CGG":"R",
 "ATT":"I","ATC":"I","ATA":"M","ATG":"M","ACT":"T","ACC":"T","ACA":"T","ACG":"T",
 "AAT":"N","AAC":"N","AAA":"K","AAG":"K","AGT":"S","AGC":"S","AGA":"*","AGG":"*",
 "GTT":"V","GTC":"V","GTA":"V","GTG":"V","GCT":"A","GCC":"A","GCA":"A","GCG":"A",
 "GAT":"D","GAC":"D","GAA":"E","GAG":"E","GGT":"G","GGC":"G","GGA":"G","GGG":"G",
}
COMP=str.maketrans("ACGTN","TGCAN"); BASES="ACGT"; PUR=set("AG")
NS_SET={"missense_variant","stop_gained","stop_lost","start_lost","initiator_codon_variant"}

def norm_gene(g):
    g=g.upper()
    return {"COB":"CYTB","CYB":"CYTB","CO1":"COX1","CO2":"COX2","CO3":"COX3",
            "COI":"COX1","COII":"COX2","COIII":"COX3"}.get(g,g)

def read_fasta(path,want=None):
    seqs={};name=None;buf=[]
    op=gzip.open if path.endswith(".gz") else open
    with op(path,"rt") as fh:
        for line in fh:
            line=line.rstrip()
            if line.startswith(">"):
                if name is not None: seqs[name]="".join(buf).upper()
                name=line[1:].split()[0];buf=[]
            else: buf.append(line)
    if name is not None: seqs[name]="".join(buf).upper()
    return seqs[want] if want else seqs

def translate(s):
    L=len(s)-len(s)%3
    return "".join(MT2.get(s[i:i+3],"X") for i in range(0,L,3))

def gff_cds_coords():
    out=[]
    for line in open(FHET_GFF):
        if line.startswith("#") or not line.strip(): continue
        p=line.rstrip("\n").split("\t")
        if len(p)<9 or p[2]!="CDS": continue
        attrs=dict(kv.split("=",1) for kv in p[8].split(";") if "=" in kv)
        out.append((norm_gene(attrs.get("gene") or attrs.get("Name","")),
                    int(p[3]),int(p[4]),p[6]))
    return out

def cds_from_genome(genome_seq,coords):
    out={}
    for gene,s,e,strand in coords:
        seq=genome_seq[s-1:e]
        if strand=="-": seq=seq.translate(COMP)[::-1]
        out[gene]=seq
    return out

def outgroup_cds():
    seqs=read_fasta(OG_FASTA)
    acc=next((k for k in seqs if k.startswith("NC_")), next(iter(seqs)))
    genome=seqs[acc]; out={}
    for line in open(OG_REPORT):
        if not line.strip(): continue
        r=json.loads(line)
        if r.get("type")!="PROTEIN_CODING": continue
        loc=r["annotations"][0]["genomicLocations"][0]["genomicRange"]
        s,e=int(loc["begin"]),int(loc["end"])
        seq=genome[s-1:e]
        if loc["orientation"]=="minus": seq=seq.translate(COMP)[::-1]
        out[norm_gene(r["symbol"])]=seq
    return out

def nw_align(a,b,match=2,mis=-1,gap=-2):
    n,m=len(a),len(b)
    F=[[0]*(m+1) for _ in range(n+1)]
    for i in range(1,n+1):F[i][0]=i*gap
    for j in range(1,m+1):F[0][j]=j*gap
    for i in range(1,n+1):
        ai=a[i-1]
        for j in range(1,m+1):
            F[i][j]=max(F[i-1][j-1]+(match if ai==b[j-1] else mis),F[i-1][j]+gap,F[i][j-1]+gap)
    i,j=n,m;A=[];B=[]
    while i>0 and j>0:
        s=F[i][j]
        if s==F[i-1][j-1]+(match if a[i-1]==b[j-1] else mis):A.append(a[i-1]);B.append(b[j-1]);i-=1;j-=1
        elif s==F[i-1][j]+gap:A.append(a[i-1]);B.append("-");i-=1
        else:A.append("-");B.append(b[j-1]);j-=1
    while i>0:A.append(a[i-1]);B.append("-");i-=1
    while j>0:A.append("-");B.append(b[j-1]);j-=1
    return "".join(reversed(A)),"".join(reversed(B))

def codon_pairs(cds1,cds2):
    p1,p2=translate(cds1),translate(cds2)
    a1,a2=nw_align(p1,p2)
    c1=[cds1[i:i+3] for i in range(0,len(cds1)-len(cds1)%3,3)]
    c2=[cds2[i:i+3] for i in range(0,len(cds2)-len(cds2)%3,3)]
    out=[];i=j=0
    for x,y in zip(a1,a2):
        k1=c1[i] if x!="-" else None; k2=c2[j] if y!="-" else None
        if x!="-":i+=1
        if y!="-":j+=1
        if k1 and k2 and "N" not in k1 and "N" not in k2 and MT2.get(k1)!="*" and MT2.get(k2)!="*":
            out.append((k1,k2))
    return out

def codon_S(codon):
    ref=MT2.get(codon)
    if ref is None or ref=="*": return None
    s=0.0
    for pos in range(3):
        for b in BASES:
            if b!=codon[pos] and MT2.get(codon[:pos]+b+codon[pos+1:],"?")==ref: s+=1/3.0
    return s

def diff_events(c1,c2):
    """yield (issyn,ists) for each substitution, averaged over NG pathways."""
    pos=[k for k in range(3) if c1[k]!=c2[k]]
    if not pos: return []
    if len(pos)==1:
        k=pos[0]
        return [(MT2[c1]==MT2[c2],(c1[k] in PUR)==(c2[k] in PUR))]
    paths=[]
    for order in permutations(pos):
        cur=c1; ev=[]; ok=True
        for k in order:
            nxt=cur[:k]+c2[k]+cur[k+1:]
            if MT2.get(nxt)=="*": ok=False; break
            ev.append((MT2[cur]==MT2[nxt],(cur[k] in PUR)==(nxt[k] in PUR))); cur=nxt
        if ok: paths.append(ev)
    if not paths: return []
    # average: return weighted events as fractional counts
    nsyn_ts=nsyn_tv=nnon_ts=nnon_tv=0.0
    for ev in paths:
        for issyn,ists in ev:
            if issyn:
                if ists: nsyn_ts+=1
                else: nsyn_tv+=1
            else:
                if ists: nnon_ts+=1
                else: nnon_tv+=1
    w=len(paths)
    return [("frac",nsyn_ts/w,nsyn_tv/w,nnon_ts/w,nnon_tv/w)]

def divergence(cds_a,cds_b,genes):
    """Return dict with S_sites,N_sites,Sd,Nd, and Ts/Tv split for syn & nonsyn."""
    S=N=0.0; syn_ts=syn_tv=non_ts=non_tv=0.0
    for g in genes:
        if g not in cds_a or g not in cds_b: continue
        for c1,c2 in codon_pairs(cds_a[g],cds_b[g]):
            s1=codon_S(c1); s2=codon_S(c2)
            if s1 is None or s2 is None: continue
            S+=(s1+s2)/2; N+=3-(s1+s2)/2
            for ev in diff_events(c1,c2):
                if ev[0]=="frac":
                    _,a,b,cc,d=ev; syn_ts+=a; syn_tv+=b; non_ts+=cc; non_tv+=d
                else:
                    issyn,ists=ev
                    if issyn:
                        if ists: syn_ts+=1
                        else: syn_tv+=1
                    else:
                        if ists: non_ts+=1
                        else: non_tv+=1
    return dict(S=S,N=N,syn_ts=syn_ts,syn_tv=syn_tv,non_ts=non_ts,non_tv=non_tv,
                Sd=syn_ts+syn_tv,Nd=non_ts+non_tv)

def k2p(P,Q):
    """Kimura 2-parameter distance from transition (P) and transversion (Q) proportions."""
    a=1-2*P-Q; b=1-2*Q
    if a<=0 or b<=0: return float("inf")
    return -0.5*math.log(a)-0.25*math.log(b)

def rates(div):
    """uncorrected and K2P-corrected dN, dS."""
    S,N=div["S"],div["N"]
    pS=div["Sd"]/S; pN=div["Nd"]/N
    dS_k=k2p(div["syn_ts"]/S,div["syn_tv"]/S)
    dN_k=k2p(div["non_ts"]/N,div["non_tv"]/N)
    return dict(pS=pS,pN=pN,dS=dS_k,dN=dN_k,
                dNdS_raw=(div["Nd"]/N)/(div["Sd"]/S),
                dNdS_k=(dN_k/dS_k) if dS_k not in (0,float("inf")) else float("nan"))

# ── Fisher exact (two-sided), pure python ─────────────────────────────────────
def fisher(a,b,c,d):
    from math import comb
    r1,r2,c1,c2=a+b,c+d,a+c,b+d; n=a+b+c+d
    def p(x):  # P(a=x) hypergeometric
        return comb(c1,x)*comb(c2,r1-x)/comb(n,r1)
    p0=p(a)*(1+1e-9)
    lo=max(0,r1-c2); hi=min(r1,c1)
    return sum(p(x) for x in range(lo,hi+1) if p(x)<=p0)

def mk(label,Pn,Ps,div,rt):
    Sd,Nd=div["Sd"],div["Nd"]
    # corrected divergence substitution numbers
    Ds_c=rt["dS"]*div["S"]; Dn_c=rt["dN"]*div["N"]
    PnPs=Pn/Ps
    NI_raw=PnPs/((Nd/Sd)); a_raw=1-NI_raw
    NI_k =PnPs/((Dn_c/Ds_c)); a_k=1-NI_k
    wa=a_k*rt["dNdS_k"]
    DoS=(Nd/(Nd+Sd))-(Pn/(Pn+Ps))
    pval=fisher(round(Ps),round(Sd),round(Pn),round(Nd))
    return dict(label=label,Pn=Pn,Ps=Ps,Sd=round(Sd,1),Nd=round(Nd,1),
                dS=round(rt["dS"],3),dN=round(rt["dN"],3),dNdS=round(rt["dNdS_k"],3),
                PnPs=round(PnPs,3),NI=round(NI_k,3),
                a_raw=round(a_raw,3),alpha=round(a_k,3),omega_a=round(wa,4),
                DoS=round(DoS,3),fisher_p=round(pval,4))

# ── VCF: clade assignment, polymorphism counts, clade consensus genome ────────
def parse_vcf():
    samples=[];var=[]   # var: (pos,ref,alt,cls,carriers[list 0/1])
    with gzip.open(VCF,"rt") as fh:
        for line in fh:
            if line.startswith("##"): continue
            if line.startswith("#CHROM"): samples=line.rstrip().split("\t")[9:]; continue
            c=line.rstrip().split("\t")
            ann=""
            for kv in c[7].split(";"):
                if kv.startswith("ANN="): ann=kv[4:]
            cls="Other"
            for ent in ann.split(","):
                f=ent.split("|")
                if len(f)<2: continue
                effs=f[1].split("&")
                if any(e in NS_SET for e in effs): cls="NS"; break
                if "synonymous_variant" in effs: cls="SYN"; break
            gi=c[8].split(":").index("GT")
            car=[1 if s.split(":")[gi]=="1" else 0 for s in c[9:]]
            var.append((int(c[1]),c[3],c[4],cls,car))
    n=len(samples)
    tot=[0]*n
    for _,_,_,_,car in var:
        for i,v in enumerate(car): tot[i]+=v
    clade=["N" if tot[i]<50 else ("S" if tot[i]>200 else "C") for i in range(n)]
    return samples,var,clade

def poly_counts(var,idx,within_clade):
    """Pn,Ps. within_clade: count sites segregating within idx; else panel-wide (idx=all)."""
    Pn=Ps=0
    for pos,ref,alt,cls,car in var:
        cc=sum(car[i] for i in idx); k=len(idx)
        seg = (0<cc<k) if within_clade else (cc>0)
        if not seg: continue
        if cls=="NS": Pn+=1
        elif cls=="SYN": Ps+=1
    return Pn,Ps

def clade_consensus_genome(ref_seq,var,idx):
    g=list(ref_seq); k=len(idx)
    for pos,ref,alt,cls,car in var:
        if sum(car[i] for i in idx) > k/2:           # clade-majority ALT
            if len(ref)==1 and len(alt)==1:           # SNP
                g[pos-1]=alt
    return "".join(g)


def main():
    coords=gff_cds_coords()
    ref_genome=read_fasta(FHET_FA); rchrom=next(iter(ref_genome)); ref_seq=ref_genome[rchrom]
    fhet_ref_cds=cds_from_genome(ref_seq,coords)
    og=outgroup_cds()
    allgenes=[g for g,_,_,_ in coords]
    samples,var,clade=parse_vcf()
    Ni=[i for i,c in enumerate(clade) if c=="N"]
    Si=[i for i,c in enumerate(clade) if c=="S"]
    allidx=list(range(len(samples)))
    print(f"clades: N={len(Ni)} S={len(Si)} C={clade.count('C')}",file=sys.stderr)

    results=[]
    # whole species: D = reference vs catenatus; P = all panel variants
    div=divergence(fhet_ref_cds,og,allgenes); rt=rates(div)
    Pn,Ps=poly_counts(var,allidx,within_clade=False)
    results.append(mk("Whole species",Pn,Ps,div,rt))

    # North clade
    ng=clade_consensus_genome(ref_seq,var,Ni); ncds=cds_from_genome(ng,coords)
    dN_=divergence(ncds,og,allgenes); rtN=rates(dN_)
    Pn,Ps=poly_counts(var,Ni,within_clade=True)
    results.append(mk("North clade",Pn,Ps,dN_,rtN))

    # South clade
    sg=clade_consensus_genome(ref_seq,var,Si); scds=cds_from_genome(sg,coords)
    dS_=divergence(scds,og,allgenes); rtS=rates(dS_)
    Pn,Ps=poly_counts(var,Si,within_clade=True)
    results.append(mk("South clade",Pn,Ps,dS_,rtS))

    # sensitivity: whole species, low-saturation genes only
    divL=divergence(fhet_ref_cds,og,[g for g in allgenes if g in LOWSAT]); rtL=rates(divL)
    PnL=PsL=0
    for pos,ref,alt,cls,car in var:
        pass
    # P restricted to low-sat genes: need gene of each variant — map by CDS coords
    gene_of={}
    for gene,s,e,strand in coords:
        for p in range(s,e+1): gene_of[p]=gene
    Pn=Ps=0
    for pos,ref,alt,cls,car in var:
        if gene_of.get(pos) in LOWSAT and sum(car)>0:
            if cls=="NS":Pn+=1
            elif cls=="SYN":Ps+=1
    results.append(mk("Whole sp. (low-sat genes)",Pn,Ps,divL,rtL))

    cols=["label","Pn","Ps","Sd","Nd","dS","dN","dNdS","PnPs","NI","a_raw","alpha","omega_a","DoS","fisher_p"]
    out=f"{MTG}/vcf/mk_{OG_NAME}.tsv"
    with open(out,"w") as fh:
        fh.write("\t".join(cols)+"\n")
        for r in results: fh.write("\t".join(str(r[c]) for c in cols)+"\n")

    print("\n===== McDonald-Kreitman vs F. catenatus (K2P-corrected divergence) =====")
    h=["test","Pn","Ps","Dn","Ds","dN/dS","Pn/Ps","NI","a_raw","alpha","omega_a","DoS","Fisher_p"]
    print(f"{h[0]:26}"+"".join(f"{x:>9}" for x in h[1:]))
    for r in results:
        print(f"{r['label']:26}"+"".join(f"{x:>9}" for x in
            [r['Pn'],r['Ps'],r['Nd'],r['Sd'],r['dNdS'],r['PnPs'],r['NI'],r['a_raw'],r['alpha'],r['omega_a'],r['DoS'],r['fisher_p']]))
    print(f"\nalpha = K2P-corrected (1 - NI);  a_raw = uncorrected Nei-Gojobori counts;")
    print(f"omega_a = alpha * (dN/dS).  Dn/Ds shown are raw NG fixed-difference counts.")
    print(f"Wrote {out}")

if __name__=="__main__":
    main()

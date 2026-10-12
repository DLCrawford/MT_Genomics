#!/usr/bin/env python3
"""
30_outgroup_dS_gate.py  —  SATURATION GATE for an outgroup-based MK test.

Before running a McDonald-Kreitman test of F. heteroclitus mt polymorphism
against an outgroup, we must verify that synonymous divergence to the outgroup
is NOT saturated (multiple-hit).  Saturated dS inflates the synonymous
divergence, biases dN/dS, and makes count-based alpha / omega_a unreliable.

This script computes Fhet-vs-outgroup divergence per CDS gene by:
  1. extracting the 13 ortholog CDS from each genome (coding sense),
  2. translating (NCBI table 2) and globally aligning the proteins
     (self-contained Needleman-Wunsch — no Biopython needed),
  3. back-translating to a codon alignment,
  4. Nei-Gojobori counting of synonymous / nonsynonymous SITES and DIFFERENCES
     (multi-difference codons averaged over mutational pathways),
  5. reporting uncorrected pN/pS and Jukes-Cantor dN/dS, per gene and overall,
  6. saturation diagnostics: transition/transversion ratio at synonymous
     differences, and the proportion of 4-fold-degenerate sites that differ.

Default outgroup: F. catenatus (NC_084079.1).  Decision rule printed at end.

Outputs: vcf/dS_<outgroup>_gate.tsv  + stdout summary.
"""
import json, math, os, sys

BASE = os.environ.get("MTG_BASE", os.path.expanduser("~/Projects/MT_Genomics_Cl_Ap2026"))
MTG  = f"{BASE}/MT_Genomics2"

OUTGROUP   = os.environ.get("MTG_OG", "catenatus")   # set MTG_OG=grandis to switch
OG_FASTA   = f"{MTG}/refs_outgroup/F_{OUTGROUP}/ncbi_dataset/data/genome.fna"
OG_REPORT  = f"{MTG}/refs_outgroup/F_{OUTGROUP}/ncbi_dataset/data/organelle_gene_report.jsonl"
FHET_FASTA = f"{MTG}/Missing_Files/SSM_MT_ref/Fhet_MT.fasta"
FHET_GFF   = f"{MTG}/Missing_Files/SSM_MT_ref/Fhet_MT.gff"

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
COMP = str.maketrans("ACGTN", "TGCAN")
BASES = "ACGT"
PURINES = set("AG")

# Fhet GFF gene name -> canonical; catenatus uses same symbols (COX*, CYTB, ND*)
def norm_gene(g):
    g = g.upper().replace("CO", "COX") if g.upper() in ("CO1","CO2","CO3") else g.upper()
    return {"COB":"CYTB","CYB":"CYTB","COI":"COX1","COII":"COX2","COIII":"COX3"}.get(g, g)


def read_fasta(path, want=None):
    seqs={}; name=None; buf=[]
    with open(path) as fh:
        for line in fh:
            line=line.rstrip()
            if line.startswith(">"):
                if name is not None: seqs[name]="".join(buf).upper()
                name=line[1:].split()[0]; buf=[]
            else: buf.append(line)
    if name is not None: seqs[name]="".join(buf).upper()
    if want: return seqs[want]
    return seqs


def translate(seq):
    L=len(seq)-(len(seq)%3)
    return "".join(MT2.get(seq[i:i+3],"X") for i in range(0,L,3))


def fhet_cds():
    g=read_fasta(FHET_FASTA); chrom=next(iter(g)); genome=g[chrom]
    out={}
    with open(FHET_GFF) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip(): continue
            p=line.rstrip("\n").split("\t")
            if len(p)<9 or p[2]!="CDS": continue
            attrs=dict(kv.split("=",1) for kv in p[8].split(";") if "=" in kv)
            gene=norm_gene(attrs.get("gene") or attrs.get("Name",""))
            s,e,strand=int(p[3]),int(p[4]),p[6]
            seq=genome[s-1:e]
            if strand=="-": seq=seq.translate(COMP)[::-1]
            out[gene]=seq
    return out


def outgroup_cds():
    seqs=read_fasta(OG_FASTA)
    acc=next((k for k in seqs if k.startswith("NC_")), next(iter(seqs)))
    genome=seqs[acc]
    out={}
    for line in open(OG_REPORT):
        if not line.strip(): continue
        r=json.loads(line)
        if r.get("type")!="PROTEIN_CODING": continue
        loc=r["annotations"][0]["genomicLocations"][0]["genomicRange"]
        s,e=int(loc["begin"]),int(loc["end"]); strand=loc["orientation"]
        seq=genome[s-1:e]
        if strand=="minus": seq=seq.translate(COMP)[::-1]
        out[norm_gene(r["symbol"])]=seq
    return out


# ── Needleman-Wunsch global protein alignment (identity scoring) ──────────────
def nw_align(a, b, match=2, mis=-1, gap=-2):
    n,m=len(a),len(b)
    F=[[0]*(m+1) for _ in range(n+1)]
    for i in range(1,n+1): F[i][0]=i*gap
    for j in range(1,m+1): F[0][j]=j*gap
    for i in range(1,n+1):
        ai=a[i-1]
        for j in range(1,m+1):
            d=F[i-1][j-1]+(match if ai==b[j-1] else mis)
            F[i][j]=max(d, F[i-1][j]+gap, F[i][j-1]+gap)
    # traceback
    i,j=n,m; A=[]; B=[]
    while i>0 and j>0:
        s=F[i][j]
        if s==F[i-1][j-1]+(match if a[i-1]==b[j-1] else mis):
            A.append(a[i-1]); B.append(b[j-1]); i-=1; j-=1
        elif s==F[i-1][j]+gap:
            A.append(a[i-1]); B.append("-"); i-=1
        else:
            A.append("-"); B.append(b[j-1]); j-=1
    while i>0: A.append(a[i-1]); B.append("-"); i-=1
    while j>0: A.append("-"); B.append(b[j-1]); j-=1
    return "".join(reversed(A)), "".join(reversed(B))


def codon_alignment(cds1, cds2):
    """Return list of (codon1, codon2) for ungapped aligned positions."""
    p1, p2 = translate(cds1), translate(cds2)
    a1, a2 = nw_align(p1, p2)
    cods1=[cds1[i:i+3] for i in range(0, len(cds1)-len(cds1)%3, 3)]
    cods2=[cds2[i:i+3] for i in range(0, len(cds2)-len(cds2)%3, 3)]
    out=[]; i=j=0
    for x,y in zip(a1,a2):
        c1 = cods1[i] if x!="-" else None
        c2 = cods2[j] if y!="-" else None
        if x!="-": i+=1
        if y!="-": j+=1
        if c1 and c2 and "N" not in c1 and "N" not in c2:
            if MT2.get(c1)=="*" or MT2.get(c2)=="*":   # skip stop codons
                continue
            out.append((c1,c2))
    return out


def codon_S_sites(codon):
    ref=MT2.get(codon)
    if ref is None or ref=="*": return None
    s=0.0
    for pos in range(3):
        nsyn=0
        for b in BASES:
            if b==codon[pos]: continue
            if MT2.get(codon[:pos]+b+codon[pos+1:],"?")==ref: nsyn+=1
        s+=nsyn/3.0
    return s


def path_diffs(c1, c2):
    """Nei-Gojobori syn/nonsyn differences between two codons (avg over paths)."""
    pos=[k for k in range(3) if c1[k]!=c2[k]]
    if not pos: return 0.0,0.0,[]
    if len(pos)==1:
        k=pos[0]
        syn = 1.0 if MT2[c1]==MT2[c2] else 0.0
        ts = (c1[k] in PURINES)==(c2[k] in PURINES)   # transition?
        return syn, 1.0-syn, [(syn>0, ts)]
    # multi-position: average over orderings, skip stop-containing paths
    import itertools
    sds=[]; nds=[]; events=[]
    for order in itertools.permutations(pos):
        cur=c1; ok=True; sd=nd=0.0; ev=[]
        for k in order:
            nxt=cur[:k]+c2[k]+cur[k+1:]
            if MT2.get(nxt)=="*": ok=False; break
            syn = 1.0 if MT2[cur]==MT2[nxt] else 0.0
            sd+=syn; nd+=(1.0-syn)
            ev.append((syn>0,(cur[k] in PURINES)==(nxt[k] in PURINES)))
            cur=nxt
        if ok: sds.append(sd); nds.append(nd); events.append(ev)
    if not sds: return 0.0,float(len(pos)),[]
    # average
    sd=sum(sds)/len(sds); nd=sum(nds)/len(nds)
    flat=[e for ev in events for e in ev]   # for Ts/Tv tally (syn diffs only)
    return sd, nd, flat


def jc(p):
    if p<=0: return 0.0
    if p>=0.75: return float("inf")
    return -0.75*math.log(1-4*p/3)


def main():
    f=fhet_cds(); o=outgroup_cds()
    genes=[g for g in ("ND1","ND2","COX1","COX2","ATP8","ATP6","COX3","ND3","ND4L","ND4","ND5","ND6","CYTB") if g in f and g in o]
    rows=[]
    tot_Ss=tot_Ns=tot_Sd=tot_Nd=0.0
    syn_ts=syn_tv=0
    ff_sites=ff_diff=0   # 4-fold degenerate sites
    for g in genes:
        pairs=codon_alignment(f[g],o[g])
        Ss=Ns=Sd=Nd=0.0
        for c1,c2 in pairs:
            s1=codon_S_sites(c1); s2=codon_S_sites(c2)
            if s1 is None or s2 is None: continue
            Ss+=(s1+s2)/2.0; Ns+=(3-(s1+s2)/2.0)
            sd,nd,ev=path_diffs(c1,c2)
            Sd+=sd; Nd+=nd
            for issyn,ists in ev:
                if issyn:
                    if ists: syn_ts+=1
                    else: syn_tv+=1
            # 4-fold site (3rd position) saturation probe
            stem=c1[:2]
            if len({MT2.get(stem+b) for b in BASES})==1:   # all 4 same aa => 4fold
                ff_sites+=1
                if c1[2]!=c2[2]: ff_diff+=1
        pS=Sd/Ss if Ss else float("nan"); pN=Nd/Ns if Ns else float("nan")
        rows.append((g,len(pairs),round(Ss,1),round(Ns,1),round(Sd,1),round(Nd,1),
                     round(pS,3),round(pN,3),round(jc(pS),3),round(jc(pN),3),
                     round((jc(pN)/jc(pS)) if jc(pS) not in (0,float("inf")) else float("nan"),3)))
        tot_Ss+=Ss; tot_Ns+=Ns; tot_Sd+=Sd; tot_Nd+=Nd

    pS=tot_Sd/tot_Ss; pN=tot_Nd/tot_Ns
    dS=jc(pS); dN=jc(pN)
    tstv=syn_ts/syn_tv if syn_tv else float("inf")
    ff_prop=ff_diff/ff_sites if ff_sites else float("nan")

    # write per-gene tsv
    out_tsv=f"{MTG}/vcf/dS_{OUTGROUP}_gate.tsv"
    with open(out_tsv,"w") as fh:
        fh.write("Gene\tcodons\tS_sites\tN_sites\tSd\tNd\tpS\tpN\tdS_JC\tdN_JC\tdNdS_JC\n")
        for r in rows: fh.write("\t".join(str(x) for x in r)+"\n")
        fh.write(f"Overall\t{sum(r[1] for r in rows)}\t{tot_Ss:.1f}\t{tot_Ns:.1f}\t{tot_Sd:.1f}\t{tot_Nd:.1f}\t{pS:.3f}\t{pN:.3f}\t{dS:.3f}\t{dN:.3f}\t{(dN/dS) if dS not in (0,float('inf')) else 'NA'}\n")

    print(f"\n===== Fhet vs F. {OUTGROUP} divergence (NCBI table 2) =====")
    print(f"{'Gene':6}{'codons':>7}{'pS':>8}{'dS_JC':>8}{'pN':>8}{'dN_JC':>8}{'dN/dS':>8}")
    for r in rows:
        print(f"{r[0]:6}{r[1]:>7}{r[6]:>8}{r[8]:>8}{r[7]:>8}{r[9]:>8}{r[10]:>8}")
    print("-"*53)
    print(f"{'ALL':6}{sum(r[1] for r in rows):>7}{pS:>8.3f}{dS:>8.3f}{pN:>8.3f}{dN:>8.3f}{(dN/dS):>8.3f}")
    print(f"\n── Saturation diagnostics ──")
    print(f"  Uncorrected synonymous divergence   pS = {pS:.3f}")
    print(f"  JC-corrected synonymous divergence  dS = {dS:.3f}   (inflation dS/pS = {dS/pS:.2f}x)")
    print(f"  Ts/Tv among synonymous differences     = {tstv:.2f}   (→0.5 at saturation)")
    print(f"  4-fold-degenerate sites differing      = {ff_prop:.3f}  ({ff_diff}/{ff_sites}; →0.75 saturated)")
    print(f"\n── Gate decision ──")
    if pS>=0.75:
        verdict="FAIL: synonymous sites SATURATED (pS>=0.75, dS undefined). Count-based MK invalid."
    elif pS>0.5 or ff_prop>0.55:
        verdict="MARGINAL: high synonymous divergence. Use a codon model (PAML codeml) for dN/dS, or a closer non-introgressed outgroup."
    else:
        verdict="PASS: synonymous divergence not saturated. Count-based MK (alpha, omega_a) is defensible."
    print("  "+verdict)
    print(f"\nWrote {out_tsv}")


if __name__=="__main__":
    main()

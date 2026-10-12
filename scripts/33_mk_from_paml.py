#!/usr/bin/env python3
"""
33_mk_from_paml.py  —  MK alpha / omega_a from PAML divergence + Fhet polymorphism.

Reads codeml/yn00 pairwise dN,dS matrices (2ML.*, 2YN.*, 2NG.*) for the
Fhet-vs-F.catenatus comparisons, combines with the 141-panel polymorphism
(whole species and per clade), and reports:

    alpha   = 1 - NI,   NI = (pN/pS_polymorphism)/(dN/dS_divergence)
    omega_a = alpha * (dN/dS)        (rate of adaptive evolution)

pN/pS_polymorphism = (Pn/Ps) * (L_S/L_N), using the project's Nei-Gojobori site
counts (L_S=2693.67, L_N=8694.33 from dnds_per_gene.tsv) so polymorphism and the
manuscript dN/dS share one site convention.

dN/dS_divergence is reported THREE ways — ML (codeml, recommended), YN
(Yang-Nielsen), NG (Nei-Gojobori) — because at this high divergence (dS~1) the
estimators straddle the alpha=0 line; showing all three makes the (lack of)
robustness explicit.

Output: vcf/mk_outgroup_paml.tsv
"""
import importlib.util, os, sys

HERE=os.path.dirname(os.path.abspath(__file__))
spec=importlib.util.spec_from_file_location("mk",os.path.join(HERE,"31_mk_outgroup.py"))
mk=importlib.util.module_from_spec(spec); spec.loader.exec_module(mk)
PAML=f"{mk.MTG}/paml" if mk.OG_NAME=="catenatus" else f"{mk.MTG}/paml_{mk.OG_NAME}"

L_S, L_N = 2693.67, 8694.33          # from vcf/dnds_per_gene.tsv (NCBI table 2)
LSLN = L_S/L_N

def read_matrix(path):
    """PHYLIP lower-tri; return F_caten row = [vs Fhet_ref, vs Fhet_N, vs Fhet_S]."""
    for line in open(path):
        if line.startswith("F_caten"):
            return [float(x) for x in line.split()[1:]]
    raise RuntimeError(f"F_caten row not found in {path}")

def main():
    div={}
    for tag in ("ML","YN","NG"):
        dN=read_matrix(f"{PAML}/2{tag}.dN"); dS=read_matrix(f"{PAML}/2{tag}.dS")
        div[tag]={"Fhet_ref":(dN[0],dS[0]),"Fhet_N":(dN[1],dS[1]),"Fhet_S":(dN[2],dS[2])}

    _,var,clade=mk.parse_vcf()
    Ni=[i for i,c in enumerate(clade) if c=="N"]
    Si=[i for i,c in enumerate(clade) if c=="S"]
    allidx=list(range(len(clade)))

    tests=[("Whole species","Fhet_ref",allidx,False),
           ("North clade","Fhet_N",Ni,True),
           ("South clade","Fhet_S",Si,True)]

    rows=[]
    for label,key,idx,within in tests:
        Pn,Ps=mk.poly_counts(var,idx,within_clade=within)
        pnps=(Pn/Ps)*LSLN
        for tag in ("ML","YN","NG"):
            dN,dS=div[tag][key]; dnds=dN/dS
            NI=pnps/dnds; alpha=1-NI; omega_a=alpha*dnds
            rows.append(dict(test=label,Pn=Pn,Ps=Ps,pN_pS=round(pnps,4),
                method=tag,dN=round(dN,4),dS=round(dS,4),dNdS=round(dnds,4),
                NI=round(NI,3),alpha=round(alpha,3),omega_a=round(omega_a,4)))

    cols=["test","Pn","Ps","pN_pS","method","dN","dS","dNdS","NI","alpha","omega_a"]
    out=f"{mk.MTG}/vcf/mk_outgroup_paml_{mk.OG_NAME}.tsv"
    with open(out,"w") as fh:
        fh.write("\t".join(cols)+"\n")
        for r in rows: fh.write("\t".join(str(r[c]) for c in cols)+"\n")

    print(f"{'test':14}{'meth':5}{'Pn':>5}{'Ps':>5}{'pN/pS':>8}{'dN':>8}{'dS':>8}{'dN/dS':>8}{'NI':>7}{'alpha':>8}{'omega_a':>9}")
    for r in rows:
        print(f"{r['test']:14}{r['method']:5}{r['Pn']:>5}{r['Ps']:>5}{r['pN_pS']:>8}{r['dN']:>8}{r['dS']:>8}{r['dNdS']:>8}{r['NI']:>7}{r['alpha']:>8}{r['omega_a']:>9}")
    print(f"\nWrote {out}")

if __name__=="__main__":
    main()

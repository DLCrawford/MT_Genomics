#!/usr/bin/env python3
"""
34_mk_test_catenatus.py  —  full MK contingency test vs F. catenatus,
whole species + each clade, with significance.

Reports the MK 2x2 two ways, because synonymous divergence to catenatus is
high (dS~1):

  OBSERVED  : Dn,Ds = raw counts of nonsyn/syn fixed differences (classic MK).
              At dS~1 the synonymous count is SATURATED (undercounted), which
              inflates observed Dn/Ds and biases alpha UP (toward false
              'adaptation').
  CORRECTED : Dn,Ds = multiple-hit-corrected substitution numbers from PAML
              codeml ML (dN*N_sites, dS*S_sites).  Removes the saturation bias.

For each: NI, alpha = 1-NI, DoS, and a two-sided Fisher exact p.
(omega_a = alpha_corrected * dN/dS.)  Corrected Fisher uses model-estimated
substitution numbers in the divergence cell, so its p is approximate.

Output: vcf/mk_test_catenatus.tsv
"""
import importlib.util, math, os

HERE=os.path.dirname(os.path.abspath(__file__))
spec=importlib.util.spec_from_file_location("mk",os.path.join(HERE,"31_mk_outgroup.py"))
mk=importlib.util.module_from_spec(spec); spec.loader.exec_module(mk)
PAML=f"{mk.MTG}/paml" if mk.OG_NAME=="catenatus" else f"{mk.MTG}/paml_{mk.OG_NAME}"

def read_mat(path):
    for line in open(path):
        if line.startswith("F_caten"):
            return [float(x) for x in line.split()[1:]]

def gtest(a,b,c,d):
    n=a+b+c+d
    if min(a,b,c,d)==0: return float("nan")
    E=[(a+b)*(a+c)/n,(a+b)*(b+d)/n,(c+d)*(a+c)/n,(c+d)*(b+d)/n]
    O=[a,b,c,d]
    return 2*sum(o*math.log(o/e) for o,e in zip(O,E) if o>0)

def stats(Pn,Ps,Dn,Ds):
    NI=(Pn/Ps)/(Dn/Ds)
    alpha=1-NI
    DoS=Dn/(Dn+Ds)-Pn/(Pn+Ps)
    p=mk.fisher(round(Ps),round(Ds),round(Pn),round(Dn))
    return NI,alpha,DoS,p

def main():
    dN_ML=read_mat(f"{PAML}/2ML.dN"); dS_ML=read_mat(f"{PAML}/2ML.dS")  # [ref,N,S]
    coords=mk.gff_cds_coords()
    ref=mk.read_fasta(mk.FHET_FA); rs=ref[next(iter(ref))]
    _,var,clade=mk.parse_vcf()
    Ni=[i for i,c in enumerate(clade) if c=="N"]; Si=[i for i,c in enumerate(clade) if c=="S"]
    allidx=list(range(len(clade)))
    cds_ref=mk.cds_from_genome(rs,coords)
    cds_N=mk.cds_from_genome(mk.clade_consensus_genome(rs,var,Ni),coords)
    cds_S=mk.cds_from_genome(mk.clade_consensus_genome(rs,var,Si),coords)
    og=mk.outgroup_cds(); genes=[g for g,_,_,_ in coords]

    setup=[("Whole species",cds_ref,allidx,False,0),
           ("North clade",  cds_N,  Ni,    True, 1),
           ("South clade",  cds_S,  Si,    True, 2)]
    rows=[]
    for label,cds,idx,within,k in setup:
        d=mk.divergence(cds,og,genes)
        Dn_o,Ds_o=d["Nd"],d["Sd"]; N,S=d["N"],d["S"]
        Pn,Ps=mk.poly_counts(var,idx,within_clade=within)
        # observed
        NIo,ao,DoSo,po=stats(Pn,Ps,Dn_o,Ds_o)
        # ML-corrected substitution numbers
        Dn_c=dN_ML[k]*N; Ds_c=dS_ML[k]*S
        NIc,ac,DoSc,pc=stats(Pn,Ps,Dn_c,Ds_c)
        omega=ac*(dN_ML[k]/dS_ML[k])
        rows.append(dict(test=label,Pn=Pn,Ps=Ps,
            Dn_obs=round(Dn_o,1),Ds_obs=round(Ds_o,1),NI_obs=round(NIo,3),
            alpha_obs=round(ao,3),DoS_obs=round(DoSo,3),fisherP_obs=round(po,4),
            dN_ML=round(dN_ML[k],4),dS_ML=round(dS_ML[k],4),
            Dn_cor=round(Dn_c,1),Ds_cor=round(Ds_c,1),NI_cor=round(NIc,3),
            alpha_cor=round(ac,3),omega_a=round(omega,4),fisherP_cor=round(pc,4)))

    cols=["test","Pn","Ps","Dn_obs","Ds_obs","NI_obs","alpha_obs","DoS_obs","fisherP_obs",
          "dN_ML","dS_ML","Dn_cor","Ds_cor","NI_cor","alpha_cor","omega_a","fisherP_cor"]
    out=f"{mk.MTG}/vcf/mk_test_{mk.OG_NAME}.tsv"
    with open(out,"w") as fh:
        fh.write("\t".join(cols)+"\n")
        for r in rows: fh.write("\t".join(str(r[c]) for c in cols)+"\n")

    print("OBSERVED counts (classic MK; synonymous SATURATED -> alpha biased UP):")
    print(f"{'test':14}{'Pn':>5}{'Ps':>5}{'Dn':>8}{'Ds':>8}{'NI':>7}{'alpha':>8}{'DoS':>7}{'Fisher_p':>10}")
    for r in rows:
        print(f"{r['test']:14}{r['Pn']:>5}{r['Ps']:>5}{r['Dn_obs']:>8}{r['Ds_obs']:>8}{r['NI_obs']:>7}{r['alpha_obs']:>8}{r['DoS_obs']:>7}{r['fisherP_obs']:>10}")
    print("\nML-CORRECTED (PAML codeml; multiple-hit corrected -> defensible):")
    print(f"{'test':14}{'Pn':>5}{'Ps':>5}{'Dn':>8}{'Ds':>8}{'NI':>7}{'alpha':>8}{'omega_a':>9}{'Fisher_p':>10}")
    for r in rows:
        print(f"{r['test']:14}{r['Pn']:>5}{r['Ps']:>5}{r['Dn_cor']:>8}{r['Ds_cor']:>8}{r['NI_cor']:>7}{r['alpha_cor']:>8}{r['omega_a']:>9}{r['fisherP_cor']:>10}")
    print(f"\nWrote {out}")

if __name__=="__main__":
    main()

#!/usr/bin/env python3
"""
01_pnps_bootstrap.py
--------------------
Within-species pN/pS per gene and per OXPHOS complex, with 95% bootstrap CIs.
This is the PRIMARY purifying-selection statistic (all 141 individuals).

  pN     = obs_N / N_sites      (nonsynonymous polymorphisms per nonsyn site)
  pS     = obs_S / S_sites      (synonymous  polymorphisms per syn  site)
  pN/pS  = pN / pS              (<< 1 => purifying selection)

N_sites / S_sites are the mutational-opportunity counts from the codon model
(Nei-Gojobori style), read from vcf/dnds_per_gene.tsv. The CI is a percentile
bootstrap: the obs_S + obs_N segregating sites of a gene (or pooled complex)
are resampled with replacement B times; N_sites/S_sites stay fixed.

Input : ../vcf/dnds_per_gene.tsv
Output: pnps_bootstrap_by_gene.tsv, pnps_bootstrap_by_complex.tsv
"""
import csv, random
from common import (COMPLEXES, GENE2CX, GENES, B, opath, vpath, fmt)

random.seed(2)   # reproducibility (matches the delivered tables)

# ---- load per-gene site counts ----
G = {}
with open(vpath("dnds_per_gene.tsv")) as f:
    for row in csv.DictReader(f, delimiter="\t"):
        if not row["Gene"] or row["Gene"] == "Overall":
            continue
        G[row["Gene"]] = dict(
            Ssite=float(row["S_sites"]), Nsite=float(row["N_sites"]),
            oS=int(row["obs_S"]), oN=int(row["obs_N"]),
        )

def pnps(oS, oN, Ss, Ns):
    if oS == 0 or Ss == 0:
        return None
    return (oN / Ns) / (oS / Ss)

def boot_group(glist):
    """Point estimate + 95% CI for pooled genes; bootstrap resamples segregating sites."""
    Ss = sum(G[g]["Ssite"] for g in glist)
    Ns = sum(G[g]["Nsite"] for g in glist)
    oS = sum(G[g]["oS"]    for g in glist)
    oN = sum(G[g]["oN"]    for g in glist)
    pe = pnps(oS, oN, Ss, Ns)
    tokens = [1] * oN + [0] * oS          # 1 = nonsyn, 0 = syn
    if len(tokens) < 2:
        return pe, None, None, oS, oN
    vals = []
    for _ in range(B):
        r  = [tokens[random.randrange(len(tokens))] for _ in tokens]
        bN = sum(r); bS = len(r) - bN
        v  = pnps(bS, bN, Ss, Ns)
        if v is not None:
            vals.append(v)
    vals.sort()
    return pe, vals[int(0.025*len(vals))], vals[int(0.975*len(vals))], oS, oN

# ---- per gene ----
with open(opath("pnps_bootstrap_by_gene.tsv"), "w") as o:
    o.write("Gene\tComplex\tobs_S\tobs_N\tpN_pS\tCI_low\tCI_high\n")
    for g in GENES:
        pe, lo, hi, oS, oN = boot_group([g])
        o.write("%s\t%s\t%d\t%d\t%s\t%s\t%s\n" %
                (g, GENE2CX[g], oS, oN, fmt(pe,4), fmt(lo,4), fmt(hi,4)))

# ---- per complex + overall ----
with open(opath("pnps_bootstrap_by_complex.tsv"), "w") as o:
    o.write("Group\tobs_S\tobs_N\tpN_pS\tCI_low\tCI_high\n")
    for cx in ("I", "III", "IV", "V"):
        pe, lo, hi, oS, oN = boot_group(COMPLEXES[cx])
        o.write("Complex_%s\t%d\t%d\t%s\t%s\t%s\n" % (cx, oS, oN, fmt(pe,4), fmt(lo,4), fmt(hi,4)))
    pe, lo, hi, oS, oN = boot_group(GENES)
    o.write("Overall\t%d\t%d\t%s\t%s\t%s\n" % (oS, oN, fmt(pe,4), fmt(lo,4), fmt(hi,4)))

print("01_pnps_bootstrap.py: wrote pnps_bootstrap_by_gene.tsv and pnps_bootstrap_by_complex.tsv")

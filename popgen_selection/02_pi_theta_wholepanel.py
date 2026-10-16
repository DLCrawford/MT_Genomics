#!/usr/bin/env python3
"""
02_pi_theta_wholepanel.py
-------------------------
pi/theta by site class (synonymous vs nonsynonymous), whole panel (all 141),
per gene and per OXPHOS complex, with 95% bootstrap CIs.

The synonymous-vs-nonsynonymous CONTRAST here is a frequency-spectrum signal of
purifying selection, but note (see 03 and the methods) that with all 141
individuals it is dominated by the deep Northern/Southern mt-clade split: the
clade-diagnostic variants sit at intermediate frequency and are overwhelmingly
synonymous, which inflates synonymous pi relative to theta. It is a BETWEEN-clade
signal; the within-clade control is script 03.

Estimator (per site class, over its biallelic segregating sites):
  pi/theta = a_n * sum(pi_site) / S ,  a_n = harmonic(141) = 5.522
  pi_site  = 2*AC*(n-AC)/(n*(n-1))   (biallelic, unbiased)

Input : ../vcf/haplotype_matrix.csv   (one 0/1 column per <id>_MT)
Output: pi_theta_wholepanel_by_gene.tsv, pi_theta_wholepanel_by_complex.tsv
"""
import csv, random
from common import (COMPLEXES, GENE2CX, GENES, A_N, N_PANEL, B,
                    pi_theta, boot_ci, pi_site_biallelic, opath, vpath, fmt)

random.seed(3)
a_n = A_N["panel"]

# ---- collect per-site (gene, class, pi_site) over segregating biallelic sites ----
sites = []
with open(vpath("haplotype_matrix.csv")) as f:
    for r in csv.DictReader(f):
        gt = [r[k] for k in r if k.endswith("_MT")]
        n  = len(gt)
        AC = sum(1 for x in gt if x == "1")
        if AC == 0 or AC == n:          # monomorphic in the panel
            continue
        vt = r["Variant_type"]
        if "synonymous" in vt and "missense" not in vt:
            cls = "syn"
        elif "missense" in vt:
            cls = "ns"
        else:
            continue                    # e.g. stop/indel effects excluded
        sites.append((r["Gene"], cls, pi_site_biallelic(AC, n)))

def subset(gene=None, cx=None, cls=None):
    pis = [pi for g, c, pi in sites
           if (gene is None or g == gene)
           and (cx   is None or GENE2CX.get(g) == cx)
           and (cls  is None or c == cls)]
    pe = pi_theta(pis, a_n)
    lo, hi = boot_ci(pis, a_n, seed=1)
    return len(pis), pe, lo, hi

CLASSES = [("syn", "synonymous"), ("ns", "nonsynonymous"), (None, "all")]

with open(opath("pi_theta_wholepanel_by_gene.tsv"), "w") as o:
    o.write("Gene\tComplex\tClass\tS_sites\tpi_theta\tCI_low\tCI_high\n")
    for g in GENES:
        for cls, lab in CLASSES:
            S, pe, lo, hi = subset(gene=g, cls=cls)
            o.write("%s\t%s\t%s\t%d\t%s\t%s\t%s\n" % (g, GENE2CX[g], lab, S, fmt(pe), fmt(lo), fmt(hi)))

with open(opath("pi_theta_wholepanel_by_complex.tsv"), "w") as o:
    o.write("Group\tClass\tS_sites\tpi_theta\tCI_low\tCI_high\n")
    for cx in ("I", "III", "IV", "V"):
        for cls, lab in CLASSES:
            S, pe, lo, hi = subset(cx=cx, cls=cls)
            o.write("Complex_%s\t%s\t%d\t%s\t%s\t%s\n" % (cx, lab, S, fmt(pe), fmt(lo), fmt(hi)))
    for cls, lab in CLASSES:
        S, pe, lo, hi = subset(cls=cls)
        o.write("Overall\t%s\t%d\t%s\t%s\t%s\n" % (lab, S, fmt(pe), fmt(lo), fmt(hi)))

print("02_pi_theta_wholepanel.py: %d segregating sites, a_n(%d)=%.3f" % (len(sites), N_PANEL, a_n))
print("  wrote pi_theta_wholepanel_by_gene.tsv and pi_theta_wholepanel_by_complex.tsv")

#!/usr/bin/env python3
"""
03_pi_theta_withinclade.py
--------------------------
pi/theta by site class, computed SEPARATELY within the Northern (n=77) and
Southern (n=63) mt clades, per gene and per complex, with 95% bootstrap CIs.

This is the reviewer's requested control. Removing the between-clade split
collapses the whole-panel synonymous>nonsynonymous gap: within each clade both
classes are dominated by rare intra-clade variants, so syn and NS pi/theta
converge (~0.2) with overlapping CIs. The whole-panel gap (script 02) is
therefore a clade-structure effect, not a within-population SFS distortion, and
the within-species selection evidence rests on pN/pS (script 01).

A site is "segregating within a clade" iff its per-clade pi_site > 0 (i.e.
0 < allele freq < 1 in that clade). a_n uses that clade's sample size.

Input : ../vcf/pi_by_clade_persite.tsv  (per-variant north_pi_site/south_pi_site, Class)
        ../vcf/haplotype_matrix.csv      (POS/REF/ALT -> Gene)
Output: pi_theta_withinclade_by_gene.tsv, pi_theta_withinclade_by_complex.tsv
"""
import csv, random
from common import (COMPLEXES, GENE2CX, GENES, A_N, N_NORTH, N_SOUTH, B,
                    pi_theta, boot_ci, opath, vpath, fmt, read_tsv)

random.seed(4)

# ---- POS/REF/ALT -> Gene ----
gmap = {}
with open(vpath("haplotype_matrix.csv")) as f:
    for row in csv.DictReader(f):
        gmap[(row["POS"], row["REF"], row["ALT"])] = row["Gene"]

# ---- collect per-clade segregating sites (gene, class, pi_site) ----
# pi_by_clade_persite.tsv carries trailing summary columns; keep only first 15.
data = {"north": [], "south": []}
for d in read_tsv(vpath("pi_by_clade_persite.tsv"), ncols=15):
    try:
        int(d["north_n_called"])
    except (ValueError, KeyError):
        continue
    g = gmap[(d["POS"], d["REF"], d["ALT"])]
    cls = {"synonymous": "syn", "nonsynonymous": "ns"}.get(d["Class"])
    if cls is None:
        continue
    for clade, fld in (("north", "north_pi_site"), ("south", "south_pi_site")):
        try:
            pi = float(d[fld])
        except ValueError:
            pi = 0.0
        if pi > 0:                       # segregating in that clade
            data[clade].append((g, cls, pi))

def subset(clade, gene=None, cx=None, cls=None):
    a_n = A_N[clade]
    pis = [pi for g, c, pi in data[clade]
           if (gene is None or g == gene)
           and (cx   is None or GENE2CX.get(g) == cx)
           and (cls  is None or c == cls)]
    pe = pi_theta(pis, a_n)
    lo, hi = boot_ci(pis, a_n, seed=1)
    return len(pis), pe, lo, hi

CLASSES = [("syn", "synonymous"), ("ns", "nonsynonymous"), (None, "all")]

with open(opath("pi_theta_withinclade_by_gene.tsv"), "w") as o:
    o.write("Clade\tGene\tComplex\tClass\tS_sites\tpi_theta\tCI_low\tCI_high\n")
    for clade in ("north", "south"):
        for g in GENES:
            for cls, lab in CLASSES:
                S, pe, lo, hi = subset(clade, gene=g, cls=cls)
                o.write("%s\t%s\t%s\t%s\t%d\t%s\t%s\t%s\n" %
                        (clade.capitalize(), g, GENE2CX[g], lab, S, fmt(pe), fmt(lo), fmt(hi)))

with open(opath("pi_theta_withinclade_by_complex.tsv"), "w") as o:
    o.write("Clade\tGroup\tClass\tS_sites\tpi_theta\tCI_low\tCI_high\n")
    for clade in ("north", "south"):
        for cx in ("I", "III", "IV", "V"):
            for cls, lab in CLASSES:
                S, pe, lo, hi = subset(clade, cx=cx, cls=cls)
                o.write("%s\tComplex_%s\t%s\t%d\t%s\t%s\t%s\n" %
                        (clade.capitalize(), cx, lab, S, fmt(pe), fmt(lo), fmt(hi)))
        for cls, lab in CLASSES:
            S, pe, lo, hi = subset(clade, cls=cls)
            o.write("%s\tOverall\t%s\t%d\t%s\t%s\t%s\n" %
                    (clade.capitalize(), lab, S, fmt(pe), fmt(lo), fmt(hi)))

print("03_pi_theta_withinclade.py: segregating sites north=%d south=%d (a_n %.3f / %.3f)" %
      (len(data["north"]), len(data["south"]), A_N["north"], A_N["south"]))
print("  wrote pi_theta_withinclade_by_gene.tsv and pi_theta_withinclade_by_complex.tsv")

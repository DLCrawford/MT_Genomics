"""
common.py — shared constants and helpers for the pNps selection-statistics pipeline.

Fundulus heteroclitus mtDNA, 13 protein-coding genes, 141 individuals.
All scripts import from here so gene->complex mapping, sample sizes, the
pi/theta estimator, and the bootstrap are defined in exactly one place.
"""
import csv, os, random

# ---- paths (scripts live in pNps/, inputs live in ../vcf/) ----
HERE = os.path.dirname(os.path.abspath(__file__))
VCF  = os.path.normpath(os.path.join(HERE, "..", "vcf"))
OUT  = HERE
def vpath(name): return os.path.join(VCF, name)
def opath(name): return os.path.join(OUT, name)

# ---- OXPHOS complex <-> gene mapping ----
# Complex II has NO mtDNA-encoded subunits and is therefore absent by definition.
COMPLEXES = {
    "I":   ["ND1", "ND2", "ND3", "ND4", "ND4L", "ND5", "ND6"],
    "III": ["CYTB"],
    "IV":  ["COX1", "COX2", "COX3"],
    "V":   ["ATP6", "ATP8"],
}
GENE2CX = {g: c for c, gs in COMPLEXES.items() for g in gs}
GENES   = ["ND1","ND2","ND3","ND4","ND4L","ND5","ND6","CYTB","COX1","COX2","COX3","ATP6","ATP8"]

# ---- sample sizes ----
# Whole panel = 141 mt haplotypes. Clade sizes from per-site called counts in
# pi_by_clade_persite.tsv (mode of north_n_called / south_n_called).
N_PANEL = 141
N_NORTH = 77
N_SOUTH = 63

# ---- bootstrap ----
B = 5000   # bootstrap replicates (site resampling)

def harmonic(n):
    """Watterson's a_n = sum_{i=1}^{n-1} 1/i."""
    return sum(1.0 / i for i in range(1, n))

A_N = {"panel": harmonic(N_PANEL), "north": harmonic(N_NORTH), "south": harmonic(N_SOUTH)}

def pi_theta(pis, a_n):
    """
    pi/theta for a set of biallelic segregating sites.
      pi     = sum of per-site heterozygosities (pi_site)
      theta  = S / a_n            (Watterson, S = number of segregating sites)
      pi/theta = a_n * sum(pi_site) / S      (gene length L cancels)
    Returns None if S == 0.
    """
    S = len(pis)
    return (a_n * sum(pis) / S) if S else None

def boot_ci(pis, a_n, B=B, alpha=0.05, seed=None):
    """Percentile bootstrap CI for pi/theta by resampling sites with replacement."""
    S = len(pis)
    if S < 2:
        return (None, None)
    rng = random.Random(seed)
    vals = sorted(a_n * sum(pis[rng.randrange(S)] for _ in range(S)) / S for _ in range(B))
    lo = vals[int((alpha/2) * B)]
    hi = vals[int((1 - alpha/2) * B)]
    return (lo, hi)

def pi_site_biallelic(AC, n):
    """Per-site nucleotide diversity for a biallelic site (unbiased, sampling without replacement)."""
    return 2.0 * AC * (n - AC) / (n * (n - 1))

def fmt(x, nd=3):
    return "NA" if x is None else ("%.*f" % (nd, x))

def read_tsv(path, ncols=None, encoding="latin-1"):
    """Yield dict rows. ncols truncates to the first ncols columns (some inputs have trailing summary columns)."""
    with open(path, encoding=encoding) as f:
        r = csv.reader(f, delimiter="\t")
        hdr = next(r)
        cols = hdr[:ncols] if ncols else hdr
        for row in r:
            if ncols and len(row) < ncols:
                continue
            yield dict(zip(cols, row[:ncols] if ncols else row))

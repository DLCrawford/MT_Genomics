#!/usr/bin/env python3
"""
35_clade_ns_reref.py — Reference-placement demonstration for the clade NS/S debate.

Recomputes the per-individual variant burden (Total / SYN / NonSYN, with NS/Total
and NS/SYN) and the per-clade NS/S summary under FOUR reference haplotypes, to show
that the apparent between-clade NS/S difference is a reference-placement artifact
(supersedes the "adaptive North" reading — see CLAUDE.md session 22).

References
---------
  N_ref     published GenBank NC_012312.1 (North-like)  -> baseline; reproduces
            the existing vs-reference per-individual table (script 18).
  112_MT    a clean South-clade haplotype (6 SNPs from the S consensus).
  77_MT     the admixed "ambiguous" fish (190 ALTs vs N_ref) — a real, near-
            intermediate haplotype that leans South but carries North alleles
            at ~30 South-divergent sites.
  random    a synthetic "truly intermediate" reference: at each of the 927
            variable sites the reference allele is drawn 0/1 with p=0.5
            (seeded). One seeded draw is used for the per-individual table;
            the per-clade summary additionally reports the mean +/- SD over
            N_DRAWS independent draws (the expectation is what is meaningfully
            "intermediate"; a single draw is an arbitrary synthetic haplotype).

Re-referencing logic
--------------------
For a biallelic SNP, the SYN/NS impact class is symmetric (REF<->ALT both
missense, or both synonymous), so SnpEff's ANN[0] class is reference-invariant
and is reused. Only the POLARITY changes: an individual "differs" from the
reference at a site when its genotype != the reference's allele at that site.
Missing focal genotypes are skipped (as in script 18). A missing genotype in a
sample-based reference (112/77) is treated as REF (0).

Inputs
------
  vcf/141_MT_variants.vcf.gz                       canonical 141-sample VCF
  vcf/pi_by_clade_persite.tsv.membership.tsv       clade membership (CRLF)

Outputs
-------
  vcf/reref_per_individual_<REF>.tsv     (4 files) per-individual burden
  vcf/reref_per_clade_summary.tsv        per-clade mean NS/SYN/(NS/SYN) x ref
  vcf/reref_summary.txt                  human-readable summary
"""
from __future__ import annotations
import gzip, random
from pathlib import Path
from statistics import mean, pstdev

ROOT = Path(__file__).resolve().parent.parent
VCF  = ROOT / "vcf/141_MT_variants.vcf.gz"
MEMB = ROOT / "vcf/pi_by_clade_persite.tsv.membership.tsv"
OUT  = ROOT / "vcf"

SEED, N_DRAWS = 42, 1000

SYN_EFF = {"synonymous_variant", "stop_retained_variant"}
NS_EFF  = {"missense_variant", "stop_gained", "stop_lost",
           "start_lost", "initiator_codon_variant"}

def classify(eff: str) -> str:
    if not eff: return "Other"
    p = eff.split("&")
    if any(x in NS_EFF for x in p): return "NS"
    if any(x in SYN_EFF for x in p): return "SYN"
    return "Other"

def ann_effect(info: str) -> str:
    for e in info.split(";"):
        if e.startswith("ANN="):
            return e[4:].split(",", 1)[0].split("|")[1]
    return ""

# ── load clade membership (CRLF) ────────────────────────────────────────────
clade = {}
for ln in MEMB.read_text().splitlines()[1:]:
    f = ln.strip().split("\t")
    if len(f) >= 3:
        clade[f[0]] = f[2].strip()

# ── parse VCF ───────────────────────────────────────────────────────────────
samples: list[str] = []
site_cls: list[str] = []
gt: list[list[str]] = []          # gt[site][sample] in {'0','1','.'}
with gzip.open(VCF, "rt") as fh:
    for line in fh:
        if line.startswith("##"): continue
        if line.startswith("#CHROM"):
            samples = line.rstrip("\n").split("\t")[9:]; continue
        f = line.rstrip("\n").split("\t")
        site_cls.append(classify(ann_effect(f[7])))
        gt.append([v.split(":", 1)[0] for v in f[9:]])

NS_SITES = len(samples); n_sites = len(site_cls)
idx = {s: i for i, s in enumerate(samples)}
north = [i for i, s in enumerate(samples) if clade.get(s) == "north"]
south = [i for i, s in enumerate(samples) if clade.get(s) == "south"]

# ── reference allele vectors ────────────────────────────────────────────────
def ref_from_sample(name: str) -> list[str]:
    j = idx[name]
    return [gt[k][j] if gt[k][j] in ("0", "1") else "0" for k in range(n_sites)]

def ref_random(seed: int) -> list[str]:
    rng = random.Random(seed)
    return ["1" if rng.random() < 0.5 else "0" for _ in range(n_sites)]

REFS = {
    "N_ref":  ["0"] * n_sites,          # GenBank NC_012312.1 = all-REF
    "112_MT": ref_from_sample("112_MT"),
    "77_MT":  ref_from_sample("77_MT"),
    "random": ref_random(SEED),
}

# ── per-individual counts under a given reference ───────────────────────────
def counts_under(ref: list[str]):
    """return {sample: (SYN, NS, Other)} differences vs ref."""
    out = {}
    for sj, s in enumerate(samples):
        cS = cN = cO = 0
        for k in range(n_sites):
            g = gt[k][sj]
            if g not in ("0", "1"):  # missing focal call -> skip
                continue
            if g != ref[k]:
                c = site_cls[k]
                if c == "SYN": cS += 1
                elif c == "NS": cN += 1
                else: cO += 1
        out[s] = (cS, cN, cO)
    return out

def write_per_individual(refname: str, ref: list[str]):
    c = counts_under(ref)
    path = OUT / f"reref_per_individual_{refname}.tsv"
    with open(path, "w") as w:
        w.write("IND\tIND_name\tClade\tTotal\tSYN\tNonSYN\tOther\tNS_over_Total\tNS_over_SYN\n")
        # sort by Total descending (like script 18)
        for s in sorted(samples, key=lambda x: -(sum(c[x]))):
            cS, cN, cO = c[s]; tot = cS + cN + cO
            rid = s.split("_")[0]
            nst = (cN / tot) if tot else 0.0
            nss = (cN / cS) if cS else 0.0
            w.write(f"{rid}\t{s}\t{clade.get(s,'NA')}\t{tot}\t{cS}\t{cN}\t{cO}\t{nst:.4f}\t{nss:.6f}\n")
    return c, path

# ── per-clade summary (mean per individual) ─────────────────────────────────
def clade_means(counts, group):
    syn = mean(counts[samples[j]][0] for j in group)
    ns  = mean(counts[samples[j]][1] for j in group)
    return ns, syn, (ns / syn if syn else float("nan"))

summary_rows = []
all_counts = {}
for refname, ref in REFS.items():
    c, p = write_per_individual(refname, ref)
    all_counts[refname] = c
    for cl_name, grp in [("North", north), ("South", south)]:
        ns, syn, r = clade_means(c, grp)
        summary_rows.append((refname, cl_name, len(grp), ns, syn, r))

# random reference: average over N_DRAWS draws for a stable "intermediate" value
rnd_stats = {}
for cl_name, grp in [("North", north), ("South", south)]:
    rs = []
    for d in range(N_DRAWS):
        c = counts_under(ref_random(SEED + d))
        ns, syn, r = clade_means(c, grp)
        rs.append((ns, syn, r))
    rnd_stats[cl_name] = (
        mean(x[0] for x in rs), mean(x[1] for x in rs),
        mean(x[2] for x in rs), pstdev(x[2] for x in rs),
    )

# ── write per-clade summary tsv ─────────────────────────────────────────────
with open(OUT / "reref_per_clade_summary.tsv", "w") as w:
    w.write("Reference\tClade\tn\tmean_NS\tmean_SYN\tNS_over_SYN\n")
    for refname, cl, n, ns, syn, r in summary_rows:
        w.write(f"{refname}\t{cl}\t{n}\t{ns:.3f}\t{syn:.3f}\t{r:.4f}\n")

# ── human-readable summary ──────────────────────────────────────────────────
L = []
L.append("Reference-placement effect on per-clade NS/S (mean per individual)")
L.append("=" * 66)
L.append(f"Panel: 141 samples (North={len(north)}, South={len(south)}, "
         f"ambiguous=1 [77_MT]). Sites: {n_sites} (SYN={site_cls.count('SYN')}, "
         f"NS={site_cls.count('NS')}, Other={site_cls.count('Other')}).")
L.append(f"Site-class ratio NS/SYN = {site_cls.count('NS')}/{site_cls.count('SYN')} "
         f"= {site_cls.count('NS')/site_cls.count('SYN'):.4f} "
         f"(expected NS/S under an unbiased reference).")
L.append("")
hdr = f"{'Reference':<10}{'Clade':<7}{'n':>4}{'mean_NS':>10}{'mean_SYN':>10}{'NS/SYN':>9}"
L.append(hdr); L.append("-" * len(hdr))
for refname, cl, n, ns, syn, r in summary_rows:
    L.append(f"{refname:<10}{cl:<7}{n:>4}{ns:>10.2f}{syn:>10.2f}{r:>9.4f}")
L.append("")
L.append(f"random reference: mean +/- SD of NS/SYN over {N_DRAWS} seeded draws")
for cl in ("North", "South"):
    ns, syn, rmean, rsd = rnd_stats[cl]
    L.append(f"  {cl:<6} NS/SYN = {rmean:.4f} +/- {rsd:.4f}  "
             f"(mean_NS={ns:.1f}, mean_SYN={syn:.1f})")
L.append("")
L.append("Reading: N_ref makes North look high-NS/S (polymorphism regime) and")
L.append("South low (divergence regime); a South reference (112/77) flips it;")
L.append("the intermediate/random reference collapses both clades onto the")
L.append("site-class ratio ~0.215. The clade NS/S 'difference' is set by where")
L.append("the reference sits on the N-S axis, not by selection.")
txt = "\n".join(L)
(OUT / "reref_summary.txt").write_text(txt + "\n")
print(txt)
print(f"\nWrote: reref_per_individual_{{{','.join(REFS)}}}.tsv, "
      f"reref_per_clade_summary.tsv, reref_summary.txt  -> {OUT}")

#!/usr/bin/env python3
"""
scripts/14_hp_from_pileup.py

Heteroplasmy detection from a per-CDS-position pileup (stage 13), independent
of the panel variant call set. This complements scripts/09_heteroplasmy_report.py
(which runs on MT_DP_AD_141.txt — variant-only positions) by exposing the
"private ALT-Hp" category that the variant-only file cannot represent.

INPUT
-----
vcf/pileup_cds_141.vcf.gz   bcftools-mpileup -a AD,DP over docs/mito_protein_coding.bed
                            across 141 BAMs (70, 125 excluded). One record per
                            CDS position; per-sample FORMAT/AD lists counts for
                            REF + each observed alt + <*> (the catch-all).

THRESHOLDS
----------
DP_MIN     = 20    (per-cell minimum read depth to consider that cell at all)
MAJOR_FRAC = 0.70  (haplotype/major call threshold — same as stage 08)
HP_MIN     = 0.10  (heteroplasmy lower bound)
HP_AD_MIN  = 4     (minimum reads supporting the Hp allele; at DP=20 this is AF=0.20,
                   well above the noise floor for a 0.1-cutoff at low DP, and at
                   DP=200 the AF is just 0.02 so DP>20 is the binding constraint)

CLASSIFICATION
--------------
For each (Individual, POS), evaluate every observed allele's AF = AD/DP:
  - "major" if AF >= MAJOR_FRAC (this is the haplotype call)
  - "Hp"    if HP_MIN <= AF < MAJOR_FRAC  AND  AD >= HP_AD_MIN

For each Hp event we further classify the Hp allele:
  - REF_Hp           : Hp allele IS the reference base at this position
  - shared_alt_Hp    : Hp allele IS an ALT base AND is a major call in
                       >=1 other individual at the same POS
  - private_alt_Hp   : Hp allele IS an ALT base AND is NOT a major call in
                       any other individual at this POS  (the category
                       09 could not produce because MT_DP_AD_141.txt only
                       includes panel-variant ALTs)

OUTPUT
------
vcf/heteroplasmy_pileup_events.tsv          one row per Hp event
vcf/heteroplasmy_pileup_per_site.tsv        per-CDS-position Hp summary
vcf/heteroplasmy_pileup_per_individual.tsv  per-individual Hp summary
vcf/heteroplasmy_pileup_summary.txt         human-readable headline counts
"""

import os
import sys
import subprocess
from collections import defaultdict
import pandas as pd

# ---------------------------------------------------------------- paths/config
VCF_DIR = os.environ.get(
    "VCF_DIR",
    "/Users/douglas_crawford/Projects/MT_Genomics_Cl_Ap2026/MT_Genomics2/vcf",
)
IN_VCF   = f"{VCF_DIR}/pileup_cds_141.vcf.gz"

OUT_EV   = f"{VCF_DIR}/heteroplasmy_pileup_events.tsv"
OUT_SITE = f"{VCF_DIR}/heteroplasmy_pileup_per_site.tsv"
OUT_IND  = f"{VCF_DIR}/heteroplasmy_pileup_per_individual.tsv"
OUT_SUM  = f"{VCF_DIR}/heteroplasmy_pileup_summary.txt"

DP_MIN     = 20
MAJOR_FRAC = 0.70
HP_MIN     = 0.10
HP_AD_MIN  = 4

# ------------------------------------------------------- parse pileup with bcftools query
# bcftools mpileup emits alleles as REF + observed alts + <*>. We pull:
#   CHROM, POS, REF, ALT, [SAMPLE, DP, AD]
# ALT is comma-joined (e.g., "A,<*>" at an invariant T site with one read of A)
# AD is comma-joined per sample, in order [REF, alt1, alt2, ..., <*>]
print(f"[14] Streaming pileup VCF: {IN_VCF}", file=sys.stderr)
cmd = [
    "bcftools", "query",
    "-f", "%CHROM\t%POS\t%REF\t%ALT[\t%SAMPLE=%DP=%AD]\n",
    IN_VCF,
]
proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True)

major_calls_by_site = defaultdict(set)       # POS -> set of alleles called major in >=1 ind
records_per_site    = defaultdict(int)       # POS -> #cells with DP>=DP_MIN
events              = []                     # list of dicts, one per Hp event
sample_set          = set()
n_lines             = 0

for line in proc.stdout:
    n_lines += 1
    parts = line.rstrip("\n").split("\t")
    chrom, pos, ref, alt = parts[0], int(parts[1]), parts[2], parts[3]
    # alleles in AD order: [REF, *ALT.split(","), where <*> is the catch-all]
    alts = [a for a in alt.split(",")]
    alleles_in_AD_order = [ref] + alts        # length = len(AD)
    # iterate per-sample triples
    for cell in parts[4:]:
        # SAMPLE=DP=AD  (AD is comma-joined)
        try:
            sname, dp_s, ad_s = cell.split("=", 2)
        except ValueError:
            continue
        sample_set.add(sname)
        try:
            dp = int(dp_s)
        except ValueError:
            continue
        if dp < DP_MIN:
            continue
        records_per_site[pos] += 1
        ads = []
        for x in ad_s.split(","):
            try:
                ads.append(int(x))
            except ValueError:
                ads.append(0)
        # zip allele -> ad, dropping the <*> catch-all
        cell_counts = {}
        for allele, ad in zip(alleles_in_AD_order, ads):
            if allele == "<*>":
                continue
            cell_counts[allele] = cell_counts.get(allele, 0) + ad
        # major
        major = None
        for a, ad in cell_counts.items():
            if ad / dp >= MAJOR_FRAC:
                major = a
                major_calls_by_site[pos].add(a)
                break
        # Hp candidates
        for a, ad in cell_counts.items():
            if a == major:
                continue
            af = ad / dp
            if af >= HP_MIN and ad >= HP_AD_MIN and af < MAJOR_FRAC:
                events.append(dict(
                    Individual=sname, Chromo=chrom, POS=pos, DP=dp, REF=ref,
                    Major=(major if major else "NONE"),
                    Hp_allele=a, Hp_AD=ad, Hp_AF=round(af, 4),
                    Hp_is_REF=(a == ref),
                ))
    if n_lines % 2000 == 0:
        print(f"      {n_lines:,} positions; {len(events):,} events so far",
              file=sys.stderr)

proc.wait()
if proc.returncode != 0:
    sys.exit(f"bcftools query failed (rc={proc.returncode})")
print(f"[14] Done streaming. positions={n_lines:,}  samples={len(sample_set)}  "
      f"events={len(events):,}", file=sys.stderr)

# ----------------------------------------------------- classify shared vs private
def classify(row):
    if row["Hp_is_REF"]:
        return "REF_Hp"
    if row["Hp_allele"] in major_calls_by_site.get(row["POS"], set()):
        return "shared_alt_Hp"
    return "private_alt_Hp"

ev = pd.DataFrame(events)
if ev.empty:
    sys.exit("No Hp events found — check thresholds / input.")
ev["Hp_class"] = ev.apply(classify, axis=1)

# ----------------------------------------------------- aggregates
n_inds_hp        = ev["Individual"].nunique()
n_sites_hp       = ev["POS"].nunique()
ind_counts       = ev.groupby("Individual").size()
n_inds_2plus     = (ind_counts >= 2).sum()

n_ref            = int((ev["Hp_class"] == "REF_Hp").sum())
n_private_alt    = int((ev["Hp_class"] == "private_alt_Hp").sum())
n_shared_alt     = int((ev["Hp_class"] == "shared_alt_Hp").sum())

per_site = (ev.groupby(["Chromo", "POS", "REF"])
              .agg(n_carriers=("Individual", "nunique"),
                   n_events=("Individual", "size"),
                   hp_alleles=("Hp_allele", lambda s: ",".join(sorted(set(s)))),
                   hp_classes=("Hp_class", lambda s: ",".join(sorted(set(s)))))
              .reset_index()
              .sort_values(["n_carriers", "POS"], ascending=[False, True]))

per_ind = (ev.groupby("Individual")
             .agg(n_hp_sites=("POS", "nunique"),
                  n_hp_events=("POS", "size"),
                  n_ref_hp=("Hp_is_REF", "sum"))
             .reset_index()
             .sort_values("n_hp_sites", ascending=False))

# ----------------------------------------------------- write outputs
ev.to_csv(OUT_EV, sep="\t", index=False)
per_site.to_csv(OUT_SITE, sep="\t", index=False)
per_ind.to_csv(OUT_IND, sep="\t", index=False)

summary = f"""Heteroplasmy report — pileup-based (stage 14)
Source     : {IN_VCF}
Thresholds : DP >= {DP_MIN}, Hp = {HP_MIN} <= AD/DP < {MAJOR_FRAC}, Hp AD >= {HP_AD_MIN}

Coverage:
  CDS positions in pileup     : {n_lines:,}
  Samples observed            : {len(sample_set)}
  Total Hp events             : {len(ev):,}

Headline:
  Q1  individuals with >=1 Hp        : {n_inds_hp}
  Q2  sites with >=1 Hp              : {n_sites_hp}
  Q3  individuals with >=2 Hp events : {n_inds_2plus}

Hp allele classification:
  Q4  Hp == REF                      : {n_ref}
  Q5  ALT Hp NOT major in any other  : {n_private_alt}  (private)
  Q6  ALT Hp IS  major in some other : {n_shared_alt}   (shared/transmitted)
  Total                              : {n_ref + n_private_alt + n_shared_alt}

Top 10 individuals by Hp-site count:
{per_ind.head(10).to_string(index=False)}

Top 10 sites by Hp-carrier count:
{per_site.head(10).to_string(index=False)}
"""

with open(OUT_SUM, "w") as f:
    f.write(summary)

print(summary)
print(f"\n[14] Wrote:", file=sys.stderr)
for p in (OUT_EV, OUT_SITE, OUT_IND, OUT_SUM):
    print(f"      {p}", file=sys.stderr)

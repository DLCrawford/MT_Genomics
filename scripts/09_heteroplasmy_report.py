#!/usr/bin/env python3
"""
Heteroplasmy analysis from MT_DP_AD_141.txt

Definitions:
    - At each (individual, POS), the allele with AD/DP >= 0.7 is the "major"
      (haplotype) call. Any other allele with AD/DP >= 0.1 (and < 0.7 by
      construction) is a heteroplasmy (Hp) event.
    - Each input row is one (ind, POS, ALT) record after norm-split, so for
      multi-allelic sites we collapse 2+ rows back into one per (ind, POS),
      with REF count + each ALT count.

Questions answered:
  (Q1) How many individuals carry >=1 Hp event?
  (Q2) How many sites carry >=1 Hp event (any individual)?
  (Q3) How many individuals carry >=2 Hp events?
  (Q4) How many Hp events are at the REF allele (i.e., the individual's
       major call is some ALT, the minor heteroplasmy is REF)?
  (Q5) How many ALT-Hp events involve an ALT that is NOT a major call
       in any other individual (i.e., private ALT-Hp)?
  (Q6) How many ALT-Hp events involve an ALT that IS a major call in
       >=1 other individual (i.e., shared/transmitted)?
"""

import sys
import pandas as pd
from collections import defaultdict

import os
# Allow override via env so the same script runs under the sandboxed bash
# mount path AND under macOS-native Python.
VCF_DIR = os.environ.get(
    "VCF_DIR",
    "/Users/douglas_crawford/Projects/MT_Genomics_Cl_Ap2026/MT_Genomics2/vcf",
)
IN            = f"{VCF_DIR}/MT_DP_AD_141.txt"
OUT_PER_SITE  = f"{VCF_DIR}/heteroplasmy_per_site.tsv"
OUT_PER_IND   = f"{VCF_DIR}/heteroplasmy_per_individual.tsv"
OUT_EVENTS    = f"{VCF_DIR}/heteroplasmy_events.tsv"
OUT_SUMMARY   = f"{VCF_DIR}/heteroplasmy_summary.txt"

MAJOR_FRAC = 0.70
HP_MIN     = 0.10

print(f"Loading {IN} ...", file=sys.stderr)
df = pd.read_csv(
    IN, sep="\t",
    dtype={"Individual": str, "Position": str, "REF": str, "ALT1": str,
           "Chromo": str, "Gene": str, "Effect": str, "Variant_types": str},
)
print(f"Loaded {len(df):,} rows", file=sys.stderr)

# Coerce numerics; some ADaltN are "NA" strings
for col in ["DP", "ADref", "ADalt1", "ADalt2", "ADalt3", "ADalt4"]:
    df[col] = pd.to_numeric(df[col], errors="coerce")

n_individuals_in_file = df["Individual"].nunique()
n_positions_in_file   = df["Position"].nunique()
print(f"Unique individuals: {n_individuals_in_file}", file=sys.stderr)
print(f"Unique positions:   {n_positions_in_file}", file=sys.stderr)

# --------------------------------------------------------------------------
# Collapse split-ALT rows back to one record per (Individual, Position)
# Each (ind, POS) has one ADref / DP value (same on every row), and each
# row contributes the count for that row's ALT1.
# ADalt2/3/4 are NA throughout (norm-split file), so we only need ADalt1
# from each row, indexed by that row's ALT base.
# --------------------------------------------------------------------------
print("Collapsing split-ALT rows per (individual, POS)...", file=sys.stderr)

records = []  # one dict per (ind, POS)
for (ind, pos), g in df.groupby(["Individual", "Position"], sort=False):
    g0 = g.iloc[0]
    dp     = int(g0["DP"])
    adref  = int(g0["ADref"]) if pd.notna(g0["ADref"]) else 0
    ref    = g0["REF"]
    chrom  = g0["Chromo"]
    posnum = int(g0["Postion #"])
    # collect per-ALT counts; ALT1 column is the row's ALT base
    alts = {}
    for _, r in g.iterrows():
        a = r["ALT1"]
        ad = int(r["ADalt1"]) if pd.notna(r["ADalt1"]) else 0
        alts[a] = ad
    records.append(dict(
        Individual=ind, Position=pos, Chromo=chrom, POS=posnum,
        REF=ref, DP=dp, ADref=adref, ALT_counts=alts,
    ))

print(f"Collapsed to {len(records):,} (ind, POS) records", file=sys.stderr)

# --------------------------------------------------------------------------
# Classify each allele at each (ind, POS) as major / heteroplasmy / noise
# --------------------------------------------------------------------------
events = []  # one row per Hp event (ind, POS, allele)
major_calls_by_site = defaultdict(set)  # POS -> set of alleles that are MAJOR in >=1 ind

for rec in records:
    dp = rec["DP"]
    if dp <= 0:
        continue
    pos = rec["Position"]
    all_alleles = {rec["REF"]: rec["ADref"], **rec["ALT_counts"]}

    # find major: the unique allele with AD/DP >= 0.7
    major_allele = None
    for a, ad in all_alleles.items():
        if ad / dp >= MAJOR_FRAC:
            major_allele = a
            major_calls_by_site[pos].add(a)
            break

    # heteroplasmy: any other allele with AD/DP >= 0.1
    for a, ad in all_alleles.items():
        if a == major_allele:
            continue
        af = ad / dp
        if af >= HP_MIN:  # AF < 0.7 here by construction (only one allele can be >=0.7)
            events.append(dict(
                Individual=rec["Individual"],
                Position=rec["Position"],
                Chromo=rec["Chromo"],
                POS=rec["POS"],
                DP=dp,
                REF=rec["REF"],
                Major=major_allele if major_allele else "NONE",
                Hp_allele=a,
                Hp_AD=ad,
                Hp_AF=round(af, 4),
                Hp_is_REF=(a == rec["REF"]),
            ))

events_df = pd.DataFrame(events)
print(f"Hp events: {len(events_df):,}", file=sys.stderr)

# --------------------------------------------------------------------------
# Classify each ALT-Hp event as private vs shared (transmitted)
# private = Hp ALT is NOT a major call in any other individual at any site
# shared  = Hp ALT IS  a major call in >=1 other individual at the SAME site
# --------------------------------------------------------------------------
def shared_status(row):
    if row["Hp_is_REF"]:
        return "REF_Hp"
    pos = row["Position"]
    a = row["Hp_allele"]
    majors_here = major_calls_by_site.get(pos, set())
    if a in majors_here:
        return "shared_with_major"
    return "private_alt_Hp"

events_df["Hp_class"] = events_df.apply(shared_status, axis=1)

# --------------------------------------------------------------------------
# Aggregate answers
# --------------------------------------------------------------------------
n_inds_with_hp        = events_df["Individual"].nunique()
n_sites_with_hp       = events_df["Position"].nunique()
ind_event_counts      = events_df.groupby("Individual").size()
n_inds_with_2plus_hp  = (ind_event_counts >= 2).sum()

n_ref_hp              = int((events_df["Hp_class"] == "REF_Hp").sum())
n_private_alt_hp      = int((events_df["Hp_class"] == "private_alt_Hp").sum())
n_shared_alt_hp       = int((events_df["Hp_class"] == "shared_with_major").sum())

# per-site Hp summary
per_site = (events_df
            .groupby(["Position", "POS", "Chromo", "REF"])
            .agg(n_carriers=("Individual", "nunique"),
                 n_events=("Individual", "size"),
                 hp_alleles=("Hp_allele", lambda s: ",".join(sorted(set(s)))),
                 hp_classes=("Hp_class", lambda s: ",".join(sorted(set(s)))))
            .reset_index()
            .sort_values(["n_carriers", "POS"], ascending=[False, True]))

# per-individual Hp summary
per_ind = (events_df
           .groupby("Individual")
           .agg(n_hp_sites=("Position", "nunique"),
                n_hp_events=("Position", "size"),
                n_ref_hp=("Hp_is_REF", "sum"))
           .reset_index()
           .sort_values("n_hp_sites", ascending=False))

# --------------------------------------------------------------------------
# Write outputs
# --------------------------------------------------------------------------
events_df.to_csv(OUT_EVENTS, sep="\t", index=False)
per_site.to_csv(OUT_PER_SITE, sep="\t", index=False)
per_ind.to_csv(OUT_PER_IND, sep="\t", index=False)

summary = f"""Heteroplasmy report — MT_DP_AD_141.txt
Thresholds: major (haplotype) call = AD/DP >= {MAJOR_FRAC}; Hp = {HP_MIN} <= AD/DP < {MAJOR_FRAC}

Input
  individuals in file        : {n_individuals_in_file}
  unique positions in file   : {n_positions_in_file}
  total (ind,POS) records    : {len(records):,}
  total Hp events            : {len(events_df):,}

Q1  individuals with >=1 Hp event   : {n_inds_with_hp}
Q2  sites with >=1 Hp event         : {n_sites_with_hp}
Q3  individuals with >=2 Hp events  : {n_inds_with_2plus_hp}

Hp-event allele classification
  Q4  Hp == REF (major call is ALT) : {n_ref_hp}
  Q5  ALT Hp NOT a major elsewhere  : {n_private_alt_hp}  (private)
  Q6  ALT Hp IS a major elsewhere   : {n_shared_alt_hp}   (shared / transmitted)
                                       --------
                                       {n_ref_hp + n_private_alt_hp + n_shared_alt_hp}  (matches total events)

Top 10 individuals by Hp-site count:
"""
top10 = per_ind.head(10).to_string(index=False)
summary += top10 + "\n\nTop 10 sites by Hp-carrier count:\n"
summary += per_site.head(10).to_string(index=False) + "\n"

with open(OUT_SUMMARY, "w") as f:
    f.write(summary)

print(summary)
print(f"\nWrote: {OUT_EVENTS}", file=sys.stderr)
print(f"Wrote: {OUT_PER_SITE}", file=sys.stderr)
print(f"Wrote: {OUT_PER_IND}", file=sys.stderr)
print(f"Wrote: {OUT_SUMMARY}", file=sys.stderr)

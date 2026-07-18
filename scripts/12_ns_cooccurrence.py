#!/usr/bin/env python3
"""
scripts/12_ns_cooccurrence.py

Find nonsynonymous (MODERATE-impact missense) sites that share identical
per-sample call vectors — i.e., sites in perfect linkage disequilibrium
across the panel.

INPUT
-----
    vcf/haplotype_matrix_nonsyn.csv   (preferred; from scripts/11_haplotypes_nonsyn.py)
    vcf/haplotype_matrix.csv          (fallback; filter to missense rows here)

WHAT "ALWAYS OCCUR TOGETHER" MEANS
----------------------------------
For each missense row in the matrix, build a tuple of its 0/1 calls across
the N "real" sample columns (names matching the pattern `\\d+_MT`).
Rows with *identical* tuples carry the same alt in exactly the same set
of samples and absent in the same set — i.e., perfect LD across the
sampled population.

INTERPRETATION GUIDE
--------------------
The grouped rows fall into three classes; the script labels each group
in the output.

1. **Fixed reference-divergence sites** (every sample carries the alt):
   Not LD per se — these are positions where the reference differs from
   the entire panel. Flag for the methods write-up; arguably these
   shouldn't enter dN/dS because they're "polymorphism" only against an
   outgroup reference.

2. **Genuine co-segregating haplogroups** (≥2 sites, intermediate carrier
   count, ≥2 carriers): a real maternal-lineage signal. The session-15
   analysis on the *F. heteroclitus* panel found one major group: 7
   missense sites in 4 genes (ATP6, ND1, ND2, ND5) co-segregating in
   64 of 141 samples (~45 %). Also six smaller groups at 2–4 carriers.

3. **Singleton aggregation artifacts** (multiple sites all at carrier
   count = 1, signature = a single 1 at the same sample): these are
   *one fish's collection of private NS mutations*, not population-level
   LD. Filter them out by `--min-carriers 2` if uninterested.

INPUTS NOT NEEDED
-----------------
The script reads only the matrix CSV (already CDS-restricted, already
SnpEff-annotated, already split by `norm -m -any`). It does NOT touch
the VCF or the FASTA.

OUTPUTS
-------
    vcf/ns_cooccurrence_groups.tsv
        One row per multi-site group, columns:
        group_id | n_sites | n_carriers | category | site_list
        category ∈ {fixed_ref_divergence, haplogroup, singleton_artifact}

    Stdout: a pretty-printed summary listing each group.

Run
---
    cd ~/Projects/MT_Genomics_Cl_Ap2026/MT_Genomics2
    python scripts/12_ns_cooccurrence.py [--min-carriers 1]
"""

from __future__ import annotations

import argparse
import re
import sys
from collections import defaultdict
from pathlib import Path

import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parent.parent
VCF_DIR = PROJECT_ROOT / "vcf"

MATRIX_NONSYN = VCF_DIR / "haplotype_matrix_nonsyn.csv"
MATRIX_ALL    = VCF_DIR / "haplotype_matrix.csv"
OUT_TSV       = VCF_DIR / "ns_cooccurrence_groups.tsv"

META_COLS = {"CHROM", "POS", "REF", "ALT", "Gene", "Variant_type", "Impact"}
SAMPLE_RE = re.compile(r"^\d+_MT$")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--min-carriers", type=int, default=1,
                    help="Skip groups whose sites are carried by fewer than this many samples "
                         "(useful for hiding singleton-artifact groups). Default 1 = show all.")
    args = ap.parse_args()

    if MATRIX_NONSYN.exists():
        src = MATRIX_NONSYN
        df = pd.read_csv(src)
        print(f"Using nonsyn matrix: {src}  ({len(df)} rows)")
    elif MATRIX_ALL.exists():
        src = MATRIX_ALL
        df = pd.read_csv(src)
        df = df[df["Variant_type"].str.contains("missense", na=False)].reset_index(drop=True)
        print(f"Using ALL-CDS matrix filtered to missense: {src}  ({len(df)} rows)")
    else:
        sys.exit(f"ERROR: neither {MATRIX_NONSYN} nor {MATRIX_ALL} found.\n"
                 "  Run scripts/11_haplotypes_nonsyn.py or scripts/08_call_haplotypes.py first.")

    sample_cols = [c for c in df.columns if SAMPLE_RE.match(c)]
    ignored = [c for c in df.columns if c not in META_COLS and c not in sample_cols]
    if ignored:
        print(f"Ignoring non-sample columns: {ignored}")
    print(f"Real sample columns ({len(sample_cols)}): {sample_cols[0]} ... {sample_cols[-1]}")

    # Coerce sample cells to 0/1 ints (handle stray strings)
    for c in sample_cols:
        df[c] = pd.to_numeric(df[c], errors="coerce").fillna(0).astype(int)

    # Drop never-called rows
    n_before = len(df)
    df = df[df[sample_cols].sum(axis=1) > 0].reset_index(drop=True)
    print(f"Rows with ≥1 carrier: {len(df)}  (dropped {n_before - len(df)} never-called)")

    # Build per-row 141-element signature
    groups: dict[tuple, list] = defaultdict(list)
    for _, row in df.iterrows():
        sig = tuple(int(row[c]) for c in sample_cols)
        groups[sig].append({
            "POS":   int(row["POS"]),
            "REF":   row["REF"],
            "ALT":   row["ALT"],
            "Gene":  row["Gene"],
            "n_car": int(sum(sig)),
        })

    multi = {sig: sites for sig, sites in groups.items() if len(sites) > 1}
    print(f"\nUnique call patterns:               {len(groups)}")
    print(f"Patterns shared by ≥2 sites:        {len(multi)}")
    print(f"NS sites in shared-pattern groups:  {sum(len(v) for v in multi.values())}")
    print(f"NS sites with unique patterns:      {sum(1 for v in groups.values() if len(v) == 1)}")

    # Classify groups
    def classify(n_carriers: int, n_samples: int) -> str:
        if n_carriers == n_samples:
            return "fixed_ref_divergence"
        if n_carriers == 1:
            return "singleton_artifact"
        return "haplogroup"

    # Write TSV + stdout report
    rows = []
    sorted_groups = sorted(multi.values(), key=lambda v: (-v[0]["n_car"], -len(v)))
    for gid, sites in enumerate(sorted_groups, 1):
        nc = sites[0]["n_car"]
        if nc < args.min_carriers:
            continue
        category = classify(nc, len(sample_cols))
        site_list = "; ".join(
            f"{s['Gene']}:{s['POS']}({s['REF']}>{s['ALT']})"
            for s in sorted(sites, key=lambda x: (x["Gene"], x["POS"]))
        )
        rows.append({
            "group_id":   gid,
            "n_sites":    len(sites),
            "n_carriers": nc,
            "n_samples":  len(sample_cols),
            "category":   category,
            "site_list":  site_list,
        })

    if not rows:
        print("\nNo multi-site groups passed --min-carriers.")
        return

    out_df = pd.DataFrame(rows)
    out_df.to_csv(OUT_TSV, sep="\t", index=False)
    print(f"\nWrote {OUT_TSV}\n")

    # Pretty print
    print(f"{'#':>3}  {'sites':>5}  {'carriers':>8}  {'category':<23}  detail")
    print("-" * 80)
    for r in rows:
        det = r["site_list"] if len(r["site_list"]) < 80 else r["site_list"][:77] + "..."
        print(f"{r['group_id']:>3}  {r['n_sites']:>5}  "
              f"{r['n_carriers']:>3}/{r['n_samples']:<4}  "
              f"{r['category']:<23}  {det}")


if __name__ == "__main__":
    main()

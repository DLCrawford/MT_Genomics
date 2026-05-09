#!/usr/bin/env python3
"""
scripts/08_call_haplotypes.py

Haplotype calling from the canonical MT variant VCF.

INPUT
-----
vcf/Fhet_MT_CDS.snps.split.vcf.gz   (SnpEff-annotated, CDS-only, split multiallelic)

CALLING RULE
------------
For each sample at each site:
    AD_alt / DP > 0.7  → 1  (alt)
    AD_alt / DP ≤ 0.7  → 0  (ref)
    missing / no data  → '.' → imputed to 0 (see below)

SPLIT-SITE IMPUTATION
---------------------
After bcftools norm -m -any, a triallelic site (e.g. T/C/A) becomes two rows:
    Row 1:  T → C  (dominant alt)
    Row 2:  T → A  (minor alt / V2)

At the V2 row, samples that carry T→C or the ref allele get GT=./. from bcftools,
producing '.' in the calling step.  '.' is imputed to 0 (absent for that allele).

This rule applies uniformly to ALL missing values:
  - Split-site '.'  : individual has a different allele → 0 is correct.
  - Coverage-gap '.': extremely rare in high-coverage MT data → 0 is conservative.

After imputation, rows where every sample = 0 (never-called alt rows, e.g. T→G
where no individual ever had G) are dropped.

OUTPUTS
-------
vcf/haplotype_matrix.csv   — long table: one row per site, columns = metadata + samples (0/1)
vcf/haplotype_calls.csv    — one row per sample: Sample, Haplotype ("C_0110..."), N_alt_sites

Run
---
cd ~/Projects/MT_Genomics_Cl_Ap2026/MT_Genomics2
python scripts/08_call_haplotypes.py
"""

from pathlib import Path
import sys

import numpy as np
import pandas as pd
import cyvcf2

# ── CONFIG ────────────────────────────────────────────────────────────────────
PROJECT_ROOT = Path(__file__).resolve().parent.parent

VCF_IN     = PROJECT_ROOT / "vcf" / "Fhet_MT_CDS.snps.split.vcf.gz"
OUT_MATRIX = PROJECT_ROOT / "vcf" / "haplotype_matrix.csv"
OUT_HAPLO  = PROJECT_ROOT / "vcf" / "haplotype_calls.csv"

AD_THRESHOLD = 0.7   # AD_alt / DP must exceed this to call alt (1)
MIN_DP       = 3     # samples with DP < MIN_DP treated as missing (.)

# Samples to exclude (no phenotype data).
# Use the sample names exactly as they appear in the VCF header.
# Stage-05 reheader strips '_0', so '10_0' → '10' (or '10_MT' if _MT was retained).
EXCLUDE_SAMPLES = {"70", "125", "70_MT", "125_MT"}


# ── HELPERS ───────────────────────────────────────────────────────────────────

def parse_sample_call(variant: cyvcf2.Variant, sample_idx: int) -> object:
    """
    Return 1, 0, or '.' for one sample at one (already-split) site.

    '.' means: no usable data — either no coverage, or bcftools set GT=./.
    because this sample carries a different allele at a split site.
    """
    gt = variant.genotypes[sample_idx]
    # cyvcf2: gt = [allele1, allele2, phased]; -1 means missing allele
    if gt[0] == -1 or gt[1] == -1:
        return "."

    try:
        ad = variant.format("AD")[sample_idx]   # [ref_depth, alt1_depth]
        dp = variant.format("DP")[sample_idx][0]
    except Exception:
        return "."

    if dp is None or dp < MIN_DP:
        return "."
    # cyvcf2 encodes missing integer fields as a large negative sentinel
    if ad is None or ad[1] < 0:
        return "."

    return 1 if (ad[1] / dp) > AD_THRESHOLD else 0


def parse_ann(variant: cyvcf2.Variant) -> tuple:
    """
    Extract (gene, effect, impact) from the first ANN entry.

    SnpEff ANN format (pipe-delimited):
      Allele | Effect | Impact | GeneName | GeneID | FeatureType | FeatureID |
      TranscriptBiotype | Rank | HGVS.c | HGVS.p | cDNA.pos | CDS.pos |
      AA.pos | Distance | ERRORS
    """
    ann_raw = variant.INFO.get("ANN", "")
    if not ann_raw:
        return "", "", ""
    fields = ann_raw.split(",")[0].split("|")
    effect = fields[1] if len(fields) > 1 else ""
    impact = fields[2] if len(fields) > 2 else ""
    gene   = fields[3] if len(fields) > 3 else ""
    return gene, effect, impact


# ── MAIN ──────────────────────────────────────────────────────────────────────

def main() -> None:
    if not VCF_IN.exists():
        sys.exit(f"ERROR: VCF not found: {VCF_IN}\n"
                 "  Run scripts/07_cds_snps_norm_mac.sh first.")

    vcf = cyvcf2.VCF(str(VCF_IN))
    all_samples = list(vcf.samples)

    keep_idx = [i for i, s in enumerate(all_samples) if s not in EXCLUDE_SAMPLES]
    samples  = [all_samples[i] for i in keep_idx]
    n_excl   = len(all_samples) - len(samples)
    print(f"Samples : {len(samples)} kept, {n_excl} excluded")

    rows = []
    for variant in vcf:
        gene, effect, impact = parse_ann(variant)

        calls = {s: parse_sample_call(variant, idx)
                 for idx, s in zip(keep_idx, samples)}

        rows.append({
            "CHROM":        variant.CHROM,
            "POS":          variant.POS,
            "REF":          variant.REF,
            "ALT":          variant.ALT[0],
            "Gene":         gene,
            "Variant_type": effect,
            "Impact":       impact,
            **calls,
        })

    vcf.close()
    df = pd.DataFrame(rows)
    meta_cols = ["CHROM", "POS", "REF", "ALT", "Gene", "Variant_type", "Impact"]
    print(f"Sites from VCF : {len(df)}")

    # ── IMPUTE '.' → 0 ────────────────────────────────────────────────────────
    # Uniform rule: missing = absent for this allele = 0.
    # Applies to both split-site GT=./. and rare true coverage gaps.
    for col in samples:
        df[col] = df[col].replace(".", 0).astype(int)

    # ── DROP NEVER-CALLED ROWS ─────────────────────────────────────────────────
    # After imputation, a row of all-zeros was never observed in any sample.
    # These arise from bcftools norm splitting a site into more alleles than
    # were actually genotyped (e.g. T→G when no individual ever had G).
    never_called = df[samples].sum(axis=1) == 0
    n_dropped = never_called.sum()
    if n_dropped:
        print(f"Dropped never-called rows : {n_dropped}")
        df = df[~never_called].reset_index(drop=True)

    print(f"Sites after filtering       : {len(df)}")

    # ── WRITE HAPLOTYPE MATRIX ─────────────────────────────────────────────────
    df.to_csv(OUT_MATRIX, index=False)
    print(f"Matrix written  : {OUT_MATRIX}")

    # ── CALL HAPLOTYPES ────────────────────────────────────────────────────────
    # Haplotype string: "C_" + binary calls concatenated across all sites
    # in VCF order (position-sorted, CDS only).
    # Example: "C_110001000..." where each character = one site.
    haplo_rows = []
    for s in samples:
        binary_calls = df[s].astype(str).tolist()
        hap_str = "C_" + "".join(binary_calls)
        n_alt   = df[s].sum()
        haplo_rows.append({"Sample": s, "Haplotype": hap_str, "N_alt_sites": int(n_alt)})

    haplo_df = pd.DataFrame(haplo_rows)

    # Summary
    counts = (haplo_df.groupby("Haplotype")
                      .agg(N=("Sample", "count"))
                      .reset_index()
                      .sort_values("N", ascending=False))
    print(f"\nUnique haplotypes : {len(counts)}")
    print(counts.to_string(index=False))

    haplo_df.to_csv(OUT_HAPLO, index=False)
    print(f"\nHaplotype calls written : {OUT_HAPLO}")


if __name__ == "__main__":
    main()

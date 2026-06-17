#!/usr/bin/env python3
"""
scripts/11_haplotypes_nonsyn.py

Haplotype calls + matrix using **nonsynonymous (MODERATE-impact) sites only**.

Same calling rule as scripts/08_call_haplotypes.py — only the site filter
differs. Stage 08 calls haplotypes on *all* CDS SNPs; this stage restricts
to amino-acid-changing variants (missense_variant), which is a stricter
biological signal: two samples sharing a nonsyn haplotype share protein-
level identity, not just silent variation.

INPUT
-----
    vcf/Fhet_MT_CDS.snps.split.vcf.gz   (SnpEff-annotated, CDS-only, split multiallelic)

SITE FILTER
-----------
Keep a row iff its **first** ANN entry's effect is in:
    missense_variant
    splice_region_variant&missense_variant
    (i.e., MODERATE-impact, amino-acid-changing)

Stop-impact effects (stop_gained, stop_lost, start_lost) — HIGH in SnpEff —
are intentionally NOT included per the user spec "(MODERATE)". To include
them, add them to the NONSYN_EFFECTS set below.

CALLING RULE (matches stage 08)
-------------------------------
For each sample at each kept site:
    AD_alt / DP > 0.7  → 1  (alt-homoplasmic)
    AD_alt / DP ≤ 0.7  → 0  (ref-dominant)
    missing / no data  → '.' → imputed to 0

SPLIT-SITE IMPUTATION
---------------------
Same as stage 08: at a row that came from a split multi-ALT site, samples
that carry the other ALT (or REF) get GT=./. → '.' → 0.

OUTPUTS
-------
    vcf/haplotype_matrix_nonsyn.csv    long table: one row per kept site, sample cols = 0/1
    vcf/haplotype_calls_nonsyn.csv     one row per sample: Sample, Haplotype ("N_0110..."), N_alt_sites

Run
---
    cd ~/Projects/MT_Genomics_Cl_Ap2026/MT_Genomics2
    python scripts/11_haplotypes_nonsyn.py
"""

from pathlib import Path
import sys

import pandas as pd
import cyvcf2

# ── CONFIG ────────────────────────────────────────────────────────────────────
PROJECT_ROOT = Path(__file__).resolve().parent.parent

VCF_IN     = PROJECT_ROOT / "vcf" / "Fhet_MT_CDS.snps.split.vcf.gz"
OUT_MATRIX = PROJECT_ROOT / "vcf" / "haplotype_matrix_nonsyn.csv"
OUT_HAPLO  = PROJECT_ROOT / "vcf" / "haplotype_calls_nonsyn.csv"

AD_THRESHOLD = 0.7   # AD_alt / DP must exceed this to call alt (1)
MIN_DP       = 3     # samples with DP < MIN_DP treated as missing (.)

# Match stage 08's exclusions: samples with no phenotype data.
EXCLUDE_SAMPLES = {"70", "125", "70_MT", "125_MT"}

# Nonsynonymous effects (MODERATE per user spec). To include HIGH-impact
# stops, add: "stop_gained", "stop_lost", "start_lost", etc.
NONSYN_EFFECTS = {
    "missense_variant",
    "splice_region_variant&missense_variant",
}


# ── HELPERS (mirror stage 08) ────────────────────────────────────────────────
def parse_sample_call(variant: cyvcf2.Variant, sample_idx: int) -> object:
    gt = variant.genotypes[sample_idx]
    if gt[0] == -1 or gt[1] == -1:
        return "."
    try:
        ad = variant.format("AD")[sample_idx]
        dp = variant.format("DP")[sample_idx][0]
    except Exception:
        return "."
    if dp is None or dp < MIN_DP:
        return "."
    if ad is None or ad[1] < 0:
        return "."
    return 1 if (ad[1] / dp) > AD_THRESHOLD else 0


def parse_ann(variant: cyvcf2.Variant) -> tuple:
    raw = variant.INFO.get("ANN", "")
    if not raw:
        return "", "", ""
    fields = raw.split(",")[0].split("|")
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
    samples = [all_samples[i] for i in keep_idx]
    print(f"Samples : {len(samples)} kept, {len(all_samples) - len(samples)} excluded")
    print(f"Filter  : effect in {sorted(NONSYN_EFFECTS)}")

    rows = []
    n_total = 0
    n_kept = 0
    for variant in vcf:
        n_total += 1
        gene, effect, impact = parse_ann(variant)
        if effect not in NONSYN_EFFECTS:
            continue
        n_kept += 1

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
    print(f"Sites total       : {n_total}")
    print(f"Sites kept (nonsyn): {n_kept}")

    if not rows:
        sys.exit("ERROR: no rows survived the nonsyn filter. Check ANN field on input VCF.")

    df = pd.DataFrame(rows)

    # ── IMPUTE '.' → 0 (same rule as stage 08) ────────────────────────────────
    for col in samples:
        df[col] = df[col].replace(".", 0).astype(int)

    # ── DROP NEVER-CALLED ROWS ────────────────────────────────────────────────
    never_called = df[samples].sum(axis=1) == 0
    n_dropped = never_called.sum()
    if n_dropped:
        print(f"Dropped never-called rows : {n_dropped}")
        df = df[~never_called].reset_index(drop=True)
    print(f"Sites after never-called drop : {len(df)}")

    # ── WRITE MATRIX ──────────────────────────────────────────────────────────
    df.to_csv(OUT_MATRIX, index=False)
    print(f"Matrix written  : {OUT_MATRIX}")

    # ── HAPLOTYPE STRINGS ─────────────────────────────────────────────────────
    # Prefix "N_" (nonsynonymous) to distinguish from stage 08's "C_" (CDS-all).
    haplo_rows = []
    for s in samples:
        binary_calls = df[s].astype(str).tolist()
        hap_str = "N_" + "".join(binary_calls)
        n_alt = df[s].sum()
        haplo_rows.append({"Sample": s, "Haplotype": hap_str, "N_alt_sites": int(n_alt)})

    haplo_df = pd.DataFrame(haplo_rows)

    counts = (haplo_df.groupby("Haplotype")
                      .agg(N=("Sample", "count"))
                      .reset_index()
                      .sort_values("N", ascending=False))
    print(f"\nUnique nonsyn haplotypes : {len(counts)}")
    print(counts.head(20).to_string(index=False))
    if len(counts) > 20:
        print(f"  ... ({len(counts) - 20} more)")

    haplo_df.to_csv(OUT_HAPLO, index=False)
    print(f"\nHaplotype calls written : {OUT_HAPLO}")


if __name__ == "__main__":
    main()

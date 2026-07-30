#!/usr/bin/env python3
"""
23_celegans_pi.py
Calculate nucleotide diversity (π) for 540 C. elegans mitochondrial
CDS variants from C_elegansFINAL.annotated.vcf.

Key differences from Fhet (19_calc_pi.py):
  - GATK diploid calls (GT: 0/0, 0/1, 1/1) — not AD-based haploid calls
  - C. elegans is ~99.9% selfing → effectively haploid
  - Genotype calling: 0/0 = REF, 1/1 = ALT, 0/1 = skip (heterozygous/ambiguous)
  - Filter: PASS sites only, SNPs only (skip indels/MNP)
  - SnpEff ANN field same format as Fhet

Formula (haploid-equivalent, one allele per individual):
    π_site = (n / n-1) * 2 * p * (1-p)
    π      = Σ π_site / L_CDS

Usage (from C_elegans/):
    conda activate SNP_env
    python 23_celegans_pi.py

Outputs:
    celegans_pi_summary.txt
    celegans_pi_per_site.tsv
"""

import pysam
import sys
import csv
import os

# ── Config ────────────────────────────────────────────────────────────────────
IN_VCF   = "C_elegansFINAL.annotated.vcf"
OUT_SUM  = "celegans_pi_summary.txt"
OUT_SITE = "celegans_pi_per_site.tsv"
L_CDS    = 10_299    # C. elegans mt CDS length (bp) — Table 2
# ─────────────────────────────────────────────────────────────────────────────

# ── SnpEff effect classification (same sets as 19_calc_pi.py) ────────────────
SYN_EFFECTS = {
    "synonymous_variant",
    "stop_retained_variant",
}
NS_EFFECTS = {
    "missense_variant",
    "stop_gained",
    "stop_lost",
    "start_lost",
    "start_gained",
    "frameshift_variant",
    "disruptive_inframe_insertion",
    "disruptive_inframe_deletion",
    "conservative_inframe_insertion",
    "conservative_inframe_deletion",
    "splice_acceptor_variant",
    "splice_donor_variant",
}
CDS_EFFECTS = SYN_EFFECTS | NS_EFFECTS
IMPACT_RANK = {"HIGH": 0, "MODERATE": 1, "LOW": 2, "MODIFIER": 3}
# ─────────────────────────────────────────────────────────────────────────────


def best_cds_effect(info, alt_allele):
    """Return (effect, impact) for best CDS annotation; (None,None) if non-CDS."""
    ann_raw = info.get("ANN", None)
    if not ann_raw:
        return None, None

    entries = list(ann_raw) if isinstance(ann_raw, (tuple, list)) else ann_raw.split(",")
    candidates, fallback = [], []

    for entry in entries:
        parts = entry.split("|")
        if len(parts) < 4:
            continue
        allele, effect, impact, gene = parts[0], parts[1], parts[2], parts[3]
        # Handle &-joined effects (e.g. synonymous_variant&splice_region_variant)
        matched = next((e for e in effect.split("&") if e in CDS_EFFECTS), None)
        if matched:
            rec = (allele, matched, impact, gene)
            (candidates if allele == alt_allele else fallback).append(rec)

    pool = candidates if candidates else fallback
    if not pool:
        return None, None
    best = min(pool, key=lambda x: IMPACT_RANK.get(x[2], 99))
    return best[1], best[2]


def call_gt(sample):
    """
    Parse diploid GT for a selfing organism:
      0/0 → 0 (REF)
      1/1 → 1 (ALT)
      0/1 or ./. → '.' (skip)
    """
    gt = sample.get("GT", None)
    if gt is None:
        return "."
    # pysam returns GT as tuple of allele indices, e.g. (0,0), (1,1), (0,1)
    if isinstance(gt, tuple):
        alleles = [a for a in gt if a is not None]
        if not alleles:
            return "."
        if all(a == 0 for a in alleles):
            return 0
        if all(a == 1 for a in alleles):
            return 1
        return "."   # heterozygous → skip
    # String fallback
    gt_str = str(gt).replace("|", "/")
    if gt_str in ("0/0", "0"):
        return 0
    if gt_str in ("1/1", "1"):
        return 1
    return "."


def pi_site(n_called, n_alt):
    if n_called < 2:
        return 0.0
    p = n_alt / n_called
    return (n_called / (n_called - 1)) * 2 * p * (1 - p)


def main():
    if not os.path.exists(IN_VCF):
        print(f"ERROR: {IN_VCF} not found. Run from C_elegans/.", file=sys.stderr)
        sys.exit(1)

    vcf = pysam.VariantFile(IN_VCF)
    samples = list(vcf.header.samples)
    N = len(samples)
    print(f"Samples : {N}", file=sys.stderr)
    print(f"L_CDS   : {L_CDS:,} bp", file=sys.stderr)

    pi_total = pi_syn = pi_ns = 0.0
    n_snp = n_cds = n_syn = n_ns = n_skip_indel = n_skip_filter = 0
    rows = []

    for rec in vcf.fetch():

        # PASS filter only
        if rec.filter.keys() and "PASS" not in rec.filter.keys():
            n_skip_filter += 1
            continue

        if rec.alts is None:
            continue

        for alt in rec.alts:
            if alt == "*":
                continue

            # SNPs only — skip indels and MNPs
            if len(rec.ref) != 1 or len(alt) != 1:
                n_skip_indel += 1
                continue

            n_snp += 1

            # CDS classification from ANN
            effect, impact = best_cds_effect(rec.info, alt)
            if effect is None:
                continue

            if effect in SYN_EFFECTS:
                site_class = "synonymous"
            elif effect in NS_EFFECTS:
                site_class = "nonsynonymous"
            else:
                continue

            n_cds += 1

            # Count alleles across samples
            n_called = n_alt_count = 0
            for sname in samples:
                call = call_gt(rec.samples[sname])
                if call != ".":
                    n_called += 1
                    if call == 1:
                        n_alt_count += 1

            if n_called < 2:
                continue

            freq   = n_alt_count / n_called
            p_site = pi_site(n_called, n_alt_count)

            pi_total += p_site
            if site_class == "synonymous":
                pi_syn += p_site
                n_syn  += 1
            else:
                pi_ns  += p_site
                n_ns   += 1

            rows.append({
                "CHROM":     rec.chrom,
                "POS":       rec.pos,
                "REF":       rec.ref,
                "ALT":       alt,
                "Effect":    effect,
                "Impact":    impact,
                "Class":     site_class,
                "n_called":  n_called,
                "n_alt":     n_alt_count,
                "freq_alt":  round(freq, 6),
                "pi_site":   round(p_site, 8),
            })

    vcf.close()

    # Normalise
    pi_t = pi_total / L_CDS
    pi_s = pi_syn   / L_CDS
    pi_n = pi_ns    / L_CDS
    pnps = pi_n / pi_s if pi_s > 0 else float("nan")

    # Per-site TSV
    fields = ["CHROM","POS","REF","ALT","Effect","Impact","Class",
              "n_called","n_alt","freq_alt","pi_site"]
    with open(OUT_SITE, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    summary = (
        f"=== Nucleotide Diversity (π) — C. elegans mt CDS ===\n"
        f"  Samples (N)          : {N}\n"
        f"  L_CDS                : {L_CDS:,} bp\n"
        f"  Total SNPs parsed    : {n_snp}  (indels/MNP skipped: {n_skip_indel})\n"
        f"  CDS variant sites    : {n_cds}  (syn={n_syn}, NS={n_ns})\n"
        f"\n"
        f"  π_total              : {pi_t:.5f}\n"
        f"  π_syn                : {pi_s:.5f}\n"
        f"  π_ns                 : {pi_n:.5f}\n"
        f"  pN/pS  (π_ns/π_syn)  : {pnps:.3f}\n"
        f"\nNote: 0/1 heterozygous calls excluded (selfing organism).\n"
        f"Per-site table: {OUT_SITE}\n"
    )
    print("\n" + summary)
    with open(OUT_SUM, "w") as f:
        f.write(summary)


if __name__ == "__main__":
    main()

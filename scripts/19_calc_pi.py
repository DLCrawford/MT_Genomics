#!/usr/bin/env python3
"""
Calculate nucleotide diversity (π) for mitochondrial CDS variants.
Splits into total CDS, synonymous, and nonsynonymous π.

Formula (per site, haploid):
    π_site = (n / n-1) * 2 * p * (1-p)
    π_total = Σ π_site / L_CDS

Usage:
    python calc_pi.py -i Fhet_MT_Annot_AD.vcf.gz \
                      -g Fhet_MT.gff \
                      -o pi_results.tsv \
                      [--threshold 0.7] [--L_CDS 11417]

Outputs:
    - Summary to stdout
    - Per-site table to <output>
"""

import pysam
import argparse
import sys
from collections import defaultdict

# ── SnpEff effect classification ──────────────────────────────────────────────
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

# ── Sample / site exclusions ─────────────────────────────────────────────────
# None needed — 141_MT_variants.vcf.gz is already filtered to 141 samples
# with monomorphic sites (private to 70/125) removed.
# The 9 invariant fixed-ALT sites (6 SYN + 3 NS reference-divergence) remain
# in the VCF but contribute π = 0 and do not affect the result.
EXCLUDE_SAMPLES = set()
EXCLUDE_POS     = set()
# ─────────────────────────────────────────────────────────────────────────────


def parse_args():
    p = argparse.ArgumentParser(description="Calculate π (total, syn, NS) from annotated VCF.")
    p.add_argument("-i", "--input", required=True,
                   help="SnpEff-annotated VCF/BCF (bgzipped + indexed)")
    p.add_argument("-g", "--gff", required=False,
                   help="GFF3 annotation file (used to compute L_CDS if not provided)")
    p.add_argument("-o", "--output", required=True,
                   help="Output per-site TSV file")
    p.add_argument("--threshold", type=float, default=0.7,
                   help="AD[ALT]/DP threshold for confident haploid call (default: 0.7)")
    p.add_argument("--L_CDS", type=int, default=None,
                   help="Total CDS length in bp (if omitted, computed from GFF)")
    p.add_argument("--L_syn", type=int, default=None,
                   help="Synonymous sites length (optional; if omitted uses L_CDS)")
    p.add_argument("--L_ns", type=int, default=None,
                   help="Nonsynonymous sites length (optional; if omitted uses L_CDS)")
    return p.parse_args()


# ── GFF parsing ───────────────────────────────────────────────────────────────
def compute_L_CDS_from_gff(gff_path):
    """Return total non-overlapping CDS length from a GFF3 file."""
    intervals = []
    with open(gff_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            if parts[2].upper() == "CDS":
                intervals.append((int(parts[3]), int(parts[4])))

    if not intervals:
        print("WARNING: no CDS features found in GFF.", file=sys.stderr)
        return None

    # Merge overlapping intervals
    intervals.sort()
    merged = [intervals[0]]
    for start, end in intervals[1:]:
        if start <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], end))
        else:
            merged.append((start, end))

    L = sum(e - s + 1 for s, e in merged)
    print(f"L_CDS from GFF: {L:,} bp ({len(merged)} merged CDS intervals)", file=sys.stderr)
    return L


# ── ANN parsing ───────────────────────────────────────────────────────────────
def best_cds_effect(info, alt_allele):
    """
    Return (effect, impact) for the best CDS annotation matching alt_allele.
    Returns (None, None) if the site is not in a CDS.
    """
    ann_raw = info.get("ANN", None)
    if not ann_raw:
        return None, None

    # pysam returns multi-value INFO fields as tuples; each element is one annotation
    if isinstance(ann_raw, (tuple, list)):
        entries = list(ann_raw)
    else:
        entries = ann_raw.split(",")

    candidates = []
    fallback = []

    for entry in entries:
        parts = entry.split("|")
        if len(parts) < 4:
            continue
        allele, effect, impact, gene = parts[0], parts[1], parts[2], parts[3]
        # SnpEff sometimes joins effects with '&' (e.g. synonymous_variant&splice_region_variant)
        # Check if any component of the joined effect is in CDS_EFFECTS
        effect_components = effect.split("&")
        matched = next((e for e in effect_components if e in CDS_EFFECTS), None)
        if matched:
            record = (allele, matched, impact, gene)
            if allele == alt_allele:
                candidates.append(record)
            else:
                fallback.append(record)

    pool = candidates if candidates else fallback
    if not pool:
        return None, None  # not a CDS site for this allele

    best = min(pool, key=lambda x: IMPACT_RANK.get(x[2], 99))
    return best[1], best[2]   # effect, impact


# ── Haploid genotype caller (same logic as ANN_Claude_1A.py) ─────────────────
def call_haplotype(fmt, dp_info, threshold):
    """
    Returns:
        1  -> confidently ALT
        0  -> confidently REF
        .  -> ambiguous / missing
    """
    gt = fmt.get("GT", ".")
    if gt in (".", "./.", ".|."):
        return "."

    ad_raw = fmt.get("AD", None)
    if ad_raw is None:
        return "."

    try:
        ad = [int(x) if x is not None else 0 for x in ad_raw] \
             if isinstance(ad_raw, (list, tuple)) else [int(ad_raw)]
    except (ValueError, TypeError):
        return "."

    if len(ad) < 2:
        return "."

    dp_raw = fmt.get("DP", None)
    try:
        dp = int(dp_raw) if dp_raw is not None else dp_info
    except (ValueError, TypeError):
        dp = dp_info

    if not dp:
        return "."

    if ad[1] > threshold * dp:
        return 1
    elif ad[0] > threshold * dp:
        return 0
    return "."


# ── π calculation ─────────────────────────────────────────────────────────────
def pi_site(n_called, n_alt):
    """
    Haploid π at one site.
    π = (n / n-1) * 2 * p * (1-p)
    """
    if n_called < 2:
        return 0.0
    p = n_alt / n_called
    q = 1.0 - p
    return (n_called / (n_called - 1)) * 2 * p * q


# ── Main ──────────────────────────────────────────────────────────────────────
def main():
    args = parse_args()

    # ── Determine CDS lengths ─────────────────────────────────────────────────
    L_CDS = args.L_CDS
    if L_CDS is None:
        if args.gff:
            L_CDS = compute_L_CDS_from_gff(args.gff)
        else:
            print("ERROR: provide --L_CDS or --gff to determine CDS length.", file=sys.stderr)
            sys.exit(1)

    L_syn = args.L_syn or L_CDS   # user can supply; otherwise same denominator
    L_ns  = args.L_ns  or L_CDS

    print(f"L_CDS = {L_CDS:,}  |  L_syn = {L_syn:,}  |  L_ns = {L_ns:,}", file=sys.stderr)

    # ── Open VCF ─────────────────────────────────────────────────────────────
    vcf = pysam.VariantFile(args.input)
    all_samples = list(vcf.header.samples)
    samples = [s for s in all_samples if s not in EXCLUDE_SAMPLES]
    n_excl = len(all_samples) - len(samples)
    print(f"Samples in VCF: {len(all_samples)}  |  excluded: {n_excl} {sorted(EXCLUDE_SAMPLES & set(all_samples))}  |  kept: {len(samples)}", file=sys.stderr)

    format_keys = None

    pi_total = 0.0
    pi_syn   = 0.0
    pi_ns    = 0.0

    rows = []   # per-site output

    for i, rec in enumerate(vcf.fetch()):
        if format_keys is None:
            format_keys = list(rec.format.keys())

        if rec.alts is None:
            continue

        if rec.pos in EXCLUDE_POS:
            continue

        dp_info = rec.info.get("DP", None)

        # Process each ALT allele separately
        for alt in rec.alts:
            if alt == "*":          # spanning deletion placeholder
                continue

            effect, impact = best_cds_effect(rec.info, alt)
            if effect is None:
                continue            # site not in CDS for this allele

            # Classify
            if effect in SYN_EFFECTS:
                site_class = "synonymous"
            elif effect in NS_EFFECTS:
                site_class = "nonsynonymous"
            else:
                continue

            # Call haplotypes across samples
            n_called = 0
            n_alt_allele = 0
            for sname in samples:
                s = rec.samples[sname]
                fmt = {}
                for key in format_keys:
                    try:
                        fmt[key] = s[key]
                    except KeyError:
                        fmt[key] = None

                call = call_haplotype(fmt, dp_info, args.threshold)
                if call != ".":
                    n_called += 1
                    if call == 1:
                        n_alt_allele += 1

            if n_called < 2:
                continue

            p_site = pi_site(n_called, n_alt_allele)
            freq   = n_alt_allele / n_called

            pi_total += p_site
            if site_class == "synonymous":
                pi_syn += p_site
            else:
                pi_ns  += p_site

            rows.append({
                "CHROM":      rec.chrom,
                "POS":        rec.pos,
                "REF":        rec.ref,
                "ALT":        alt,
                "Effect":     effect,
                "Impact":     impact,
                "Class":      site_class,
                "n_called":   n_called,
                "n_alt":      n_alt_allele,
                "freq_alt":   round(freq, 6),
                "pi_site":    round(p_site, 8),
            })

        if i % 500 == 0:
            print(f"  Processed {i} VCF records ...", file=sys.stderr)

    vcf.close()

    # ── Normalise by length ───────────────────────────────────────────────────
    pi_total_norm = pi_total / L_CDS
    pi_syn_norm   = pi_syn   / L_syn
    pi_ns_norm    = pi_ns    / L_ns
    pn_ps         = (pi_ns_norm / pi_syn_norm) if pi_syn_norm > 0 else float("nan")

    # ── Write per-site table ──────────────────────────────────────────────────
    import csv
    fieldnames = ["CHROM","POS","REF","ALT","Effect","Impact","Class",
                  "n_called","n_alt","freq_alt","pi_site"]
    with open(args.output, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    print(f"\nPer-site table: {len(rows)} CDS variant sites → {args.output}", file=sys.stderr)

    # ── Print summary ─────────────────────────────────────────────────────────
    n_syn = sum(1 for r in rows if r["Class"] == "synonymous")
    n_ns  = sum(1 for r in rows if r["Class"] == "nonsynonymous")

    print("\n" + "="*55)
    print("  Nucleotide Diversity (π) — Fhet Mitochondrial CDS")
    print("="*55)
    print(f"  Samples (N)          : {len(samples)}  (excluded: {n_excl})")
    print(f"  CDS variant sites    : {len(rows)}  (syn={n_syn}, NS={n_ns})")
    print(f"  L_CDS                : {L_CDS:,} bp")
    print(f"  π_total              : {pi_total_norm:.5f}")
    print(f"  π_syn                : {pi_syn_norm:.5f}")
    print(f"  π_ns                 : {pi_ns_norm:.5f}")
    print(f"  pN/pS  (π_ns/π_syn)  : {pn_ps:.3f}")
    print("="*55)


if __name__ == "__main__":
    main()

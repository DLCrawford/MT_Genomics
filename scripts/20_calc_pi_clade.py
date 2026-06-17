#!/usr/bin/env python3
"""
20_calc_pi_clade.py — Nucleotide diversity (π) split by mitochondrial clade.

Adapted from 19_calc_pi.py. Instead of one π for all 141 samples, this
computes π (total / synonymous / nonsynonymous) SEPARATELY for the
northern and southern mitochondrial clades.

Why split by clade
------------------
Variants are polarized against a northern REF (NC_012312.1). Raw per-sample
ALT counts therefore measure *distance from the reference*, not within-group
polymorphism: north-clade individuals look ALT-poor by construction and
south-clade individuals ALT-rich. π is reference-independent — it is the
mean pairwise difference among individuals *within a group* — so computing
π separately within each clade gives a fair comparison of how polymorphic
each clade actually is. Clade-defining fixed differences are ~monomorphic
within a clade and contribute π ≈ 0, as they should.

Clade assignment (bimodal per-individual ALT burden)
---------------------------------------------------
Pass 1 counts confident ALT calls per sample across all VCF records using
the same haploid caller as π (AD[ALT]/DP ≥ threshold). The panel is
strongly bimodal:
    north-clade : < ~50 ALTs   (observed 15-28)
    south-clade : > ~200 ALTs  (observed 217-233)
Samples between the two cutoffs are AMBIGUOUS and excluded from π (by
default one sample, 77_MT at ~193 ALTs).

Formula (per site, haploid), computed within each clade:
    π_site  = (n / n-1) * 2 * p * (1-p)
    π_total = Σ π_site / L

Usage:
    python 20_calc_pi_clade.py -i vcf/141_MT_variants.vcf.gz \
                               -g Missing_Files/SSM_MT_ref/Fhet_MT.gff \
                               -o vcf/pi_by_clade_persite.tsv \
                               [--threshold 0.7] [--L_CDS 11417] \
                               [--north-max 50] [--south-min 200]

Outputs:
    - Per-site table to <output> (n_called/n_alt/freq/π per clade)
    - Clade-membership table to <output>.membership.tsv
    - Summary (per-clade π_total, π_syn, π_ns, pN/pS) to stdout
"""

import pysam
import argparse
import sys
import csv

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

EXCLUDE_SAMPLES = set()
EXCLUDE_POS     = set()

CLADES = ("north", "south")


def parse_args():
    p = argparse.ArgumentParser(description="Calculate π (total, syn, NS) per mt clade.")
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
    p.add_argument("--north-max", type=int, default=50,
                   help="Max per-individual ALT count to be called north clade (default: 50)")
    p.add_argument("--south-min", type=int, default=200,
                   help="Min per-individual ALT count to be called south clade (default: 200)")
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
    """Return (effect, impact) for the best CDS annotation matching alt_allele."""
    ann_raw = info.get("ANN", None)
    if not ann_raw:
        return None, None

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
        return None, None

    best = min(pool, key=lambda x: IMPACT_RANK.get(x[2], 99))
    return best[1], best[2]


# ── Haploid genotype caller (same logic as 19_calc_pi.py) ────────────────────
def call_haplotype(fmt, dp_info, threshold):
    """1 -> confidently ALT, 0 -> confidently REF, '.' -> ambiguous/missing."""
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


def get_fmt(sample_rec, format_keys):
    fmt = {}
    for key in format_keys:
        try:
            fmt[key] = sample_rec[key]
        except KeyError:
            fmt[key] = None
    return fmt


# ── π calculation ─────────────────────────────────────────────────────────────
def pi_site(n_called, n_alt):
    """Haploid π at one site: (n / n-1) * 2 * p * (1-p)."""
    if n_called < 2:
        return 0.0
    p = n_alt / n_called
    q = 1.0 - p
    return (n_called / (n_called - 1)) * 2 * p * q


# ── Pass 1: per-sample ALT burden → clade assignment ─────────────────────────
def assign_clades(vcf_path, samples, threshold, north_max, south_min):
    """Count confident ALT calls per sample across all records, classify clades."""
    vcf = pysam.VariantFile(vcf_path)
    format_keys = None
    counts = {s: 0 for s in samples}

    for rec in vcf.fetch():
        if format_keys is None:
            format_keys = list(rec.format.keys())
        if rec.alts is None:
            continue
        dp_info = rec.info.get("DP", None)
        for sname in samples:
            fmt = get_fmt(rec.samples[sname], format_keys)
            if call_haplotype(fmt, dp_info, threshold) == 1:
                counts[sname] += 1
    vcf.close()

    clade_of = {}
    for s, c in counts.items():
        if c < north_max:
            clade_of[s] = "north"
        elif c > south_min:
            clade_of[s] = "south"
        else:
            clade_of[s] = "ambiguous"
    return counts, clade_of


# ── Main ──────────────────────────────────────────────────────────────────────
def main():
    args = parse_args()

    L_CDS = args.L_CDS
    if L_CDS is None:
        if args.gff:
            L_CDS = compute_L_CDS_from_gff(args.gff)
        else:
            print("ERROR: provide --L_CDS or --gff to determine CDS length.", file=sys.stderr)
            sys.exit(1)

    L_syn = args.L_syn or L_CDS
    L_ns  = args.L_ns  or L_CDS
    L = {"total": L_CDS, "syn": L_syn, "ns": L_ns}
    print(f"L_CDS = {L_CDS:,}  |  L_syn = {L_syn:,}  |  L_ns = {L_ns:,}", file=sys.stderr)

    # Sample list
    vcf = pysam.VariantFile(args.input)
    all_samples = list(vcf.header.samples)
    vcf.close()
    samples = [s for s in all_samples if s not in EXCLUDE_SAMPLES]

    # ── Pass 1: clade assignment ──────────────────────────────────────────────
    print("Pass 1: counting per-sample ALT burden for clade assignment ...", file=sys.stderr)
    counts, clade_of = assign_clades(args.input, samples, args.threshold,
                                     args.north_max, args.south_min)
    clade_members = {c: [s for s in samples if clade_of[s] == c] for c in CLADES}
    ambiguous = [s for s in samples if clade_of[s] == "ambiguous"]

    print(f"  north (<{args.north_max} ALT) : {len(clade_members['north'])}", file=sys.stderr)
    print(f"  south (>{args.south_min} ALT) : {len(clade_members['south'])}", file=sys.stderr)
    print(f"  ambiguous (excluded)         : {len(ambiguous)} "
          f"{[(s, counts[s]) for s in ambiguous]}", file=sys.stderr)

    # Write membership table
    memb_path = args.output + ".membership.tsv"
    with open(memb_path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["Sample", "n_alt_panel", "Clade"])
        for s in sorted(samples, key=lambda x: counts[x]):
            w.writerow([s, counts[s], clade_of[s]])
    print(f"Clade membership → {memb_path}", file=sys.stderr)

    # ── Pass 2: π per clade ───────────────────────────────────────────────────
    # accumulators: pi[clade][class]
    pi = {c: {"total": 0.0, "syn": 0.0, "ns": 0.0} for c in CLADES}
    nsites = {c: {"syn": 0, "ns": 0} for c in CLADES}

    vcf = pysam.VariantFile(args.input)
    format_keys = None
    rows = []

    for i, rec in enumerate(vcf.fetch()):
        if format_keys is None:
            format_keys = list(rec.format.keys())
        if rec.alts is None or rec.pos in EXCLUDE_POS:
            continue

        dp_info = rec.info.get("DP", None)

        for alt in rec.alts:
            if alt == "*":
                continue
            effect, impact = best_cds_effect(rec.info, alt)
            if effect is None:
                continue
            if effect in SYN_EFFECTS:
                site_class = "syn"
            elif effect in NS_EFFECTS:
                site_class = "ns"
            else:
                continue

            row = {
                "CHROM": rec.chrom, "POS": rec.pos, "REF": rec.ref, "ALT": alt,
                "Effect": effect, "Impact": impact,
                "Class": "synonymous" if site_class == "syn" else "nonsynonymous",
            }

            for clade in CLADES:
                n_called = n_alt = 0
                for sname in clade_members[clade]:
                    fmt = get_fmt(rec.samples[sname], format_keys)
                    call = call_haplotype(fmt, dp_info, args.threshold)
                    if call != ".":
                        n_called += 1
                        if call == 1:
                            n_alt += 1
                p_site = pi_site(n_called, n_alt)
                freq = (n_alt / n_called) if n_called else float("nan")

                pi[clade]["total"] += p_site
                pi[clade][site_class] += p_site
                if p_site > 0:
                    nsites[clade][site_class] += 1

                row[f"{clade}_n_called"] = n_called
                row[f"{clade}_n_alt"]    = n_alt
                row[f"{clade}_freq"]     = round(freq, 6) if n_called else ""
                row[f"{clade}_pi_site"]  = round(p_site, 8)

            rows.append(row)

        if i % 200 == 0:
            print(f"  Pass 2: processed {i} VCF records ...", file=sys.stderr)

    vcf.close()

    # ── Write per-site table ──────────────────────────────────────────────────
    fieldnames = ["CHROM", "POS", "REF", "ALT", "Effect", "Impact", "Class"]
    for clade in CLADES:
        fieldnames += [f"{clade}_n_called", f"{clade}_n_alt",
                       f"{clade}_freq", f"{clade}_pi_site"]
    with open(args.output, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    print(f"\nPer-site table: {len(rows)} CDS variant rows → {args.output}", file=sys.stderr)

    # ── Summary ───────────────────────────────────────────────────────────────
    print("\n" + "=" * 64)
    print("  Nucleotide Diversity (π) by mitochondrial clade — Fhet CDS")
    print("=" * 64)
    print(f"  L_CDS = {L_CDS:,} bp   threshold = {args.threshold}")
    print(f"  cutoffs: north < {args.north_max} ALT,  south > {args.south_min} ALT")
    print("-" * 64)
    header = f"  {'metric':<18}" + "".join(f"{c:>14}" for c in CLADES)
    print(header)
    print("-" * 64)
    print(f"  {'N samples':<18}" + "".join(f"{len(clade_members[c]):>14}" for c in CLADES))
    for cls, lab in [("total", "variant sites"), ("syn", "  syn sites"), ("ns", "  NS sites")]:
        if cls == "total":
            vals = "".join(f"{nsites[c]['syn'] + nsites[c]['ns']:>14}" for c in CLADES)
        else:
            vals = "".join(f"{nsites[c][cls]:>14}" for c in CLADES)
        print(f"  {lab:<18}" + vals)
    print("-" * 64)
    for cls, lab in [("total", "π_total"), ("syn", "π_syn"), ("ns", "π_ns")]:
        vals = "".join(f"{pi[c][cls] / L[cls]:>14.6f}" for c in CLADES)
        print(f"  {lab:<18}" + vals)
    pnps = []
    for c in CLADES:
        ps = pi[c]["syn"] / L_syn
        pn = pi[c]["ns"] / L_ns
        pnps.append(pn / ps if ps > 0 else float("nan"))
    print(f"  {'pN/pS':<18}" + "".join(f"{v:>14.3f}" for v in pnps))
    print("=" * 64)
    if ambiguous:
        print(f"  Excluded (ambiguous): {[(s, counts[s]) for s in ambiguous]}")


if __name__ == "__main__":
    main()

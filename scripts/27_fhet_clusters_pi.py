#!/usr/bin/env python3
"""
27_fhet_clusters_pi.py
Calculate π and Watterson's θ for Fhet mt CDS split into two haplogroup clusters:
  - North clade : n_variants_total <  50  (~77 samples)
  - South clade : n_variants_total > 200  (~65 samples, 77_MT excluded)

Cluster membership from vcf/per_individual_burden_pileup.tsv.
Variant calling from vcf/141_MT_variants.vcf.gz (haploid, AD-based).

Usage (from MT_Genomics2/):
    conda activate SNP_env
    python scripts/27_fhet_clusters_pi.py

Outputs:
    vcf/fhet_clusters_pi_summary.txt
    vcf/fhet_clusters_pi_per_site.tsv
"""

import pysam
import csv
import sys
import os

# ── Config ────────────────────────────────────────────────────────────────────
VCF         = "vcf/141_MT_variants.vcf.gz"
BURDEN      = "vcf/per_individual_burden_pileup.tsv"
GFF         = "Missing_Files/SSM_MT_ref/Fhet_MT.gff"
OUT_SUM     = "vcf/fhet_clusters_pi_summary.txt"
OUT_SITE    = "vcf/fhet_clusters_pi_per_site.tsv"
L_CDS       = 11_417
THRESHOLD   = 0.7    # AD-based haploid call threshold
NORTH_MAX   = 50     # < 50 variants → north
SOUTH_MIN   = 200    # > 200 variants → south
EXCLUDE     = {"77_MT"}   # admixed individual — exclude from both clusters
# ─────────────────────────────────────────────────────────────────────────────

SYN_EFFECTS = {"synonymous_variant", "stop_retained_variant"}
NS_EFFECTS  = {
    "missense_variant", "stop_gained", "stop_lost", "start_lost",
    "start_gained", "frameshift_variant",
    "disruptive_inframe_insertion", "disruptive_inframe_deletion",
    "conservative_inframe_insertion", "conservative_inframe_deletion",
    "splice_acceptor_variant", "splice_donor_variant",
}
CDS_EFFECTS = SYN_EFFECTS | NS_EFFECTS
IMPACT_RANK = {"HIGH": 0, "MODERATE": 1, "LOW": 2, "MODIFIER": 3}


def best_cds_effect(info, alt_allele):
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
        matched = next((e for e in effect.split("&") if e in CDS_EFFECTS), None)
        if matched:
            rec = (allele, matched, impact, gene)
            (candidates if allele == alt_allele else fallback).append(rec)
    pool = candidates if candidates else fallback
    if not pool:
        return None, None
    best = min(pool, key=lambda x: IMPACT_RANK.get(x[2], 99))
    return best[1], best[2]


def call_haplotype(sample, threshold):
    gt = sample.get("GT", ".")
    if gt in (".", "./.", ".|."):
        return "."
    ad_raw = sample.get("AD", None)
    if ad_raw is None:
        return "."
    try:
        ad = [int(x) if x is not None else 0 for x in ad_raw] \
             if isinstance(ad_raw, (tuple, list)) else [int(ad_raw)]
    except (ValueError, TypeError):
        return "."
    if len(ad) < 2:
        return "."
    dp_raw = sample.get("DP", None)
    try:
        dp = int(dp_raw) if dp_raw is not None else sum(ad)
    except (ValueError, TypeError):
        dp = sum(ad)
    if not dp:
        return "."
    if ad[1] > threshold * dp:
        return 1
    elif ad[0] > threshold * dp:
        return 0
    return "."


def pi_site(n_called, n_alt):
    if n_called < 2:
        return 0.0
    p = n_alt / n_called
    return (n_called / (n_called - 1)) * 2 * p * (1 - p)


def harmonic(n):
    return sum(1/i for i in range(1, n))


def main():
    # ── Assign cluster membership ─────────────────────────────────────────────
    north, south = set(), set()
    with open(BURDEN) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            ind  = row["Individual"]
            nvt  = int(row["n_variants_total"])
            if ind in EXCLUDE:
                continue
            if nvt < NORTH_MAX:
                north.add(ind)
            elif nvt > SOUTH_MIN:
                south.add(ind)

    print(f"North cluster (<{NORTH_MAX} variants): {len(north)} samples", file=sys.stderr)
    print(f"South cluster (>{SOUTH_MIN} variants): {len(south)} samples", file=sys.stderr)
    print(f"Excluded: {EXCLUDE}", file=sys.stderr)

    # ── Open VCF ─────────────────────────────────────────────────────────────
    vcf      = pysam.VariantFile(VCF)
    all_samp = list(vcf.header.samples)

    # Map cluster sample sets to VCF indices
    north_idx = [i for i, s in enumerate(all_samp) if s in north]
    south_idx = [i for i, s in enumerate(all_samp) if s in south]
    print(f"North samples in VCF: {len(north_idx)}", file=sys.stderr)
    print(f"South samples in VCF: {len(south_idx)}", file=sys.stderr)

    # Accumulate per cluster
    clusters = {
        "north": {"idx": north_idx, "N": len(north_idx),
                  "pi_t": 0.0, "pi_s": 0.0, "pi_n": 0.0,
                  "S": 0, "S_syn": 0, "S_ns": 0},
        "south": {"idx": south_idx, "N": len(south_idx),
                  "pi_t": 0.0, "pi_s": 0.0, "pi_n": 0.0,
                  "S": 0, "S_syn": 0, "S_ns": 0},
    }

    fmt_keys = None
    rows_out = []

    for rec in vcf.fetch():
        if fmt_keys is None:
            fmt_keys = list(rec.format.keys())
        if rec.alts is None:
            continue

        for alt in rec.alts:
            if alt == "*":
                continue
            effect, impact = best_cds_effect(rec.info, alt)
            if effect is None:
                continue
            if effect in SYN_EFFECTS:
                site_class = "synonymous"
            elif effect in NS_EFFECTS:
                site_class = "nonsynonymous"
            else:
                continue

            # Call haplotypes for all samples once
            calls = []
            for sname in all_samp:
                s   = rec.samples[sname]
                fmt = {k: s[k] for k in fmt_keys if k in rec.format}
                calls.append(call_haplotype(fmt, THRESHOLD))

            row_entry = {
                "POS": rec.pos, "REF": rec.ref, "ALT": alt,
                "Effect": effect, "Class": site_class,
            }

            for cname, cd in clusters.items():
                n_called = n_alt = 0
                for i in cd["idx"]:
                    c = calls[i]
                    if c != ".":
                        n_called += 1
                        if c == 1:
                            n_alt += 1

                if n_called < 2:
                    row_entry[f"{cname}_n"] = n_called
                    row_entry[f"{cname}_freq"] = "."
                    row_entry[f"{cname}_pi"] = 0.0
                    continue

                freq   = n_alt / n_called
                p_site = pi_site(n_called, n_alt)
                cd["pi_t"] += p_site
                if n_alt > 0:
                    cd["S"] += 1
                if site_class == "synonymous":
                    cd["pi_s"] += p_site
                    if n_alt > 0:
                        cd["S_syn"] += 1
                else:
                    cd["pi_n"] += p_site
                    if n_alt > 0:
                        cd["S_ns"] += 1

                row_entry[f"{cname}_n"]    = n_called
                row_entry[f"{cname}_freq"] = round(freq, 6)
                row_entry[f"{cname}_pi"]   = round(p_site, 8)

            rows_out.append(row_entry)

    vcf.close()

    # ── Write per-site table ──────────────────────────────────────────────────
    fields = ["POS","REF","ALT","Effect","Class",
              "north_n","north_freq","north_pi",
              "south_n","south_freq","south_pi"]
    with open(OUT_SITE, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t",
                           extrasaction="ignore")
        w.writeheader()
        w.writerows(rows_out)

    # ── Recompute S_syn / S_ns cleanly from rows_out ──────────────────────────
    # (more reliable than in-loop accumulation)
    for cname, cd in clusters.items():
        pi_col = f"{cname}_pi"
        cd["pi_t"] = sum(float(r[pi_col]) for r in rows_out if pi_col in r)
        cd["pi_s"] = sum(float(r[pi_col]) for r in rows_out
                         if r.get("Class") == "synonymous" and pi_col in r)
        cd["pi_n"] = sum(float(r[pi_col]) for r in rows_out
                         if r.get("Class") == "nonsynonymous" and pi_col in r)
        # S: sites where this cluster has a non-zero freq (segregating in cluster)
        freq_col = f"{cname}_freq"
        cd["S"]     = sum(1 for r in rows_out
                          if r.get(freq_col) not in (".", None, "", "0.0", 0.0)
                          and float(r.get(freq_col, 0)) > 0)
        cd["S_syn"] = sum(1 for r in rows_out
                          if r.get("Class") == "synonymous"
                          and r.get(freq_col) not in (".", None, "", "0.0", 0.0)
                          and float(r.get(freq_col, 0)) > 0)
        cd["S_ns"]  = sum(1 for r in rows_out
                          if r.get("Class") == "nonsynonymous"
                          and r.get(freq_col) not in (".", None, "", "0.0", 0.0)
                          and float(r.get(freq_col, 0)) > 0)

    # ── Summary ───────────────────────────────────────────────────────────────
    lines = [f"=== π and Watterson's θ — Fhet mt CDS by haplogroup cluster ===",
             f"  VCF             : {VCF}",
             f"  L_CDS           : {L_CDS:,} bp",
             f"  Excluded        : {EXCLUDE}",
             f""]

    for cname, cd in clusters.items():
        N   = cd["N"]
        a1  = harmonic(N)
        pi_t = cd["pi_t"] / L_CDS
        pi_s = cd["pi_s"] / L_CDS
        pi_n = cd["pi_n"] / L_CDS
        th_t = cd["S"]     / (a1 * L_CDS) if cd["S"]     > 0 else float("nan")
        th_s = cd["S_syn"] / (a1 * L_CDS) if cd["S_syn"] > 0 else float("nan")
        th_n = cd["S_ns"]  / (a1 * L_CDS) if cd["S_ns"]  > 0 else float("nan")
        pnps_pi = pi_n / pi_s if pi_s > 0 else float("nan")
        pnps_th = th_n / th_s if th_s > 0 else float("nan")
        label = "North (< 50 variants)" if cname == "north" else "South (> 200 variants)"
        lines += [
            f"  {label}",
            f"    N               : {N}",
            f"    a₁              : {a1:.4f}",
            f"    S_total         : {cd['S']}",
            f"    S_syn           : {cd['S_syn']}",
            f"    S_ns            : {cd['S_ns']}",
            f"    θ_W_total       : {th_t:.5f}",
            f"    θ_W_syn         : {th_s:.5f}",
            f"    θ_W_ns          : {th_n:.5f}",
            f"    pN/pS (θ_W)     : {pnps_th:.3f}",
            f"    π_total         : {pi_t:.5f}",
            f"    π_syn           : {pi_s:.5f}",
            f"    π_ns            : {pi_n:.5f}",
            f"    pN/pS (π)       : {pnps_pi:.3f}",
            f"",
        ]

    lines.append(f"Per-site table  : {OUT_SITE}")
    summary = "\n".join(lines)
    print("\n" + summary)
    with open(OUT_SUM, "w") as fh:
        fh.write(summary + "\n")


if __name__ == "__main__":
    main()

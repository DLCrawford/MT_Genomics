#!/usr/bin/env python3
"""
26_amr_pi.py
Calculate nucleotide diversity (π) for the AMR (Admixed American) population
from gnomAD v3.1 mitochondrial VCF.

Source: gnomad.genomes.v3.1.sites.chrM.vcf.bgz
        (gnomad.broadinstitute.org/downloads)

Population order in pop_AN / pop_AC_hom / pop_AC_het:
  ['afr', 'ami', 'amr', 'asj', 'eas', 'fin', 'nfe', 'oth', 'sas', 'mid']
  AMR index = 2

Approach (haploid mt):
  - Use AC_hom (≥95% heteroplasmy) as the ALT allele count → effectively haploid
  - Exclude AC_het (heteroplasmic) — not true population-level segregation
  - Filter: PASS sites, SNPs only, CDS sites only (protein_coding VEP)
  - π_site = (n / n-1) * 2 * p * (1-p)  where n = pop_AN[amr], p = AC_hom[amr]/n
  - L_CDS = 11,395 bp (13 human mt protein-coding genes, rCRS annotation)

Usage (from Hm_Mt/):
    conda activate SNP_env
    python 26_amr_pi.py

Outputs:
    amr_pi_summary.txt
    amr_pi_per_site.tsv
"""

import gzip
import csv
import sys
import os

# ── Config ────────────────────────────────────────────────────────────────────
VCF      = "gnomad.genomes.v3.1.sites.chrM.vcf.bgz"
OUT_SUM  = "amr_pi_summary.txt"
OUT_SITE = "amr_pi_per_site.tsv"
L_CDS    = 11_395
AMR_IDX  = 2      # index in pop_AN / pop_AC_hom lists
# ─────────────────────────────────────────────────────────────────────────────

SYN_TERMS = {
    "synonymous_variant",
    "stop_retained_variant",
}
NS_TERMS = {
    "missense_variant",
    "stop_gained",
    "stop_lost",
    "start_lost",
    "start_gained",
    "frameshift_variant",
}
CDS_TERMS = SYN_TERMS | NS_TERMS

IMPACT_RANK = {"HIGH": 0, "MODERATE": 1, "LOW": 2, "MODIFIER": 3}


def parse_info(info_str):
    """Parse INFO field → dict."""
    d = {}
    for token in info_str.split(";"):
        if "=" in token:
            k, v = token.split("=", 1)
            d[k] = v
    return d


def get_pop_val(field_str, idx):
    """Extract numeric value at population index from pipe-separated string."""
    if not field_str:
        return 0
    parts = field_str.split("|")
    try:
        return int(parts[idx])
    except (IndexError, ValueError):
        return 0


def best_cds_effect(vep_str):
    """
    Parse VEP field → return (consequence, impact) for best CDS annotation.
    VEP format: Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|Biotype|...
    Returns (None, None) if no CDS annotation found.
    """
    if not vep_str:
        return None, None

    candidates = []
    for entry in vep_str.split(","):
        fields = entry.split("|")
        if len(fields) < 8:
            continue
        biotype    = fields[7]
        impact     = fields[2]
        conseq_raw = fields[1]
        if biotype != "protein_coding":
            continue
        # Handle &-joined consequences
        for c in conseq_raw.split("&"):
            if c in CDS_TERMS:
                candidates.append((c, impact))
                break

    if not candidates:
        return None, None
    best = min(candidates, key=lambda x: IMPACT_RANK.get(x[1], 99))
    return best


def pi_site(n, n_alt):
    """Haploid π at one site."""
    if n < 2 or n_alt == 0:
        return 0.0
    p = n_alt / n
    return (n / (n - 1)) * 2 * p * (1 - p)


def main():
    if not os.path.exists(VCF):
        print(f"ERROR: {VCF} not found. Run from Hm_Mt/.", file=sys.stderr)
        sys.exit(1)

    print(f"Reading {VCF} ...", file=sys.stderr)

    pi_total = pi_syn = pi_ns = 0.0
    n_snp = n_cds = n_syn = n_ns = 0
    # S_amr: segregating in AMR (AC_hom > 0) — used for θ_W
    s_amr = s_syn_amr = s_ns_amr = 0
    n_amr_max = 0
    rows = []

    opener = gzip.open if VCF.endswith(".gz") or VCF.endswith(".bgz") else open

    with opener(VCF, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")
            chrom, pos, _, ref, alt = fields[0], int(fields[1]), fields[2], fields[3], fields[4]
            filt = fields[6]
            info_str = fields[7]

            # PASS only
            if filt != "PASS":
                continue

            # SNPs only
            if len(ref) != 1 or len(alt) != 1 or alt == "*":
                continue

            n_snp += 1
            info = parse_info(info_str)

            # AMR allele counts
            n_amr  = get_pop_val(info.get("pop_AN", ""),      AMR_IDX)
            ac_hom = get_pop_val(info.get("pop_AC_hom", ""),  AMR_IDX)

            if n_amr < 2:
                continue
            n_amr_max = max(n_amr_max, n_amr)

            # CDS classification from VEP
            conseq, impact = best_cds_effect(info.get("vep", ""))
            if conseq is None:
                continue

            if conseq in SYN_TERMS:
                site_class = "synonymous"
            elif conseq in NS_TERMS:
                site_class = "nonsynonymous"
            else:
                continue

            n_cds += 1
            p_site = pi_site(n_amr, ac_hom)
            freq   = ac_hom / n_amr

            pi_total += p_site
            if site_class == "synonymous":
                pi_syn += p_site
                n_syn  += 1
                if ac_hom > 0:
                    s_syn_amr += 1
            else:
                pi_ns  += p_site
                n_ns   += 1
                if ac_hom > 0:
                    s_ns_amr += 1
            if ac_hom > 0:
                s_amr += 1

            rows.append({
                "CHROM":      chrom,
                "POS":        pos,
                "REF":        ref,
                "ALT":        alt,
                "Consequence":conseq,
                "Impact":     impact,
                "Class":      site_class,
                "n_amr":      n_amr,
                "AC_hom_amr": ac_hom,
                "freq_alt":   round(freq, 6),
                "pi_site":    round(p_site, 8),
            })

    # Normalise
    pi_t = pi_total / L_CDS
    pi_s = pi_syn   / L_CDS
    pi_n = pi_ns    / L_CDS
    pnps = pi_n / pi_s if pi_s > 0 else float("nan")

    # Watterson's θ — use only sites segregating in AMR (AC_hom > 0)
    N    = n_amr_max
    a1   = sum(1/i for i in range(1, N))
    th_t = s_amr     / (a1 * L_CDS)
    th_s = s_syn_amr / (a1 * L_CDS)
    th_n = s_ns_amr  / (a1 * L_CDS)
    pnps_th = th_n / th_s if th_s > 0 else float("nan")

    # Write per-site TSV
    fields_out = ["CHROM","POS","REF","ALT","Consequence","Impact","Class",
                  "n_amr","AC_hom_amr","freq_alt","pi_site"]
    with open(OUT_SITE, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields_out, delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    summary = (
        f"=== π and Watterson's θ — gnomAD v3.1 AMR mt CDS ===\n"
        f"  Source          : gnomAD v3.1 chrM (pop_AC_hom, AMR index={AMR_IDX})\n"
        f"  N (max pop_AN)  : {N:,}\n"
        f"  a₁              : {a1:.4f}\n"
        f"  L_CDS           : {L_CDS:,} bp\n"
        f"  Total PASS SNPs : {n_snp:,}\n"
        f"  CDS SNPs (any pop)   : {n_cds}  (syn={n_syn}, NS={n_ns})\n"
        f"  S_AMR (AC_hom>0)     : {s_amr}  (syn={s_syn_amr}, NS={s_ns_amr})\n"
        f"\n"
        f"  π_total         : {pi_t:.5f}\n"
        f"  π_syn           : {pi_s:.5f}\n"
        f"  π_ns            : {pi_n:.5f}\n"
        f"  pN/pS (π)       : {pnps:.3f}\n"
        f"\n"
        f"  θ_W_total       : {th_t:.5f}\n"
        f"  θ_W_syn         : {th_s:.5f}\n"
        f"  θ_W_ns          : {th_n:.5f}\n"
        f"  pN/pS (θ_W)     : {pnps_th:.3f}\n"
        f"\nNote: AC_het (heteroplasmic) excluded; AC_hom (≥95%) only.\n"
        f"Per-site table  : {OUT_SITE}\n"
    )
    print("\n" + summary)
    with open(OUT_SUM, "w") as fh:
        fh.write(summary)


if __name__ == "__main__":
    main()

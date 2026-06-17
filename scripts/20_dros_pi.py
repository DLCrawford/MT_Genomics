#!/usr/bin/env python3
"""
20_dros_pi.py
Calculate nucleotide diversity (π) for Drosophila DGRP mitochondrial
coding variants from Dros_Mt_coding.csv.

Data format:
  - Rows: variant sites (synonymous, missense, MNP, or unclassified NaN)
  - Columns 12+: 169 DGRP line genotype calls (0=REF, >=1=ALT, NaN=missing)
  - Values >1 present in 43 rows (likely heteroplasmy / data encoding);
    treated as ALT (binary presence/absence).

Formula (haploid):
    π_site = (n / n-1) * 2 * p * (1-p)
    π      = Σ π_site / L_CDS

Usage (from MT_Genomics2/):
    conda activate SNP_env
    python scripts/20_dros_pi.py

Reference:
    Mackay et al. 2012 Nature — Drosophila Genetic Reference Panel (DGRP)
    L_CDS = 11,173 bp (Drosophila melanogaster mitochondrial coding sequence)
"""

import pandas as pd
import numpy as np
import sys
import os

# ── Config ────────────────────────────────────────────────────────────────────
DATA_FILE = "data_files_May/Dros_Mt_coding.csv"
OUT_FILE  = "vcf/dros_pi_results.tsv"
L_CDS     = 11_173   # bp — Drosophila mt CDS length (Table 2)
# ─────────────────────────────────────────────────────────────────────────────


def pi_site(n_called, n_alt):
    """Haploid π at one site: (n/n-1) * 2 * p * (1-p)"""
    if n_called < 2:
        return 0.0
    p = n_alt / n_called
    return (n_called / (n_called - 1)) * 2 * p * (1 - p)


def main():
    if not os.path.exists(DATA_FILE):
        print(f"ERROR: {DATA_FILE} not found. Run from MT_Genomics2/.", file=sys.stderr)
        sys.exit(1)

    df = pd.read_csv(DATA_FILE)
    sample_cols = df.columns[11:]
    n_samples   = len(sample_cols)

    print(f"Input          : {DATA_FILE}", file=sys.stderr)
    print(f"Total rows     : {len(df)}", file=sys.stderr)
    print(f"Sample columns : {n_samples}", file=sys.stderr)
    print(f"FuncType counts:\n{df['FuncType'].value_counts(dropna=False)}\n", file=sys.stderr)

    # ── Filter: exclude MNP rows ──────────────────────────────────────────────
    df_snp = df[df['FuncType'] != 'MNP'].copy()
    print(f"After removing MNP: {len(df_snp)} SNP rows", file=sys.stderr)

    # ── Per-site π ────────────────────────────────────────────────────────────
    rows_out = []
    pi_acc   = {'total': 0.0, 'syn': 0.0, 'ns': 0.0, 'unclassified': 0.0}
    n_sites  = {'total': 0,   'syn': 0,   'ns': 0,   'unclassified': 0}

    for _, row in df_snp.iterrows():
        gts    = pd.to_numeric(row[sample_cols], errors='coerce')
        called = gts.dropna()
        n      = len(called)
        if n < 2:
            continue

        n_alt   = int((called >= 1).sum())   # treat >=1 as ALT
        freq    = n_alt / n
        pi_val  = pi_site(n, n_alt)
        ftype   = row['FuncType']

        if pd.isna(ftype):
            site_class = 'unclassified'
        elif ftype == 'synonymous':
            site_class = 'synonymous'
        elif ftype == 'missense':
            site_class = 'nonsynonymous'
        else:
            site_class = 'unclassified'

        pi_acc['total']  += pi_val
        n_sites['total'] += 1
        if site_class == 'synonymous':
            pi_acc['syn']  += pi_val
            n_sites['syn'] += 1
        elif site_class == 'nonsynonymous':
            pi_acc['ns']   += pi_val
            n_sites['ns']  += 1
        else:
            pi_acc['unclassified']  += pi_val
            n_sites['unclassified'] += 1

        rows_out.append({
            'Chrom':     row.get('Chromosome', '.'),
            'Pos':       row.get('Pos', '.'),
            'Ref':       row.get('RefAllele', '.'),
            'Alt':       row.get('AltAllele', '.'),
            'Gene':      row.get('Gene', '.'),
            'FuncType':  ftype,
            'Class':     site_class,
            'n_called':  n,
            'n_alt':     n_alt,
            'freq_alt':  round(freq, 6),
            'pi_site':   round(pi_val, 8),
        })

    # ── Normalise ─────────────────────────────────────────────────────────────
    # Re-do accumulation cleanly
    pi_total = sum(r['pi_site'] for r in rows_out)
    pi_syn   = sum(r['pi_site'] for r in rows_out if r['Class'] == 'synonymous')
    pi_ns    = sum(r['pi_site'] for r in rows_out if r['Class'] == 'nonsynonymous')
    pi_unc   = sum(r['pi_site'] for r in rows_out if r['Class'] == 'unclassified')

    n_total = len(rows_out)
    n_syn   = sum(1 for r in rows_out if r['Class'] == 'synonymous')
    n_ns    = sum(1 for r in rows_out if r['Class'] == 'nonsynonymous')
    n_unc   = sum(1 for r in rows_out if r['Class'] == 'unclassified')

    pi_t_norm = pi_total / L_CDS
    pi_s_norm = pi_syn   / L_CDS
    pi_n_norm = pi_ns    / L_CDS
    pn_ps     = (pi_n_norm / pi_s_norm) if pi_s_norm > 0 else float('nan')

    # ── Watterson's θ ─────────────────────────────────────────────────────────
    a1    = sum(1/i for i in range(1, n_samples))
    th_t  = n_total / (a1 * L_CDS)
    th_s  = n_syn   / (a1 * L_CDS)
    th_n  = n_ns    / (a1 * L_CDS)
    pnps_th = th_n / th_s if th_s > 0 else float('nan')

    # ── Write per-site table ──────────────────────────────────────────────────
    os.makedirs(os.path.dirname(OUT_FILE), exist_ok=True)
    pd.DataFrame(rows_out).to_csv(OUT_FILE, sep='\t', index=False)
    print(f"Per-site table : {OUT_FILE}  ({len(rows_out)} rows)", file=sys.stderr)

    # ── Summary ───────────────────────────────────────────────────────────────
    print("\n" + "="*55)
    print("  π and Watterson's θ — Drosophila DGRP mt CDS")
    print("="*55)
    print(f"  Samples (N)          : {n_samples}")
    print(f"  a₁                   : {a1:.4f}")
    print(f"  L_CDS                : {L_CDS:,} bp")
    print(f"  CDS variant sites (S): {n_total}  (syn={n_syn}, NS={n_ns}, unclassified={n_unc})")
    print(f"")
    print(f"  θ_W_total            : {th_t:.5f}")
    print(f"  θ_W_syn              : {th_s:.5f}")
    print(f"  θ_W_ns               : {th_n:.5f}")
    print(f"  pN/pS (θ_W)          : {pnps_th:.3f}")
    print(f"")
    print(f"  π_total              : {pi_t_norm:.5f}")
    print(f"  π_syn                : {pi_s_norm:.5f}")
    print(f"  π_ns                 : {pi_n_norm:.5f}")
    print(f"  pN/pS (π)            : {pn_ps:.3f}")
    print("="*55)
    if n_unc > 0:
        print(f"\nNote: {n_unc} unclassified rows (NaN FuncType) contribute to")
        print(f"  π_total/θ_W_total but not syn or NS.")
    print(f"  Values >1 in genotype matrix treated as ALT (binary).")


if __name__ == "__main__":
    main()

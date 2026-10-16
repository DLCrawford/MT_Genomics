#!/usr/bin/env python3
"""
25_comparison_table.py
Build unified π + Watterson's θ comparison table across all datasets.

Reads per-site TSV outputs from scripts 19–26 and computes:
  - S (segregating sites total, syn, NS)
  - a₁ (harmonic number for each dataset's N)
  - θ_W = S / (a₁ × L_CDS)
  - π   = Σ π_site / L_CDS  (from per-site TSVs)
  - pN/pS from both θ_W and π

For gnomAD AMR: S is counted only where AC_hom_amr > 0 (truly segregating
in AMR); all-population sites with AMR AC=0 are excluded from θ_W but
contribute 0 to π anyway.

Usage (from MT_Genomics2/):
    conda activate SNP_env
    python scripts/25_comparison_table.py

Output:
    vcf/comparison_table.tsv   — machine-readable
    vcf/comparison_table.txt   — formatted for publication
"""

import csv
import sys
import os

BASE = os.path.expanduser("~/Projects/MT_Genomics_Cl_Ap2026")

# ── Dataset definitions ───────────────────────────────────────────────────────
# s_filter: optional (col, min_val) — row only counts toward S if col > min_val
DATASETS = [
    {
        "name":      "Fhet",
        "N":         141,
        "L_CDS":     11_417,
        "tsv":       f"{BASE}/MT_Genomics2/vcf/pi_results.tsv",
        "class_col": "Class",
        "pi_col":    "pi_site",
        "syn_val":   "synonymous",
        "ns_val":    "nonsynonymous",
        "s_filter":  None,
    },
    {
        "name":      "Drosophila",
        "N":         169,
        "L_CDS":     11_173,
        "tsv":       f"{BASE}/MT_Genomics2/vcf/dros_pi_results.tsv",
        "class_col": "Class",
        "pi_col":    "pi_site",
        "syn_val":   "synonymous",
        "ns_val":    "nonsynonymous",
        "s_filter":  None,
    },
    {
        "name":      "C. elegans",
        "N":         540,
        "L_CDS":     10_299,
        "tsv":       f"{BASE}/C_elegans/celegans_pi_per_site.tsv",
        "class_col": "Class",
        "pi_col":    "pi_site",
        "syn_val":   "synonymous",
        "ns_val":    "nonsynonymous",
        "s_filter":  None,
    },
    {
        "name":      "Human African (Lankheet 2026)",
        "N":         1_176,
        "L_CDS":     11_395,
        "tsv":       f"{BASE}/Human_mt/human_mt_cds_pi_per_site.tsv",
        "class_col": "Effect",
        "pi_col":    "pi_site",
        "syn_val":   "SYN",
        "ns_val":    "NS",
        "s_filter":  None,
    },
    {
        "name":      "Human AMR (gnomAD v3.1)",
        "N":         5_718,
        "L_CDS":     11_395,
        "tsv":       f"{BASE}/Hm_Mt/amr_pi_per_site.tsv",
        "class_col": "Class",
        "pi_col":    "pi_site",
        "syn_val":   "synonymous",
        "ns_val":    "nonsynonymous",
        "s_filter":  ("AC_hom_amr", 0),   # only count AMR-segregating sites
    },
    {
        "name":      "Yeast (lit. θ_W only)",
        "N":         1_011,
        "L_CDS":     6_684,
        "tsv":       None,               # π unreliable; literature values used
        "class_col": None,
        "pi_col":    "pi_site",
        "syn_val":   None,
        "ns_val":    None,
        "s_filter":  None,
        # Hardcoded from De Chiara et al. 2020 (Table 2)
        "lit_S":     384,
        "lit_theta": 0.00766,
    },
]
# ─────────────────────────────────────────────────────────────────────────────


def harmonic(n):
    return sum(1/i for i in range(1, n))


def load_tsv(path):
    if not path or not os.path.exists(path):
        return None
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def compute_stats(ds):
    L    = ds["L_CDS"]
    N    = ds["N"]
    a1   = harmonic(N)
    cc   = ds["class_col"]
    pc   = ds["pi_col"]
    sf   = ds["s_filter"]    # (col, min_val) or None

    # ── Yeast: use literature values ──────────────────────────────────────────
    if ds.get("lit_theta"):
        S_lit = ds["lit_S"]
        th_t  = ds["lit_theta"]
        return {
            "Dataset":     ds["name"],
            "N":           N,
            "a1":          round(a1, 4),
            "L_CDS":       L,
            "S":           S_lit,
            "S_syn":       "N/A",
            "S_ns":        "N/A",
            "theta_total": th_t,
            "theta_syn":   "N/A",
            "theta_ns":    "N/A",
            "pNpS_theta":  "N/A",
            "pi_total":    "N/A",
            "pi_syn":      "N/A",
            "pi_ns":       "N/A",
            "pNpS_pi":     "N/A",
        }

    rows = load_tsv(ds["tsv"])
    if rows is None:
        print(f"  WARNING: TSV not found — {ds['tsv']}", file=sys.stderr)
        return None

    pi_total = pi_syn = pi_ns = 0.0
    S = S_syn = S_ns = 0

    for r in rows:
        try:
            pv = float(r[pc])
        except (KeyError, ValueError):
            continue

        # Determine if this row counts toward S
        counts_as_seg = True
        if sf:
            try:
                counts_as_seg = float(r[sf[0]]) > sf[1]
            except (KeyError, ValueError):
                counts_as_seg = False

        pi_total += pv

        if counts_as_seg:
            S += 1

        if cc and r.get(cc) == ds["syn_val"]:
            pi_syn += pv
            if counts_as_seg:
                S_syn += 1
        elif cc and r.get(cc) == ds["ns_val"]:
            pi_ns += pv
            if counts_as_seg:
                S_ns += 1

    pi_t = pi_total / L
    pi_s = pi_syn   / L if S_syn > 0 else float("nan")
    pi_n = pi_ns    / L if S_ns  > 0 else float("nan")

    theta_t = S     / (a1 * L) if S     > 0 else float("nan")
    theta_s = S_syn / (a1 * L) if S_syn > 0 else float("nan")
    theta_n = S_ns  / (a1 * L) if S_ns  > 0 else float("nan")

    pnps_pi    = pi_n    / pi_s    if pi_s    > 0 else float("nan")
    pnps_theta = theta_n / theta_s if theta_s > 0 else float("nan")

    def fmt5(v): return round(v, 5) if isinstance(v, float) and not (v != v) else "N/A"
    def fmt3(v): return round(v, 3) if isinstance(v, float) and not (v != v) else "N/A"

    return {
        "Dataset":     ds["name"],
        "N":           N,
        "a1":          round(a1, 4),
        "L_CDS":       L,
        "S":           S,
        "S_syn":       S_syn if cc else "N/A",
        "S_ns":        S_ns  if cc else "N/A",
        "theta_total": fmt5(theta_t),
        "theta_syn":   fmt5(theta_s) if cc else "N/A",
        "theta_ns":    fmt5(theta_n) if cc else "N/A",
        "pNpS_theta":  fmt3(pnps_theta) if cc else "N/A",
        "pi_total":    fmt5(pi_t),
        "pi_syn":      fmt5(pi_s) if cc else "N/A",
        "pi_ns":       fmt5(pi_n) if cc else "N/A",
        "pNpS_pi":     fmt3(pnps_pi) if cc else "N/A",
    }


def main():
    results = []
    for ds in DATASETS:
        print(f"Processing {ds['name']} ...", file=sys.stderr)
        stat = compute_stats(ds)
        if stat is None:
            print(f"  Skipping.", file=sys.stderr)
        else:
            results.append(stat)

    if not results:
        print("No results.", file=sys.stderr)
        sys.exit(1)

    # ── Machine-readable TSV ──────────────────────────────────────────────────
    out_tsv = "vcf/comparison_table.tsv"
    fields  = ["Dataset","N","a1","L_CDS",
               "S","S_syn","S_ns",
               "theta_total","theta_syn","theta_ns","pNpS_theta",
               "pi_total","pi_syn","pi_ns","pNpS_pi"]
    with open(out_tsv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        w.writeheader()
        w.writerows(results)

    # ── Formatted text table ──────────────────────────────────────────────────
    out_txt = "vcf/comparison_table.txt"
    col_w = [28, 6, 7, 6, 6, 6, 6, 8, 8, 8, 8, 8, 8, 8, 8]
    hdr = ["Dataset","N","a₁","L_CDS","S","S_syn","S_ns",
           "θ_tot","θ_syn","θ_ns","pN/pS_θ",
           "π_tot","π_syn","π_ns","pN/pS_π"]

    def row_str(vals):
        return "  ".join(str(v).rjust(w) for v, w in zip(vals, col_w))

    sep   = "-" * (sum(col_w) + 2 * (len(col_w) - 1))
    lines = [sep, row_str(hdr), sep]
    for r in results:
        lines.append(row_str([
            r["Dataset"], r["N"], r["a1"], r["L_CDS"],
            r["S"], r["S_syn"], r["S_ns"],
            r["theta_total"], r["theta_syn"], r["theta_ns"], r["pNpS_theta"],
            r["pi_total"], r["pi_syn"], r["pi_ns"], r["pNpS_pi"],
        ]))
    lines.append(sep)
    lines.append("\nNote: Yeast π unreliable (Biopython alignment artefact); θ_W from De Chiara et al. 2020.")
    lines.append("      Human African θ_W from Lankheet et al. 2026 (N=1,176); AMR θ_W from gnomAD v3.1 (N=5,718).")

    table_str = "\n".join(lines)
    print("\n" + table_str)
    with open(out_txt, "w") as fh:
        fh.write(table_str + "\n")

    print(f"\nTSV : {out_tsv}", file=sys.stderr)
    print(f"TXT : {out_txt}", file=sys.stderr)


if __name__ == "__main__":
    main()

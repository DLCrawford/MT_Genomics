#!/usr/bin/env python3
"""
28_mk_test.py
McDonald-Kreitman (MK) test for Fhet mt CDS using the reference genome
(NC_012312.1 / Fhet_MT.fasta) as the outgroup.

Fixed differences (D): sites where ALL called samples carry the ALT allele
    (freq_alt = 1.0) — divergence from reference.
Polymorphisms (P): sites where ALT is segregating (0 < freq_alt < 1).

MK table (per class):
          Polymorphic (P)   Fixed (D)
  Syn          Ps               Ds
  NS           Pn               Dn

Statistics:
  NI    = (Pn/Ps) / (Dn/Ds)        Neutrality Index (NI=1 → neutral)
  α     = 1 − NI                    Proportion adaptive (α>0 → + selection)
  DoS   = Dn/(Dn+Ds) − Pn/(Pn+Ps)  Direction of Selection
  Fisher's exact test (2×2 table)

Analyses:
  1. Overall Fhet (141 samples) from vcf/pi_results.tsv
  2. North cluster (<50 variants) from vcf/fhet_clusters_pi_per_site.tsv
  3. South cluster (>200 variants) from vcf/fhet_clusters_pi_per_site.tsv

CAVEAT: The reference (NC_012312.1) is conspecific (F. heteroclitus), not
a true outgroup. Fixed differences may include intraspecific divergence
from the sampled population. Ideally a congener (e.g. F. similis) would
serve as outgroup. Results should be interpreted accordingly.

Usage (from MT_Genomics2/):
    conda activate SNP_env
    python scripts/28_mk_test.py

Output:
    vcf/mk_test_results.txt
"""

import csv
import sys
import os
import random
from scipy.stats import fisher_exact

# ── Input files ───────────────────────────────────────────────────────────────
PI_TSV      = "vcf/pi_results.tsv"
CLUSTER_TSV = "vcf/fhet_clusters_pi_per_site.tsv"
OUT_FILE    = "vcf/mk_test_results.txt"
# ─────────────────────────────────────────────────────────────────────────────


def dos(Ps, Pn, Ds, Dn):
    """Direction of Selection: Dn/(Dn+Ds) - Pn/(Pn+Ps)."""
    if (Dn + Ds) == 0 or (Pn + Ps) == 0:
        return float('nan')
    return Dn / (Dn + Ds) - Pn / (Pn + Ps)


def bootstrap_dos(Ps, Pn, Ds, Dn, n_boot=10_000, seed=42, ci=95):
    """
    Bootstrap 95% CI for DoS by resampling sites within each of the
    4 MK categories (Ps, Pn, Ds, Dn) independently with replacement.
    Returns (DoS_observed, ci_lo, ci_hi, p_value).
    p_value = fraction of bootstrap replicates where DoS has the
    opposite sign to observed (two-tailed permutation-style).
    """
    rng = random.Random(seed)
    obs = dos(Ps, Pn, Ds, Dn)
    if obs != obs:   # nan
        return obs, float('nan'), float('nan'), float('nan')

    boot_vals = []
    for _ in range(n_boot):
        # Resample counts within each category using Poisson approximation
        # (equivalent to multinomial resampling for large N)
        ps_b = sum(rng.random() < Ps / (Ps + Pn) for _ in range(Ps + Pn))
        pn_b = (Ps + Pn) - ps_b
        ds_b = sum(rng.random() < Ds / max(Ds + Dn, 1) for _ in range(Ds + Dn))
        dn_b = (Ds + Dn) - ds_b
        v = dos(ps_b, pn_b, ds_b, dn_b)
        if v == v:   # not nan
            boot_vals.append(v)

    if not boot_vals:
        return obs, float('nan'), float('nan'), float('nan')

    boot_vals.sort()
    lo = (100 - ci) / 2
    hi = 100 - lo
    ci_lo = boot_vals[int(lo / 100 * len(boot_vals))]
    ci_hi = boot_vals[int(hi / 100 * len(boot_vals))]

    # Two-tailed p: fraction of bootstraps with DoS opposite sign to observed
    if obs > 0:
        p = sum(1 for v in boot_vals if v <= 0) / len(boot_vals) * 2
    else:
        p = sum(1 for v in boot_vals if v >= 0) / len(boot_vals) * 2
    p = min(p, 1.0)

    return obs, ci_lo, ci_hi, p


def mk_stats(Ps, Pn, Ds, Dn, n_boot=10_000):
    """Compute NI, α, DoS + bootstrap CI, Fisher p-value."""
    odds, pval = fisher_exact([[Ps, Ds], [Pn, Dn]], alternative='two-sided')
    NI    = (Pn / Ps) / (Dn / Ds) if (Ps > 0 and Ds > 0) else float('nan')
    alpha = 1 - NI if NI == NI else float('nan')
    DoS_obs, dos_lo, dos_hi, dos_p = bootstrap_dos(Ps, Pn, Ds, Dn, n_boot)
    return NI, alpha, DoS_obs, dos_lo, dos_hi, dos_p, pval, odds


def format_mk(label, Ps, Pn, Ds, Dn, N, n_boot=10_000, caveat=""):
    NI, alpha, DoS_obs, dos_lo, dos_hi, dos_p, pval, odds = \
        mk_stats(Ps, Pn, Ds, Dn, n_boot)

    def f4(v): return f"{v:.4f}" if v == v else "N/A"

    lines = [
        f"\n{'='*60}",
        f"  {label}  (N = {N})",
        f"{'='*60}",
        f"  MK table:",
        f"                      Polymorphic (P)   Fixed (D)",
        f"    Synonymous    :        {Ps:>6}           {Ds:>4}",
        f"    Nonsynonymous :        {Pn:>6}           {Dn:>4}",
        f"",
        f"  Pn/Ps           : {f4(Pn/Ps if Ps>0 else float('nan'))}",
        f"  Dn/Ds           : {f4(Dn/Ds if Ds>0 else float('nan'))}",
        f"  NI              : {f4(NI)}",
        f"  α (adaptive)    : {f4(alpha)}",
        f"  DoS             : {f4(DoS_obs)}  "
            f"(95% CI: {f4(dos_lo)} – {f4(dos_hi)},  bootstrap p = {f4(dos_p)})",
        f"  Fisher p-val    : {f4(pval)}",
        f"  Odds ratio      : {f4(odds)}",
        f"  Bootstrap reps  : {n_boot:,}",
    ]
    if caveat:
        lines.append(f"\n  NOTE: {caveat}")
    return "\n".join(lines)


def load_tsv(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter='\t'))


def classify_overall(rows):
    """Return Ps, Pn, Ds, Dn from overall pi_results.tsv."""
    Ps = Pn = Ds = Dn = 0
    for r in rows:
        freq  = float(r['freq_alt'])
        cls   = r['Class']
        if freq == 1.0:
            if cls == 'synonymous':    Ds += 1
            elif cls == 'nonsynonymous': Dn += 1
        elif 0 < freq < 1.0:
            if cls == 'synonymous':    Ps += 1
            elif cls == 'nonsynonymous': Pn += 1
    return Ps, Pn, Ds, Dn


def classify_cluster(rows, freq_col):
    """Return Ps, Pn, Ds, Dn for a cluster from fhet_clusters_pi_per_site.tsv."""
    Ps = Pn = Ds = Dn = 0
    for r in rows:
        freq_str = r.get(freq_col, ".")
        if freq_str in (".", "", None):
            continue
        try:
            freq = float(freq_str)
        except ValueError:
            continue
        cls = r['Class']
        if freq == 1.0:
            if cls == 'synonymous':      Ds += 1
            elif cls == 'nonsynonymous': Dn += 1
        elif 0 < freq < 1.0:
            if cls == 'synonymous':      Ps += 1
            elif cls == 'nonsynonymous': Pn += 1
    return Ps, Pn, Ds, Dn


def main():
    conspecific_caveat = (
        "Reference NC_012312.1 is conspecific (F. heteroclitus), not a true "
        "outgroup. Fixed differences reflect divergence from reference individual, "
        "not interspecific divergence. Use congener (e.g. F. similis) for a "
        "standard MK test."
    )

    output = ["McDonald-Kreitman Test — Fhet mt CDS",
              "Reference outgroup: NC_012312.1 (F. heteroclitus)",
              f"\n{conspecific_caveat}"]

    # ── 1. Overall Fhet (141 samples) ─────────────────────────────────────────
    if not os.path.exists(PI_TSV):
        print(f"ERROR: {PI_TSV} not found.", file=sys.stderr)
        sys.exit(1)

    rows_all = load_tsv(PI_TSV)
    Ps, Pn, Ds, Dn = classify_overall(rows_all)
    output.append(format_mk("Overall Fhet — 141 samples", Ps, Pn, Ds, Dn, 141, n_boot=10_000))

    # ── 2. North and South clusters ───────────────────────────────────────────
    if os.path.exists(CLUSTER_TSV):
        rows_cl = load_tsv(CLUSTER_TSV)

        # North
        Ps, Pn, Ds, Dn = classify_cluster(rows_cl, "north_freq")
        output.append(format_mk("North cluster — < 50 variants", Ps, Pn, Ds, Dn, "~77", n_boot=10_000))

        # South
        Ps, Pn, Ds, Dn = classify_cluster(rows_cl, "south_freq")
        output.append(format_mk("South cluster — > 200 variants (excl. 77_MT)", Ps, Pn, Ds, Dn, "~65", n_boot=10_000))
    else:
        output.append(f"\nCluster TSV not found ({CLUSTER_TSV}) — run script 27 first.")

    # ── Interpretation guide ──────────────────────────────────────────────────
    output.append("""
Interpretation:
  NI = 1     → neutral evolution
  NI < 1, α > 0 → excess fixed NS → positive/adaptive selection
  NI > 1, α < 0 → deficit fixed NS → slightly deleterious NS (purifying selection)
  DoS > 0    → NS substitutions favoured relative to polymorphism
  DoS < 0    → NS substitutions disfavoured (purifying selection)
  Fisher p < 0.05 → significant departure from neutrality
""")

    result_str = "\n".join(output)
    print(result_str)
    with open(OUT_FILE, "w") as fh:
        fh.write(result_str + "\n")
    print(f"\nResults written to {OUT_FILE}", file=sys.stderr)


if __name__ == "__main__":
    main()

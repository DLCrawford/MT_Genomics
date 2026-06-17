#!/usr/bin/env python3
"""
18_variant_burden_per_individual.py — Per-sample panel-variant counts
                                       (Total / SYN / NS / Other), joined
                                       with per-individual heteroplasmy
                                       burden. Enables the test of
                                       "do low-variant samples carry
                                       more heteroplasmy?"

Background
----------
The *F. heteroclitus* panel shows a bimodal per-individual variant
count: north-clade individuals carry < ~50 panel ALTs vs the GenBank
reference (NC_012312.1, north), while south-clade individuals carry
> ~190. We want to test whether the low-variant-count (north)
samples carry significantly more heteroplasmy than high-variant-count
(south) samples — or vice versa.

To do that we need per-sample:

  (a) `n_variants_*`  — how many panel ALT calls (AF ≥ 0.7) does this
                        sample carry, broken down by Effect_class
                        from SnpEff `ANN[0]` (the same SYN / NS /
                        Other buckets used by scripts 10, 11, 16, 17).
  (b) `n_hp_*`        — how many heteroplasmy events this sample has,
                        from stage 09 (variant-only) and/or stage 14
                        (CDS pileup). Zero-filled for samples that
                        don't appear in the per_individual file.

Inputs
------
    vcf/Fhet_MT_CDS.snps.split.vcf.gz   canonical, SnpEff-annotated, ploidy=1
    vcf/heteroplasmy_per_individual.tsv             (stage 09 — optional)
    vcf/heteroplasmy_pileup_per_individual_all.tsv  (stage 14 — optional)

Outputs
-------
    vcf/per_individual_burden_pileup.tsv    all canonical samples × {variant
                                            counts, stage-14 Hp counts}
    vcf/per_individual_burden_variant.tsv   all canonical samples × {variant
                                            counts, stage-09 Hp counts}
    vcf/per_individual_burden_summary.txt   bimodality check + 2x2-style
                                            (low/high variant × Hp burden)

Per-sample variant-count rule (mtDNA, ploidy = 1)
-------------------------------------------------
For each panel variant row, the per-sample GT is one of:
    `0`  REF (not counted toward this sample)
    `1`  ALT (counted toward this sample)
    `.`  missing
    `0:255,255:...` (with allele depth) — same convention
The script reads only the leading GT character before `:`, so it's
robust to either `1` or `1:255,0:1185:47,1134` formats.

Effect_class buckets
--------------------
Same as scripts 16/17:
    SYN     synonymous_variant
    NS      missense_variant, stop_gained, stop_lost, start_lost,
            initiator_codon_variant
    Other   anything else SnpEff returned (expected ~0 on the
            CDS-restricted canonical)

Sample-name normalization
-------------------------
Canonical sample columns are `{ID}_MT` (e.g., `33_MT`).
Stage 09 per_individual uses the same form (`{ID}_MT`).
Stage 14 per_individual uses BAM-path form
    (`MT_only_bams/{ID}_0_MT_only.bam`).
The script extracts the leading numeric ID from each and matches on
that, so the join is robust to whichever convention each table uses.

Run
---
    cd ~/Projects/MT_Genomics_Cl_Ap2026/MT_Genomics2
    conda activate SNP_env
    python scripts/18_variant_burden_per_individual.py
"""

from __future__ import annotations

import gzip
import re
import sys
from pathlib import Path

import pandas as pd

# ── PATHS ─────────────────────────────────────────────────────────────────────
PROJECT_ROOT = Path(__file__).resolve().parent.parent
VCF_DIR = PROJECT_ROOT / "vcf"

CANONICAL_VCF = VCF_DIR / "Fhet_MT_CDS.snps.split.vcf.gz"
HP_09_PER_IND = VCF_DIR / "heteroplasmy_per_individual.tsv"
HP_14_PER_IND = VCF_DIR / "heteroplasmy_pileup_per_individual_all.tsv"

OUT_14   = VCF_DIR / "per_individual_burden_pileup.tsv"
OUT_09   = VCF_DIR / "per_individual_burden_variant.tsv"
OUT_SUM  = VCF_DIR / "per_individual_burden_summary.txt"

# ── EFFECT CLASSIFICATION (same as scripts 16/17) ───────────────────────────
SYN_EFFECTS = {"synonymous_variant"}
NS_EFFECTS = {
    "missense_variant",
    "stop_gained",
    "stop_lost",
    "start_lost",
    "initiator_codon_variant",
}


def classify_effect(eff: str) -> str:
    if not eff:
        return "Other"
    parts = eff.split("&")
    if any(p in NS_EFFECTS for p in parts):
        return "NS"
    if any(p in SYN_EFFECTS for p in parts):
        return "SYN"
    return "Other"


def first_ann_effect(info: str) -> str:
    for entry in info.split(";"):
        if entry.startswith("ANN="):
            first = entry[4:].split(",", 1)[0]
            cols = first.split("|")
            return cols[1] if len(cols) > 1 else ""
    return ""


# ── ID EXTRACTION ────────────────────────────────────────────────────────────
ID_RE = re.compile(r"(\d+)")


def to_numeric_id(name: str) -> str | None:
    """Extract the leading numeric ID from a sample / BAM-path name.

    `33_MT`                                 -> '33'
    `MT_only_bams/33_0_MT_only.bam`         -> '33'
    `103`                                   -> '103'
    """
    m = ID_RE.search(name)
    return m.group(1) if m else None


# ── CANONICAL VCF PARSE: per-sample variant counts ──────────────────────────
def count_variants_per_sample(vcf_path: Path
                              ) -> tuple[list[str], dict[str, dict]]:
    """
    Parse the canonical VCF and return:
      samples: list of sample names in column order
      counts:  {sample_name: {"Total": n, "SYN": n, "NS": n, "Other": n}}
    """
    samples: list[str] = []
    counts: dict[str, dict] = {}
    opener = gzip.open if vcf_path.suffix == ".gz" else open

    with opener(vcf_path, "rt") as fh:
        for line in fh:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                header = line.rstrip("\n").split("\t")
                samples = header[9:]
                counts = {s: {"Total": 0, "SYN": 0, "NS": 0, "Other": 0}
                          for s in samples}
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) < 10:
                continue
            info     = fields[7]
            sample_v = fields[9:]
            if len(sample_v) != len(samples):
                continue  # shouldn't happen, but bail safely

            eff_class = classify_effect(first_ann_effect(info))

            for s, v in zip(samples, sample_v):
                gt = v.split(":", 1)[0]
                if gt == "1":
                    c = counts[s]
                    c["Total"] += 1
                    c[eff_class] += 1

    return samples, counts


# ── LOAD HP PER-INDIVIDUAL TABLES (zero-fill missing) ───────────────────────
HP_COLS = ["n_hp_sites", "n_hp_events", "n_ref_hp"]


def load_hp_per_individual(path: Path, sample_ids: list[str]
                           ) -> pd.DataFrame:
    """
    Returns a DataFrame indexed by numeric sample ID, with HP_COLS
    columns, zero-filled for any sample_id not present in `path`.
    """
    base = pd.DataFrame({"id": sample_ids})
    base["id_key"] = base["id"].astype(str)

    if not path.exists():
        for c in HP_COLS:
            base[c] = 0
        return base.drop(columns=["id_key"]).set_index("id")

    src = pd.read_csv(path, sep="\t")
    # Extract numeric ID from Individual column
    src["id_key"] = src["Individual"].astype(str).map(
        lambda n: to_numeric_id(n) or n
    )
    # Keep only HP_COLS + id_key
    src = src[["id_key"] + HP_COLS].copy()
    # Some tables may have multiple rows per individual — sum them
    src = src.groupby("id_key", as_index=False)[HP_COLS].sum()

    merged = base.merge(src, on="id_key", how="left")
    for c in HP_COLS:
        merged[c] = merged[c].fillna(0).astype(int)
    return merged.drop(columns=["id_key"]).set_index("id")


# ── BUILD MASTER TABLE ──────────────────────────────────────────────────────
def build_burden_table(samples: list[str], counts: dict[str, dict],
                       hp_path: Path) -> pd.DataFrame:
    rows = []
    for s in samples:
        rid = to_numeric_id(s) or s
        c = counts[s]
        rows.append({
            "Individual":       s,
            "id":               rid,
            "n_variants_total": c["Total"],
            "n_variants_SYN":   c["SYN"],
            "n_variants_NS":    c["NS"],
            "n_variants_Other": c["Other"],
        })
    df = pd.DataFrame(rows).set_index("id")

    sample_ids = df.index.tolist()
    hp = load_hp_per_individual(hp_path, sample_ids)
    out = df.join(hp, how="left")
    for c in HP_COLS:
        out[c] = out[c].fillna(0).astype(int)
    return out.reset_index().sort_values("n_variants_total", ascending=False)


def summarize_bimodality(df: pd.DataFrame, label: str) -> str:
    """Quick low/high split on n_variants_total and Hp burden by group."""
    n = len(df)
    low_cut, high_cut = 50, 190

    low  = df[df["n_variants_total"] < low_cut]
    high = df[df["n_variants_total"] > high_cut]
    mid  = df[(df["n_variants_total"] >= low_cut)
              & (df["n_variants_total"] <= high_cut)]

    def grp_stats(g: pd.DataFrame, name: str) -> str:
        if len(g) == 0:
            return f"  {name}: 0 samples"
        return (f"  {name}: n={len(g):>3}  "
                f"n_variants mean={g['n_variants_total'].mean():.1f} "
                f"median={g['n_variants_total'].median():.0f}  |  "
                f"n_hp_sites mean={g['n_hp_sites'].mean():.1f} "
                f"median={g['n_hp_sites'].median():.0f}  "
                f"max={g['n_hp_sites'].max()}  "
                f"n_with_hp={int((g['n_hp_sites'] > 0).sum())}")

    return "\n".join([
        f"== {label} ({n} samples) ==",
        f"  variant bimodality at <{low_cut} vs >{high_cut}:",
        grp_stats(low,  f"low  (<{low_cut} variants)"),
        grp_stats(mid,  f"mid  ({low_cut}-{high_cut} variants)"),
        grp_stats(high, f"high (>{high_cut} variants)"),
        "",
    ])


# ── MAIN ────────────────────────────────────────────────────────────────────
def main() -> int:
    if not CANONICAL_VCF.exists():
        print(f"ERROR: canonical VCF not found: {CANONICAL_VCF}", file=sys.stderr)
        return 2

    print(f"Parsing canonical VCF: {CANONICAL_VCF}")
    samples, counts = count_variants_per_sample(CANONICAL_VCF)
    print(f"  samples: {len(samples)}")
    total_records = sum(c["Total"] for c in counts.values())
    print(f"  sum across samples of GT=1 calls: {total_records:,}")
    print()

    lines = []
    for label, hp_path, out_path in [
        ("stage 14 (CDS pileup)",     HP_14_PER_IND, OUT_14),
        ("stage 09 (variant-only)",   HP_09_PER_IND, OUT_09),
    ]:
        if not hp_path.exists():
            print(f"  WARN: {hp_path.name} not found; "
                  f"writing variant counts only, Hp columns zero-filled.")
        df = build_burden_table(samples, counts, hp_path)
        df.to_csv(out_path, sep="\t", index=False)
        print(f"Wrote {out_path}  ({len(df)} rows)")
        lines.append(summarize_bimodality(df, label))

    OUT_SUM.write_text("\n".join(lines))
    print()
    print("\n".join(lines))
    print(f"Summary written to: {OUT_SUM}")
    return 0


if __name__ == "__main__":
    sys.exit(main())

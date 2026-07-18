#!/usr/bin/env python3
"""
16_annotate_hp.py — Annotate heteroplasmy event tables with SYN / NS / Other
                    coding effect, derived from the canonical VCF's SnpEff ANN.

Inputs  (in vcf/):
    Fhet_MT_CDS.snps.split.vcf.gz   canonical, SnpEff-annotated, one row per ALT
    heteroplasmy_events.tsv         stage 09 events (variant-only input)
    heteroplasmy_pileup_events_all.tsv
                                    stage 14 events (CDS pileup input)

Outputs (in vcf/):
    heteroplasmy_events_annot.tsv               stage 09 + Effect_class
    heteroplasmy_pileup_events_all_annot.tsv    stage 14 + Effect_class
    heteroplasmy_annot_summary.txt              Hp_class × Effect_class counts

Effect classification (from the first ANN entry per row, per SnpEff convention
that ANN[0] is the highest-impact gene-internal call):

    SYN          synonymous_variant (with any & variants thereof)
    NS           missense_variant, stop_gained, stop_lost, start_lost,
                 initiator_codon_variant (with any & variants)
                 — matches the impact set used by scripts 10 and 11.
    Other        anything else SnpEff returned (e.g., upstream_gene_variant,
                 intergenic_region). Should be 0 rows for events from
                 Fhet_MT_CDS.snps.split.vcf.gz (CDS-restricted), but kept
                 for safety.
    Unannotated  (POS, non-REF allele) not present in the canonical VCF.
                 Expected to be a small number for stage 14's
                 private_alt_Hp events that fell at non-canonical
                 (POS, ALT) combinations not produced by panel variant
                 calling. To compute SYN/NS for these, use the codon /
                 translation-table-2 machinery from scripts/10_dnds_per_gene.py.

Per-row rule for "which allele to annotate" — exactly one non-REF allele
exists per Hp event:

    Hp_is_REF == True   ->  non-REF allele = Major
    Hp_is_REF == False  ->  non-REF allele = Hp_allele

Run:
    cd ~/Projects/MT_Genomics_Cl_Ap2026/MT_Genomics2
    conda activate SNP_env
    python scripts/16_annotate_hp.py

Standalone, no external bcftools / cyvcf2 dependency — parses the gzipped
VCF text directly, since the canonical is small (~1k records).
"""

from __future__ import annotations

import gzip
import sys
from pathlib import Path

import pandas as pd

# --- paths ---------------------------------------------------------------

REPO = Path(__file__).resolve().parent.parent
VCF_DIR = REPO / "vcf"

CANONICAL_VCF = VCF_DIR / "Fhet_MT_CDS.snps.split.vcf.gz"
EVENTS_09     = VCF_DIR / "heteroplasmy_events.tsv"
EVENTS_14     = VCF_DIR / "heteroplasmy_pileup_events_all.tsv"

OUT_09        = VCF_DIR / "heteroplasmy_events_annot.tsv"
OUT_14        = VCF_DIR / "heteroplasmy_pileup_events_all_annot.tsv"
OUT_SUMMARY   = VCF_DIR / "heteroplasmy_annot_summary.txt"

# --- effect classification ----------------------------------------------

SYN_EFFECTS = {
    "synonymous_variant",
}
NS_EFFECTS = {
    "missense_variant",
    "stop_gained",
    "stop_lost",
    "start_lost",
    "initiator_codon_variant",
}


def classify_effect(eff: str) -> str:
    """SYN / NS / Other for a SnpEff effect string (may contain '&' joins)."""
    if not eff:
        return "Other"
    parts = eff.split("&")
    if any(p in NS_EFFECTS for p in parts):
        return "NS"
    if any(p in SYN_EFFECTS for p in parts):
        return "SYN"
    return "Other"


# --- canonical VCF parser -----------------------------------------------

def build_ann_map(vcf_path: Path) -> dict[tuple[int, str], dict]:
    """
    Parse the canonical (one-row-per-ALT, SnpEff-annotated) VCF and return
    a dict keyed by (POS, ALT) -> {effect, effect_class, gene, hgvs_p}.

    SnpEff ANN format (pipe-separated, comma-separated across annotations):
        Allele | Effect | Impact | GeneName | GeneID | FeatureType | FeatureID
        | TranscriptBiotype | Rank | HGVS.c | HGVS.p | cDNA.pos | CDS.pos
        | AA.pos | Distance | ERRORS

    We take ANN[0] (highest-impact gene-internal call) per SnpEff convention.
    """
    ann_map: dict[tuple[int, str], dict] = {}
    opener = gzip.open if vcf_path.suffix == ".gz" else open
    with opener(vcf_path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                continue
            try:
                pos = int(fields[1])
            except ValueError:
                continue
            alt = fields[4]
            info = fields[7]
            ann_field = None
            for entry in info.split(";"):
                if entry.startswith("ANN="):
                    ann_field = entry[4:]
                    break
            if ann_field is None:
                continue
            first = ann_field.split(",", 1)[0]
            cols = first.split("|")
            # Guard against malformed ANN strings.
            effect    = cols[1] if len(cols) > 1  else ""
            gene      = cols[3] if len(cols) > 3  else ""
            hgvs_p    = cols[10] if len(cols) > 10 else ""
            ann_map[(pos, alt)] = {
                "effect":       effect,
                "effect_class": classify_effect(effect),
                "gene":         gene,
                "hgvs_p":       hgvs_p,
            }
    return ann_map


# --- table annotation ---------------------------------------------------

REQUIRED_COLS = ["POS", "REF", "Major", "Hp_allele", "Hp_is_REF", "Hp_class"]


def annotate_events(df: pd.DataFrame, ann_map: dict[tuple[int, str], dict]
                    ) -> pd.DataFrame:
    missing = [c for c in REQUIRED_COLS if c not in df.columns]
    if missing:
        raise ValueError(f"Input table is missing required columns: {missing}")

    pos = df["POS"].astype(int).to_numpy()

    # Normalize Hp_is_REF: file may store it as True/False, "True"/"False",
    # or 1/0. Coerce to boolean.
    raw = df["Hp_is_REF"]
    if raw.dtype == bool:
        is_ref = raw.to_numpy()
    else:
        is_ref = raw.astype(str).str.strip().str.lower().isin(
            {"true", "t", "1", "yes"}
        ).to_numpy()

    major     = df["Major"].astype(str).to_numpy()
    hp_allele = df["Hp_allele"].astype(str).to_numpy()

    non_ref = [m if r else h for m, r, h in zip(major, is_ref, hp_allele)]

    eff_class = []
    effect    = []
    gene      = []
    hgvs_p    = []
    for p, a in zip(pos, non_ref):
        info = ann_map.get((int(p), a))
        if info is None:
            eff_class.append("Unannotated")
            effect.append("")
            gene.append("")
            hgvs_p.append("")
        else:
            eff_class.append(info["effect_class"])
            effect.append(info["effect"])
            gene.append(info["gene"])
            hgvs_p.append(info["hgvs_p"])

    out = df.copy()
    out["non_REF_allele"] = non_ref
    out["Effect_class"]   = eff_class
    out["Effect"]         = effect
    out["Gene"]           = gene
    out["HGVS_p"]         = hgvs_p
    return out


def summarize(label: str, df: pd.DataFrame) -> str:
    lines = [f"== {label} ({len(df)} events) ==", ""]
    ct = (
        df.groupby(["Hp_class", "Effect_class"], dropna=False)
          .size()
          .unstack(fill_value=0)
    )
    # Preferred column order
    for c in ["SYN", "NS", "Other", "Unannotated"]:
        if c not in ct.columns:
            ct[c] = 0
    ct = ct[["SYN", "NS", "Other", "Unannotated"]]
    ct["Total"] = ct.sum(axis=1)
    ct.loc["Total"] = ct.sum(axis=0)
    lines.append(ct.to_string())
    lines.append("")
    return "\n".join(lines)


# --- main ---------------------------------------------------------------

def main() -> int:
    if not CANONICAL_VCF.exists():
        print(f"ERROR: canonical VCF not found: {CANONICAL_VCF}", file=sys.stderr)
        return 2

    print(f"Parsing canonical VCF: {CANONICAL_VCF}")
    ann_map = build_ann_map(CANONICAL_VCF)
    syn = sum(1 for v in ann_map.values() if v["effect_class"] == "SYN")
    ns  = sum(1 for v in ann_map.values() if v["effect_class"] == "NS")
    oth = sum(1 for v in ann_map.values() if v["effect_class"] == "Other")
    print(f"  ANN entries: {len(ann_map):,} unique (POS, ALT)  "
          f"[SYN: {syn:,}  NS: {ns:,}  Other: {oth:,}]")

    summary_parts = []

    for label, in_path, out_path in [
        ("stage 09 (variant-only)", EVENTS_09, OUT_09),
        ("stage 14 (CDS pileup)",   EVENTS_14, OUT_14),
    ]:
        if not in_path.exists():
            print(f"WARN: {in_path} not found, skipping {label}")
            continue
        print(f"Annotating {label}: {in_path}")
        df = pd.read_csv(in_path, sep="\t")
        ann = annotate_events(df, ann_map)
        ann.to_csv(out_path, sep="\t", index=False)
        print(f"  wrote {out_path}  ({len(ann):,} rows)")
        summary_parts.append(summarize(label, ann))

    summary = "\n".join(summary_parts)
    OUT_SUMMARY.write_text(summary)
    print()
    print(summary)
    print(f"Summary also written to: {OUT_SUMMARY}")
    return 0


if __name__ == "__main__":
    sys.exit(main())

#!/usr/bin/env python3
"""
17_annotate_hp_codon.py — Fill in SYN / NS for heteroplasmy events that
                          script 16 marked `Unannotated`, by computing the
                          coding effect from first principles (reference
                          FASTA + GFF + vertebrate mt translation table 2).

Why this script exists
----------------------
Script 16 joins SnpEff `ANN` from the canonical VCF onto the Hp event
tables. Rows whose (POS, non-REF allele) combo is not in the canonical
get `Effect_class = Unannotated`. The canonical VCF only contains
panel-called variants, so by construction:

  - All 27 stage-14 `private_alt_Hp` events are `Unannotated` (their
    minor allele was never called as a panel variant).
  - Some stage-14 `REF_Hp` events are `Unannotated` when the
    pileup-major call is `NONE` (insufficient depth to make a major
    call) or is an ALT the panel caller never produced at that POS.
  - Stage-09 events are derived from `MT_DP_AD_141.txt`, which may
    have a small number of (POS, ALT) combos that don't survive into
    the canonical CDS-restricted VCF; those land here too.

For any Unannotated row whose POS sits inside a CDS, we *can* still
classify SYN vs NS — we just need the codon containing that position
and the substituted base. That's what this script does.

Method
------
For each CDS in the GFF (using `transl_table=2`, the vertebrate
mitochondrial code):

  1. Extract the coding-sense CDS sequence (reverse-complement for the
     − strand gene ND6); same logic as `scripts/10_dnds_per_gene.py`.
  2. For each Hp event with `Effect_class == "Unannotated"` and a POS
     falling inside this CDS:
        - Translate POS to a 1-based offset within the coding-sense
          CDS:
              + strand:  offset = POS − gene_start + 1
              − strand:  offset = gene_end − POS  + 1
        - Determine the codon index and the position within the codon
          (1, 2, or 3).
        - Build `new_base` = the non-REF allele (already chosen by 16
          as `Major` if `Hp_is_REF` else `Hp_allele`). For − strand,
          complement it.
        - Build the mutated codon, translate REF and mutated codons
          via NCBI table 2. SYN if the AA is unchanged, NS if it
          changed (including to/from stop).
        - Sanity-check: the REF codon's base at the substituted
          position should equal REF (for + strand) or complement(REF)
          (for −). Rows that fail this check are marked
          `Codon_REF_mismatch` so they don't silently corrupt the
          output.

Rows whose POS doesn't sit inside any CDS are left as `Other` (with
`Effect = "non_CDS"`). Rows with `NONE` or non-ACGT non-REF allele
are left as `Other_codon_invalid_allele`.

Inputs
------
    Missing_Files/SSM_MT_ref/Fhet_MT.fasta   (reference)
    Missing_Files/SSM_MT_ref/Fhet_MT.gff     (CDS coords + strand)
    vcf/heteroplasmy_events_annot.tsv        (stage 09, from script 16)
    vcf/heteroplasmy_pileup_events_all_annot.tsv
                                             (stage 14, from script 16)

Outputs
-------
    vcf/heteroplasmy_events_annot_codon.tsv               (stage 09 + codon-filled)
    vcf/heteroplasmy_pileup_events_all_annot_codon.tsv    (stage 14 + codon-filled)
    vcf/heteroplasmy_annot_codon_summary.txt              (before/after counts)

New columns added to the output tables on top of script 16's columns:
    Annotation_source   "snpeff" | "codon" | "snpeff_then_unannotated"
                        Tells you whether the row's SYN/NS came from
                        the canonical's ANN field (script 16) or from
                        codon computation (this script). If a row was
                        `Unannotated` after 16 AND its POS doesn't fall
                        in a CDS (rare), the field stays "snpeff" and
                        Effect_class stays "Other" — flagged with
                        Effect = "non_CDS".
    Codon_ref           REF codon at the (POS, codon-position) — only
                        populated for rows annotated by this script.
    Codon_alt           Mutated codon under the non-REF allele.
    AA_ref / AA_alt     Translated amino acids under MT_CODE table 2.

Run
---
    cd ~/Projects/MT_Genomics_Cl_Ap2026/MT_Genomics2
    conda activate SNP_env
    python scripts/17_annotate_hp_codon.py

Reuses MT_CODE + read_fasta + read_gff_cds + extract_cds_sequence from
the same source-of-truth in `scripts/10_dnds_per_gene.py`. No
dependency on cyvcf2 or bcftools — this script operates on the TSVs
produced by script 16.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

# ── PATHS ─────────────────────────────────────────────────────────────────────
PROJECT_ROOT = Path(__file__).resolve().parent.parent
VCF_DIR      = PROJECT_ROOT / "vcf"

REF_DIRS = [
    PROJECT_ROOT / "Missing_Files" / "SSM_MT_ref",
    PROJECT_ROOT / "refs",
]


def _find(name: str) -> Path:
    for d in REF_DIRS:
        p = d / name
        if p.exists():
            return p
    sys.exit(f"ERROR: {name} not found under any of {REF_DIRS}")


REF_FASTA = _find("Fhet_MT.fasta")
REF_GFF   = _find("Fhet_MT.gff")

IN_09  = VCF_DIR / "heteroplasmy_events_annot.tsv"
IN_14  = VCF_DIR / "heteroplasmy_pileup_events_all_annot.tsv"
OUT_09 = VCF_DIR / "heteroplasmy_events_annot_codon.tsv"
OUT_14 = VCF_DIR / "heteroplasmy_pileup_events_all_annot_codon.tsv"
OUT_SUM = VCF_DIR / "heteroplasmy_annot_codon_summary.txt"

# ── GENETIC CODE (NCBI table 2: Vertebrate Mitochondrial) ─────────────────────
# Identical to scripts/10_dnds_per_gene.py — single source of truth.
MT_CODE: dict[str, str] = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "TGT": "C", "TGC": "C", "TGA": "W", "TGG": "W",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "ATT": "I", "ATC": "I", "ATA": "M", "ATG": "M",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "AGT": "S", "AGC": "S", "AGA": "*", "AGG": "*",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

COMPLEMENT = str.maketrans("ACGTN", "TGCAN")
VALID_BASES = set("ACGT")


# ── REFERENCE / GFF READERS (lifted from script 10) ──────────────────────────
def read_fasta(path: Path) -> dict[str, str]:
    seqs: dict[str, str] = {}
    name = None
    buf: list[str] = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(buf).upper()
                name = line[1:].split()[0]
                buf = []
            else:
                buf.append(line)
    if name is not None:
        seqs[name] = "".join(buf).upper()
    return seqs


def read_gff_cds(path: Path) -> list[dict]:
    cdss: list[dict] = []
    with open(path) as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] != "CDS":
                continue
            attrs = dict(kv.split("=", 1) for kv in parts[8].split(";") if "=" in kv)
            gene = attrs.get("gene") or attrs.get("Name", "")
            cdss.append({
                "gene":   gene,
                "chrom":  parts[0],
                "start":  int(parts[3]),
                "end":    int(parts[4]),
                "strand": parts[6],
            })
    return cdss


def extract_cds_sequence(genome: dict[str, str], cds: dict) -> str:
    seq = genome[cds["chrom"]][cds["start"] - 1 : cds["end"]]
    if cds["strand"] == "-":
        seq = seq.translate(COMPLEMENT)[::-1]
    return seq


# ── CODON-LEVEL EFFECT FOR A SINGLE Hp ROW ───────────────────────────────────
def annotate_one(pos: int, ref_base: str, new_base: str,
                 cds_for_pos: dict | None,
                 cds_seqs: dict[str, str]) -> dict:
    """
    Return a dict with codon-level annotation for one Hp event:
        Effect_class    SYN / NS / Other / Other_codon_invalid_allele / Codon_REF_mismatch
        Effect          short string ("synonymous_codon", "missense_codon", ...)
        Gene            gene name (or "" if outside any CDS)
        HGVS_p          synthetic p.XnnnY where nn is 1-based AA position
                        (calculated; not from SnpEff)
        Codon_ref       AGT-style REF codon
        Codon_alt       AGT-style mutated codon
        AA_ref / AA_alt single-letter (or '*') under MT_CODE
    """
    blank = {
        "Effect_class": "Other",
        "Effect":       "non_CDS",
        "Gene":         "",
        "HGVS_p":       "",
        "Codon_ref":    "",
        "Codon_alt":    "",
        "AA_ref":       "",
        "AA_alt":       "",
    }
    if cds_for_pos is None:
        return blank

    new_base = new_base.upper()
    if new_base not in VALID_BASES:
        out = dict(blank)
        out["Effect_class"] = "Other"
        out["Effect"]       = "invalid_alt_allele"
        out["Gene"]         = cds_for_pos["gene"]
        return out

    gene  = cds_for_pos["gene"]
    s_pos = cds_for_pos["start"]
    e_pos = cds_for_pos["end"]
    strand = cds_for_pos["strand"]
    cds_seq = cds_seqs[gene]

    if strand == "+":
        cds_offset0 = pos - s_pos              # 0-based within coding-sense CDS
        coded_new   = new_base
        coded_ref   = ref_base.upper()
    else:  # − strand
        cds_offset0 = e_pos - pos              # 0-based within reversed CDS
        coded_new   = new_base.translate(COMPLEMENT)
        coded_ref   = ref_base.upper().translate(COMPLEMENT)

    codon_idx     = cds_offset0 // 3
    pos_in_codon  = cds_offset0 % 3            # 0, 1, or 2
    codon_start   = codon_idx * 3
    ref_codon     = cds_seq[codon_start : codon_start + 3]

    if len(ref_codon) != 3 or "N" in ref_codon:
        out = dict(blank)
        out["Effect_class"] = "Other"
        out["Effect"]       = "partial_or_N_codon"
        out["Gene"]         = gene
        return out

    # Sanity check: the codon's nucleotide at pos_in_codon should equal
    # the (possibly complemented) REF base. If not, the GFF / FASTA /
    # event-table positioning don't agree and the result would be garbage.
    if ref_codon[pos_in_codon] != coded_ref:
        out = dict(blank)
        out["Effect_class"] = "Codon_REF_mismatch"
        out["Effect"]       = (f"expected_{coded_ref}_at_codon_pos_"
                               f"{pos_in_codon + 1}_got_{ref_codon[pos_in_codon]}")
        out["Gene"]         = gene
        out["Codon_ref"]    = ref_codon
        return out

    alt_codon = ref_codon[:pos_in_codon] + coded_new + ref_codon[pos_in_codon + 1 :]
    aa_ref = MT_CODE.get(ref_codon, "?")
    aa_alt = MT_CODE.get(alt_codon, "?")

    if aa_ref == "?" or aa_alt == "?":
        eff_class = "Other"
        effect    = "unknown_codon"
    elif aa_ref == aa_alt:
        eff_class = "SYN"
        effect    = "synonymous_codon"
    elif aa_ref == "*" and aa_alt != "*":
        eff_class = "NS"
        effect    = "stop_lost_codon"
    elif aa_alt == "*" and aa_ref != "*":
        eff_class = "NS"
        effect    = "stop_gained_codon"
    else:
        eff_class = "NS"
        effect    = "missense_codon"

    aa_pos = codon_idx + 1
    hgvs_p = f"p.{aa_ref}{aa_pos}{aa_alt}"

    return {
        "Effect_class": eff_class,
        "Effect":       effect,
        "Gene":         gene,
        "HGVS_p":       hgvs_p,
        "Codon_ref":    ref_codon,
        "Codon_alt":    alt_codon,
        "AA_ref":       aa_ref,
        "AA_alt":       aa_alt,
    }


# ── DRIVER ───────────────────────────────────────────────────────────────────
def find_cds_for_pos(pos: int, cdss: list[dict]) -> dict | None:
    for cds in cdss:
        if cds["start"] <= pos <= cds["end"]:
            return cds
    return None


REQUIRED_COLS = ["POS", "REF", "Major", "Hp_allele", "Hp_is_REF",
                 "Hp_class", "non_REF_allele", "Effect_class",
                 "Effect", "Gene", "HGVS_p"]


def process_table(in_path: Path, out_path: Path,
                  cdss: list[dict], cds_seqs: dict[str, str]) -> dict:
    df = pd.read_csv(in_path, sep="\t")
    missing = [c for c in REQUIRED_COLS if c not in df.columns]
    if missing:
        raise ValueError(f"{in_path} missing columns from script 16 output: {missing}")

    n_total = len(df)
    n_unann_before = int((df["Effect_class"] == "Unannotated").sum())

    # New columns to be populated row-wise.
    new_codon_ref = [""] * n_total
    new_codon_alt = [""] * n_total
    new_aa_ref    = [""] * n_total
    new_aa_alt    = [""] * n_total
    new_source    = ["snpeff"] * n_total

    new_eff_class = df["Effect_class"].tolist()
    new_effect    = df["Effect"].fillna("").tolist()
    new_gene      = df["Gene"].fillna("").tolist()
    new_hgvs_p    = df["HGVS_p"].fillna("").tolist()

    counts = {
        "SYN": 0, "NS": 0, "Other": 0,
        "Codon_REF_mismatch": 0,
    }

    for i, row in df.iterrows():
        if row["Effect_class"] != "Unannotated":
            continue
        pos      = int(row["POS"])
        ref_base = str(row["REF"]).upper()
        new_base = str(row["non_REF_allele"]).upper()

        cds = find_cds_for_pos(pos, cdss)
        ann = annotate_one(pos, ref_base, new_base, cds, cds_seqs)

        new_eff_class[i] = ann["Effect_class"]
        new_effect[i]    = ann["Effect"]
        new_gene[i]      = ann["Gene"] or new_gene[i]
        new_hgvs_p[i]    = ann["HGVS_p"]
        new_codon_ref[i] = ann["Codon_ref"]
        new_codon_alt[i] = ann["Codon_alt"]
        new_aa_ref[i]    = ann["AA_ref"]
        new_aa_alt[i]    = ann["AA_alt"]
        new_source[i]    = "codon"
        counts[ann["Effect_class"]] = counts.get(ann["Effect_class"], 0) + 1

    out = df.copy()
    out["Effect_class"]      = new_eff_class
    out["Effect"]            = new_effect
    out["Gene"]              = new_gene
    out["HGVS_p"]            = new_hgvs_p
    out["Codon_ref"]         = new_codon_ref
    out["Codon_alt"]         = new_codon_alt
    out["AA_ref"]            = new_aa_ref
    out["AA_alt"]            = new_aa_alt
    out["Annotation_source"] = new_source
    out.to_csv(out_path, sep="\t", index=False)

    n_unann_after = int((out["Effect_class"] == "Unannotated").sum())
    return {
        "n_total":          n_total,
        "n_unann_before":   n_unann_before,
        "n_unann_after":    n_unann_after,
        "n_filled":         n_unann_before - n_unann_after,
        "codon_breakdown":  counts,
        "out_path":         out_path,
    }


def main() -> int:
    if not IN_09.exists() and not IN_14.exists():
        print(f"ERROR: neither {IN_09.name} nor {IN_14.name} found. "
              f"Run scripts/16_annotate_hp.py first.", file=sys.stderr)
        return 2

    print(f"Reference FASTA : {REF_FASTA}")
    print(f"Reference GFF   : {REF_GFF}")
    genome   = read_fasta(REF_FASTA)
    cdss     = read_gff_cds(REF_GFF)
    cds_seqs = {c["gene"]: extract_cds_sequence(genome, c) for c in cdss}
    print(f"CDS records read: {len(cdss)} "
          f"({sum(1 for c in cdss if c['strand'] == '-')} on − strand)")
    print()

    results = []
    for label, in_path, out_path in [
        ("stage 09 (variant-only)", IN_09, OUT_09),
        ("stage 14 (CDS pileup)",   IN_14, OUT_14),
    ]:
        if not in_path.exists():
            print(f"WARN: {in_path} not found, skipping {label}")
            continue
        print(f"Processing {label}: {in_path}")
        res = process_table(in_path, out_path, cdss, cds_seqs)
        results.append((label, res))
        print(f"  rows: {res['n_total']:,}  "
              f"Unannotated before: {res['n_unann_before']}  "
              f"after: {res['n_unann_after']}  "
              f"filled by codon: {res['n_filled']}")
        print(f"  codon-classified: {res['codon_breakdown']}")
        print(f"  wrote {res['out_path']}")
        print()

    # Summary file
    lines = ["heteroplasmy_annot_codon_summary.txt",
             "----------------------------------------",
             ""]
    for label, res in results:
        lines.append(f"== {label} ==")
        lines.append(f"  total rows:             {res['n_total']:,}")
        lines.append(f"  Unannotated (post-16):  {res['n_unann_before']}")
        lines.append(f"  Filled by codon (17):   {res['n_filled']}")
        lines.append(f"    SYN:                  {res['codon_breakdown'].get('SYN', 0)}")
        lines.append(f"    NS:                   {res['codon_breakdown'].get('NS', 0)}")
        lines.append(f"    Other (non-CDS / N):  {res['codon_breakdown'].get('Other', 0)}")
        lines.append(f"    REF mismatch:         {res['codon_breakdown'].get('Codon_REF_mismatch', 0)}")
        lines.append(f"  Still Unannotated:      {res['n_unann_after']}")
        lines.append("")
    OUT_SUM.write_text("\n".join(lines))
    print("\n".join(lines))
    print(f"Summary written to: {OUT_SUM}")
    return 0


if __name__ == "__main__":
    sys.exit(main())

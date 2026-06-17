#!/usr/bin/env python3
"""
scripts/10_dnds_per_gene.py

Per-gene + whole-mtDNA dN, dS, dN/dS using the simple Nei-Gojobori method
on the **vertebrate mitochondrial genetic code** (NCBI translation table 2).

DEFINITIONS
-----------
    dN  = (# observed nonsynonymous variants in gene) / (# possible nonsynonymous sites in gene)
    dS  = (# observed synonymous variants in gene)    / (# possible synonymous sites in gene)
    dN/dS  = ratio

"Sites" are counted by Nei-Gojobori: for each codon position, count what
fraction of the 3 possible single-nucleotide substitutions are synonymous
(s_i / 3) vs nonsynonymous ((3 - s_i) / 3) under the *vertebrate mt code*.
Stop-codon-creating substitutions are classified as nonsynonymous.
Per-gene S sites = sum over codons of (s_1 + s_2 + s_3) / 3.
Per-gene N sites = (3 * n_codons) − S sites.

POLYA-COMPLETED STOP CODONS
---------------------------
Five of the 13 *F. heteroclitus* mt protein-coding genes (ND2, COX2,
COX3, ND3, ND4) end one base short of a complete terminal stop codon
in the genomic sequence — the final TAA is created post-transcriptionally
by 3'-polyadenylation of the mRNA (the GFF records this as
`Note=TAA stop codon is completed by the addition of 3' A residues`).
Their CDS lengths are therefore (3*n + 1), not 3*n.

This script handles that correctly: `extract_cds_sequence()` returns the
genomic CDS exactly as annotated, and `count_sites_nei_gojobori()` then
truncates to the nearest codon multiple
(`L = len - (len % 3)`), dropping the trailing partial codon. The
incomplete final codon is also conceptually a stop and would be skipped
anyway. For complete-stop genes (the other 8), the final TAA / TAG /
AGA / AGG codon is recognized as `*` and skipped by the same logic.
Either way the terminal stop contributes nothing to S or N counts —
matching the standard treatment. Previous dN/dS attempts on this
genome failed because the partial last codon was passed unfiltered
into the translator and produced garbage codons; this script avoids
that by construction.

The Jukes-Cantor-corrected dN and dS are also reported:
    dN_JC = -3/4 * ln(1 - (4/3) * p_N),  where p_N = obs_N / N_sites
    dS_JC = -3/4 * ln(1 - (4/3) * p_S),  where p_S = obs_S / S_sites
NaN when the JC log argument goes ≤ 0 (very high divergence, doesn't
happen at the polymorphism level we see in this panel).

OBSERVED COUNTS
---------------
The numerator comes from the canonical SnpEff-annotated VCF, using the
**first** ANN entry per row (the highest-priority transcript):
    synonymous_variant                            → S
    splice_region_variant&synonymous_variant      → S
    stop_retained_variant                         → S  (encodes same Stop)
    splice_region_variant&stop_retained_variant   → S
    missense_variant                              → N
    stop_gained                                   → N
    stop_lost                                     → N
    start_lost                                    → N
Other ANN effects (upstream_gene_variant, etc.) are ignored — they
shouldn't appear in the CDS-restricted canonical, and if they do they
don't change a CDS amino acid.

INPUTS
------
    refs/Fhet_MT.gff                                 (gene + CDS coords, strand)
    refs/Fhet_MT.fasta                               (reference sequence)
    vcf/Fhet_MT_CDS.snps.split.vcf.gz                (canonical, SnpEff-annotated)

OUTPUTS
-------
    vcf/dnds_per_gene.tsv                            (one row per gene + Overall)
    stdout: short summary table

Run
---
    cd ~/Projects/MT_Genomics_Cl_Ap2026/MT_Genomics2
    python scripts/10_dnds_per_gene.py
"""

from __future__ import annotations

import math
from pathlib import Path
import sys

import cyvcf2

# ── PATHS ─────────────────────────────────────────────────────────────────────
PROJECT_ROOT = Path(__file__).resolve().parent.parent

# The Fhet_MT reference + GFF live under Missing_Files/SSM_MT_ref (the Mac
# canonical location). refs/ is the T2 location and may not exist locally.
REF_DIRS = [
    PROJECT_ROOT / "Missing_Files" / "SSM_MT_ref",
    PROJECT_ROOT / "refs",
]


def _find(filename: str) -> Path:
    for d in REF_DIRS:
        p = d / filename
        if p.exists():
            return p
    sys.exit(f"ERROR: {filename} not found under any of {REF_DIRS}")


REF_FASTA = _find("Fhet_MT.fasta")
REF_GFF   = _find("Fhet_MT.gff")
VCF_IN    = PROJECT_ROOT / "vcf" / "Fhet_MT_CDS.snps.split.vcf.gz"
OUT_TSV   = PROJECT_ROOT / "vcf" / "dnds_per_gene.tsv"


# ── GENETIC CODE (NCBI table 2: Vertebrate Mitochondrial) ─────────────────────
# Diff from standard code:
#   TGA = Trp (W) instead of Stop
#   ATA = Met (M) instead of Ile
#   AGA = Stop instead of Arg
#   AGG = Stop instead of Arg
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

# ANN effects → S or N
SYN_EFFECTS = {
    "synonymous_variant",
    "splice_region_variant&synonymous_variant",
    "stop_retained_variant",
    "splice_region_variant&stop_retained_variant",
}
NONSYN_EFFECTS = {
    "missense_variant",
    "splice_region_variant&missense_variant",
    "stop_gained",
    "splice_region_variant&stop_gained",
    "stop_lost",
    "splice_region_variant&stop_lost",
    "start_lost",
}


# ── FASTA + GFF READERS (minimal — no Biopython dependency) ───────────────────
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
    """Return list of CDS records with name, chrom, start (1-based incl), end, strand."""
    cdss: list[dict] = []
    with open(path) as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] != "CDS":
                continue
            attrs = dict(
                kv.split("=", 1) for kv in parts[8].split(";") if "=" in kv
            )
            gene = attrs.get("gene") or attrs.get("Name", "")
            cdss.append({
                "gene":   gene,
                "chrom":  parts[0],
                "start":  int(parts[3]),   # GFF is 1-based, inclusive
                "end":    int(parts[4]),
                "strand": parts[6],
            })
    return cdss


def extract_cds_sequence(genome: dict[str, str], cds: dict) -> str:
    """Return the *coding-sense* sequence for a CDS (reverse-complemented if strand=-)."""
    seq = genome[cds["chrom"]][cds["start"] - 1 : cds["end"]]
    if cds["strand"] == "-":
        seq = seq.translate(COMPLEMENT)[::-1]
    return seq


# ── NEI-GOJOBORI SITE COUNTING ────────────────────────────────────────────────
def count_sites_nei_gojobori(cds_seq: str) -> tuple[float, float, int]:
    """
    Count synonymous and nonsynonymous sites in a CDS using NCBI table 2.

    For each codon position i (1, 2, 3):
        s_i = number of the 3 possible single-nucleotide substitutions
              at that position that leave the encoded amino acid unchanged
              (under the vertebrate mt code).
    Per-codon S sites = (s_1 + s_2 + s_3) / 3
    Per-codon N sites = 3 − S sites
    Sum across all codons (skipping the stop codon if present at the end,
    and skipping any codon containing N).

    Returns (S_sites, N_sites, n_codons_scored).
    """
    bases = "ACGT"
    S_total = 0.0
    n_codons = 0

    L = len(cds_seq) - (len(cds_seq) % 3)  # truncate to codon multiple
    for i in range(0, L, 3):
        codon = cds_seq[i : i + 3]
        if "N" in codon or codon not in MT_CODE:
            continue
        ref_aa = MT_CODE[codon]
        # Skip the terminal stop codon (not subject to syn/nonsyn pressure
        # in the normal sense — substitutions there are mostly nonsense).
        # Internal stops (would-be premature stops in REF) are unusual and
        # also skipped.
        if ref_aa == "*":
            continue
        # For each position, count synonymous alternatives
        codon_S = 0.0
        for pos in range(3):
            n_syn = 0
            for b in bases:
                if b == codon[pos]:
                    continue
                new_codon = codon[:pos] + b + codon[pos + 1 :]
                new_aa = MT_CODE.get(new_codon, "?")
                # Stop-creating substitutions count as nonsynonymous.
                if new_aa == ref_aa:
                    n_syn += 1
            codon_S += n_syn / 3.0
        S_total += codon_S
        n_codons += 1

    N_total = 3 * n_codons - S_total
    return S_total, N_total, n_codons


# ── ANN PARSING ───────────────────────────────────────────────────────────────
def parse_first_ann(variant: cyvcf2.Variant) -> tuple[str, str]:
    """Return (gene, effect) from the first ANN entry; ('','') if no ANN."""
    raw = variant.INFO.get("ANN", "")
    if not raw:
        return "", ""
    fields = raw.split(",")[0].split("|")
    effect = fields[1] if len(fields) > 1 else ""
    gene   = fields[3] if len(fields) > 3 else ""
    return gene, effect


# ── JC CORRECTION ─────────────────────────────────────────────────────────────
def jukes_cantor(p: float) -> float:
    """Jukes-Cantor distance from raw proportion. NaN if math undefined."""
    arg = 1.0 - (4.0 / 3.0) * p
    if arg <= 0:
        return float("nan")
    return -0.75 * math.log(arg)


# ── MAIN ──────────────────────────────────────────────────────────────────────
def main() -> None:
    if not VCF_IN.exists():
        sys.exit(f"ERROR: VCF not found: {VCF_IN}\n"
                 "  Run scripts/07_cds_snps_norm_mac.sh first.")

    print(f"Reference FASTA : {REF_FASTA}")
    print(f"Reference GFF   : {REF_GFF}")
    print(f"VCF             : {VCF_IN}")

    genome = read_fasta(REF_FASTA)
    cdss = read_gff_cds(REF_GFF)
    print(f"CDS records in GFF: {len(cdss)}")

    # ── Denominator: S and N sites per gene ──────────────────────────────────
    gene_sites: dict[str, dict] = {}
    for cds in cdss:
        seq = extract_cds_sequence(genome, cds)
        S, N, n_codons = count_sites_nei_gojobori(seq)
        gene_sites[cds["gene"]] = {
            "strand":   cds["strand"],
            "length":   cds["end"] - cds["start"] + 1,
            "n_codons": n_codons,
            "S_sites":  S,
            "N_sites":  N,
        }

    # ── Numerator: observed S and N variants per gene from canonical VCF ────
    gene_obs: dict[str, dict] = {g: {"obs_S": 0, "obs_N": 0, "other": 0}
                                  for g in gene_sites}
    skipped_no_ann = 0
    skipped_no_gene = 0

    for variant in cyvcf2.VCF(str(VCF_IN)):
        gene, effect = parse_first_ann(variant)
        if not effect:
            skipped_no_ann += 1
            continue
        if gene not in gene_obs:
            skipped_no_gene += 1
            continue
        if effect in SYN_EFFECTS:
            gene_obs[gene]["obs_S"] += 1
        elif effect in NONSYN_EFFECTS:
            gene_obs[gene]["obs_N"] += 1
        else:
            gene_obs[gene]["other"] += 1

    # ── Compute dN, dS, ratio, JC-corrected ──────────────────────────────────
    rows: list[dict] = []
    for gene, sites in gene_sites.items():
        obs = gene_obs[gene]
        S_sites = sites["S_sites"]
        N_sites = sites["N_sites"]
        obs_S = obs["obs_S"]
        obs_N = obs["obs_N"]

        p_S = obs_S / S_sites if S_sites > 0 else float("nan")
        p_N = obs_N / N_sites if N_sites > 0 else float("nan")
        ratio = (p_N / p_S) if p_S > 0 else float("nan")

        dN_JC = jukes_cantor(p_N) if not math.isnan(p_N) else float("nan")
        dS_JC = jukes_cantor(p_S) if not math.isnan(p_S) else float("nan")
        ratio_JC = (dN_JC / dS_JC) if dS_JC and dS_JC > 0 else float("nan")

        rows.append({
            "Gene":       gene,
            "Strand":     sites["strand"],
            "CDS_len":    sites["length"],
            "n_codons":   sites["n_codons"],
            "S_sites":    f"{S_sites:.2f}",
            "N_sites":    f"{N_sites:.2f}",
            "obs_S":      obs_S,
            "obs_N":      obs_N,
            "other_eff":  obs["other"],
            "pN":         f"{p_N:.5f}" if not math.isnan(p_N) else "NA",
            "pS":         f"{p_S:.5f}" if not math.isnan(p_S) else "NA",
            "pN_over_pS": f"{ratio:.4f}" if not math.isnan(ratio) else "NA",
            "dN_JC":      f"{dN_JC:.5f}" if not math.isnan(dN_JC) else "NA",
            "dS_JC":      f"{dS_JC:.5f}" if not math.isnan(dS_JC) else "NA",
            "dN_dS_JC":   f"{ratio_JC:.4f}" if not math.isnan(ratio_JC) else "NA",
        })

    # ── Overall (pooled across all 13 genes) ─────────────────────────────────
    total_S_sites = sum(g["S_sites"] for g in gene_sites.values())
    total_N_sites = sum(g["N_sites"] for g in gene_sites.values())
    total_obs_S = sum(g["obs_S"] for g in gene_obs.values())
    total_obs_N = sum(g["obs_N"] for g in gene_obs.values())
    total_other = sum(g["other"] for g in gene_obs.values())
    total_codons = sum(g["n_codons"] for g in gene_sites.values())
    total_len = sum(g["length"] for g in gene_sites.values())

    p_S_all = total_obs_S / total_S_sites if total_S_sites > 0 else float("nan")
    p_N_all = total_obs_N / total_N_sites if total_N_sites > 0 else float("nan")
    ratio_all = (p_N_all / p_S_all) if p_S_all > 0 else float("nan")
    dN_JC_all = jukes_cantor(p_N_all) if not math.isnan(p_N_all) else float("nan")
    dS_JC_all = jukes_cantor(p_S_all) if not math.isnan(p_S_all) else float("nan")
    ratio_JC_all = (dN_JC_all / dS_JC_all) if dS_JC_all and dS_JC_all > 0 else float("nan")

    rows.append({
        "Gene":       "Overall",
        "Strand":     "",
        "CDS_len":    total_len,
        "n_codons":   total_codons,
        "S_sites":    f"{total_S_sites:.2f}",
        "N_sites":    f"{total_N_sites:.2f}",
        "obs_S":      total_obs_S,
        "obs_N":      total_obs_N,
        "other_eff":  total_other,
        "pN":         f"{p_N_all:.5f}" if not math.isnan(p_N_all) else "NA",
        "pS":         f"{p_S_all:.5f}" if not math.isnan(p_S_all) else "NA",
        "pN_over_pS": f"{ratio_all:.4f}" if not math.isnan(ratio_all) else "NA",
        "dN_JC":      f"{dN_JC_all:.5f}" if not math.isnan(dN_JC_all) else "NA",
        "dS_JC":      f"{dS_JC_all:.5f}" if not math.isnan(dS_JC_all) else "NA",
        "dN_dS_JC":   f"{ratio_JC_all:.4f}" if not math.isnan(ratio_JC_all) else "NA",
    })

    # ── Write TSV ────────────────────────────────────────────────────────────
    header = list(rows[0].keys())
    with open(OUT_TSV, "w") as out:
        out.write("\t".join(header) + "\n")
        for r in rows:
            out.write("\t".join(str(r[k]) for k in header) + "\n")

    # ── Stdout summary ───────────────────────────────────────────────────────
    print(f"\nVariants with no ANN     : {skipped_no_ann}")
    print(f"Variants outside named genes: {skipped_no_gene}")
    print(f"\nWrote {OUT_TSV}\n")

    # Pretty print
    col_widths = {k: max(len(k), max(len(str(r[k])) for r in rows)) for k in header}
    print(" | ".join(f"{k:<{col_widths[k]}}" for k in header))
    print("-+-".join("-" * col_widths[k] for k in header))
    for r in rows:
        print(" | ".join(f"{str(r[k]):<{col_widths[k]}}" for k in header))


if __name__ == "__main__":
    main()

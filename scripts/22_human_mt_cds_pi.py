#!/usr/bin/env python3
"""
22_human_mt_cds_pi.py
Calculate CDS-specific nucleotide diversity (π_syn, π_ns) for the
Lankheet et al. 2026 human mt dataset, using:
  - per-site π values from 21_human_mt_pi.py  (human_mt_pi_per_site.tsv)
  - rCRS (NC_012920) sequence and hardcoded CDS coordinates
  - Vertebrate mitochondrial genetic code (NCBI translation table 2)

Approach:
  For each variable CDS site, substitute the minor allele into the rCRS
  codon context → translate → classify SYN / NS.
  π_syn and π_ns are summed over classified sites, normalised by L_CDS.

Usage (from Human_mt/):
    conda activate SNP_env
    python 22_human_mt_cds_pi.py

Output:
    human_mt_cds_pi_summary.txt
    human_mt_cds_pi_per_site.tsv
"""

import sys
import csv
from Bio import SeqIO

# ── Input files ───────────────────────────────────────────────────────────────
SITES_TSV  = "human_mt_pi_per_site.tsv"
RCRS_FILE  = "rCRS_NC_012920.fasta"
OUT_SUM    = "human_mt_cds_pi_summary.txt"
OUT_SITES  = "human_mt_cds_pi_per_site.tsv"

# ── rCRS CDS coordinates (1-based, inclusive) ─────────────────────────────────
# Source: NC_012920.1 GenBank annotation
# strand: +1 or -1 (ND6 is the only minus-strand gene)
CDS_GENES = [
    ("MT-ND1",  3307,  4262,  1),
    ("MT-ND2",  4470,  5511,  1),
    ("MT-CO1",  5904,  7445,  1),
    ("MT-CO2",  7586,  8269,  1),
    ("MT-ATP8", 8366,  8572,  1),
    ("MT-ATP6", 8527,  9207,  1),
    ("MT-CO3",  9207,  9990,  1),
    ("MT-ND3", 10059, 10404,  1),
    ("MT-ND4L",10470, 10766,  1),
    ("MT-ND4", 10760, 12137,  1),
    ("MT-ND5", 12337, 14148,  1),
    ("MT-ND6", 14149, 14673, -1),   # minus strand
    ("MT-CYB", 14747, 15887,  1),
]

L_CDS = sum(end - start + 1 for _, start, end, _ in CDS_GENES)   # 11,395 bp

# ── Vertebrate mitochondrial genetic code (NCBI table 2) ─────────────────────
MT_CODE = {
    'TTT':'Phe','TTC':'Phe','TTA':'Leu','TTG':'Leu',
    'CTT':'Leu','CTC':'Leu','CTA':'Leu','CTG':'Leu',
    'ATT':'Ile','ATC':'Ile','ATA':'Met','ATG':'Met',   # ATA=Met (mt)
    'GTT':'Val','GTC':'Val','GTA':'Val','GTG':'Val',
    'TCT':'Ser','TCC':'Ser','TCA':'Ser','TCG':'Ser',
    'CCT':'Pro','CCC':'Pro','CCA':'Pro','CCG':'Pro',
    'ACT':'Thr','ACC':'Thr','ACA':'Thr','ACG':'Thr',
    'GCT':'Ala','GCC':'Ala','GCA':'Ala','GCG':'Ala',
    'TAT':'Tyr','TAC':'Tyr','TAA':'Stop','TAG':'Stop',
    'CAT':'His','CAC':'His','CAA':'Gln','CAG':'Gln',
    'AAT':'Asn','AAC':'Asn','AAA':'Lys','AAG':'Lys',
    'GAT':'Asp','GAC':'Asp','GAA':'Glu','GAG':'Glu',
    'TGT':'Cys','TGC':'Cys','TGA':'Trp','TGG':'Trp', # TGA=Trp (mt)
    'CGT':'Arg','CGC':'Arg','CGA':'Arg','CGG':'Arg',
    'AGT':'Ser','AGC':'Ser','AGA':'Stop','AGG':'Stop',# AGA/AGG=Stop (mt)
    'GGT':'Gly','GGC':'Gly','GGA':'Gly','GGG':'Gly',
}

RC = str.maketrans('ACGT','TGCA')

def revcomp(seq):
    return seq.translate(RC)[::-1]

def translate_mt(codon):
    return MT_CODE.get(codon.upper(), '?')


def pos_to_cds(pos_1based):
    """
    Return (gene, strand, cds_pos_0based) if pos is in a CDS, else None.
    cds_pos_0based is the position within the gene's coding sequence (0-based),
    accounting for strand.
    """
    for gene, start, end, strand in CDS_GENES:
        if start <= pos_1based <= end:
            if strand == 1:
                cds_pos = pos_1based - start      # 0-based within CDS
            else:
                cds_pos = end - pos_1based         # 0-based, 5'→3' on minus strand
            return gene, strand, cds_pos
    return None


def classify_substitution(pos_1based, ref_allele, alt_allele, ref_seq):
    """
    Return 'SYN', 'NS', or 'non_CDS'.
    Uses rCRS codon context; substitutes alt_allele at the variable position.
    For minus-strand genes, reverse-complements the codon before translation.
    """
    info = pos_to_cds(pos_1based)
    if info is None:
        return 'non_CDS'

    gene, strand, cds_pos = info
    codon_idx = cds_pos // 3        # which codon (0-based)
    codon_off = cds_pos % 3         # position within codon (0-based)

    # Get CDS start/end for this gene
    for g, start, end, s in CDS_GENES:
        if g == gene:
            gene_start, gene_end = start, end
            break

    if strand == 1:
        codon_start_genome = gene_start + codon_idx * 3   # 1-based
        codon_seq = ref_seq[codon_start_genome - 1 : codon_start_genome + 2].upper()
        if len(codon_seq) < 3:
            return 'non_CDS'
        # Substitute alt allele at codon position
        codon_alt = codon_seq[:codon_off] + alt_allele.upper() + codon_seq[codon_off+1:]
    else:
        # Minus strand: codon runs from genome right to left
        # cds_pos = end - pos (0-based), codon_idx, codon_off as above
        codon_start_genome = gene_end - codon_idx * 3      # rightmost genome coord of codon
        codon_seq_fwd = ref_seq[codon_start_genome - 3 : codon_start_genome].upper()
        codon_seq = revcomp(codon_seq_fwd)                 # 5'→3' on minus strand
        if len(codon_seq) < 3:
            return 'non_CDS'
        alt_on_minus = revcomp(alt_allele.upper())         # complement of genome ALT
        codon_alt = codon_seq[:codon_off] + alt_on_minus + codon_seq[codon_off+1:]

    aa_ref = translate_mt(codon_seq)
    aa_alt = translate_mt(codon_alt)

    if aa_ref == '?' or aa_alt == '?':
        return 'non_CDS'
    return 'SYN' if aa_ref == aa_alt else 'NS'


def main():
    # ── Load rCRS ─────────────────────────────────────────────────────────────
    ref_record = next(SeqIO.parse(RCRS_FILE, "fasta"))
    ref_seq    = str(ref_record.seq)
    print(f"rCRS: {ref_record.id}  ({len(ref_seq):,} bp)", file=sys.stderr)
    print(f"L_CDS (13 genes): {L_CDS:,} bp", file=sys.stderr)

    # ── Load per-site π table ─────────────────────────────────────────────────
    with open(SITES_TSV) as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        sites  = list(reader)
    print(f"Variable sites (whole genome): {len(sites)}", file=sys.stderr)

    # ── Classify and accumulate ───────────────────────────────────────────────
    pi_total = pi_syn = pi_ns = 0.0
    n_cds = n_syn = n_ns = n_noncds = 0
    rows_out = []

    for s in sites:
        pos       = int(s['POS'])
        ref_base  = s['ref_base']
        minor     = s['minor']
        pi_val    = float(s['pi_site'])

        if minor == '.' or minor not in 'ACGT':
            continue

        info = pos_to_cds(pos)
        if info is None:
            n_noncds += 1
            continue

        gene, strand, cds_pos = info
        effect = classify_substitution(pos, ref_base, minor, ref_seq)

        pi_total += pi_val
        n_cds    += 1

        if effect == 'SYN':
            pi_syn += pi_val
            n_syn  += 1
        elif effect == 'NS':
            pi_ns += pi_val
            n_ns  += 1

        rows_out.append({
            'POS':       pos,
            'Gene':      gene,
            'Strand':    strand,
            'ref_base':  ref_base,
            'minor':     minor,
            'n_called':  s['n_called'],
            'freq_minor':s['freq_minor'],
            'Effect':    effect,
            'pi_site':   s['pi_site'],
        })

    # ── Normalise ─────────────────────────────────────────────────────────────
    pi_t = pi_total / L_CDS
    pi_s = pi_syn   / L_CDS
    pi_n = pi_ns    / L_CDS
    pnps = pi_n / pi_s if pi_s > 0 else float('nan')

    # ── Watterson's θ ─────────────────────────────────────────────────────────
    N_seq = 1_176   # Lankheet 2026 sequences
    a1    = sum(1/i for i in range(1, N_seq))
    th_t  = n_cds / (a1 * L_CDS)
    th_s  = n_syn / (a1 * L_CDS)
    th_n  = n_ns  / (a1 * L_CDS)
    pnps_th = th_n / th_s if th_s > 0 else float('nan')

    # ── Write per-site CDS table ──────────────────────────────────────────────
    fields = ['POS','Gene','Strand','ref_base','minor','n_called','freq_minor','Effect','pi_site']
    with open(OUT_SITES, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter='\t')
        w.writeheader()
        w.writerows(rows_out)

    # ── Summary ───────────────────────────────────────────────────────────────
    summary = (
        f"=== π and Watterson's θ — Lankheet et al. 2026 Human mt CDS ===\n"
        f"  Reference       : rCRS NC_012920  ({len(ref_seq):,} bp)\n"
        f"  Sequences (N)   : {N_seq:,}\n"
        f"  a₁              : {a1:.4f}\n"
        f"  L_CDS (13 genes): {L_CDS:,} bp\n"
        f"  CDS variable sites (S): {n_cds}  (syn={n_syn}, NS={n_ns})\n"
        f"  Non-CDS sites skipped : {n_noncds}\n"
        f"\n"
        f"  θ_W_total       : {th_t:.5f}\n"
        f"  θ_W_syn         : {th_s:.5f}\n"
        f"  θ_W_ns          : {th_n:.5f}\n"
        f"  pN/pS (θ_W)     : {pnps_th:.3f}\n"
        f"\n"
        f"  π_total         : {pi_t:.5f}\n"
        f"  π_syn           : {pi_s:.5f}\n"
        f"  π_ns            : {pi_n:.5f}\n"
        f"  pN/pS (π)       : {pnps:.3f}\n"
        f"\nPer-site table   : {OUT_SITES}\n"
    )
    print('\n' + summary)
    with open(OUT_SUM, 'w') as f:
        f.write(summary)


if __name__ == "__main__":
    main()

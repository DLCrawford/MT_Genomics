#!/usr/bin/env python3
"""
21_human_mt_pi.py
Calculate whole-genome nucleotide diversity (π) for 1,176 complete human
mitochondrial genomes (Lankheet et al. 2026, Commun Biol).

Approach: reference-based (no MSA tool required).
  1. Download rCRS (NC_012920) via NCBI efetch.
  2. Align each sequence to rCRS using Biopython PairwiseAligner.
  3. Record the called base at each rCRS position across all 1,176 sequences.
  4. Compute π = Σ_sites (n/n-1) * 2*p*(1-p) / L

Usage (from Human_mt/):
    conda activate SNP_env
    python 21_human_mt_pi.py

Output:
    human_mt_pi_summary.txt   — π summary
    human_mt_pi_per_site.tsv  — per-site allele frequencies and π contribution
"""

import sys
import os
import urllib.request
import urllib.parse
import time
from collections import Counter

from Bio import SeqIO
from Bio.Align import PairwiseAligner

# ── Config ────────────────────────────────────────────────────────────────────
IN_FASTA   = "lankheet_2026_mt.fasta"
RCRS_FILE  = "rCRS_NC_012920.fasta"
OUT_SUM    = "human_mt_pi_summary.txt"
OUT_SITES  = "human_mt_pi_per_site.tsv"
API_KEY    = "bdc92badbf2cdedaa04b11712a3af3b8a109"
RCRS_ACC   = "NC_012920"
# ─────────────────────────────────────────────────────────────────────────────


def download_rcrs(api_key, out_file):
    """Download rCRS (NC_012920) from NCBI if not already present."""
    if os.path.exists(out_file):
        print(f"rCRS already present: {out_file}", file=sys.stderr)
        return
    print("Downloading rCRS (NC_012920) ...", file=sys.stderr)
    data = urllib.parse.urlencode({
        "db": "nuccore", "id": RCRS_ACC,
        "rettype": "fasta", "retmode": "text",
        "api_key": api_key,
    }).encode("utf-8")
    req = urllib.request.Request(
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi",
        data=data, method="POST"
    )
    with urllib.request.urlopen(req, timeout=60) as resp:
        text = resp.read().decode("utf-8")
    with open(out_file, "w") as f:
        f.write(text)
    print(f"  Saved: {out_file}", file=sys.stderr)


def make_aligner():
    """Biopython PairwiseAligner configured for near-identical mt sequences."""
    aligner = PairwiseAligner()
    aligner.mode            = "global"
    aligner.match_score     =  2
    aligner.mismatch_score  = -1
    aligner.open_gap_score  = -5
    aligner.extend_gap_score = -0.5
    return aligner


def align_to_ref(query_seq, ref_seq, aligner):
    """
    Return a list of length len(ref_seq) with the query base at each
    reference position (or '-' for deletion, None for insertion gaps).
    """
    alignment = aligner.align(ref_seq, query_seq)[0]
    ref_aln, qry_aln = alignment[0], alignment[1]

    ref_bases = []
    ref_pos   = 0
    for r, q in zip(ref_aln, qry_aln):
        if r != '-':                  # reference base (not insertion in ref)
            ref_bases.append((ref_pos, q))
            ref_pos += 1

    # Build position array indexed to reference
    result = [None] * len(ref_seq)
    for rpos, qbase in ref_bases:
        result[rpos] = qbase
    return result


def pi_site(counts, n):
    """π at one site given a Counter of bases and total called n."""
    if n < 2:
        return 0.0
    pi = 0.0
    for base, cnt in counts.items():
        if base in "ACGT":
            p = cnt / n
            pi += p * (1 - p)
    return pi * (n / (n - 1)) * 2   # multiply by 2 accounts for all pairs


def main():
    # ── Download rCRS ─────────────────────────────────────────────────────────
    download_rcrs(API_KEY, RCRS_FILE)
    ref_record = next(SeqIO.parse(RCRS_FILE, "fasta"))
    ref_seq    = str(ref_record.seq).upper()
    L          = len(ref_seq)
    print(f"rCRS length: {L:,} bp", file=sys.stderr)

    # ── Read all query sequences ──────────────────────────────────────────────
    records = list(SeqIO.parse(IN_FASTA, "fasta"))
    N       = len(records)
    print(f"Query sequences: {N}", file=sys.stderr)

    aligner = make_aligner()

    # site_counts[pos] = Counter of bases across all sequences
    site_counts = [Counter() for _ in range(L)]

    for i, rec in enumerate(records):
        if i % 100 == 0:
            print(f"  Aligning {i+1}/{N} ...", file=sys.stderr, flush=True)
        query = str(rec.seq).upper()
        try:
            pos_bases = align_to_ref(query, ref_seq, aligner)
            for pos, base in enumerate(pos_bases):
                if base and base in "ACGT":
                    site_counts[pos][base] += 1
        except Exception as e:
            print(f"  WARNING: failed for {rec.id}: {e}", file=sys.stderr)

    # ── Calculate π ──────────────────────────────────────────────────────────
    print("\nCalculating π ...", file=sys.stderr)
    pi_total = 0.0
    n_variable = 0

    site_rows = []
    for pos, counts in enumerate(site_counts):
        n_called = sum(counts.values())
        if n_called < 2:
            continue
        n_alleles = len([b for b in counts if b in "ACGT"])
        p_site = pi_site(counts, n_called)
        pi_total += p_site

        if n_alleles > 1:
            n_variable += 1
            major = counts.most_common(1)[0][0]
            minor_counts = {b: c for b, c in counts.items() if b != major and b in "ACGT"}
            minor_base = max(minor_counts, key=minor_counts.get) if minor_counts else "."
            minor_freq  = minor_counts.get(minor_base, 0) / n_called if minor_counts else 0
            site_rows.append({
                "POS":        pos + 1,
                "ref_base":   ref_seq[pos],
                "major":      major,
                "minor":      minor_base,
                "n_called":   n_called,
                "n_minor":    minor_counts.get(minor_base, 0),
                "freq_minor": round(minor_freq, 6),
                "pi_site":    round(p_site, 8),
            })

    pi_per_site = pi_total / L

    # ── Write per-site table ──────────────────────────────────────────────────
    import csv
    with open(OUT_SITES, "w", newline="") as fh:
        fields = ["POS","ref_base","major","minor","n_called","n_minor","freq_minor","pi_site"]
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        w.writeheader()
        w.writerows(site_rows)

    # ── Summary ───────────────────────────────────────────────────────────────
    summary = (
        f"=== Nucleotide Diversity (π) — Lankheet et al. 2026 Human mt ===\n"
        f"  Reference      : rCRS NC_012920 ({L:,} bp)\n"
        f"  Sequences (N)  : {N}\n"
        f"  Variable sites : {n_variable}\n"
        f"  π (whole mt)   : {pi_per_site:.5f}\n"
        f"  π (whole mt)   : {pi_per_site:.2e}\n"
        f"\nPer-site table  : {OUT_SITES} ({len(site_rows)} variable sites)\n"
    )
    print("\n" + summary)
    with open(OUT_SUM, "w") as f:
        f.write(summary)
    print(f"Summary written : {OUT_SUM}", file=sys.stderr)


if __name__ == "__main__":
    main()

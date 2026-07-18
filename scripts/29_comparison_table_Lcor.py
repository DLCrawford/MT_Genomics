#!/usr/bin/env python3
"""
29_comparison_table_Lcor.py  —  SITE-LENGTH-CORRECTED cross-species comparison.

Background
----------
`25_comparison_table.py` normalises BOTH the synonymous and the nonsynonymous
diversity by the *total* CDS length (L_CDS).  Because the same L appears in
numerator and denominator, the reported pN/pS columns collapse to the raw
*count ratio*  S_ns / S_syn  (for theta) and the frequency ratio
piN_raw / piS_raw  (for pi).  That is NOT the conventional, dN/dS-comparable
pN/pS, and it is not strictly comparable across species that use different
genetic codes / base compositions (vertebrate-mt table 2 vs invertebrate-mt
table 5).

This script recomputes the per-class rates with the proper denominators:

    theta_S = S_syn / (a1 * L_S)        theta_N = S_ns / (a1 * L_N)
    pi_S    = piS_raw / L_S             pi_N    = piN_raw / L_N
    pN/pS   = (S_ns/L_N)/(S_syn/L_S) = (S_ns/S_syn) * (L_S/L_N)        [theta]
              (piN_raw/L_N)/(piS_raw/L_S) = (piN/piS_raw)*(L_S/L_N)    [pi]

L_S and L_N are Nei-Gojobori site counts:
  * Fhet                 — exact, from Fhet_MT.fasta + Fhet_MT.gff  (NCBI table 2)
  * Human Afr & AMR      — exact, from rCRS + hardcoded CDS coords  (NCBI table 2)
  * Drosophila, C.elegans— code-level estimate under invertebrate-mt (table 5);
                           their reference CDS are not in the repo.  Flagged
                           "code-est" in the L_source column.  Validated:
                           the table-2 *code-level* synonymous fraction is
                           printed next to the Fhet/Human *actual-CDS*
                           fraction so the size of the approximation is visible.
  * Yeast                — N/A (literature theta only; no syn/NS split).

Data reading (TSV paths, class columns, s_filter for AMR) is reproduced from
25_comparison_table.py so the segregating-site counts and pi sums are identical
to the original table; only the denominators of the syn/NS terms change.

Outputs (in vcf/):
    comparison_table_L_cor.tsv     machine-readable
    comparison_table_L_cor.txt     formatted
"""
import csv, os, sys

BASE = os.environ.get("MTG_BASE", os.path.expanduser("~/Projects/MT_Genomics_Cl_Ap2026"))
MTG  = f"{BASE}/MT_Genomics2"

# ── Genetic codes ─────────────────────────────────────────────────────────────
# NCBI table 2 (vertebrate mt): TGA=W, ATA=M, AGA=*, AGG=*
MT2 = {
    "TTT":"F","TTC":"F","TTA":"L","TTG":"L","TCT":"S","TCC":"S","TCA":"S","TCG":"S",
    "TAT":"Y","TAC":"Y","TAA":"*","TAG":"*","TGT":"C","TGC":"C","TGA":"W","TGG":"W",
    "CTT":"L","CTC":"L","CTA":"L","CTG":"L","CCT":"P","CCC":"P","CCA":"P","CCG":"P",
    "CAT":"H","CAC":"H","CAA":"Q","CAG":"Q","CGT":"R","CGC":"R","CGA":"R","CGG":"R",
    "ATT":"I","ATC":"I","ATA":"M","ATG":"M","ACT":"T","ACC":"T","ACA":"T","ACG":"T",
    "AAT":"N","AAC":"N","AAA":"K","AAG":"K","AGT":"S","AGC":"S","AGA":"*","AGG":"*",
    "GTT":"V","GTC":"V","GTA":"V","GTG":"V","GCT":"A","GCC":"A","GCA":"A","GCG":"A",
    "GAT":"D","GAC":"D","GAA":"E","GAG":"E","GGT":"G","GGC":"G","GGA":"G","GGG":"G",
}
# NCBI table 5 (invertebrate mt) = table 2 but AGA=S, AGG=S (Ser, not Stop)
MT5 = dict(MT2); MT5["AGA"] = "S"; MT5["AGG"] = "S"

COMP = str.maketrans("ACGTN", "TGCAN")
BASES = "ACGT"


def codon_S_sites(codon, code):
    """Nei-Gojobori synonymous-site count for one codon (0..3)."""
    ref = code.get(codon)
    if ref is None or ref == "*":
        return None
    s = 0.0
    for pos in range(3):
        nsyn = 0
        for b in BASES:
            if b == codon[pos]:
                continue
            alt = codon[:pos] + b + codon[pos+1:]
            if code.get(alt, "?") == ref:
                nsyn += 1
        s += nsyn / 3.0
    return s


def count_sites(cds_seq, code):
    """Return (L_S, L_N, n_codons) over a coding-sense CDS (truncated to codons)."""
    S = 0.0; n = 0
    L = len(cds_seq) - (len(cds_seq) % 3)
    for i in range(0, L, 3):
        cs = codon_S_sites(cds_seq[i:i+3], code)
        if cs is None:
            continue
        S += cs; n += 1
    return S, 3*n - S, n


def code_level_fraction(code):
    """Equal-weight mean synonymous fraction over all sense codons (S_sites/3)."""
    tot = 0.0; n = 0
    for c in (a+b+d for a in BASES for b in BASES for d in BASES):
        cs = codon_S_sites(c, code)
        if cs is None:
            continue
        tot += cs; n += 1
    return (tot / n) / 3.0       # fraction of sites that are synonymous


# ── FASTA / GFF readers ───────────────────────────────────────────────────────
def read_fasta(path):
    seqs = {}; name = None; buf = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(buf).upper()
                name = line[1:].split()[0]; buf = []
            else:
                buf.append(line)
    if name is not None:
        seqs[name] = "".join(buf).upper()
    return seqs


def read_gff_cds(path):
    out = []
    with open(path) as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 9 or p[2] != "CDS":
                continue
            out.append({"chrom": p[0], "start": int(p[3]), "end": int(p[4]), "strand": p[6]})
    return out


def cds_seq(genome, chrom, start, end, strand):
    seq = genome[chrom][start-1:end]
    if strand in ("-", -1):
        seq = seq.translate(COMP)[::-1]
    return seq


# ── L_S / L_N per dataset ─────────────────────────────────────────────────────
def L_fhet():
    g = read_fasta(f"{MTG}/Missing_Files/SSM_MT_ref/Fhet_MT.fasta")
    cdss = read_gff_cds(f"{MTG}/Missing_Files/SSM_MT_ref/Fhet_MT.gff")
    S = N = 0.0
    for c in cdss:
        s, n, _ = count_sites(cds_seq(g, c["chrom"], c["start"], c["end"], c["strand"]), MT2)
        S += s; N += n
    return S, N, len(cdss)


# rCRS CDS coordinates (1-based incl.) + strand — from 22_human_mt_cds_pi.py
HUMAN_CDS = [
    ("MT-ND1",3307,4262,1),("MT-ND2",4470,5511,1),("MT-CO1",5904,7445,1),
    ("MT-CO2",7586,8269,1),("MT-ATP8",8366,8572,1),("MT-ATP6",8527,9207,1),
    ("MT-CO3",9207,9990,1),("MT-ND3",10059,10404,1),("MT-ND4L",10470,10766,1),
    ("MT-ND4",10760,12137,1),("MT-ND5",12337,14148,1),("MT-ND6",14149,14673,-1),
    ("MT-CYB",14747,15887,1),
]
def L_human():
    g = read_fasta(f"{BASE}/Human_mt/rCRS_NC_012920.fasta")
    chrom = next(iter(g))
    S = N = 0.0
    for _, st, en, strand in HUMAN_CDS:
        s, n, _ = count_sites(cds_seq(g, chrom, st, en, strand), MT2)
        S += s; N += n
    return S, N, len(HUMAN_CDS)


# ── Dataset table (paths + class cols reproduced from script 25) ──────────────
DATASETS = [
    {"name":"Fhet","N":141,"L_CDS":11417,"tsv":f"{MTG}/vcf/pi_results.tsv",
     "cc":"Class","pc":"pi_site","syn":"synonymous","ns":"nonsynonymous","sf":None,"code":"table2"},
    {"name":"Drosophila","N":169,"L_CDS":11173,"tsv":f"{MTG}/vcf/dros_pi_results.tsv",
     "cc":"Class","pc":"pi_site","syn":"synonymous","ns":"nonsynonymous","sf":None,"code":"table5"},
    {"name":"C. elegans","N":540,"L_CDS":10299,"tsv":f"{BASE}/C_elegans/celegans_pi_per_site.tsv",
     "cc":"Class","pc":"pi_site","syn":"synonymous","ns":"nonsynonymous","sf":None,"code":"table5"},
    {"name":"Human African (Lankheet 2026)","N":1176,"L_CDS":11395,
     "tsv":f"{BASE}/Human_mt/human_mt_cds_pi_per_site.tsv",
     "cc":"Effect","pc":"pi_site","syn":"SYN","ns":"NS","sf":None,"code":"table2"},
    {"name":"Human AMR (gnomAD v3.1)","N":5718,"L_CDS":11395,
     "tsv":f"{BASE}/Hm_Mt/amr_pi_per_site.tsv",
     "cc":"Class","pc":"pi_site","syn":"synonymous","ns":"nonsynonymous",
     "sf":("AC_hom_amr",0),"code":"table2"},
    {"name":"Yeast (lit. theta only)","N":1011,"L_CDS":6684,"tsv":None,
     "cc":None,"pc":"pi_site","syn":None,"ns":None,"sf":None,"code":None,
     "lit_S":384,"lit_theta":0.00766},
]


def harmonic(n):
    return sum(1.0/i for i in range(1, n))


def read_counts(ds):
    """Return S, S_syn, S_ns, pi_total, piS_raw, piN_raw (script-25 logic)."""
    with open(ds["tsv"]) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    cc, pc, sf = ds["cc"], ds["pc"], ds["sf"]
    pt = ps = pn = 0.0; S = Ss = Sn = 0
    for r in rows:
        try:
            pv = float(r[pc])
        except (KeyError, ValueError):
            continue
        seg = True
        if sf:
            try:
                seg = float(r[sf[0]]) > sf[1]
            except (KeyError, ValueError):
                seg = False
        pt += pv
        if seg: S += 1
        if cc and r.get(cc) == ds["syn"]:
            ps += pv
            if seg: Ss += 1
        elif cc and r.get(cc) == ds["ns"]:
            pn += pv
            if seg: Sn += 1
    return S, Ss, Sn, pt, ps, pn


def main():
    # ---- L_S / L_N ----
    fS, fN, fn = L_fhet()
    hS, hN, hn = L_human()
    frac2_code = code_level_fraction(MT2)
    frac5_code = code_level_fraction(MT5)
    frac_fhet  = fS / (fS + fN)
    frac_human = hS / (hS + hN)

    print("── L_S / L_N provenance ──", file=sys.stderr)
    print(f"Fhet  (actual CDS, table2): L_S={fS:.1f} L_N={fN:.1f}  syn-frac={frac_fhet:.4f} ({fn} CDS)", file=sys.stderr)
    print(f"Human (actual rCRS, table2): L_S={hS:.1f} L_N={hN:.1f}  syn-frac={frac_human:.4f} ({hn} CDS)", file=sys.stderr)
    print(f"code-level syn-frac: table2={frac2_code:.4f}  table5={frac5_code:.4f}", file=sys.stderr)
    print(f"  (Fhet/Human actual vs table2 code-level: {frac_fhet:.4f}/{frac_human:.4f} vs {frac2_code:.4f}"
          f"  → code-level proxy error ≈ {abs(frac_fhet-frac2_code):.4f}/{abs(frac_human-frac2_code):.4f})",
          file=sys.stderr)

    def LS_LN(ds):
        if ds["name"] == "Fhet":
            return fS, fN, "exact(t2)"
        if ds["name"].startswith("Human"):
            return hS, hN, "exact(t2)"
        if ds["code"] == "table5":
            LS = frac5_code * ds["L_CDS"]
            return LS, ds["L_CDS"] - LS, "code-est(t5)"
        return None, None, "N/A"

    rows_out = []
    for ds in DATASETS:
        if ds.get("lit_theta"):
            rows_out.append({"Dataset":ds["name"],"N":ds["N"],"a1":round(harmonic(ds["N"]),4),
                "L_CDS":ds["L_CDS"],"L_S":"N/A","L_N":"N/A","L_source":"N/A",
                "S":ds["lit_S"],"S_syn":"N/A","S_ns":"N/A",
                "theta_total":ds["lit_theta"],"theta_syn":"N/A","theta_ns":"N/A","pNpS_theta":"N/A",
                "pi_total":"N/A","pi_syn":"N/A","pi_ns":"N/A","pNpS_pi":"N/A"})
            continue
        a1 = harmonic(ds["N"])
        S, Ss, Sn, pt, ps, pn = read_counts(ds)
        LS, LN, src = LS_LN(ds)
        th_t = S  / (a1 * ds["L_CDS"])
        th_s = Ss / (a1 * LS)
        th_n = Sn / (a1 * LN)
        pnps_th = th_n / th_s
        pi_t = pt / ds["L_CDS"]
        pi_s = ps / LS
        pi_n = pn / LN
        pnps_pi = pi_n / pi_s
        rows_out.append({"Dataset":ds["name"],"N":ds["N"],"a1":round(a1,4),
            "L_CDS":ds["L_CDS"],"L_S":round(LS,1),"L_N":round(LN,1),"L_source":src,
            "S":S,"S_syn":Ss,"S_ns":Sn,
            "theta_total":round(th_t,5),"theta_syn":round(th_s,5),"theta_ns":round(th_n,5),
            "pNpS_theta":round(pnps_th,3),
            "pi_total":round(pi_t,5),"pi_syn":round(pi_s,5),"pi_ns":round(pi_n,5),
            "pNpS_pi":round(pnps_pi,3)})

    fields = ["Dataset","N","a1","L_CDS","L_S","L_N","L_source","S","S_syn","S_ns",
              "theta_total","theta_syn","theta_ns","pNpS_theta",
              "pi_total","pi_syn","pi_ns","pNpS_pi"]
    out_tsv = f"{MTG}/vcf/comparison_table_L_cor.tsv"
    with open(out_tsv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        w.writeheader()
        for r in rows_out:
            w.writerow(r)
    print(f"\nWrote {out_tsv}", file=sys.stderr)

    out_txt = f"{MTG}/vcf/comparison_table_L_cor.txt"
    with open(out_txt, "w") as fh:
        fh.write("Site-length-corrected cross-species mt-CDS diversity comparison\n")
        fh.write("pN/pS = (S_ns/L_N)/(S_syn/L_S) [theta] and (piN/L_N)/(piS/L_S) [pi]\n")
        fh.write(f"L_S/L_N: Fhet & Human exact (Nei-Gojobori, NCBI table 2); "
                 f"Drosophila & C.elegans code-level estimate (table 5, "
                 f"syn-frac={frac5_code:.4f}); Yeast N/A.\n")
        fh.write(f"Validation: table-2 code-level syn-frac={frac2_code:.4f} vs "
                 f"actual Fhet {frac_fhet:.4f} / Human {frac_human:.4f}.\n\n")
        hdr = ["Dataset","N","L_CDS","L_S","L_N","src","S","Ssyn","Sns",
               "th_tot","th_syn","th_ns","pN/pS_th","pi_tot","pi_syn","pi_ns","pN/pS_pi"]
        fh.write("  ".join(f"{h:>9}" if i else f"{h:<30}" for i,h in enumerate(hdr)) + "\n")
        for r in rows_out:
            vals = [r["Dataset"],r["N"],r["L_CDS"],r["L_S"],r["L_N"],r["L_source"],
                    r["S"],r["S_syn"],r["S_ns"],r["theta_total"],r["theta_syn"],
                    r["theta_ns"],r["pNpS_theta"],r["pi_total"],r["pi_syn"],
                    r["pi_ns"],r["pNpS_pi"]]
            fh.write("  ".join(f"{str(v):>9}" if i else f"{str(v):<30}" for i,v in enumerate(vals)) + "\n")
    print(f"Wrote {out_txt}", file=sys.stderr)

    # echo to stdout
    print("\n".join("\t".join(str(r[f]) for f in fields) for r in [dict(zip(fields,fields))]+rows_out))


if __name__ == "__main__":
    main()

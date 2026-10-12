#!/usr/bin/env python3
"""
32_make_paml_input.py  —  build PAML-ready codon alignment + control files for
the outgroup dN/dS step (Fhet clades vs F. catenatus).

Produces, in MT_Genomics2/paml/:
  fhet_caten.phy   4-sequence codon alignment (PHYLIP sequential):
                     Fhet_ref  (reference NC_012312.1, ~northern)
                     Fhet_N    (North-clade consensus)
                     Fhet_S    (South-clade consensus)
                     F_caten   (F. catenatus NC_084079.1)
                   13 mt CDS concatenated, codon-aligned, equal length.
  codeml.ctl       pairwise ML dN/dS (runmode=-2), vertebrate-mt code (icode=1)
  yn00.ctl         Yang-Nielsen counting estimate (cross-check)

The three Fhet sequences share the reference's coordinates (consensus differs
from REF only by SNP substitutions, no indels), so they are mutually aligned;
only the Fhet-vs-catenatus frame needs alignment, done per gene at the protein
level and applied identically to ref/N/S.  Codon columns containing a stop or N
in any sequence are dropped (keeps frame, avoids PAML stop-codon errors).
"""
import importlib.util, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("mk", os.path.join(HERE, "31_mk_outgroup.py"))
mk = importlib.util.module_from_spec(spec); spec.loader.exec_module(mk)

MTG = mk.MTG
OUTDIR = f"{MTG}/paml" if mk.OG_NAME == "catenatus" else f"{MTG}/paml_{mk.OG_NAME}"
os.makedirs(OUTDIR, exist_ok=True)

def build():
    coords = mk.gff_cds_coords()
    ref = mk.read_fasta(mk.FHET_FA); rs = ref[next(iter(ref))]
    _, var, clade = mk.parse_vcf()
    Ni = [i for i, c in enumerate(clade) if c == "N"]
    Si = [i for i, c in enumerate(clade) if c == "S"]
    gN = mk.clade_consensus_genome(rs, var, Ni)
    gS = mk.clade_consensus_genome(rs, var, Si)
    cds_ref = mk.cds_from_genome(rs, coords)
    cds_N   = mk.cds_from_genome(gN, coords)
    cds_S   = mk.cds_from_genome(gS, coords)
    cds_O   = mk.outgroup_cds()
    genes = [g for g, _, _, _ in coords]

    seqs = {"Fhet_ref": [], "Fhet_N": [], "Fhet_S": [], "F_caten": []}
    for g in genes:
        if g not in cds_O:
            continue
        ref_cds = cds_ref[g]
        # codon lists (truncate to codon multiple); ref/N/S share length & frame
        def codons(s): return [s[i:i+3] for i in range(0, len(s)-len(s)%3, 3)]
        cr, cn, cs_, co = codons(ref_cds), codons(cds_N[g]), codons(cds_S[g]), codons(cds_O[g])
        pr, po = mk.translate(ref_cds), mk.translate(cds_O[g])
        ar, ao = mk.nw_align(pr, po)
        i = j = 0
        for x, y in zip(ar, ao):
            kr = cr[i] if x != "-" else "---"
            kn = cn[i] if x != "-" else "---"
            ks = cs_[i] if x != "-" else "---"
            ko = co[j] if y != "-" else "---"
            if x != "-": i += 1
            if y != "-": j += 1
            quad = [kr, kn, ks, ko]
            # drop columns with N or stop in any sequence
            bad = any(("N" in c) or (mk.MT2.get(c) == "*") for c in quad if c != "---")
            if bad:
                continue
            seqs["Fhet_ref"].append(kr); seqs["Fhet_N"].append(kn)
            seqs["Fhet_S"].append(ks);  seqs["F_caten"].append(ko)

    aln = {k: "".join(v) for k, v in seqs.items()}
    L = len(aln["Fhet_ref"])
    assert all(len(v) == L for v in aln.values()), "unequal lengths"
    assert L % 3 == 0, "not codon multiple"

    phy = f"{OUTDIR}/fhet_caten.phy"
    with open(phy, "w") as fh:
        fh.write(f" 4 {L}\n")
        for name in ("Fhet_ref", "Fhet_N", "Fhet_S", "F_caten"):
            fh.write(f"{name:<12}{aln[name]}\n")
    print(f"Wrote {phy}  ({L} bp, {L//3} codons)")

    with open(f"{OUTDIR}/codeml.ctl", "w") as fh:
        fh.write(
            "      seqfile = fhet_caten.phy\n"
            "     treefile = fhet_caten.tree\n"
            "      outfile = codeml_pairwise.out\n"
            "        noisy = 3\n      verbose = 1\n"
            "      runmode = -2\n"      # pairwise
            "      seqtype = 1\n"       # codons
            "    CodonFreq = 2\n"       # F3x4
            "        icode = 1\n"       # vertebrate/mammalian mt = NCBI table 2
            "        model = 0\n      NSsites = 0\n"
            "    fix_kappa = 0\n        kappa = 2\n"
            "    fix_omega = 0\n        omega = 0.4\n"
            "    cleandata = 0\n")
    # dummy tree (runmode=-2 ignores topology but codeml wants the file present)
    with open(f"{OUTDIR}/fhet_caten.tree", "w") as fh:
        fh.write("4 1\n(Fhet_ref,Fhet_N,Fhet_S,F_caten);\n")
    with open(f"{OUTDIR}/yn00.ctl", "w") as fh:
        fh.write(
            "      seqfile = fhet_caten.phy\n"
            "      outfile = yn00.out\n"
            "      verbose = 1\n"
            "        icode = 1\n"       # vertebrate mt
            "    weighting = 0\n"
            "   commonf3x4 = 0\n")
    print(f"Wrote {OUTDIR}/codeml.ctl, yn00.ctl, fhet_caten.tree")

if __name__ == "__main__":
    build()

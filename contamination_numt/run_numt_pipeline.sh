#!/usr/bin/env bash
# SwarmGenomics NUMT identification (Kylmanen et al. 2026, MER, doi:10.1111/1755-0998.70119)
# for Fundulus heteroclitus MU-UCD_Fhet_4.1 (GCF_011125445.2), plus a heteroplasmy
# cross-check that tests whether NUMTs could explain the observed heteroplasmy sites.
#
# Usage:  ./run_numt_pipeline.sh  <nuclear_genome.fna>  <mito.fasta>  <heteroplasmy_events.tsv>
# Requires: blastn, makeblastdb, bedtools, python3, awk
set -euo pipefail

GENOME="${1:?nuclear genome FASTA}"
MITO="${2:?mito query FASTA}"
HPEVENTS="${3:-}"          # optional: stage-14 pileup events tsv for the cross-check
OUT="numt_out"; mkdir -p "$OUT"
MINLEN=200                 # SwarmGenomics: drop fragments < 200 bp
WORDSIZE=20                # SwarmGenomics / Lammers 2017 / Vendrami 2022

echo "[1/6] makeblastdb"
makeblastdb -in "$GENOME" -parse_seqids -dbtype nucl >/dev/null

echo "[2/6] blastn (mito query vs nuclear genome, word_size=$WORDSIZE)"
blastn -db "$GENOME" -query "$MITO" -outfmt \
  '7 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore' \
  -word_size "$WORDSIZE" -out "$OUT/blast_out.tsv"

echo "[3/6] clean table (drop comment lines)"
grep -v '^#' "$OUT/blast_out.tsv" > "$OUT/blast_clean.tsv" || true

echo "[4/6] BED of NUMT loci on the nuclear genome (length >= ${MINLEN} bp)"
# columns: qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
awk -v OFS='\t' -v ML=$MINLEN '$4>=ML{
  s=($9<$10)?$9:$10; e=($9<$10)?$10:$9;             # subject (nuclear) start/end, 1-based
  name=$1"|mt:"$7"-"$8"|pid"$3"|len"$4;
  print $2, s-1, e, name, $3, ($9<$10?"+":"-")
}' "$OUT/blast_clean.tsv" | sort -k1,1 -k2,2n > "$OUT/numt.bed"

echo "[5/6] extract NUMT sequences"
bedtools getfasta -fi "$GENOME" -bed "$OUT/numt.bed" -name -fo "$OUT/numt.fasta" 2>/dev/null || \
bedtools getfasta -fi "$GENOME" -bed "$OUT/numt.bed" -fo "$OUT/numt.fasta"

echo "[6/6] summaries"
python3 - "$OUT" "$HPEVENTS" <<'PY'
import sys, csv, os
outdir, hpevents = sys.argv[1], (sys.argv[2] if len(sys.argv)>2 else "")
rows=[]
with open(os.path.join(outdir,"blast_clean.tsv")) as f:
    for L in f:
        p=L.rstrip("\n").split("\t")
        if len(p)<12: continue
        rows.append(dict(qseqid=p[0],sseqid=p[1],pident=float(p[2]),length=int(p[3]),
                         qstart=int(p[6]),qend=int(p[7]),sstart=int(p[8]),send=int(p[9])))
hits=[r for r in rows if r["length"]>=200]
scaf={}
for r in hits: scaf.setdefault(r["sseqid"],[]).append(r)
with open(os.path.join(outdir,"numt_summary.txt"),"w") as o:
    o.write("NUMT BLASTN summary (SwarmGenomics; length>=200bp)\n")
    o.write(f"  total HSPs (all)      : {len(rows)}\n")
    o.write(f"  NUMT HSPs (>=200bp)   : {len(hits)}\n")
    o.write(f"  nuclear scaffolds hit : {len(scaf)}\n")
    o.write(f"  total NUMT bp         : {sum(r['length'] for r in hits)}\n")
    if hits:
        pid=sorted(r['pident'] for r in hits)
        o.write(f"  %identity min/median/max: {pid[0]:.1f} / {pid[len(pid)//2]:.1f} / {pid[-1]:.1f}\n")
    o.write("\n  Top scaffolds by summed NUMT length:\n")
    for s,rs in sorted(scaf.items(),key=lambda kv:-sum(r['length'] for r in kv[1]))[:15]:
        o.write(f"    {s:20s} n_hsp={len(rs):3d}  bp={sum(r['length'] for r in rs):6d}  "
                f"maxpid={max(r['pident'] for r in rs):.1f}\n")
    # mito-coordinate coverage: which mito positions are covered by ANY NUMT
    covered=set()
    for r in hits:
        a,b=sorted((r["qstart"],r["qend"]))
        covered.update(range(a,b+1))
    o.write(f"\n  Mitochondrial positions covered by >=1 NUMT HSP: {len(covered)}\n")
    # cross-check vs heteroplasmy sites
    if hpevents and os.path.exists(hpevents):
        hp=[]; priv=[]
        with open(hpevents) as f:
            rd=csv.DictReader(f,delimiter="\t")
            for row in rd:
                pos=int(row["POS"]); hp.append(pos)
                if row.get("Hp_class")=="private_alt_Hp": priv.append(pos)
        hpset=set(hp); privset=set(priv)
        o.write("\n  === Heteroplasmy cross-check ===\n")
        o.write(f"  distinct Hp sites            : {len(hpset)}\n")
        o.write(f"    ...falling inside a NUMT    : {len(hpset & covered)}\n")
        o.write(f"    ...NOT explainable by NUMT  : {len(hpset - covered)}\n")
        o.write(f"  distinct private Hp sites     : {len(privset)}\n")
        o.write(f"    ...falling inside a NUMT    : {len(privset & covered)}\n")
        o.write(f"    ...NOT explainable by NUMT  : {len(privset - covered)}\n")
        o.write("  (An Hp site OUTSIDE every NUMT-aligned mito region cannot be a NUMT artifact.)\n")
print(open(os.path.join(outdir,"numt_summary.txt")).read())
PY
echo "DONE -> $OUT/{blast_out.tsv,numt.bed,numt.fasta,numt_summary.txt}"

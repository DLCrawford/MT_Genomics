#!/usr/bin/env bash
#BSUB -J simpleVCF_contam
#BSUB -P fun_gen_1
#BSUB -q normal
#BSUB -n 12
#BSUB -R "rusage[mem=24000M]"
#BSUB -W 72:00
#BSUB -o /projectnb/dcrawford/MT_Genomics2/NUMT_nuclear/simpleVCF.out
#BSUB -e /projectnb/dcrawford/MT_Genomics2/NUMT_nuclear/simpleVCF.err
#
# "Simple VCF" contamination test from the BAMs (NOT the GATK hard.final VCF, which can
# truncate low-frequency 3rd alleles via calling thresholds / MAF filters).
# Raw bcftools mpileup -a AD keeps every allele it observes -> preserves the ~13% signal.
# Metric = heteroplasmy-style (0.10<=AD/DP<0.35); contrast = 53 hp vs 88 no_hp (groups_hp.tsv).
# Submit as a BATCH job (it is long):   bsub < run_simpleVCF_allfish.sh
set -euo pipefail

### ─── EDIT THESE ──────────────────────────────────────────────────────────
WORKDIR=/projectnb/dcrawford/MT_Genomics2/NUMT_nuclear
GENOME=/projectnb/dcrawford/Fhet_genome_4.1/GCF_011125445.2/GCF_011125445.2_MU-UCD_Fhet_4.1_genomic.fna
BAMLIST=$WORKDIR/bamlist_panel.txt   # all 141 panel BAMs (full paths). Subset to speed up.
FOCALS=$WORKDIR/focals.txt
GROUPSFILE=$WORKDIR/groups_hp.tsv
SCRIPTDIR=$WORKDIR                    # holds contam_allfish.py
REGION=NC_046361.1:1-10000000        # 10 Mb of chr1 (~>100k sites). Widen once it works.
THREADS=12
### ────────────────────────────────────────────────────────────────────────

source /projectnb/dcrawford/MT_Genomics2/jobs/config.sh 2>/dev/null || true
command -v bcftools >/dev/null || { echo "ERROR: bcftools not on PATH"; exit 1; }
command -v python3  >/dev/null || { echo "ERROR: python3 not on PATH"; exit 1; }
echo "bcftools: $(bcftools --version | head -1)"
ver=$(bcftools --version | head -1 | grep -oE '[0-9]+\.[0-9]+' | head -1)
awk -v v="$ver" 'BEGIN{split(v,a,"."); exit !(a[1]>1 || (a[1]==1 && a[2]>=15))}' || { echo "bcftools $ver too old"; exit 1; }
[ -f "$GENOME.fai" ] || { echo "ERROR: $GENOME.fai missing (samtools faidx $GENOME)"; exit 1; }

cd "$WORKDIR"; OUT=simpleVCF_allfish_out; mkdir -p "$OUT"; VCF="$OUT/chr_snps.vcf.gz"
if [ -s "$VCF" ]; then
  echo "[1/4] reusing $VCF (delete to force re-pileup)"
else
  echo "[1/4] mpileup $(wc -l < "$BAMLIST") BAMs over $REGION (raw AD; MAPQ>=30 BQ>=20) ..."
  bcftools mpileup -f "$GENOME" -b "$BAMLIST" -r "$REGION" \
      -a AD,DP -q 30 -Q 20 --max-depth 500 --threads "$THREADS" -Ou \
    | bcftools norm -m -any -f "$GENOME" -Ou \
    | bcftools view -m2 -M2 -v snps -Oz -o "$VCF"
fi
echo "[2/4] sample names (cleaned) -> $OUT/samples.txt"
bcftools query -l "$VCF" | sed 's#.*/##; s/\.merged\.bam$//; s/\.bam$//' > "$OUT/samples.txt"
while read -r f; do [ -z "$f" ] && continue
  grep -qxF "$f" "$OUT/samples.txt" || { echo "focal '$f' not in samples (see $OUT/samples.txt)"; exit 1; }
done < "$FOCALS"
echo "[3/4] extract per-sample AD"
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%AD]\n' "$VCF" > "$OUT/ad.tsv"
echo "    sites: $(wc -l < "$OUT/ad.tsv")"
echo "[4/4] metric + hp-vs-no_hp contrast"
python3 "$SCRIPTDIR/contam_allfish.py" "$OUT/ad.tsv" "$OUT/samples.txt" "$FOCALS" "$OUT" "$GROUPSFILE"
echo "DONE -> $OUT/{per_sample_contam_allfish.tsv, verdict.txt}"
cat "$OUT/verdict.txt"

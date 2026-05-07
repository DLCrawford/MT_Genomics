#!/bin/bash
###############################################################################
# 04_bwa_align_mt.sh
# Stage : MT alignment (bwa mem | samtools view | samtools sort), then index
# Input : $TRIM_DIR/{sample}_1_val_1.fq.gz, {sample}_2_val_2.fq.gz
#         (Trim Galore paired output naming — NOT Trimmomatic _p/_up; see
#         CLAUDE.md "Sample / sample-list conventions" and CHANGELOG 2026-05-07.)
#         $REF (Fhet_MT.fasta — must be bwa+samtools-indexed first; see docs)
# Output: $BAMS_DIR/{sample}_MT.bam + {sample}_MT.bam.bai
# Notes : Reference indexing is a one-time login-node step:
#             samtools faidx $REF
#             bwa index $REF
# Submit: bsub < jobs/04_bwa_align_mt.sh
###############################################################################

#Job Name and array call.
#BSUB -J fhet_align_mt[1-144]%30
#BSUB -P fun_gen_1
#BSUB -q normal
#BSUB -n 12
#BSUB -R "rusage[mem=12000M]"
#BSUB -W 72:00
#BSUB -o /projectnb/dcrawford/MT_Genomics2/logs/04_align_mt_%I.out
#BSUB -e /projectnb/dcrawford/MT_Genomics2/logs/04_align_mt_%I.err

set -euo pipefail

### ─── LOAD ENVIRONMENT + CONFIG ───────────────────────────────────────────
module load anaconda3
source /projectnb/dcrawford/MT_Genomics2/jobs/config.sh
conda activate "$CONDA_ENV"
# $CONDA_ENV must provide: bwa, samtools
### ─── PER-TASK SETUP ──────────────────────────────────────────────────────
THREADS=12
mkdir -p "$BAMS_DIR"

sample=$(sed -n "${LSB_JOBINDEX}p" "$SAMPLE_LIST")

r1=$TRIM_DIR/${sample}_1_val_1.fq.gz
r2=$TRIM_DIR/${sample}_2_val_2.fq.gz
OUT_BAM=$BAMS_DIR/${sample}_MT.bam

if [[ ! -f "$r1" || ! -f "$r2" ]]; then
    echo "[ERROR] Missing trimmed input for ${sample}: $r1 or $r2 not found" >&2
    exit 1
fi

### ─── RUN: ALIGN -> SORT -> INDEX ─────────────────────────────────────────
echo "[$(date)] Aligning $sample to $(basename "$REF")"

bwa mem -t "$THREADS" -M "$REF" "$r1" "$r2" \
  | samtools view -bS - \
  | samtools sort -@ "$THREADS" -o "$OUT_BAM"

samtools index "$OUT_BAM"

echo "[$(date)] Done: $sample"

#!/bin/bash
###############################################################################
# 02_trim_pe.sh
# Stage : Adapter+quality trimming with Trimmomatic PE
# Input : $RAW_DIR/{sample}_1.fq.gz, {sample}_2.fq.gz
# Output: $TRIM_DIR/{sample}_1_p.fq.gz   (paired R1)
#         $TRIM_DIR/{sample}_1_up.fq.gz  (unpaired R1)
#         $TRIM_DIR/{sample}_2_p.fq.gz   (paired R2)
#         $TRIM_DIR/{sample}_2_up.fq.gz  (unpaired R2)
# Notes : $ADAPTERS is TruSeq3-PE.fa + NexteraPE-PE.fa concatenated:
#             cat TruSeq3-PE.fa NexteraPE-PE.fa > CombinedAdapters.fa
# Submit: bsub < jobs/02_trim_pe.sh
###############################################################################

#Job Name and array call.
#BSUB -J fhet_trim_pe[1-287]%50
#BSUB -P fun_gen_1
#BSUB -q normal
#BSUB -n 15
#BSUB -R "rusage[mem=12000M]"
#BSUB -W 72:00
#BSUB -o /projectnb/dcrawford/MT_Genomics2/logs/02_trim_pe_%I.out
#BSUB -e /projectnb/dcrawford/MT_Genomics2/logs/02_trim_pe_%I.err

set -euo pipefail

### ─── LOAD ENVIRONMENT + CONFIG ───────────────────────────────────────────
# module load + conda shell hook are handled inside config.sh
source /projectnb/dcrawford/MT_Genomics2/jobs/config.sh
conda activate "$CONDA_ENV"
# $CONDA_ENV must provide: java (for Trimmomatic jar) — or install
# trimmomatic via conda and call `trimmomatic PE` instead of `java -jar`.

### ─── PER-TASK SETUP ──────────────────────────────────────────────────────
THREADS=15
mkdir -p "$TRIM_DIR"

SAMPLE=$(sed -n "${LSB_JOBINDEX}p" "$SAMPLE_LIST")

R1="$RAW_DIR/${SAMPLE}_1.fq.gz"
R2="$RAW_DIR/${SAMPLE}_2.fq.gz"
OUT1PAIR="$TRIM_DIR/${SAMPLE}_1_p.fq.gz"
OUT1UNPAIR="$TRIM_DIR/${SAMPLE}_1_up.fq.gz"
OUT2PAIR="$TRIM_DIR/${SAMPLE}_2_p.fq.gz"
OUT2UNPAIR="$TRIM_DIR/${SAMPLE}_2_up.fq.gz"

if [[ ! -f "$R1" || ! -f "$R2" ]]; then
    echo "[ERROR] Missing input for ${SAMPLE}: $R1 or $R2 not found" >&2
    exit 1
fi

### ─── RUN ─────────────────────────────────────────────────────────────────
cd "$TRIM_DIR"

echo "[$(date)] Trimmomatic PE on $SAMPLE"
java -jar "$TRIMJAR" PE \
    -threads "$THREADS" \
    "$R1" "$R2" \
    "$OUT1PAIR" "$OUT1UNPAIR" \
    "$OUT2PAIR" "$OUT2UNPAIR" \
    ILLUMINACLIP:"$ADAPTERS":2:30:10:2:True \
    LEADING:3 TRAILING:3 MINLEN:36
echo "[$(date)] Done: $SAMPLE"

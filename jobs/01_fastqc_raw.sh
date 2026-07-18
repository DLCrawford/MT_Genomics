#!/bin/bash
###############################################################################
# 01_fastqc_raw.sh
# Stage : Raw-read QC (FastQC) on paired raw fastqs
# Input : $RAW_DIR/{sample}_1.fq.gz, {sample}_2.fq.gz   (one per array task)
# Output: $RAW_QC_DIR/{sample}_*_fastqc.{html,zip}
# Notes : Run multiqc separately AFTER all array tasks finish:
#             conda activate mito_genomics
#             multiqc $RAW_QC_DIR -o $RAW_QC_DIR
# Submit: bsub < jobs/01_fastqc_raw.sh
###############################################################################

#Job Name and array call. Square brackets = array size, % = max concurrent.
#BSUB -J fhet_fastqc_raw[1-287]%20
#BSUB -P fun_gen_1
#BSUB -q normal
#BSUB -n 12
#BSUB -R "rusage[mem=12000M]"
#BSUB -W 12:00
#BSUB -o /projectnb/dcrawford/MT_Genomics2/logs/01_fastqc_raw_%I.out
#BSUB -e /projectnb/dcrawford/MT_Genomics2/logs/01_fastqc_raw_%I.err

set -euo pipefail

### ─── LOAD ENVIRONMENT + CONFIG ───────────────────────────────────────────
# module load + conda shell hook are handled inside config.sh
source /projectnb/dcrawford/MT_Genomics2/jobs/config.sh
conda activate "$CONDA_ENV"
# $CONDA_ENV must provide: fastqc

### ─── PER-TASK SETUP ──────────────────────────────────────────────────────
THREADS=12
mkdir -p "$RAW_QC_DIR"

# LSB_JOBINDEX runs from 1..N. With no header in $SAMPLE_LIST use NR==jindex.
sample=$(sed -n "${LSB_JOBINDEX}p" "$SAMPLE_LIST")

R1=$RAW_DIR/${sample}_1.fq.gz
R2=$RAW_DIR/${sample}_2.fq.gz

if [[ ! -f "$R1" || ! -f "$R2" ]]; then
    echo "[ERROR] Missing input for ${sample}: $R1 or $R2 not found" >&2
    exit 1
fi

### ─── RUN ─────────────────────────────────────────────────────────────────
echo "[$(date)] FastQC on $sample"
fastqc --threads "$THREADS" "$R1" "$R2" -o "$RAW_QC_DIR"
echo "[$(date)] Done: $sample"

#!/bin/bash
###############################################################################
# 03_fastqc_trimmed.sh
# Stage : Post-trim QC (FastQC) on Trimmomatic paired outputs
# Input : $TRIM_DIR/{sample}_1_p.fq.gz, {sample}_2_p.fq.gz
# Output: $TRIM_QC_DIR/{sample}_*_p_fastqc.{html,zip}
# Notes : Roll up with multiqc separately after all array tasks finish:
#             conda activate mito_genomics
#             multiqc $TRIM_QC_DIR -o $TRIM_QC_DIR
# Submit: bsub < jobs/03_fastqc_trimmed.sh
###############################################################################

#Job Name and array call.
#BSUB -J fhet_fastqc_trimmed[1-287]%20
#BSUB -P fun_gen_1
#BSUB -q normal
#BSUB -n 12
#BSUB -R "rusage[mem=12000M]"
#BSUB -W 12:00
#BSUB -o /projectnb/dcrawford/MT_Genomics2/logs/03_fastqc_trimmed_%I.out
#BSUB -e /projectnb/dcrawford/MT_Genomics2/logs/03_fastqc_trimmed_%I.err

set -euo pipefail

### ─── LOAD ENVIRONMENT + CONFIG ───────────────────────────────────────────
module load anaconda3
source /projectnb/dcrawford/MT_Genomics2/jobs/config.sh
conda activate "$CONDA_ENV"
# $CONDA_ENV must provide: fastqc, multiqc

### ─── PER-TASK SETUP ──────────────────────────────────────────────────────
THREADS=12
mkdir -p "$TRIM_QC_DIR"

sample=$(sed -n "${LSB_JOBINDEX}p" "$SAMPLE_LIST")

R1=$TRIM_DIR/${sample}_1_p.fq.gz
R2=$TRIM_DIR/${sample}_2_p.fq.gz

if [[ ! -f "$R1" || ! -f "$R2" ]]; then
    echo "[ERROR] Missing trimmed input for ${sample}: $R1 or $R2 not found" >&2
    exit 1
fi

### ─── RUN ─────────────────────────────────────────────────────────────────
echo "[$(date)] FastQC (post-trim) on $sample"
fastqc --threads "$THREADS" "$R1" "$R2" -o "$TRIM_QC_DIR"
echo "[$(date)] Done: $sample"

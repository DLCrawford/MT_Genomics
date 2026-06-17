#!/bin/bash
###############################################################################
# 05d_persample_call.sh
# Stage 05 RE-RUN v4: per-sample call (replicating the historical recipe).
#
# Why this exists:
#   The joint-call runs in stages 05 (v1), 05b (v2), 05c (v3) all converged
#   to ~150 SNPs regardless of -Q/-q/--ploidy. The historical pipeline that
#   produced 1133 SNPs (merged_144.vcf.gz, Jul 2025; archived recipe in
#   archive/Notes_dlcs/Inital_call_wo_AD.txt) used per-sample call followed
#   by bcftools merge — a different architecture, not a different parameter
#   set. This stage replicates that architecture, but adds -a AD,DP so the
#   merged output carries the per-sample allele depths needed by stage 08
#   (haplotype caller).
#
# Recipe (per sample, parallel array):
#   bcftools mpileup -f REF BAM -a AD,DP --max-depth 10000 -Ou \
#     | bcftools call -mv --ploidy 1 -Oz -o ${SAMPLE}.vcf.gz
#   bcftools index ${SAMPLE}.vcf.gz
#   bcftools norm -m -any -Oz -o ${SAMPLE}_norm.vcf.gz ${SAMPLE}.vcf.gz
#   bcftools index ${SAMPLE}_norm.vcf.gz
#
# Notes vs current 05/05b/05c (joint call):
#   - No -A in `bcftools call` (the joint scripts have it; it's what produced
#     the 96%-multiallelic phantom-alt soup).
#   - No -Q / -q overrides — uses bcftools defaults (-Q 13, -q 0), matching
#     the historical recipe.
#   - --ploidy 1 (correct for haploid mtDNA).
#   - Per-sample VCFs land under ${VCF_DIR}/persample/; merge step is 05e.
#
# Submit:
#   bsub < jobs/05d_persample_call.sh
#   # then, AFTER the array completes successfully:
#   bsub < jobs/05e_merge_persample.sh
###############################################################################

#BSUB -J fhet_persample[1-143]%24
#BSUB -P fun_gen_1
#BSUB -q normal
#BSUB -n 4
#BSUB -R "rusage[mem=8000M] span[hosts=1]"
#BSUB -W 02:00
#BSUB -o /projectnb/dcrawford/MT_Genomics2/logs/05d_persample_%J_%I.out
#BSUB -e /projectnb/dcrawford/MT_Genomics2/logs/05d_persample_%J_%I.err

set -euo pipefail
cd /projectnb/dcrawford/MT_Genomics2

# === Load environment + paths ===
source /projectnb/dcrawford/MT_Genomics2/jobs/config.sh
#conda activate "$CONDA_ENV"

PERSAMPLE_DIR="${VCF_DIR}/persample"
mkdir -p "$PERSAMPLE_DIR"

# === Resolve sample for this array task ===
# BAM_LIST has 143 entries like "10_0_MT.bam" (1_0 excluded — no paired trimmed reads).
BAM_FILE_NAME=$(sed -n "${LSB_JOBINDEX}p" "$BAM_LIST")
if [[ -z "$BAM_FILE_NAME" ]]; then
    echo "ERROR: empty BAM filename at index ${LSB_JOBINDEX} of ${BAM_LIST}" >&2
    exit 1
fi
SAMPLE="${BAM_FILE_NAME%_MT.bam}"
BAM="${BAMS_DIR}/${BAM_FILE_NAME}"

if [[ ! -s "$BAM" ]]; then
    echo "ERROR: BAM not found or empty: ${BAM}" >&2
    exit 1
fi

OUT_RAW="${PERSAMPLE_DIR}/${SAMPLE}.vcf.gz"
OUT_NORM="${PERSAMPLE_DIR}/${SAMPLE}_norm.vcf.gz"

echo "[$(date)] task=${LSB_JOBINDEX}  sample=${SAMPLE}  bam=${BAM}"
echo "  ref:    ${REF}"
echo "  out:    ${OUT_NORM}"
echo "  bcftools: $(bcftools --version | head -1)"

# === Per-sample call (haploid, AD/DP recorded, no -A) ===
bcftools mpileup \
    -f "$REF" \
    "$BAM" \
    -a AD,DP \
    --max-depth 10000 \
    -Ou \
  | bcftools call \
        -mv \
        --ploidy 1 \
        -Oz \
        -o "$OUT_RAW"
bcftools index -f "$OUT_RAW"

# === Normalize: split multiallelics into one ALT per row ===
# Done per-sample so the post-merge VCF has 1 row per (POS, ALT).
bcftools norm \
    -m -any \
    -Oz \
    -o "$OUT_NORM" \
    "$OUT_RAW"
bcftools index -f "$OUT_NORM"

echo "[$(date)] DONE  sample=${SAMPLE}  norm=${OUT_NORM}"

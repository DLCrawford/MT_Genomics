#!/bin/bash
###############################################################################
# 05_bcftools_mpileup_call_AD.sh
# Stage : Joint variant calling across all MT BAMs with allele depth (AD)
# Input : $REF (Fhet_MT.fasta)
#         $BAMS_DIR/*.bam     (sorted+indexed per-sample MT BAMs)
#         $BAM_LIST           (one BAM filename per line, e.g. 10_0_MT.bam)
# Output: ${OUTPREFIX}_fullAD.vcf.gz       (all sites, with AD/DP)
#         ${OUTPREFIX}_variantsAD.vcf.gz   (variant-only subset)
#         + .csi indices and bcftools stats files
# Notes : This supersedes the older per-sample SNP caller (now archived).
#         Sample names in the VCF are reheadered to strip the trailing "_0".
# Submit: bsub < jobs/05_bcftools_mpileup_call_AD.sh
###############################################################################

#Job Name (NOT an array — single joint call across all samples).
#BSUB -J fhet_mpileup_AD
#BSUB -P fun_gen_1
#BSUB -q normal
#BSUB -n 8
#BSUB -R "rusage[mem=16000M] span[hosts=1]"
#BSUB -W 72:00
#BSUB -o /projectnb/dcrawford/MT_Genomics2/logs/05_mpileup_AD_%J.out
#BSUB -e /projectnb/dcrawford/MT_Genomics2/logs/05_mpileup_AD_%J.err

set -euo pipefail

### ─── LOAD ENVIRONMENT + CONFIG ───────────────────────────────────────────
# module load + conda shell hook are handled inside config.sh
source /projectnb/dcrawford/MT_Genomics2/jobs/config.sh
conda activate "$CONDA_ENV"
# $CONDA_ENV must provide: bcftools

### ─── PER-RUN SETUP ───────────────────────────────────────────────────────
mkdir -p "$VCF_DIR"
OUTPREFIX="${VCF_DIR}/Fhet_mt"
BAM_LIST_FULL="${VCF_DIR}/bam_full_paths.txt"

### ─── PREPARE FULL-PATH BAM LIST ──────────────────────────────────────────
echo "[$(date)] Preparing BAM file list..."
awk -v dir="$BAMS_DIR" '{print dir "/" $1}' "$BAM_LIST" > "$BAM_LIST_FULL"

echo "=== BAM files to be used ==="
head "$BAM_LIST_FULL"
echo "  ..."
tail -n 3 "$BAM_LIST_FULL"
echo "============================"

### ─── MPILEUP + CALL (ALL SITES, WITH AD/DP) ──────────────────────────────
# -A in `bcftools call` means: keep all positions (not just variants).
# -a AD,DP records per-allele depth and total depth (needed for haplotype work).
# -Q 30 -q 30 filters low-quality bases/mappings; -d 100000 raises depth cap.
echo "[$(date)] Running bcftools mpileup + call across all sites..."
bcftools mpileup \
    -f "$REF" \
    -b "$BAM_LIST_FULL" \
    -a AD,DP \
    -Q 30 -q 30 -d 100000 \
    -Ou \
  | bcftools call \
        -mv -A -Oz \
        -o "${OUTPREFIX}_fullAD.vcf.gz"

bcftools index -f "${OUTPREFIX}_fullAD.vcf.gz"

### ─── RENAME SAMPLES (STRIP trailing "_0") ────────────────────────────────
echo "[$(date)] Renaming samples to drop '_0' suffix..."
bcftools reheader \
    -s <(bcftools query -l "${OUTPREFIX}_fullAD.vcf.gz" | sed 's/_0//') \
    -o "${OUTPREFIX}_fullAD.renamed.vcf.gz" \
    "${OUTPREFIX}_fullAD.vcf.gz"

mv "${OUTPREFIX}_fullAD.renamed.vcf.gz" "${OUTPREFIX}_fullAD.vcf.gz"
bcftools index -f "${OUTPREFIX}_fullAD.vcf.gz"

### ─── VARIANT-ONLY SUBSET ─────────────────────────────────────────────────
echo "[$(date)] Extracting variant-only sites..."
bcftools view \
    -v snps,indels \
    -Oz -o "${OUTPREFIX}_variantsAD.vcf.gz" \
    "${OUTPREFIX}_fullAD.vcf.gz"

bcftools index -f "${OUTPREFIX}_variantsAD.vcf.gz"

# Reheader the variant VCF too, so sample names match the full VCF.
bcftools reheader \
    -s <(bcftools query -l "${OUTPREFIX}_variantsAD.vcf.gz" | sed 's/_0//') \
    -o "${OUTPREFIX}_variantsAD.renamed.vcf.gz" \
    "${OUTPREFIX}_variantsAD.vcf.gz"

mv "${OUTPREFIX}_variantsAD.renamed.vcf.gz" "${OUTPREFIX}_variantsAD.vcf.gz"
bcftools index -f "${OUTPREFIX}_variantsAD.vcf.gz"

### ─── STATS ───────────────────────────────────────────────────────────────
bcftools stats "${OUTPREFIX}_fullAD.vcf.gz"     > "${OUTPREFIX}_fullAD_stats.txt"
bcftools stats "${OUTPREFIX}_variantsAD.vcf.gz" > "${OUTPREFIX}_variantsAD_stats.txt"

echo "=== DONE ==="
echo "Full VCF (all sites + AD): ${OUTPREFIX}_fullAD.vcf.gz"
echo "Variant-only VCF:          ${OUTPREFIX}_variantsAD.vcf.gz"
echo "Stats:"
echo "  ${OUTPREFIX}_fullAD_stats.txt"
echo "  ${OUTPREFIX}_variantsAD_stats.txt"

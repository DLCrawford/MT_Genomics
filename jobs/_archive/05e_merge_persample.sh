#!/bin/bash
###############################################################################
# 05e_merge_persample.sh
# Stage 05 RE-RUN v4 (continued): merge per-sample VCFs into a single
# population VCF, replicating the historical merged_144.vcf.gz step.
#
# Why this exists:
#   05d emits one normalized VCF per sample under ${VCF_DIR}/persample/.
#   This script unions them with `bcftools merge -m none` (one row per
#   POS+ALT), reheaders to strip the trailing "_0" from sample names, and
#   writes stats + a self-documenting manifest.
#
# Expected SNP count (target):
#   ~1133 SNPs / ~1140 records / 144 samples in the historical merged_144.vcf.gz.
#   We use 143 BAMs (1_0 excluded), so a slightly lower count is expected.
#
# Pre-flight:
#   Run only AFTER 05d's array has completed and all 143 *_norm.vcf.gz
#   files exist. A simple check is included below.
#
# Submit:
#   bsub < jobs/05e_merge_persample.sh
###############################################################################

#BSUB -J fhet_persample_merge
#BSUB -P fun_gen_1
#BSUB -q normal
#BSUB -n 8
#BSUB -R "rusage[mem=8000M] span[hosts=1]"
#BSUB -W 01:00
#BSUB -o /projectnb/dcrawford/MT_Genomics2/logs/05e_merge_%J.out
#BSUB -e /projectnb/dcrawford/MT_Genomics2/logs/05e_merge_%J.err

set -euo pipefail
cd /projectnb/dcrawford/MT_Genomics2

# === Load environment + paths ===
source /projectnb/dcrawford/MT_Genomics2/jobs/config.sh
conda activate "$CONDA_ENV"

PERSAMPLE_DIR="${VCF_DIR}/persample"
RUN_TAG="persample_merged"
OUT="${VCF_DIR}/Fhet_mt_${RUN_TAG}.vcf.gz"
LIST="${VCF_DIR}/Fhet_mt_${RUN_TAG}_input_list.txt"
MANIFEST="${VCF_DIR}/Fhet_mt_${RUN_TAG}_run_manifest.txt"

# === Pre-flight: every BAM in BAM_LIST has a *_norm.vcf.gz ===
echo "[$(date)] Pre-flight: checking per-sample VCFs exist..."
MISSING=0
while read -r BAM_FILE_NAME; do
    SAMPLE="${BAM_FILE_NAME%_MT.bam}"
    NORM="${PERSAMPLE_DIR}/${SAMPLE}_norm.vcf.gz"
    if [[ ! -s "$NORM" || ! -s "${NORM}.csi" ]]; then
        echo "  MISSING: ${NORM}"
        MISSING=$((MISSING + 1))
    fi
done < "$BAM_LIST"

if [[ $MISSING -gt 0 ]]; then
    echo "ERROR: ${MISSING} per-sample VCFs missing. Re-run 05d for the failed array tasks." >&2
    exit 1
fi
echo "  all $(wc -l < "$BAM_LIST") per-sample VCFs present."

# === Build the merge input list (sorted to match BAM_LIST order) ===
awk -v dir="$PERSAMPLE_DIR" '{ s = $1; sub(/_MT\.bam$/, "", s); print dir "/" s "_norm.vcf.gz" }' \
    "$BAM_LIST" > "$LIST"
N_IN=$(wc -l < "$LIST")
echo "[$(date)] Merging ${N_IN} per-sample VCFs..."

# === Manifest BEFORE the long step ===
{
    echo "=== Stage-05d/05e per-sample-then-merge manifest ==="
    echo "timestamp:        $(date -Iseconds)"
    echo "host:             $(hostname)"
    echo "lsb_jobid:        ${LSB_JOBID:-not_set}"
    echo
    echo "RUN_TAG:          ${RUN_TAG}"
    echo "architecture:     per-sample call (05d) → bcftools merge -m none (05e)"
    echo "rationale:        replicate historical merged_144 recipe; recover the"
    echo "                  ~1133-SNP signal that joint -mv suppresses on mtDNA"
    echo "                  with a divergent reference."
    echo
    echo "REF:              ${REF}"
    echo "BAMS_DIR:         ${BAMS_DIR}"
    echo "BAM_LIST:         ${BAM_LIST}"
    echo "n_inputs:         ${N_IN}"
    echo "PERSAMPLE_DIR:    ${PERSAMPLE_DIR}"
    echo "OUT:              ${OUT}"
    echo
    echo "bcftools version: $(bcftools --version | head -1)"
    echo
    echo "Per-sample command (in 05d):"
    echo "  bcftools mpileup -f \$REF \$BAM -a AD,DP --max-depth 10000 -Ou \\"
    echo "    | bcftools call -mv --ploidy 1 -Oz -o \${SAMPLE}.vcf.gz"
    echo "  bcftools norm -m -any -Oz -o \${SAMPLE}_norm.vcf.gz \${SAMPLE}.vcf.gz"
    echo
    echo "Merge command (this script):"
    echo "  bcftools merge -m none --threads 8 -l \$LIST -Oz -o \$OUT"
} > "$MANIFEST"

# === Merge ===
bcftools merge \
    -m none \
    --threads 8 \
    -l "$LIST" \
    -Oz \
    -o "$OUT"
bcftools index -f "$OUT"

# === Strip trailing "_0" from sample names ===
echo "[$(date)] Reheadering (strip trailing '_0')..."
bcftools reheader \
    -s <(bcftools query -l "$OUT" | sed 's/_0$//') \
    -o "${OUT}.renamed" \
    "$OUT"
mv "${OUT}.renamed" "$OUT"
bcftools index -f "$OUT"

# === Stats ===
STATS="${VCF_DIR}/Fhet_mt_${RUN_TAG}_stats.txt"
bcftools stats "$OUT" > "$STATS"

# === Append summary to manifest ===
N_RECORDS=$(awk -F'\t' '/^SN.*number of records:/{print $NF; exit}' "$STATS")
N_SNPS=$(awk -F'\t' '/^SN.*number of SNPs:/{print $NF; exit}' "$STATS")
N_INDELS=$(awk -F'\t' '/^SN.*number of indels:/{print $NF; exit}' "$STATS")
N_MULTI=$(awk -F'\t' '/^SN.*number of multiallelic sites:/{print $NF; exit}' "$STATS")
TS_TV=$(awk '/^TSTV/{print $5; exit}' "$STATS")
TS_TV_1ST=$(awk '/^TSTV/{print $8; exit}' "$STATS")
N_SAMPLES=$(awk -F'\t' '/^SN.*number of samples:/{print $NF; exit}' "$STATS")

{
    echo
    echo "=== Output summary ==="
    echo "n_samples:                ${N_SAMPLES}"
    echo "n_records:                ${N_RECORDS}"
    echo "n_SNPs:                   ${N_SNPS}"
    echo "n_indels:                 ${N_INDELS}"
    echo "n_multiallelic_sites:     ${N_MULTI}"
    echo "ts/tv (all alts):         ${TS_TV}"
    echo "ts/tv (1st alt only):     ${TS_TV_1ST}"
    echo
    echo "Historical target (merged_144.vcf.gz, 144 samples):"
    echo "  records=1140  SNPs=1133  multiallelic=29  ts/tv=7.92  ts/tv(1st)=9.12"
    echo
    echo "VCF:    ${OUT}"
    echo "Stats:  ${STATS}"
} >> "$MANIFEST"

echo "=== DONE ==="
echo "RUN_TAG:    ${RUN_TAG}"
echo "n_samples:  ${N_SAMPLES}    n_records:  ${N_RECORDS}    n_SNPs: ${N_SNPS}"
echo "ts/tv (all alts): ${TS_TV}    ts/tv (1st alt): ${TS_TV_1ST}"
echo "Historical target: 1133 SNPs across 144 samples (target slightly higher because we use 143)."
echo "Manifest:   ${MANIFEST}"

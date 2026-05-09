#!/bin/bash
###############################################################################
# scripts/run_stage05_core.sh
# Parametrized mpileup+call core for stage 05.
#
# Usage:
#   bash scripts/run_stage05_core.sh <RUN_TAG> <MIN_BQ> <MIN_MQ> <PLOIDY>
#
#   RUN_TAG  - short identifier baked into output filenames
#              (e.g. "v2_Q13_q20_p1")
#   MIN_BQ   - bcftools mpileup -Q (per-base quality threshold)
#   MIN_MQ   - bcftools mpileup -q (mapping quality threshold)
#   PLOIDY   - bcftools call --ploidy (1 for haploid mtDNA)
#
# Writes (under $VCF_DIR):
#   Fhet_mt_${RUN_TAG}_fullAD.vcf.gz       (variant-only; -A keeps all alts)
#   Fhet_mt_${RUN_TAG}_variantsAD.vcf.gz   (alias of fullAD; downstream compat)
#   Fhet_mt_${RUN_TAG}_*_stats.txt
#   Fhet_mt_${RUN_TAG}_run_manifest.txt    (parameters + provenance)
#   + .csi indices on each VCF
#
# Why this exists:
#   The original strict run (-Q 30 -q 30 ploidy=2) produced only 152 SNPs vs
#   ~950 expected from a comparable historical pipeline. To diagnose, we
#   re-run with relaxed filters and the biologically-correct haploid model:
#     v2: -Q 13 -q 20 --ploidy 1   (relaxed; primary candidate)
#     v3: -Q 13 -q  0 --ploidy 1   (no MAPQ filter; matches archived
#                                   per-sample caller for comparison)
#   All three results (v1 strict, v2, v3) are kept side-by-side under
#   distinct RUN_TAGs so the parameter sweep is fully reproducible.
###############################################################################

set -euo pipefail

if [[ $# -ne 4 ]]; then
    echo "Usage: $0 <RUN_TAG> <MIN_BQ> <MIN_MQ> <PLOIDY>" >&2
    echo "Example: $0 v2_Q13_q20_p1 13 20 1" >&2
    exit 2
fi

RUN_TAG="$1"
MIN_BQ="$2"
MIN_MQ="$3"
PLOIDY="$4"

# === Load environment + paths ===
source /projectnb/dcrawford/MT_Genomics2/jobs/config.sh
conda activate "$CONDA_ENV"

mkdir -p "$VCF_DIR"
OUTPREFIX="${VCF_DIR}/Fhet_mt_${RUN_TAG}"
BAM_LIST_FULL="${VCF_DIR}/bam_full_paths.txt"
MANIFEST="${OUTPREFIX}_run_manifest.txt"

# === Build BAM list (idempotent; same content across runs) ===
echo "[$(date)] Preparing BAM file list..."
awk -v dir="$BAMS_DIR" '{print dir "/" $1}' "$BAM_LIST" > "$BAM_LIST_FULL"

N_BAMS=$(wc -l < "$BAM_LIST_FULL")
echo "=== ${N_BAMS} BAM files ==="
head -n 3 "$BAM_LIST_FULL"
echo "  ..."
tail -n 3 "$BAM_LIST_FULL"
echo "================================"

# === Write run manifest BEFORE the long step ===
# Captures all parameters + provenance so the run is fully explainable later,
# even if the stats files get separated from the script.
{
    echo "=== Stage-05 run manifest ==="
    echo "timestamp:        $(date -Iseconds)"
    echo "host:             $(hostname)"
    echo "lsb_jobid:        ${LSB_JOBID:-not_set}"
    echo "lsb_hosts:        ${LSB_HOSTS:-not_set}"
    echo
    echo "RUN_TAG:          ${RUN_TAG}"
    echo "mpileup -Q (BQ):  ${MIN_BQ}"
    echo "mpileup -q (MQ):  ${MIN_MQ}"
    echo "call --ploidy:    ${PLOIDY}"
    echo
    echo "REF:              ${REF}"
    echo "BAMS_DIR:         ${BAMS_DIR}"
    echo "BAM_LIST:         ${BAM_LIST}"
    echo "n_bams:           ${N_BAMS}"
    echo "VCF_DIR:          ${VCF_DIR}"
    echo "OUTPREFIX:        ${OUTPREFIX}"
    echo
    echo "bcftools version: $(bcftools --version | head -1)"
    echo
    echo "Pipeline command:"
    echo "  bcftools mpileup -f \$REF -b \$BAM_LIST_FULL -a AD,DP \\"
    echo "      -Q ${MIN_BQ} -q ${MIN_MQ} -d 100000 -Ou \\"
    echo "    | bcftools call --ploidy ${PLOIDY} -mv -A -Oz \\"
    echo "        -o \${OUTPREFIX}_fullAD.vcf.gz"
} > "$MANIFEST"

# === MPILEUP + CALL ===
# NOTE on flags:
#   -A in `bcftools call`: keeps all alternate alleles at variant sites.
#                          Does NOT emit non-variant sites.
#   -mv:                   multiallelic caller, variants only.
#   --ploidy 1:            haploid model (correct for mtDNA).
# The "_fullAD" filename is historical; output is variant-only because of -mv.
echo "[$(date)] Running bcftools mpileup + call (-Q ${MIN_BQ} -q ${MIN_MQ} --ploidy ${PLOIDY})..."
bcftools mpileup \
    -f "$REF" \
    -b "$BAM_LIST_FULL" \
    -a AD,DP \
    -Q "${MIN_BQ}" -q "${MIN_MQ}" -d 100000 \
    -Ou \
  | bcftools call \
        --ploidy "${PLOIDY}" \
        -mv -A -Oz \
        -o "${OUTPREFIX}_fullAD.vcf.gz"

bcftools index -f "${OUTPREFIX}_fullAD.vcf.gz"

# === RENAME SAMPLES (strip trailing "_0") ===
echo "[$(date)] Renaming samples (strip trailing '_0')..."
bcftools reheader \
    -s <(bcftools query -l "${OUTPREFIX}_fullAD.vcf.gz" | sed 's/_0$//') \
    -o "${OUTPREFIX}_fullAD.renamed.vcf.gz" \
    "${OUTPREFIX}_fullAD.vcf.gz"
mv "${OUTPREFIX}_fullAD.renamed.vcf.gz" "${OUTPREFIX}_fullAD.vcf.gz"
bcftools index -f "${OUTPREFIX}_fullAD.vcf.gz"

# === VARIANT-ONLY ALIAS ===
# fullAD is already variant-only because of -mv; this step is effectively a
# copy with the legacy filename so existing downstream scripts (verify_stage05,
# 06_snpeff, 07_cds_snps_norm) can pick up what they expect via _variantsAD.
bcftools view -v snps,indels -Oz \
    -o "${OUTPREFIX}_variantsAD.vcf.gz" \
    "${OUTPREFIX}_fullAD.vcf.gz"
bcftools index -f "${OUTPREFIX}_variantsAD.vcf.gz"

# === STATS ===
bcftools stats "${OUTPREFIX}_fullAD.vcf.gz"     > "${OUTPREFIX}_fullAD_stats.txt"
bcftools stats "${OUTPREFIX}_variantsAD.vcf.gz" > "${OUTPREFIX}_variantsAD_stats.txt"

# === Append summary to manifest ===
N_RECORDS=$(awk -F'\t' '/^SN.*number of records:/{print $NF; exit}' "${OUTPREFIX}_fullAD_stats.txt")
N_SNPS=$(awk -F'\t' '/^SN.*number of SNPs:/{print $NF; exit}' "${OUTPREFIX}_fullAD_stats.txt")
N_INDELS=$(awk -F'\t' '/^SN.*number of indels:/{print $NF; exit}' "${OUTPREFIX}_fullAD_stats.txt")
N_MULTI=$(awk -F'\t' '/^SN.*number of multiallelic sites:/{print $NF; exit}' "${OUTPREFIX}_fullAD_stats.txt")
TS_TV=$(awk '/^TSTV/{print $5; exit}' "${OUTPREFIX}_fullAD_stats.txt")
TS_TV_1ST=$(awk '/^TSTV/{print $8; exit}' "${OUTPREFIX}_fullAD_stats.txt")

{
    echo
    echo "=== Output summary ==="
    echo "n_records:                ${N_RECORDS}"
    echo "n_SNPs:                   ${N_SNPS}"
    echo "n_indels:                 ${N_INDELS}"
    echo "n_multiallelic_sites:     ${N_MULTI}"
    echo "ts/tv (all alts):         ${TS_TV}"
    echo "ts/tv (1st alt only):     ${TS_TV_1ST}"
    echo
    echo "VCFs:"
    echo "  ${OUTPREFIX}_fullAD.vcf.gz"
    echo "  ${OUTPREFIX}_variantsAD.vcf.gz"
    echo "Stats:"
    echo "  ${OUTPREFIX}_fullAD_stats.txt"
    echo "  ${OUTPREFIX}_variantsAD_stats.txt"
} >> "$MANIFEST"

echo "=== DONE ==="
echo "RUN_TAG:    ${RUN_TAG}"
echo "n_records:  ${N_RECORDS}    n_SNPs: ${N_SNPS}    ts/tv(1st alt): ${TS_TV_1ST}"
echo "Manifest:   ${MANIFEST}"

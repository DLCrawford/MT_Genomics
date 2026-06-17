#!/bin/bash
###############################################################################
# 05f_joint_call.sh
# Stage 05 canonical caller: JOINT mpileup + call over all 143 slim mt-only
# BAMs. One pass. Emits one VCF row per variant position with per-cell DP/AD
# filled in for every (POS × sample) by construction — no .:.:.:. cells, no
# backfill needed.
#
# Replaces 05d2 + 05e2 (per-sample call → bcftools merge), which were the
# workaround for the bcftools-1.6 high-depth-pileup collapse bug documented
# in CHANGELOG 2026-05-15 (sessions 6–10). With bcftools 1.23.1 from
# ${HOME}/software/local/bin (PATH-injected by config.sh), the joint -mv
# pipeline produces the full variant set directly. The per-sample workaround
# is no longer necessary AND it had a structural defect for downstream
# heteroplasmy work: per-sample VCFs only contained variant rows, so REF
# samples appeared as .:.:.:. after merge — exactly the problem stage 07's
# backfill was trying to plug.
#
# Why this is also the right architecture for heteroplasmy:
#   Joint mpileup pileups every sample at every position. At a variant
#   position, every sample gets a DP and an AD = (REF_count, ALT_count).
#   Heteroplasmy at a canonical ALT in any sample shows up as nonzero
#   AD_alt with AD_alt/DP between 0.1 and 0.7 — visible without any further
#   force-calling step. Stage 08 (haplotype @ AD_alt/DP > 0.7) and stage 09
#   (heteroplasmy @ 0.1 ≤ AD_alt/DP < 0.7, AD_alt ≥ 50) both read this
#   file's per-cell AD directly.
#
# Recipe:
#   bcftools mpileup -f REF -b slim_bamlist -a AD,DP --max-depth 100000 -Ou
#     | bcftools call -mv --ploidy 1 -Oz -o Fhet_mt_joint.vcf.gz
#
# Notes on flag choices:
#   - No -A on `bcftools call`. CHANGELOG 2026-05-09 documents that the
#     joint-call experiments (05/05b/05c) had -A and produced "96 %
#     multiallelic phantom-alt soup" — phantom ALTs at every position from
#     keep-all-possible-alleles. Without -A, only ALT alleles with
#     meaningful evidence enter the call, which is what we want for a
#     clean variant set.
#   - No -Q / -q overrides on mpileup. Uses bcftools defaults
#     (-Q 13, -q 0), matching 05d2's choice (see jobs/_archive/
#     05d2_persample_call.sh header line 26).
#   - --ploidy 1 (correct for haploid mtDNA). Does not affect AD/DP — AD
#     is computed from read counts and is independent of ploidy. Only the
#     emitted GT collapses to majority allele, which is biologically
#     correct for a haploid genome.
#   - --max-depth 100000. mt depth on slim BAMs runs 5,000–20,000 ×
#     per sample; default 8,000 would truncate. 100k is well above
#     observed maximum.
#   - --threads 8 on mpileup (matches BSUB -n 8); call is single-threaded.
#
# Inputs:
#   Slim mt-only BAMs in ${PROJECT_ROOT}/MT_only_bams/${SAMPLE}_MT_only.bam
#   (+ .bai), produced by BSUB_Slim_BAM_mt.sh. Run that first if not done.
#   Reference: ${REF} (Fhet_MT.fasta + .fai)
#
# Output:
#   ${VCF_DIR}/Fhet_mt_joint.vcf.gz      (+ .csi)
#   ${VCF_DIR}/Fhet_mt_joint_stats.txt
#   ${VCF_DIR}/Fhet_mt_joint_run_manifest.txt
#   ${VCF_DIR}/slim_bamlist.txt          (full paths, used by mpileup -b)
#
# Submit:
#   bsub < jobs/BSUB_Slim_BAM_mt.sh      # if slim BAMs not yet produced
#   bsub < jobs/05f_joint_call.sh        # this script
#
# Walltime expectation: 30–60 min wallclock for 143 slim BAMs on T2.
###############################################################################

#BSUB -J fhet_joint_call
#BSUB -P fun_gen_1
#BSUB -q normal
#BSUB -n 8
#BSUB -R "rusage[mem=16000M] span[hosts=1]"
#BSUB -W 04:00
#BSUB -o /projectnb/dcrawford/MT_Genomics2/logs/05f_joint_%J.out
#BSUB -e /projectnb/dcrawford/MT_Genomics2/logs/05f_joint_%J.err

set -euo pipefail
cd /projectnb/dcrawford/MT_Genomics2

# === Load environment + paths (sources local-bin export for bcftools 1.23.1) ===
source /projectnb/dcrawford/MT_Genomics2/jobs/config.sh

SLIM_BAMS_DIR="${PROJECT_ROOT}/MT_only_bams"
SLIM_BAMLIST="${VCF_DIR}/slim_bamlist.txt"
RUN_TAG="joint"
OUT="${VCF_DIR}/Fhet_mt_${RUN_TAG}.vcf.gz"
STATS="${VCF_DIR}/Fhet_mt_${RUN_TAG}_stats.txt"
MANIFEST="${VCF_DIR}/Fhet_mt_${RUN_TAG}_run_manifest.txt"

# === Pre-flight ===
echo "[$(date)] Pre-flight checks..."

[[ -f "$REF"      ]] || { echo "ERROR: REF not found: $REF"; exit 1; }
[[ -f "${REF}.fai" ]] || { echo "ERROR: REF FASTA index not found: ${REF}.fai"; exit 1; }
[[ -f "$BAM_LIST" ]] || { echo "ERROR: bam_list not found: $BAM_LIST"; exit 1; }
[[ -d "$SLIM_BAMS_DIR" ]] || {
    echo "ERROR: slim BAM dir not found: $SLIM_BAMS_DIR"
    echo "  → run jobs/BSUB_Slim_BAM_mt.sh first"
    exit 1
}

# Confirm every BAM_LIST entry has its corresponding slim BAM + .bai
echo "[$(date)] Verifying slim BAMs..."
MISSING=0
while read -r BAM_FILE_NAME; do
    SAMPLE="${BAM_FILE_NAME%_MT.bam}"
    SLIM="${SLIM_BAMS_DIR}/${SAMPLE}_MT_only.bam"
    if [[ ! -s "$SLIM" || ! -s "${SLIM}.bai" ]]; then
        echo "  MISSING: ${SLIM} (or .bai)"
        MISSING=$((MISSING + 1))
    fi
done < "$BAM_LIST"

if [[ $MISSING -gt 0 ]]; then
    echo "ERROR: ${MISSING} slim BAMs or indexes missing. Re-run BSUB_Slim_BAM_mt.sh." >&2
    exit 1
fi
echo "  all $(wc -l < "$BAM_LIST") slim BAMs + indexes present."

echo "  bcftools: $(bcftools --version | head -1)"

# === Build slim BAM list (full paths) from BAM_LIST, in BAM_LIST order ===
awk -v dir="$SLIM_BAMS_DIR" '{ s = $1; sub(/_MT\.bam$/, "_MT_only.bam", s); print dir "/" s }' \
    "$BAM_LIST" > "$SLIM_BAMLIST"
N_BAMS=$(wc -l < "$SLIM_BAMLIST")
echo "Slim BAMs to pile up: $N_BAMS"

# === Manifest BEFORE the long step (records intent even if call dies) ===
{
    echo "=== Stage-05f joint-call manifest ==="
    echo "timestamp:        $(date -Iseconds)"
    echo "host:             $(hostname)"
    echo "lsb_jobid:        ${LSB_JOBID:-not_set}"
    echo
    echo "RUN_TAG:          ${RUN_TAG}"
    echo "architecture:     slim mt-only BAMs (BSUB_Slim_BAM_mt.sh) →"
    echo "                  JOINT bcftools mpileup + call -mv --ploidy 1"
    echo "rationale:        single-pass joint call. Variant positions only;"
    echo "                  per-cell DP/AD on every sample by construction"
    echo "                  (no .:.:.:. REF cells, no backfill step needed)."
    echo "                  Replaces 05d2 + 05e2 (per-sample → merge), which"
    echo "                  was a 1.6-bug workaround that left REF cells empty."
    echo
    echo "REF:              ${REF}"
    echo "SLIM_BAMS_DIR:    ${SLIM_BAMS_DIR}"
    echo "BAM_LIST:         ${BAM_LIST}"
    echo "n_inputs:         ${N_BAMS}"
    echo "OUT:              ${OUT}"
    echo
    echo "bcftools version: $(bcftools --version | head -1)"
    echo
    echo "Joint call command (this script):"
    echo "  bcftools mpileup -f \$REF -b \$SLIM_BAMLIST -a AD,DP \\"
    echo "    --max-depth 100000 --threads 8 -Ou \\"
    echo "    | bcftools call -mv --ploidy 1 -Oz -o \$OUT"
} > "$MANIFEST"

# === Joint mpileup → call ===
echo "[$(date)] Joint mpileup + call across ${N_BAMS} slim BAMs..."

bcftools mpileup \
    -f "$REF" \
    -b "$SLIM_BAMLIST" \
    -a AD,DP \
    --max-depth 100000 \
    --threads 8 \
    -Ou \
  | bcftools call \
        -mv \
        --ploidy 1 \
        -Oz \
        -o "$OUT"
bcftools index -f "$OUT"

# === Rename samples (slim BAMs have no @RG SM:, fall back to BAM path) ===
# Convert .../10_0_MT_only.bam → 10_MT, etc.
echo "[$(date)] Renaming samples..."
bcftools reheader \
    -s <(bcftools query -l "$OUT" | sed -e 's|.*/||' -e 's/_0_MT_only\.bam$/_MT/') \
    -o "${OUT}.renamed" \
    "$OUT"
mv "${OUT}.renamed" "$OUT"
bcftools index -f "$OUT"

# === Stats ===
bcftools stats "$OUT" > "$STATS"

# === Baseline comparison: 05f vs per-sample → merge baseline =================
# Joint -mv and per-sample -mv → merge can legitimately produce different
# variant sets at marginal sites. The session-10 baseline
# Fhet_mt_persample_merged.vcf.gz (1128 SNPs, the proven number after
# closing the bcftools 1.6 collapse diagnosis) is the floor: 05f must
# discover every position the per-sample chain found. Sites in 05f that
# are NOT in the baseline are fine (joint can recover singletons that
# per-sample emitted at low confidence). Sites in the baseline that are
# NOT in 05f are a regression — surface them.
BASELINE="${VCF_DIR}/Fhet_mt_persample_merged.vcf.gz"
if [[ -f "$BASELINE" ]]; then
    echo ""
    echo "[$(date)] Comparing 05f to per-sample baseline (${BASELINE})..."

    JOINT_POSITIONS=$(mktemp)
    BASE_POSITIONS=$(mktemp)
    LOST_POSITIONS="${VCF_DIR}/Fhet_mt_${RUN_TAG}_vs_baseline_lost.tsv"
    GAINED_POSITIONS="${VCF_DIR}/Fhet_mt_${RUN_TAG}_vs_baseline_gained.tsv"

    # Compare by unique (CHROM, POS) — splitting status doesn't matter here.
    # We're asking "did 05f see this site at all", not "did it preserve all ALTs".
    bcftools view -H "$OUT"      | awk 'BEGIN{OFS="\t"}{print $1, $2}' | sort -u > "$JOINT_POSITIONS"
    bcftools view -H "$BASELINE" | awk 'BEGIN{OFS="\t"}{print $1, $2}' | sort -u > "$BASE_POSITIONS"

    comm -23 "$BASE_POSITIONS"  "$JOINT_POSITIONS" > "$LOST_POSITIONS"     # baseline \ 05f
    comm -13 "$BASE_POSITIONS"  "$JOINT_POSITIONS" > "$GAINED_POSITIONS"   # 05f \ baseline

    N_BASE=$(wc -l < "$BASE_POSITIONS")
    N_JOINT=$(wc -l < "$JOINT_POSITIONS")
    N_LOST=$(wc -l < "$LOST_POSITIONS")
    N_GAINED=$(wc -l < "$GAINED_POSITIONS")

    echo "  baseline unique positions:   ${N_BASE}"
    echo "  05f unique positions:        ${N_JOINT}"
    echo "  positions in baseline NOT in 05f (lost): ${N_LOST}"
    echo "  positions in 05f NOT in baseline (gained): ${N_GAINED}"

    if [[ "$N_LOST" -gt 0 ]]; then
        echo ""
        echo "============================================================"
        echo "  REGRESSION: 05f is missing ${N_LOST} positions present in"
        echo "  the per-sample baseline. Do NOT flow into stages 06/07."
        echo "  Lost positions written to: ${LOST_POSITIONS}"
        echo "============================================================"
        # Don't exit non-zero here — the VCF is still produced and writing
        # the diff files is useful for diagnosis. But the manifest will
        # record the regression so it's not silently passed downstream.
    fi

    rm -f "$JOINT_POSITIONS" "$BASE_POSITIONS"
else
    echo ""
    echo "[$(date)] WARNING: baseline ${BASELINE} not found — skipping site-set"
    echo "  comparison. (If 05d2/05e2 output isn't on T2, this is expected for"
    echo "  a clean rebuild. Verify the SNP count is in the 1100–1200 range and"
    echo "  proceed only with explicit confirmation.)"
    N_LOST=""
    N_GAINED=""
fi

# === Sanity checks ===
N_RECORDS=$(awk -F'\t' '/^SN.*number of records:/{print $NF; exit}' "$STATS")
N_SNPS=$(awk -F'\t' '/^SN.*number of SNPs:/{print $NF; exit}' "$STATS")
N_INDELS=$(awk -F'\t' '/^SN.*number of indels:/{print $NF; exit}' "$STATS")
N_MULTI=$(awk -F'\t' '/^SN.*number of multiallelic sites:/{print $NF; exit}' "$STATS")
TS_TV=$(awk '/^TSTV/{print $5; exit}' "$STATS")
N_SAMPLES=$(awk -F'\t' '/^SN.*number of samples:/{print $NF; exit}' "$STATS")
N_OUT=$(bcftools view -H "$OUT" | wc -l)
N_MISSING_DP=$(bcftools query -f '[%DP\n]' "$OUT" | awk '$1=="."' | wc -l)
N_CELLS=$(( N_OUT * N_SAMPLES ))

# === Append summary to manifest ===
{
    echo
    echo "=== Output summary ==="
    echo "n_samples:                ${N_SAMPLES}"
    echo "n_records:                ${N_RECORDS}"
    echo "n_SNPs:                   ${N_SNPS}"
    echo "n_indels:                 ${N_INDELS}"
    echo "n_multiallelic_sites:     ${N_MULTI}"
    echo "ts/tv (all alts):         ${TS_TV}"
    echo "total cells:              ${N_CELLS}"
    echo "cells with DP=.:          ${N_MISSING_DP}  (target: 0 — joint call should fill every cell)"
    echo
    echo "=== Baseline comparison (vs Fhet_mt_persample_merged.vcf.gz) ==="
    if [[ -n "${N_LOST}" ]]; then
        echo "baseline unique positions:      ${N_BASE}"
        echo "05f unique positions:           ${N_JOINT}"
        echo "positions lost vs baseline:     ${N_LOST}  (target: 0 — any value > 0 is a regression)"
        echo "positions gained vs baseline:   ${N_GAINED}  (these are fine; joint can find sites per-sample missed)"
        if [[ "$N_LOST" -gt 0 ]]; then
            echo "lost positions written to:      ${LOST_POSITIONS}"
            echo "DO NOT FLOW INTO STAGES 06/07 UNTIL REGRESSION IS DIAGNOSED."
        fi
    else
        echo "(baseline not present on T2 — comparison skipped)"
    fi
    echo
    echo "Historical reference points:"
    echo "  Fhet_mt_fullAD (Oct 2025, joint -mv -A, bcftools 1.22): 1142 records, 152 SNPs"
    echo "  merged_144 (Jul 2025, per-sample then merge, 144 samples): 1140 records, 1133 SNPs"
    echo "  Fhet_mt_persample_merged (session 10, 05d2/05e2): 1128 SNPs / 143 samples / ts/tv 8.17"
    echo
    echo "VCF:    ${OUT}"
    echo "Stats:  ${STATS}"
} >> "$MANIFEST"

echo ""
echo "=== DONE ==="
echo "RUN_TAG:           ${RUN_TAG}"
echo "n_samples:         ${N_SAMPLES}    n_records: ${N_RECORDS}    n_SNPs: ${N_SNPS}"
echo "ts/tv:             ${TS_TV}"
echo "cells with DP=.:   ${N_MISSING_DP}  (target: 0)"
if [[ -n "${N_LOST}" ]]; then
    echo "vs baseline:       lost=${N_LOST}  gained=${N_GAINED}  (lost should be 0)"
fi
echo "VCF:               ${OUT}"
echo "Stats:             ${STATS}"
echo "Manifest:          ${MANIFEST}"
echo ""
echo "rsync to Mac:"
echo "  rsync -avP \\"
echo "    dcrawford@t2.idsc.miami.edu:${OUT} \\"
echo "    dcrawford@t2.idsc.miami.edu:${OUT}.csi \\"
echo "    ~/Projects/MT_Genomics_Cl_Ap2026/MT_Genomics2/vcf/"
echo ""
echo "Next: scripts/06_snpeff_mac.sh (input is now Fhet_mt_joint.vcf.gz)"

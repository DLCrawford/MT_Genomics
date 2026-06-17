#!/bin/bash
###############################################################################
# 05_1_mpileup_merge.sh
# Stage 05 canonical caller: joint bcftools mpileup over all 143 slim mt-only
# BAMs, then call -mv, then norm -m -any (split multi-ALT into one row per
# ALT). One pass; emits one VCF row per (POS, ALT) with per-cell DP/AD on
# every (POS × sample) by construction.
#
# Background:
#   - The sessions 6-10 diagnosis closed with PER-SAMPLE call + merge (05d2 +
#     05e2) on bcftools 1.23.1 producing 1128 SNPs / 143 samples / ts/tv 8.17
#     (Fhet_mt_persample_merged.vcf.gz, May 15). That established the
#     ground-truth variant set for this panel.
#   - The per-sample chain's structural defect: each per-sample VCF only
#     contained that sample's variant rows, so after bcftools merge -m none
#     unioned them, samples that didn't call a variant at a position had no
#     row to contribute and appeared as ".:.:.:." in the merged VCF. The
#     downstream stage 07 backfill was designed to plug that hole but had
#     its own bug (multi-row positions, secondary-ALT loss; see CHANGELOG
#     2026-05-17).
#   - This stage replaces the per-sample chain with a JOINT mpileup that
#     pileups every sample at every position in one pass. Same underlying
#     read evidence, no .:.:.:. cells, no backfill, no targets-file
#     plumbing. norm -m -any is done HERE (matching where 05d2 did the
#     split), so the output is already in the one-row-per-ALT representation
#     stages 06 (SnpEff) and 07 (CDS) want.
#
# Validation gate:
#   Joint -mv and per-sample -mv can legitimately produce different variant
#   sets at marginal positions. This script diffs its output's unique-
#   position set against the May-15 baseline (Fhet_mt_persample_merged.vcf.gz)
#   and surfaces the count of positions lost / gained in the manifest.
#   Positions LOST vs baseline are a regression; positions GAINED are fine
#   (joint can recover singletons per-sample emitted at low confidence).
#
# Flag choices (all match the sessions 6-10 settled-on per-sample recipe):
#   - bcftools 1.23.1 from $HOME/software/local/bin (PATH-injected by
#     jobs/config.sh; the bioconda 1.6 pin on linux-ppc64le silently
#     collapses high-depth mtDNA pileups).
#   - bcftools mpileup defaults for -Q / -q (-Q 13, -q 0; 05d2 used the
#     same). --max-depth 100000 (mt depth runs 5-20k× per sample; default
#     8000 truncates).
#   - bcftools call -mv --ploidy 1.  No -A (the 05/05b/05c joint
#     experiments with -A produced 96% multiallelic "phantom-alt soup";
#     without -A only ALT alleles with meaningful evidence enter the call).
#     --ploidy 1 is correct for haploid mtDNA and does NOT lose
#     heteroplasmic evidence — AD reflects read counts and is independent
#     of ploidy, only GT collapses to majority allele.
#   - bcftools norm -m -any -f REF: split multi-ALT rows into one row per
#     ALT (matches 05d2's per-sample norm), anchor against REF.
#
# Inputs:
#   Slim mt-only BAMs in ${PROJECT_ROOT}/MT_only_bams/${SAMPLE}_MT_only.bam
#   (+ .bai), produced by jobs/BSUB_Slim_BAM_mt.sh.
#   Reference: ${REF} (Fhet_MT.fasta + .fai)
#   Optional: ${VCF_DIR}/Fhet_mt_persample_merged.vcf.gz — the May-15
#     baseline. If present, used for the regression check. If absent, the
#     comparison is skipped with a warning.
#
# Output:
#   ${VCF_DIR}/05_1_Fhet_mt_persample_merged.vcf.gz            (+ .csi)
#   ${VCF_DIR}/05_1_Fhet_mt_persample_merged_stats.txt
#   ${VCF_DIR}/05_1_Fhet_mt_persample_merged_run_manifest.txt
#   ${VCF_DIR}/05_1_Fhet_mt_persample_merged_vs_baseline_lost.tsv
#   ${VCF_DIR}/05_1_Fhet_mt_persample_merged_vs_baseline_gained.tsv
#   ${VCF_DIR}/slim_bamlist.txt                                 (mpileup -b input)
#
# Submit:
#   bsub < jobs/BSUB_Slim_BAM_mt.sh      # if slim BAMs not yet produced
#   bsub < jobs/05_1_mpileup_merge.sh    # this script
#
# Walltime expectation: 30-60 min wallclock for 143 slim BAMs on T2.
###############################################################################

#BSUB -J fhet_mpileup_merge
#BSUB -P fun_gen_1
#BSUB -q normal
#BSUB -n 8
#BSUB -R "rusage[mem=16000M] span[hosts=1]"
#BSUB -W 04:00
#BSUB -o /projectnb/dcrawford/MT_Genomics2/logs/05_1_mpileup_merge_%J.out
#BSUB -e /projectnb/dcrawford/MT_Genomics2/logs/05_1_mpileup_merge_%J.err

set -euo pipefail
cd /projectnb/dcrawford/MT_Genomics2

# === Load environment + paths (sources local-bin export for bcftools 1.23.1) ===
source /projectnb/dcrawford/MT_Genomics2/jobs/config.sh

SLIM_BAMS_DIR="${PROJECT_ROOT}/MT_only_bams"
SLIM_BAMLIST="${VCF_DIR}/slim_bamlist.txt"
RUN_TAG="05_1_Fhet_mt_persample_merged"
OUT="${VCF_DIR}/${RUN_TAG}.vcf.gz"
STATS="${VCF_DIR}/${RUN_TAG}_stats.txt"
MANIFEST="${VCF_DIR}/${RUN_TAG}_run_manifest.txt"

# === Pre-flight ===
echo "[$(date)] Pre-flight checks..."

[[ -f "$REF"       ]] || { echo "ERROR: REF not found: $REF"; exit 1; }
[[ -f "${REF}.fai" ]] || { echo "ERROR: REF FASTA index not found: ${REF}.fai"; exit 1; }
[[ -f "$BAM_LIST"  ]] || { echo "ERROR: bam_list not found: $BAM_LIST"; exit 1; }
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

# === Build slim BAM list (full paths) in BAM_LIST order ===
awk -v dir="$SLIM_BAMS_DIR" '{ s = $1; sub(/_MT\.bam$/, "_MT_only.bam", s); print dir "/" s }' \
    "$BAM_LIST" > "$SLIM_BAMLIST"
N_BAMS=$(wc -l < "$SLIM_BAMLIST")
echo "Slim BAMs to pile up: $N_BAMS"

# === Manifest BEFORE the long step ===
{
    echo "=== Stage 05_1 (mpileup_merge) manifest ==="
    echo "timestamp:        $(date -Iseconds)"
    echo "host:             $(hostname)"
    echo "lsb_jobid:        ${LSB_JOBID:-not_set}"
    echo
    echo "RUN_TAG:          ${RUN_TAG}"
    echo "architecture:     slim mt-only BAMs (BSUB_Slim_BAM_mt.sh) →"
    echo "                  joint bcftools mpileup + call -mv --ploidy 1 →"
    echo "                  bcftools norm -m -any -f REF (split multi-ALT)"
    echo "rationale:        single-pass joint mpileup fills per-cell DP/AD"
    echo "                  on every (POS × sample) at every variant"
    echo "                  position by construction. No .:.:.:. REF cells;"
    echo "                  no backfill required. norm -m -any happens"
    echo "                  here (matching where 05d2 did the per-sample"
    echo "                  split) so SnpEff (stage 06) annotates each ALT"
    echo "                  independently."
    echo
    echo "REF:              ${REF}"
    echo "SLIM_BAMS_DIR:    ${SLIM_BAMS_DIR}"
    echo "BAM_LIST:         ${BAM_LIST}"
    echo "n_inputs:         ${N_BAMS}"
    echo "OUT:              ${OUT}"
    echo
    echo "bcftools version: $(bcftools --version | head -1)"
    echo
    echo "Command:"
    echo "  bcftools mpileup -f \$REF -b \$SLIM_BAMLIST -a AD,DP \\"
    echo "      --max-depth 100000 --threads 8 -Ou \\"
    echo "    | bcftools call -mv --ploidy 1 -Ou \\"
    echo "    | bcftools norm -m -any -f \$REF -Oz -o \$OUT"
} > "$MANIFEST"

# === Joint mpileup → call → norm-split (one pipeline, no intermediate writes) ===
echo "[$(date)] Joint mpileup + call + norm-split across ${N_BAMS} slim BAMs..."

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
        -Ou \
  | bcftools norm \
        -m -any \
        -f "$REF" \
        -Oz \
        -o "$OUT"
bcftools index -f "$OUT"

# === Rename samples (slim BAMs have no @RG SM:, fall back to BAM path) ===
# Convert .../10_0_MT_only.bam → 10_MT
echo "[$(date)] Renaming samples..."
bcftools reheader \
    -s <(bcftools query -l "$OUT" | sed -e 's|.*/||' -e 's/_0_MT_only\.bam$/_MT/') \
    -o "${OUT}.renamed" \
    "$OUT"
mv "${OUT}.renamed" "$OUT"
bcftools index -f "$OUT"

# === Stats ===
bcftools stats "$OUT" > "$STATS"

# === Baseline comparison ===
# May-15 baseline is the per-sample → merge canonical (Fhet_mt_persample_merged.vcf.gz,
# 1128 SNPs / 143 samples / ts/tv 8.17). 05_1 must not lose any of its
# unique positions. Sites gained are fine.
BASELINE="${VCF_DIR}/Fhet_mt_persample_merged.vcf.gz"
LOST_POSITIONS="${VCF_DIR}/${RUN_TAG}_vs_baseline_lost.tsv"
GAINED_POSITIONS="${VCF_DIR}/${RUN_TAG}_vs_baseline_gained.tsv"

if [[ -f "$BASELINE" ]]; then
    echo ""
    echo "[$(date)] Comparing 05_1 to per-sample baseline (${BASELINE})..."

    NEW_POS=$(mktemp)
    BASE_POS=$(mktemp)

    # Compare on unique (CHROM, POS) only. We're asking "did 05_1 see this
    # site at all?", not "did it match every ALT" (norm-split handles ALTs).
    bcftools view -H "$OUT"      | awk 'BEGIN{OFS="\t"}{print $1, $2}' | sort -u > "$NEW_POS"
    bcftools view -H "$BASELINE" | awk 'BEGIN{OFS="\t"}{print $1, $2}' | sort -u > "$BASE_POS"

    comm -23 "$BASE_POS" "$NEW_POS" > "$LOST_POSITIONS"     # in baseline but not in 05_1
    comm -13 "$BASE_POS" "$NEW_POS" > "$GAINED_POSITIONS"   # in 05_1 but not in baseline

    N_BASE=$(wc -l < "$BASE_POS")
    N_NEW=$(wc -l < "$NEW_POS")
    N_LOST=$(wc -l < "$LOST_POSITIONS")
    N_GAINED=$(wc -l < "$GAINED_POSITIONS")

    echo "  baseline unique positions:      ${N_BASE}"
    echo "  05_1 unique positions:          ${N_NEW}"
    echo "  positions lost vs baseline:     ${N_LOST}  (target: 0)"
    echo "  positions gained vs baseline:   ${N_GAINED}  (fine; joint may find singletons per-sample missed)"

    if [[ "$N_LOST" -gt 0 ]]; then
        echo ""
        echo "============================================================"
        echo "  REGRESSION: 05_1 is missing ${N_LOST} positions present"
        echo "  in the per-sample baseline. DO NOT FLOW INTO STAGES 06/07."
        echo "  Lost positions written to: ${LOST_POSITIONS}"
        echo "============================================================"
    fi

    rm -f "$NEW_POS" "$BASE_POS"
else
    echo ""
    echo "[$(date)] WARNING: baseline ${BASELINE} not present — skipping site-set"
    echo "  comparison. Verify the SNP count is in the 1100–1200 range manually."
    N_BASE=""
    N_NEW=""
    N_LOST=""
    N_GAINED=""
fi

# === Sanity numbers ===
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
    echo "n_multiallelic_sites:     ${N_MULTI}  (expected 0 — norm-split done)"
    echo "ts/tv (all alts):         ${TS_TV}"
    echo "total cells:              ${N_CELLS}"
    echo "cells with DP=.:          ${N_MISSING_DP}  (target: 0)"
    echo
    echo "=== Baseline comparison (vs Fhet_mt_persample_merged.vcf.gz) ==="
    if [[ -n "${N_LOST}" ]]; then
        echo "baseline unique positions:      ${N_BASE}"
        echo "05_1 unique positions:          ${N_NEW}"
        echo "positions lost vs baseline:     ${N_LOST}  (target: 0 — REGRESSION if > 0)"
        echo "positions gained vs baseline:   ${N_GAINED}  (fine)"
        if [[ "$N_LOST" -gt 0 ]]; then
            echo "lost positions written to:      ${LOST_POSITIONS}"
            echo "DO NOT FLOW INTO STAGES 06/07 UNTIL DIAGNOSED."
        fi
    else
        echo "(baseline not present on T2 — comparison skipped)"
    fi
    echo
    echo "Reference points:"
    echo "  Fhet_mt_persample_merged (per-sample → merge, May 15, the floor): 1128 SNPs / 143 samples / ts/tv 8.17"
    echo "  merged_144 (Jul 2025, historical):                                1133 SNPs / 144 samples / ts/tv 7.92"
    echo
    echo "VCF:    ${OUT}"
    echo "Stats:  ${STATS}"
} >> "$MANIFEST"

echo ""
echo "=== DONE ==="
echo "RUN_TAG:           ${RUN_TAG}"
echo "n_samples:         ${N_SAMPLES}    n_records: ${N_RECORDS}    n_SNPs: ${N_SNPS}"
echo "ts/tv:             ${TS_TV}"
echo "n_multiallelic:    ${N_MULTI}  (expected 0)"
echo "cells with DP=.:   ${N_MISSING_DP}  (target: 0)"
if [[ -n "${N_LOST}" ]]; then
    echo "vs baseline:       lost=${N_LOST}  gained=${N_GAINED}  (lost MUST be 0)"
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
echo "Next: scripts/06_snpeff_mac.sh (input is ${RUN_TAG}.vcf.gz)"

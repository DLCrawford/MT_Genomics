#!/bin/bash
###############################################################################
# 07b_backfill_AD.sh
# Joint pileup at the canonical SNP positions to fill in per-sample AD/DP
# for samples that the per-sample caller (05d2) emitted as REF (no record).
#
# Why this exists:
#   05d2's `bcftools call -mv` only emits variant rows per sample. After
#   05e2 merge, samples that didn't pass the per-sample variant threshold
#   appear as ".:.:.:." at the merged variant positions — no DP/AD info.
#   Because mtDNA is heteroplasmic, a "REF" sample can still carry ≥10 %
#   ALT reads at canonical positions; the heteroplasmy-aware stage 08 rule
#   needs explicit per-cell AD/DP at every (POS × sample) to detect that.
#
# Recipe:
#   - Take the canonical SNP positions + alleles from Mac-side stage 07
#     (Fhet_MT_CDS.snps.split.vcf.gz, rsync'd to T2 before submission).
#   - Build a tabix-indexed targets TSV (CHROM, POS, REF, ALT) from it.
#   - Joint mpileup over all 143 slim BAMs at exactly those positions
#     (-T <targets>), with -a AD,DP -A.
#   - Force-call the canonical alleles via `bcftools call -m -A -C alleles
#     -T <targets>` so the output has one row per canonical (POS, REF,
#     ALT) and per-sample AD reports REF and canonical-ALT counts —
#     including for samples whose original per-sample call was REF.
#
# Inputs:
#   ${VCF_DIR}/Fhet_MT_CDS.snps.split.vcf.gz       (rsync from Mac stage 07)
#   ${VCF_DIR}/Fhet_MT_CDS.snps.split.vcf.gz.csi
#   ${REF}                                         (Fhet_MT.fasta)
#   slim BAMs in ${PROJECT_ROOT}/MT_only_bams/     (from BSUB_Slim_BAM_mt.sh)
#
# Output:
#   ${VCF_DIR}/Fhet_MT_CDS_backfilled.vcf.gz       (+ .csi)
#   ${VCF_DIR}/Fhet_MT_CDS_backfilled_stats.txt
#   ${VCF_DIR}/Fhet_MT_CDS_targets.tsv.gz          (+ .tbi, the targets file)
#
# Submit (3 steps — runs on Mac → T2 → Mac):
#   # 1) On Mac: rsync the canonical to T2
#   rsync -avP \
#     ~/Projects/MT_Genomics_Cl_Ap2026/MT_Genomics2/vcf/Fhet_MT_CDS.snps.split.vcf.gz \
#     ~/Projects/MT_Genomics_Cl_Ap2026/MT_Genomics2/vcf/Fhet_MT_CDS.snps.split.vcf.gz.csi \
#     dcrawford@t2.idsc.miami.edu:/projectnb/dcrawford/MT_Genomics2/vcf/
#
#   # 2) On T2: submit
#   bsub < jobs/07b_backfill_AD.sh
#
#   # 3) On Mac: rsync the backfilled VCF back
#   rsync -avP \
#     dcrawford@t2.idsc.miami.edu:/projectnb/dcrawford/MT_Genomics2/vcf/Fhet_MT_CDS_backfilled.vcf.gz \
#     dcrawford@t2.idsc.miami.edu:/projectnb/dcrawford/MT_Genomics2/vcf/Fhet_MT_CDS_backfilled.vcf.gz.csi \
#     ~/Projects/MT_Genomics_Cl_Ap2026/MT_Genomics2/vcf/
###############################################################################

#BSUB -J fhet_backfill
#BSUB -P fun_gen_1
#BSUB -q normal
#BSUB -n 4
#BSUB -R "rusage[mem=8000M] span[hosts=1]"
#BSUB -W 01:00
#BSUB -o /projectnb/dcrawford/MT_Genomics2/logs/07b_backfill_%J.out
#BSUB -e /projectnb/dcrawford/MT_Genomics2/logs/07b_backfill_%J.err

set -euo pipefail
cd /projectnb/dcrawford/MT_Genomics2

# === Load environment + paths (sources local-bin export for bcftools 1.23.1) ===
source /projectnb/dcrawford/MT_Genomics2/jobs/config.sh

SLIM_BAMS_DIR="${PROJECT_ROOT}/MT_only_bams"
CANONICAL="${VCF_DIR}/Fhet_MT_CDS.snps.split.vcf.gz"
TARGETS="${VCF_DIR}/Fhet_MT_CDS_targets.tsv.gz"
SLIM_BAMLIST="${VCF_DIR}/slim_bamlist.txt"
BACKFILLED="${VCF_DIR}/Fhet_MT_CDS_backfilled.vcf.gz"

# === Pre-flight ===
echo "[$(date)] Pre-flight checks..."
[[ -f "$CANONICAL" ]] || {
    echo "ERROR: canonical not found: $CANONICAL"
    echo "  → rsync from Mac first (see header)"
    exit 1
}
[[ -f "$REF"       ]] || { echo "ERROR: REF not found: $REF"; exit 1; }
[[ -f "$BAM_LIST"  ]] || { echo "ERROR: bam_list not found: $BAM_LIST"; exit 1; }
[[ -d "$SLIM_BAMS_DIR" ]] || {
    echo "ERROR: slim BAM dir not found: $SLIM_BAMS_DIR"
    echo "  → run jobs/BSUB_Slim_BAM_mt.sh first"
    exit 1
}

echo "  bcftools: $(bcftools --version | head -1)"

# === Build slim BAM list (full paths) from BAM_LIST ===
# BAM_LIST has 143 entries like "10_0_MT.bam"; convert to slim equivalent
# at $SLIM_BAMS_DIR/${SAMPLE}_MT_only.bam.
awk -v dir="$SLIM_BAMS_DIR" '{ s = $1; sub(/_MT\.bam$/, "_MT_only.bam", s); print dir "/" s }' \
    "$BAM_LIST" > "$SLIM_BAMLIST"
N_BAMS=$(wc -l < "$SLIM_BAMLIST")
echo "Slim BAMs to pile up: $N_BAMS"

FIRST_BAM=$(head -1 "$SLIM_BAMLIST")
[[ -s "$FIRST_BAM" ]] || {
    echo "ERROR: first slim BAM missing: $FIRST_BAM"
    echo "  → run jobs/BSUB_Slim_BAM_mt.sh first"
    exit 1
}

# === Build the targets file (CHROM\tPOS\tREF\tALT) and tabix-index it ===
# Used by both bcftools mpileup -T (position restriction) and bcftools call
# -C alleles -T (force the canonical alleles into the output).
echo "[$(date)] Building targets file from canonical..."
TARGETS_PLAIN="${TARGETS%.gz}"
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' "$CANONICAL" > "$TARGETS_PLAIN"
N_POS=$(wc -l < "$TARGETS_PLAIN")
echo "Targets: $N_POS positions"
bgzip -f "$TARGETS_PLAIN"
tabix -f -s1 -b2 -e2 "$TARGETS"

# === Joint mpileup → call -C alleles (force canonical alleles) ===
# Output: one row per canonical position. Per-sample FORMAT carries
# AD = (REF_count, canonical_ALT_count) and DP = total depth, even for
# samples whose original per-sample call was REF. Heteroplasmy at the
# canonical alts surfaces as nonzero ALT_count in those samples' AD.
echo "[$(date)] Joint mpileup over $N_BAMS slim BAMs at $N_POS positions..."

bcftools mpileup \
    -f "$REF" \
    -b "$SLIM_BAMLIST" \
    -T "$TARGETS" \
    -a AD,DP \
    --max-depth 100000 \
    -Q 30 -q 30 \
    -Ou \
  | bcftools call \
        -m -A \
        -C alleles \
        -T "$TARGETS" \
        --ploidy 1 \
        -Oz \
        -o "$BACKFILLED"
bcftools index -f "$BACKFILLED"

# === Rename samples (same recipe as 05e2 / 07) ===
# Slim BAMs have no @RG SM: tag, so sample names came in as full BAM paths.
# Convert .../10_0_MT_only.bam → 10_MT.
echo "[$(date)] Renaming samples..."
bcftools reheader \
    -s <(bcftools query -l "$BACKFILLED" | sed -e 's|.*/||' -e 's/_0_MT_only\.bam$/_MT/') \
    -o "${BACKFILLED}.renamed" \
    "$BACKFILLED"
mv "${BACKFILLED}.renamed" "$BACKFILLED"
bcftools index -f "$BACKFILLED"

# === Stats ===
STATS="${VCF_DIR}/Fhet_MT_CDS_backfilled_stats.txt"
bcftools stats "$BACKFILLED" > "$STATS"

# === Sanity checks ===
N_OUT=$(bcftools view -H "$BACKFILLED" | wc -l)
N_SAMPLES=$(bcftools query -l "$BACKFILLED" | wc -l)
N_MISSING=$(bcftools query -f '[%DP\n]' "$BACKFILLED" | awk '$1=="."' | wc -l)
N_CELLS=$(( N_OUT * N_SAMPLES ))

echo ""
echo "=== DONE ==="
echo "Backfilled VCF : $BACKFILLED"
echo "Stats          : $STATS"
echo "Records        : $N_OUT  (expected: $N_POS)"
echo "Samples        : $N_SAMPLES"
echo "Total cells    : $N_CELLS"
echo "Cells with DP=.: $N_MISSING  (target: 0 — every cell should have a DP)"
echo ""
echo "Sample names (first 3 + last 3):"
bcftools query -l "$BACKFILLED" | head -3
echo "  ..."
bcftools query -l "$BACKFILLED" | tail -3
echo ""
echo "rsync to Mac:"
echo "  rsync -avP \\"
echo "    dcrawford@t2.idsc.miami.edu:${BACKFILLED} \\"
echo "    dcrawford@t2.idsc.miami.edu:${BACKFILLED}.csi \\"
echo "    ~/Projects/MT_Genomics_Cl_Ap2026/MT_Genomics2/vcf/"

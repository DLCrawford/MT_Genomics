#!/usr/bin/env bash
###############################################################################
# scripts/13_pileup_cds_AD.sh
# Stage : Per-CDS-position pileup of DP + AD across all 141 panel BAMs
#         (excluding individuals 70 and 125). Output is a VCF whose AD/DP
#         per cell is derived directly from the pileup at every CDS base,
#         INDEPENDENT of which positions/ALTs were panel-variant called.
#         This lets stage 14 detect ALT-Hp at bases that never reached the
#         panel variant threshold (the "truly private" category that
#         MT_DP_AD_141.txt cannot expose).
# Input : MT_only_bams/${ID}_0_MT_only.bam        (141 of 143; 70, 125 excluded)
#         Missing_Files/SSM_MT_ref/Fhet_MT.fasta  (+ .fai)
#         docs/mito_protein_coding.bed            (13 CDS intervals)
# Output: vcf/slim_bamlist_141.txt                (BAM paths feeding mpileup)
#         vcf/pileup_cds_141.vcf.gz               (per-cell DP/AD at every CDS pos)
#         vcf/pileup_cds_141.vcf.gz.tbi
# Run   : conda activate SNP_env   (bcftools)
#         bash scripts/13_pileup_cds_AD.sh
# Notes : --no-BAQ matches the canonical caller (05_1) so AD numbers here are
#         on the same scale as the canonical AD field. -d 100000 is the
#         mtDNA-depth ceiling used elsewhere in the pipeline.
###############################################################################
set -euo pipefail
cd "$(dirname "$0")/.."

REF=Missing_Files/SSM_MT_ref/Fhet_MT.fasta
BED=docs/mito_protein_coding.bed
BAMDIR=MT_only_bams
OUT_BAMLIST=vcf/slim_bamlist_141.txt
OUT_VCF=vcf/pileup_cds_141.vcf.gz

# Sanity
command -v bcftools >/dev/null || { echo "ERROR: bcftools not on PATH (conda activate SNP_env?)"; exit 1; }
[[ -f "$REF"     ]] || { echo "ERROR: missing reference $REF"; exit 1; }
[[ -f "$REF.fai" ]] || { echo "ERROR: missing $REF.fai (run: samtools faidx $REF)"; exit 1; }
[[ -f "$BED"     ]] || { echo "ERROR: missing CDS BED $BED"; exit 1; }
[[ -d "$BAMDIR"  ]] || { echo "ERROR: missing BAM dir $BAMDIR"; exit 1; }

# 1) Build 141-BAM list, excluding 70_0 and 125_0
echo "[13] Building 141-sample BAM list (excluding 70, 125)..."
ls "$BAMDIR"/*_MT_only.bam \
    | grep -Ev "/(70|125)_0_MT_only\.bam$" \
    > "$OUT_BAMLIST"
N=$(wc -l < "$OUT_BAMLIST")
echo "      $N BAMs in $OUT_BAMLIST"
[[ "$N" -eq 141 ]] || { echo "WARN: expected 141 BAMs, got $N"; }

# 2) bcftools mpileup at every CDS position, per-sample DP + AD
echo "[13] Running bcftools mpileup over $(awk '{s+=$3-$2} END{print s}' $BED) CDS positions x $N samples..."
echo "      (a few minutes; --no-BAQ matches the canonical caller)"
bcftools mpileup \
    -f "$REF" \
    -b "$OUT_BAMLIST" \
    -R "$BED" \
    -a AD,DP \
    -d 100000 \
    --no-BAQ \
    -Oz -o "$OUT_VCF"

bcftools index -f -t "$OUT_VCF"

echo "[13] DONE."
echo "      $OUT_VCF"
echo "      Next: python scripts/14_hp_from_pileup.py"

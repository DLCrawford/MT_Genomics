#!/usr/bin/env bash
# make_141_vcf.sh
# Create 141-sample canonical VCF from the 143-sample CDS canonical by
# removing samples 70 and 125 only.
#
# Position exclusions (invariant / fixed-ALT sites) are NOT applied here —
# those sites contribute π = 0 and are handled at the analysis level.
# The 9 invariant sites (6 SYN + 3 NS fixed reference-divergence) remain
# in the VCF so the record count matches haplotype_matrixDP10.csv (927 rows).
#
# Input:  vcf/Fhet_MT_CDS.snps.split.vcf.gz  (143 samples, 927 records)
# Output: vcf/141_MT_variants.vcf.gz          (141 samples, 927 records)
#                                              bgzipped + tabix indexed
#
# Run from anywhere:
#   conda activate SNP_env
#   bash scripts/make_141_vcf.sh

set -euo pipefail

# Resolve paths relative to the script's own directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"

IN="$PROJECT_DIR/vcf/Fhet_MT_CDS.snps.split.vcf.gz"
OUT="$PROJECT_DIR/vcf/141_MT_variants.vcf.gz"

echo "Input:  $IN"
echo "Output: $OUT"
echo "Excluding samples: 70_MT, 125_MT"

bcftools view \
    --samples-file <(bcftools query -l "$IN" | grep -vE "^(70_MT|125_MT)$") \
    "$IN" \
  | bcftools view -e 'AC=0' \
    -Oz -o "$OUT"

bcftools index --tbi "$OUT"

echo ""
echo "=== Done ==="
bcftools stats "$OUT" | grep -E "^SN" | grep -E "samples|SNPs|records"
echo "Output: $OUT  ($(ls -lh "$OUT" | awk '{print $5}'))"

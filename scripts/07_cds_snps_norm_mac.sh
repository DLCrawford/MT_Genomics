#!/usr/bin/env bash
###############################################################################
# scripts/07_cds_snps_norm_mac.sh
# Stage : CDS restriction + SNP filter + split multiallelic → canonical output
# Input : vcf/Fhet_mt_variantsAD_ann.vcf.gz   (from scripts/06_snpeff_mac.sh)
#         Missing_Files/SSM_MT_ref/Fhet_MT.gff (CDS coordinates)
#         Missing_Files/SSM_MT_ref/Fhet_MT.fasta (for bcftools norm -f)
# Output: vcf/Fhet_MT_CDS.snps.split.vcf.gz   ← CANONICAL (frozen once produced)
#         vcf/MT_CDS.regions.gz + .tbi
#         vcf/Fhet_MT_CDS.snps.split_stats.txt
# Run   : bash scripts/07_cds_snps_norm_mac.sh
#         (from the project root: ~/Projects/MT_Genomics_Cl_Ap2026/MT_Genomics2/)
# Needs : bcftools (brew), bgzip + tabix (come with bcftools/htslib on Mac)
###############################################################################

set -euo pipefail

### ─── PATHS ──────────────────────────────────────────────────────────────────
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"

REF="${PROJECT_ROOT}/Missing_Files/SSM_MT_ref/Fhet_MT.fasta"
GFF="${PROJECT_ROOT}/Missing_Files/SSM_MT_ref/Fhet_MT.gff"
VCF_DIR="${PROJECT_ROOT}/vcf"

INPUT="${VCF_DIR}/Fhet_mt_variantsAD_ann.vcf.gz"
CDS_REGIONS="${VCF_DIR}/MT_CDS.regions"
CANONICAL="${VCF_DIR}/Fhet_MT_CDS.snps.split.vcf.gz"

### ─── PRE-FLIGHT ─────────────────────────────────────────────────────────────
echo "[$(date)] Pre-flight checks..."
[[ -f "$INPUT" ]] || { echo "ERROR: stage-06 annotated VCF not found: $INPUT"; exit 1; }
[[ -f "$GFF"   ]] || { echo "ERROR: GFF not found: $GFF"; exit 1; }
[[ -f "$REF"   ]] || { echo "ERROR: REF FASTA not found: $REF"; exit 1; }

# bcftools norm -f requires a samtools FASTA index (.fai); create if absent.
if [[ ! -f "${REF}.fai" ]]; then
    echo "[$(date)] Creating FASTA index (samtools faidx)..."
    samtools faidx "$REF"
fi

### ─── DERIVE CDS REGIONS FROM GFF ───────────────────────────────────────────
# bcftools -R expects tab-delimited CHROM / FROM / TO (1-based inclusive).
# GFF is already 1-based inclusive, so no coordinate shift needed.
echo "[$(date)] Extracting CDS regions from GFF..."

awk 'BEGIN{OFS="\t"} !/^#/ && $3 == "CDS" {print $1, $4, $5}' "$GFF" \
    | sort -k1,1 -k2,2n \
    > "$CDS_REGIONS"

echo "CDS intervals: $(wc -l < "$CDS_REGIONS")"
cat "$CDS_REGIONS"    # MT genome is tiny — safe to print all rows

bgzip -f "$CDS_REGIONS"
tabix -s1 -b2 -e3 "${CDS_REGIONS}.gz"

### ─── CDS RESTRICT → SNPs ONLY → SPLIT MULTIALLELIC ────────────────────────
echo "[$(date)] CDS restrict + SNPs only + bcftools norm..."

bcftools view \
    -R "${CDS_REGIONS}.gz" \
    "$INPUT" \
  | bcftools view \
        -v snps \
  | bcftools norm \
        -m -any \
        -f "$REF" \
        -Oz -o "$CANONICAL"

bcftools index -f "$CANONICAL"

### ─── STATS + SUMMARY ────────────────────────────────────────────────────────
STATS="${VCF_DIR}/Fhet_MT_CDS.snps.split_stats.txt"
bcftools stats "$CANONICAL" > "$STATS"

echo ""
echo "=== DONE — CANONICAL OUTPUT PRODUCED ==="
echo "  $CANONICAL"
echo ""
echo "=== Summary counts ==="
grep "^SN" "$STATS"
echo ""
echo "Next: Mac-side Python haplotype parsing"

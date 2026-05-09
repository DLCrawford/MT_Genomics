#!/usr/bin/env bash
###############################################################################
# scripts/06_snpeff_mac.sh
# Stage : SnpEff annotation — runs on Mac (ppc64le makes Triton 2 unusable)
# Input : vcf/Fhet_mt_variantsAD.vcf.gz        (rsync'd from Triton 2 after stage 05)
# Output: vcf/Fhet_mt_variantsAD_ann.vcf.gz
#         vcf/snpEff_summary.html
#         vcf/snpEff_genes.txt
# Run   : bash scripts/06_snpeff_mac.sh
#         (from the project root: ~/Projects/MT_Genomics_Cl_Ap2026/MT_Genomics2/)
###############################################################################
#
# PRE-REQUISITE — download stage-05 VCF from Triton 2 once bjobs shows it done:
#
#   rsync -avP \
#     dcrawford@scc1.bu.edu:/projectnb/dcrawford/MT_Genomics2/vcf/Fhet_mt_variantsAD.vcf.gz \
#     ~/Projects/MT_Genomics_Cl_Ap2026/MT_Genomics2/vcf/
#   rsync -avP \
#     dcrawford@scc1.bu.edu:/projectnb/dcrawford/MT_Genomics2/vcf/Fhet_mt_variantsAD.vcf.gz.csi \
#     ~/Projects/MT_Genomics_Cl_Ap2026/MT_Genomics2/vcf/
#
# SNPEFF DATABASE NOTES:
#   Database : Fhet_MT  (custom-built from ~/snpEff/data/Fhet_MT/genes.gff +
#                        sequences.fa; chromosome = NC_012312.1)
#   Config   : ~/snpEff/snpEff.config
#              Entry:  Fhet_MT.genome : Fundulus heteroclitus mitochondrion
#   Built with: java -jar snpEff.jar build -gff3 -v -noCheckCds -noCheckProtein \
#                  -c ~/snpEff/snpEff.config Fhet_MT
#   (-noCheckCds -noCheckProtein required due to GFF Parent= tag issues)
#
# The VCF chromosome (NC_012312.1) matches the database sequences — no rename needed.
#
###############################################################################

set -euo pipefail

### ─── PATHS ──────────────────────────────────────────────────────────────────
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"

# SnpEff install: custom database lives under ~/snpEff/
# JAR is inside the micromamba SNP_env (snpeff 5.2).
# If you ever reinstall, update SNPEFF_JAR to match the new path — check with:
#   find ~/micromamba/envs/SNP_env -name "snpEff.jar" 2>/dev/null
SNPEFF_DIR="${HOME}/snpEff"
SNPEFF_JAR="${HOME}/micromamba/envs/SNP_env/share/snpeff-5.2-1/snpEff.jar"
SNPEFF_CFG="${SNPEFF_DIR}/snpEff.config"
DB_NAME="Fhet_MT"           # custom database built from genes.gff + sequences.fa

VCF_DIR="${PROJECT_ROOT}/vcf"
INPUT="${VCF_DIR}/Fhet_mt_variantsAD.vcf.gz"
OUTPUT="${VCF_DIR}/Fhet_mt_variantsAD_ann.vcf.gz"
STATS_HTML="${VCF_DIR}/snpEff_summary.html"

mkdir -p "$VCF_DIR"

### ─── PRE-FLIGHT ─────────────────────────────────────────────────────────────
echo "[$(date)] Pre-flight checks..."

[[ -f "$SNPEFF_JAR" ]] \
    || { echo "ERROR: snpEff.jar not found at $SNPEFF_JAR"
         echo "  → Check: find ~/micromamba/envs/SNP_env -name 'snpEff.jar'"
         exit 1; }

[[ -f "$SNPEFF_CFG" ]] \
    || { echo "ERROR: snpEff.config not found at $SNPEFF_CFG"; exit 1; }

[[ -f "${SNPEFF_DIR}/data/${DB_NAME}/snpEffectPredictor.bin" ]] \
    || { echo "ERROR: database '${DB_NAME}' not built — snpEffectPredictor.bin missing"
         echo "  Expected: ${SNPEFF_DIR}/data/${DB_NAME}/snpEffectPredictor.bin"
         echo "  Rebuild:  cd ~/snpEff && java -Xmx4g -jar \$SNPEFF_JAR build -gff3 -v"
         echo "              -noCheckCds -noCheckProtein -c snpEff.config Fhet_MT"
         exit 1; }

[[ -f "$INPUT" ]] \
    || { echo "ERROR: Input VCF not found: $INPUT"
         echo "  → rsync from Triton 2 first (see script header)"; exit 1; }

echo "Java version:"
java -version 2>&1 | head -1

echo "VCF contigs (must contain NC_012312.1):"
bcftools view -h "$INPUT" | grep "^##contig"

echo "Sample count:"
bcftools query -l "$INPUT" | wc -l

### ─── ANNOTATE ───────────────────────────────────────────────────────────────
echo "[$(date)] Running SnpEff annotation (db=${DB_NAME})..."

# Run from SNPEFF_DIR so snpEff resolves relative data/ paths correctly,
# then write outputs to absolute VCF_DIR paths.
cd "$SNPEFF_DIR"

java -Xmx4g -jar "$SNPEFF_JAR" ann \
    -c "$SNPEFF_CFG" \
    -v \
    -stats "$STATS_HTML" \
    "$DB_NAME" \
    "$INPUT" \
  | bcftools view -Oz -o "$OUTPUT"

cd "$PROJECT_ROOT"

bcftools index -f "$OUTPUT"

### ─── QUICK SANITY CHECK ────────────────────────────────────────────────────
echo ""
echo "=== Annotated variant count ==="
bcftools stats "$OUTPUT" | grep "^SN"

echo ""
echo "=== Variants WITH ANN field ==="
bcftools view "$OUTPUT" | grep -v "^#" | grep -c "ANN=" || true

echo ""
echo "=== DONE ==="
echo "Annotated VCF : $OUTPUT"
echo "SnpEff stats  : $STATS_HTML"
echo "Gene counts   : ${STATS_HTML%.html}.genes.txt"
echo ""
echo "Next: bash scripts/07_cds_snps_norm_mac.sh"

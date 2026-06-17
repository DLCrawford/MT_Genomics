#!/usr/bin/env bash
###############################################################################
# scripts/06_snpeff_mac.sh
# Stage : SnpEff annotation — runs on Mac (ppc64le makes Triton 2 unusable for SnpEff)
# Input : vcf/05_1_Fhet_mt_persample_merged.vcf.gz   (rsync'd from T2 after 05_1)
# Output: vcf/05_1_Fhet_mt_persample_merged_ann.vcf.gz
#         vcf/snpEff_summary.html
#         vcf/snpEff_genes.txt
# Run   : bash scripts/06_snpeff_mac.sh
#         (from the project root: ~/Projects/MT_Genomics_Cl_Ap2026/MT_Genomics2/)
###############################################################################
#
# PRE-REQUISITE — download 05_1 output from Triton 2 once 05_1_mpileup_merge.sh
# finishes (AND its baseline-comparison check shows "positions lost = 0"):
 ### DLC I believe it is "05_1_Fhet_mt_persample_merged.vcf.gz"
#
#   rsync -avP \
#     dcrawford@t2.idsc.miami.edu:/projectnb/dcrawford/MT_Genomics2/vcf/05_1_Fhet_mt_persample_merged.vcf.gz \
#     dcrawford@t2.idsc.miami.edu:/projectnb/dcrawford/MT_Genomics2/vcf/05_1_Fhet_mt_persample_merged.vcf.gz.csi \
#     ~/Projects/MT_Genomics_Cl_Ap2026/MT_Genomics2/vcf/
#
# SNPEFF DATABASE NOTES (updated 2026-05-15):
#   Database : NC_012312.1  (NCBI-built via SnpEff's buildDbNcbi.sh)
#   Built with: cd ~/snpEff && ./scripts/buildDbNcbi.sh NC_012312.1
#               (downloads NCBI GenBank record + builds snpEffectPredictor.bin
#                + sequence.bin under ~/snpEff/data/NC_012312.1/)
#   Config   : ~/snpEff/snpEff.config — required entry:
#               NC_012312.1.codonTable : Vertebrate_Mitochondrial
#              (without this, SnpEff translates with the standard nuclear code
#               and mis-calls coding-effect annotations on mtDNA.)
#
#   Note: the older Fhet_MT custom build (~/snpEff/data/Fhet_MT/) has only the
#   raw genes.gff + sequences.fa; snpEffectPredictor.bin was never produced.
#   Do not point this script at Fhet_MT — it'll fail at the pre-flight check.
#
# The VCF chromosome (NC_012312.1) matches the database — no rename needed.
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
SNPEFF_JAR="${SNPEFF_DIR}/snpEff.jar"      # standalone install, version-matched to the data/ dir
SNPEFF_CFG="${SNPEFF_DIR}/snpEff.config"
DB_NAME="NC_012312.1"       # NCBI-built via ~/snpEff/scripts/buildDbNcbi.sh (SnpEff 5.4)

# JAR / DB version compatibility:
#   Use the snpEff.jar that lives in ~/snpEff/ (NOT the conda env's
#   snpeff-5.2-1/snpEff.jar). buildDbNcbi.sh builds a 5.4-format database;
#   the conda env's 5.2 JAR will refuse it with:
#     Database version : '5.4'  Program version : '5.2'  Compatible : '[5.2, 5.0, 5.1]'
#   The standalone install under ~/snpEff/ is the one that built the database
#   and the one that can read it.

# SnpEff 5.2+ is built against Java 11+ (class file 55). The macOS system java
# is 1.8 (class file 52) and produces UnsupportedClassVersionError. Use the
# Java that ships inside the SNP_env conda env — modern (currently 23.0.2),
# runs both the 5.2 and 5.4 SnpEff JARs.
JAVA_BIN="${HOME}/micromamba/envs/SNP_env/lib/jvm/bin/java"

VCF_DIR="${PROJECT_ROOT}/vcf"
INPUT="${VCF_DIR}/05_1_Fhet_mt_persample_merged.vcf.gz"
OUTPUT="${VCF_DIR}/05_1_Fhet_mt_persample_merged_ann.vcf.gz"
STATS_HTML="${VCF_DIR}/snpEff_summary.html"

mkdir -p "$VCF_DIR"

### ─── PRE-FLIGHT ─────────────────────────────────────────────────────────────
echo "[$(date)] Pre-flight checks..."

[[ -f "$SNPEFF_JAR" ]] \
    || { echo "ERROR: snpEff.jar not found at $SNPEFF_JAR"
         echo "  → Expected the standalone install at ~/snpEff/snpEff.jar"
         echo "    (NOT the conda env's snpeff-5.2-1/snpEff.jar — version mismatch with the 5.4 DB.)"
         exit 1; }

[[ -x "$JAVA_BIN" ]] \
    || { echo "ERROR: env Java not found at $JAVA_BIN"
         echo "  → Check: find ~/micromamba/envs/SNP_env -name 'java' -type f"
         echo "  (system java is 1.8 and won't run SnpEff 5.2 — must use env Java.)"
         exit 1; }

[[ -f "$SNPEFF_CFG" ]] \
    || { echo "ERROR: snpEff.config not found at $SNPEFF_CFG"; exit 1; }

[[ -f "${SNPEFF_DIR}/data/${DB_NAME}/snpEffectPredictor.bin" ]] \
    || { echo "ERROR: database '${DB_NAME}' not built — snpEffectPredictor.bin missing"
         echo "  Expected: ${SNPEFF_DIR}/data/${DB_NAME}/snpEffectPredictor.bin"
         echo "  Rebuild:  cd ~/snpEff && ./scripts/buildDbNcbi.sh ${DB_NAME}"
         echo "  And confirm ~/snpEff/snpEff.config contains:"
         echo "    ${DB_NAME}.codonTable : Vertebrate_Mitochondrial"
         exit 1; }

[[ -f "$INPUT" ]] \
    || { echo "ERROR: Input VCF not found: $INPUT"
         echo "  → rsync from Triton 2 first (see script header)"; exit 1; }

echo "Java version (env Java, NOT system java):"
"$JAVA_BIN" -version 2>&1 | head -1

echo "VCF contigs (must contain NC_012312.1):"
bcftools view -h "$INPUT" | grep "^##contig"

echo "Sample count:"
bcftools query -l "$INPUT" | wc -l

### ─── ANNOTATE ───────────────────────────────────────────────────────────────
echo "[$(date)] Running SnpEff annotation (db=${DB_NAME})..."

# Run from SNPEFF_DIR so snpEff resolves relative data/ paths correctly,
# then write outputs to absolute VCF_DIR paths.
cd "$SNPEFF_DIR"

"$JAVA_BIN" -Xmx4g -jar "$SNPEFF_JAR" ann \
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

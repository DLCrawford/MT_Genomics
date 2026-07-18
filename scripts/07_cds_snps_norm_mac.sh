#!/usr/bin/env bash
###############################################################################
# scripts/07_cds_snps_norm_mac.sh
# Stage : CDS restriction + SNP filter + sample rename → canonical output.
#         Runs on Mac after stage 06 (SnpEff annotation).
# Input : vcf/05_1_Fhet_mt_persample_merged_ann.vcf.gz   (from stage 06)
#         docs/mito_protein_coding.bed                   (13 mt protein-coding genes)
# Output: vcf/Fhet_MT_CDS.snps.split.vcf.gz              ← CANONICAL (FROZEN
#                                                          once validated),
#                                                          with per-cell AD/DP + ANN
#         vcf/Fhet_MT_CDS.snps.split_stats.txt
# Run   : bash scripts/07_cds_snps_norm_mac.sh
# Needs : bcftools (SNP_env conda env)
###############################################################################
#
# Why this stage is short (architecture as of 2026-05-18, session 14):
#   - Stage 05 (jobs/05_1_mpileup_merge.sh) already does the joint mpileup +
#     call + `bcftools norm -m -any -f REF` (split multi-ALT into one row
#     per ALT, anchored to REF). The 05_1 output is therefore already in
#     the one-row-per-ALT representation that stages 06 and 07 expect.
#   - Stage 06 (scripts/06_snpeff_mac.sh) annotates each split row
#     independently — one ANN per (POS, ALT) — which is what we want for
#     per-ALT effect prediction.
#   - Stage 07 therefore only needs to: restrict to CDS via BED, filter to
#     SNPs only, and clean up sample names. No norm-split here; no
#     backfill; no targets-file machinery. The split multiallelic
#     representation from 05_1 + ANN from 06 both pass through transparently.
#
# Bugs the new architecture sidesteps (see CHANGELOG 2026-05-17/18):
#   - The .:.:.:. REF cells from the per-sample → merge architecture are
#     gone (05_1's joint mpileup fills DP/AD on every cell).
#   - The targets-file-format bug that lost 17 secondary ALTs at multi-row
#     positions in the session-11/12 backfill is structurally impossible
#     now (no `call -C alleles -T` step against a post-norm VCF anywhere
#     in the pipeline).
###############################################################################

set -euo pipefail

### ─── PATHS ──────────────────────────────────────────────────────────────────
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"

BED="${PROJECT_ROOT}/docs/mito_protein_coding.bed"     # 13 mt protein-coding genes
VCF_DIR="${PROJECT_ROOT}/vcf"

INPUT="${VCF_DIR}/05_1_Fhet_mt_persample_merged_ann.vcf.gz"
CANONICAL="${VCF_DIR}/Fhet_MT_CDS.snps.split.vcf.gz"
STATS="${VCF_DIR}/Fhet_MT_CDS.snps.split_stats.txt"

### ─── PRE-FLIGHT ─────────────────────────────────────────────────────────────
echo "[$(date)] Pre-flight checks..."
[[ -f "$INPUT" ]] || { echo "ERROR: stage-06 annotated VCF not found: $INPUT"; exit 1; }
[[ -f "$BED"   ]] || { echo "ERROR: CDS BED not found: $BED"; exit 1; }

echo "bcftools: $(bcftools --version | head -1)"
echo "Input:    $INPUT"
echo "CDS BED:"
cat "$BED"

### ─── PIPELINE ─────────────────────────────────────────────────────────────
# 1) CDS restrict via BED
# 2) SNPs only
# 3) reheader to strip BAM-path artifacts from sample names
#
# norm -m -any was already done in 05_1 (joint mpileup + call + norm split,
# all in one pipeline). The input here is already in one-row-per-ALT form.
# ANN (SnpEff effect prediction) from stage 06 lives in the INFO field and
# passes through bcftools view transparently — no annotate step needed.
echo "[$(date)] CDS restrict + SNPs only..."

bcftools view \
    -R "$BED" \
    "$INPUT" \
  | bcftools view \
        -v snps \
        -Oz -o "${CANONICAL}.tmp"
bcftools index -f "${CANONICAL}.tmp"

### ─── SAMPLE RENAME ────────────────────────────────────────────────────────
# Idempotent: 05_1's own reheader step should have already cleaned the
# names (.../10_0_MT_only.bam → 10_MT). If they're already clean this
# sed is a no-op; if not, it fixes them here.
echo "[$(date)] Renaming samples (idempotent)..."
bcftools reheader \
    -s <(bcftools query -l "${CANONICAL}.tmp" | sed -e 's|.*/||' -e 's/_0_MT_only\.bam$/_MT/') \
    -o "$CANONICAL" \
    "${CANONICAL}.tmp"
bcftools index -f "$CANONICAL"
rm -f "${CANONICAL}.tmp" "${CANONICAL}.tmp.csi"

### ─── STATS + VERIFICATION ───────────────────────────────────────────────────
bcftools stats "$CANONICAL" > "$STATS"

# Sanity numbers
N_OUT=$(bcftools view -H "$CANONICAL" | wc -l | tr -d ' ')
N_SAMPLES=$(bcftools query -l "$CANONICAL" | wc -l | tr -d ' ')
N_MISSING_DP=$(bcftools query -f '[%DP\n]' "$CANONICAL" | awk '$1=="."' | wc -l | tr -d ' ')
N_CELLS=$(( N_OUT * N_SAMPLES ))
N_WITH_ANN=$(bcftools view "$CANONICAL" | grep -v "^#" | grep -c "ANN=" || true)

echo ""
echo "=== Canonical output ==="
echo "  Canonical:    $CANONICAL"
echo "  Stats:        $STATS"
echo ""
echo "Records:         $N_OUT"
echo "Samples:         $N_SAMPLES"
echo "Total cells:     $N_CELLS"
echo "Cells with DP=.: $N_MISSING_DP  (target: 0 — must be 0 to freeze)"
echo "Records w/ ANN:  $N_WITH_ANN   (must equal record count to freeze)"
echo ""
echo "FREEZE GATE: do NOT mark canonical frozen until both 'cells with DP=.'"
echo "and 'records w/ ANN'  pass (= 0 and = record count, respectively),"
echo "AND 05_1's manifest showed 'positions lost vs baseline: 0'."
echo ""
echo "Sample names (first 5 + last 3):"
bcftools query -l "$CANONICAL" | head -5
echo "  ..."
bcftools query -l "$CANONICAL" | tail -3
echo ""
echo "=== Summary counts ==="
grep "^SN" "$STATS"
echo ""
echo "Next: python scripts/08_call_haplotypes.py"

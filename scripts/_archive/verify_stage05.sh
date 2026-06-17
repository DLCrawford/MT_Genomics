#!/usr/bin/env bash
###############################################################################
# scripts/verify_stage05.sh
# Purpose : Sanity-check the stage-05 VCF on Triton 2 before rsyncing to Mac.
# Run     : bash scripts/verify_stage05.sh
#           (from project root: /projectnb/dcrawford/MT_Genomics2/)
# Expects : vcf/Fhet_mt_variantsAD.vcf.gz  (+  .csi index)
###############################################################################

set -euo pipefail

source jobs/config.sh      # loads CONDA_MODULE, VCF_DIR, etc.
module load "$CONDA_MODULE"
eval "$(conda shell.bash hook)"
conda activate "$CONDA_ENV"

VCF="${VCF_DIR}/Fhet_mt_variantsAD.vcf.gz"
PASS=0
FAIL=0

check_pass() { echo "  [PASS] $1"; ((PASS++)) || true; }
check_fail() { echo "  [FAIL] $1"; ((FAIL++)) || true; }

echo "================================================="
echo " Stage-05 VCF verification"
echo " $(date)"
echo " VCF: $VCF"
echo "================================================="
echo ""

### 1. File exists and is non-empty
echo "--- File existence ---"
if [[ -f "$VCF" ]]; then
    SIZE=$(du -sh "$VCF" | cut -f1)
    check_pass "VCF exists (size: $SIZE)"
else
    check_fail "VCF not found: $VCF"
    echo "Stage 05 may still be running. Exiting."
    exit 1
fi

### 2. Index exists
echo ""
echo "--- Index ---"
if [[ -f "${VCF}.csi" ]]; then
    check_pass "CSI index exists"
else
    echo "  [WARN] .csi index missing — creating it now..."
    bcftools index -f "$VCF"
    check_pass "CSI index created"
fi

### 3. VCF is not truncated (bcftools can read it without error)
echo ""
echo "--- File integrity ---"
if bcftools view "$VCF" > /dev/null 2>&1; then
    check_pass "bcftools can read VCF without errors"
else
    check_fail "bcftools reports errors reading VCF — file may be truncated"
fi

### 4. Sample count = 143
echo ""
echo "--- Sample count ---"
N_SAMPLES=$(bcftools query -l "$VCF" | wc -l | tr -d ' ')
echo "  Sample count: ${N_SAMPLES}"
if [[ "$N_SAMPLES" -eq 143 ]]; then
    check_pass "143 samples present"
else
    check_fail "Expected 143 samples, found ${N_SAMPLES}"
fi

### 5. Chromosome = NC_012312.1
echo ""
echo "--- Chromosome name ---"
CHROMS=$(bcftools view -h "$VCF" | grep "^##contig" | grep -o 'ID=[^,>]*' | cut -d= -f2 | tr '\n' ' ')
echo "  Contigs: ${CHROMS}"
if echo "$CHROMS" | grep -q "NC_012312.1"; then
    check_pass "Chromosome NC_012312.1 present (matches SnpEff Fhet_MT database)"
else
    check_fail "NC_012312.1 not found in contig headers — SnpEff annotation will fail"
fi

### 6. AD and DP FORMAT fields present
echo ""
echo "--- FORMAT fields ---"
if bcftools view -h "$VCF" | grep -q '##FORMAT=<ID=AD,'; then
    check_pass "FORMAT/AD field present"
else
    check_fail "FORMAT/AD field MISSING — haplotype caller needs AD"
fi

if bcftools view -h "$VCF" | grep -q '##FORMAT=<ID=DP,'; then
    check_pass "FORMAT/DP field present"
else
    check_fail "FORMAT/DP field MISSING — haplotype caller needs DP"
fi

### 7. SNP count in a plausible range (50–5000 for 16 kb MT genome, 143 samples)
echo ""
echo "--- Variant counts ---"
bcftools stats "$VCF" | grep "^SN"
N_SNPS=$(bcftools stats "$VCF" | grep "^SN.*number of SNPs" | awk '{print $NF}')
echo "  SNPs: ${N_SNPS}"
if [[ "$N_SNPS" -gt 50 && "$N_SNPS" -lt 5000 ]]; then
    check_pass "SNP count ${N_SNPS} is in the expected range (50–5000)"
elif [[ "$N_SNPS" -eq 0 ]]; then
    check_fail "0 SNPs called — mpileup/call likely failed or produced only indels"
else
    echo "  [WARN] SNP count ${N_SNPS} is outside the expected range — review manually"
fi

### 8. Spot-check a data line (GT and AD both populated)
echo ""
echo "--- Spot check (first 2 variant lines) ---"
bcftools view "$VCF" | grep -v "^#" | head -2

### 9. Summary
echo ""
echo "================================================="
echo " RESULT: ${PASS} passed, ${FAIL} failed"
echo "================================================="

if [[ "$FAIL" -eq 0 ]]; then
    echo ""
    echo "All checks passed. Ready to rsync to Mac:"
    echo ""
    echo "  rsync -avP \\"
    echo "    dcrawford@scc1.bu.edu:${VCF} \\"
    echo "    ~/Projects/MT_Genomics_Cl_Ap2026/MT_Genomics2/vcf/"
    echo ""
    echo "  rsync -avP \\"
    echo "    dcrawford@scc1.bu.edu:${VCF}.csi \\"
    echo "    ~/Projects/MT_Genomics_Cl_Ap2026/MT_Genomics2/vcf/"
    echo ""
    echo "Then on Mac:"
    echo "  bash scripts/06_snpeff_mac.sh"
else
    echo ""
    echo "${FAIL} check(s) failed — review output above before rsyncing."
    exit 1
fi

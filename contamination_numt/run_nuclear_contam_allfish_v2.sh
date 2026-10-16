#!/usr/bin/env bash
# All-fish nuclear contamination check using an existing multi-sample nuclear VCF with FORMAT/AD.
# Heteroplasmy-style metric (non-major allele at 0.10<=AD/DP<0.35), each fish scored vs the cohort,
# with a GROUP CONTRAST: focals (suspected) vs private_mt / single_NS (genuine heteroplasmy) vs cohort.
# Run:  bsub -Is -n 4 -P fun_gen_1 /bin/bash run_nuclear_contam_allfish.sh
set -euo pipefail

### ─── EDIT THESE ──────────────────────────────────────────────────────────
VCF=/projectnb/dcrawford/SSM_WGS/ALLFISH_NUCLEAR.vcf.gz   # <-- set to your all-fish nuclear VCF (bgzipped)
SCRIPTDIR=/projectnb/dcrawford/MT_Genomics2/NUMT_nuclear  # holds contam_allfish.py, focals.txt, groups.txt
FOCALS=$SCRIPTDIR/focals.txt
GROUPSFILE=$SCRIPTDIR/groups_hp.tsv   # hp vs no_hp contrast (53 vs 88)
REGION=NC_046361.1        # chr1. Set REGION="" for the whole genome (slower). Add more: NC_046361.1,NC_046362.1
PASS_ONLY=yes             # GATK hard-filtered VCF: keep only FILTER=PASS sites (recommended)
OUT=$SCRIPTDIR/nuclear_contam_allfish_out
### ────────────────────────────────────────────────────────────────────────

source /projectnb/dcrawford/MT_Genomics2/jobs/config.sh 2>/dev/null || true
command -v bcftools >/dev/null || { echo "ERROR: bcftools not on PATH"; exit 1; }
command -v python3  >/dev/null || { echo "ERROR: python3 not on PATH"; exit 1; }
mkdir -p "$OUT"

# 1) VCF sanity: AD present, and indexed if a region is requested
bcftools view -h "$VCF" | grep -q 'ID=AD' || { echo "ERROR: no FORMAT/AD declared in $VCF"; exit 1; }
if [ -n "$REGION" ] && ! { [ -f "$VCF.tbi" ] || [ -f "$VCF.csi" ]; }; then
  echo "ERROR: $VCF not indexed (needed for -r). Run: bcftools index -t $VCF"; exit 1; fi

# 2) sample names — must match focals.txt / groups.txt IDs
bcftools query -l "$VCF" > "$OUT/samples.txt"
echo "samples in VCF: $(wc -l < "$OUT/samples.txt")  (first few: $(head -3 "$OUT/samples.txt" | tr '\n' ' '))"
miss=0
while read -r f; do [ -z "$f" ] && continue
  grep -qxF "$f" "$OUT/samples.txt" || { echo "  focal '$f' NOT in VCF sample names"; miss=1; }
done < "$FOCALS"
[ "$miss" = 0 ] || { echo "ERROR: fix focals.txt (and groups.txt) IDs to match $OUT/samples.txt"; exit 1; }

# 3) extract per-sample AD
RG=(); [ -n "$REGION" ] && RG=(-r "$REGION")
FILT=(); [ "$PASS_ONLY" = yes ] && FILT=(-i 'FILTER="PASS"')
echo "extracting AD ${REGION:+over $REGION} ..."
bcftools query "${RG[@]}" "${FILT[@]}" -f '%CHROM\t%POS\t%REF\t%ALT[\t%AD]\n' "$VCF" > "$OUT/ad.tsv"
echo "sites: $(wc -l < "$OUT/ad.tsv")"

# 4) metric + group contrast
python3 "$SCRIPTDIR/contam_allfish.py" "$OUT/ad.tsv" "$OUT/samples.txt" "$FOCALS" "$OUT" "$GROUPSFILE"
echo "DONE -> $OUT/{per_sample_contam_allfish.tsv, verdict.txt}"

#!/bin/bash
###############################################################################
# BSUB_Slim_BAM_mt.sh
# Pre-extract mitochondrial reads from each WGS BAM into a small mt-only BAM.
#
# Why this exists:
#   The input BAMs in /projectnb/dcrawford/MT_Genomics2/bams/ are 15-19 GB
#   each because stage 04 retained unmapped reads. Downstream
#   `bcftools mpileup` walks the entire BAM, which is unnecessary when we
#   only call mt variants. Extracting just the NC_012312.1-mapped reads cuts
#   each BAM to ~50-100 MB, so per-sample calling (05d2) and any future
#   re-runs are minutes instead of hours.
#
# Output:
#   /projectnb/dcrawford/MT_Genomics2/MT_only_bams/${SAMPLE}_MT_only.bam (+ .bai)
#
# Submit:
#   bsub < jobs/BSUB_Slim_BAM_mt.sh
#   # then, AFTER the array completes successfully:
#   bsub < jobs/05d2_persample_call.sh
###############################################################################

#Job Name and array call. Square brackets = array size, % = max number of jobs to run at a time.
#BSUB -J "fhet_slim[1-143]%30"

#Project
#BSUB -P fun_gen_1

#Queue
#BSUB -q normal

#Cores
#BSUB -n 4

#RAM per node/host (use M and G for Mb and Gb). span[hosts=1] forces all cores on the same host.
#BSUB -R "rusage[mem=4000M] span[hosts=1]"

#Walltime (HH:MM)
#BSUB -W 01:00

#Output File (Absolute path!)
#BSUB -o /projectnb/dcrawford/MT_Genomics2/logs/Slim_%J_%I.out
#BSUB -e /projectnb/dcrawford/MT_Genomics2/logs/Slim_%J_%I.err

###################################################################
#SETUP

#LOAD ENVIRONMENTS for Triton 2
# Use locally-built samtools 1.23.1 (NOT the mito_genomics env's pinned 1.6).
export PATH="$HOME/software/local/bin:$PATH"
export LD_LIBRARY_PATH="$HOME/software/local/lib:${LD_LIBRARY_PATH:-}"

#CHANGE INTO THE APPROPRIATE DIRECTORY
cd /projectnb/dcrawford/MT_Genomics2

# pipefail stops on the first error so you don't silently continue with bad data
set -euo pipefail

### === VARIABLES ===
BAMS_DIR=/projectnb/dcrawford/MT_Genomics2/bams
BAM_LIST=/projectnb/dcrawford/MT_Genomics2/bam_list.txt
OUT_DIR=/projectnb/dcrawford/MT_Genomics2/MT_only_bams
MT_CHROM=NC_012312.1
mkdir -p "$OUT_DIR"

### === RESOLVE SAMPLE FOR THIS ARRAY TASK ===
# BAM_LIST has 143 entries like "10_0_MT.bam" (1_0 excluded — no paired trimmed reads).
BAM_FILE_NAME=$(sed -n "${LSB_JOBINDEX}p" "$BAM_LIST")
if [[ -z "$BAM_FILE_NAME" ]]; then
    echo "ERROR: empty BAM filename at index ${LSB_JOBINDEX} of ${BAM_LIST}" >&2
    exit 1
fi
SAMPLE="${BAM_FILE_NAME%_MT.bam}"
IN_BAM="${BAMS_DIR}/${BAM_FILE_NAME}"
OUT_BAM="${OUT_DIR}/${SAMPLE}_MT_only.bam"

if [[ ! -s "$IN_BAM" ]]; then
    echo "ERROR: input BAM not found or empty: ${IN_BAM}" >&2
    exit 1
fi
if [[ ! -s "${IN_BAM}.bai" ]]; then
    echo "ERROR: BAM index not found: ${IN_BAM}.bai (samtools view needs it for region seek)" >&2
    exit 1
fi

echo "[$(date)] task=${LSB_JOBINDEX}  sample=${SAMPLE}"
echo "  in:       ${IN_BAM}"
echo "  out:      ${OUT_BAM}"
echo "  region:   ${MT_CHROM}"
echo "  samtools: $(samtools --version | head -1)"

### === EXTRACT MT-MAPPED READS ===
# -F 2308 = drop unmapped (4) + secondary (256) + supplementary (2048) reads,
#           leaving primary mt-mapped alignments only.
# Region argument uses the BAI to seek directly to NC_012312.1 — fast even
# on a 19 GB input.
samtools view \
    -@ 4 \
    -b \
    -F 2308 \
    -o "$OUT_BAM" \
    "$IN_BAM" \
    "$MT_CHROM"

### === INDEX EXTRACTED BAM ===
samtools index -@ 4 "$OUT_BAM"

echo "[$(date)] DONE  sample=${SAMPLE}  size=$(du -h "$OUT_BAM" | cut -f1)"

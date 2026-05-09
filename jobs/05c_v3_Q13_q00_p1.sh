#!/bin/bash
###############################################################################
# 05c_v3_Q13_q00_p1.sh
# Stage 05 RE-RUN v3: no mapping-quality filter + haploid model.
#   -Q 13       (per-base quality, default)
#   -q 0        (no MAPQ filter)
#   --ploidy 1  (correct for mtDNA)
#
# Why:
#   This is the apples-to-apples comparison vs the archived per-sample
#   caller (archive/Previous_jobs/BSUB_1_MT_SNPcalls.sh), which used no
#   -Q/-q filters. If v2 (-q 20) still under-counts but v3 (-q 0) hits
#   ~950 SNPs, the gap is localized to MAPQ filtering. If both v2 and v3
#   undershoot ~950, the gap is upstream — most likely BAM content from
#   stage 04 (e.g., differences in duplicate handling or read filtering
#   vs the historical pipeline).
#
# Outputs use the RUN_TAG "v3_Q13_q00_p1" so the strict canonical run's
# output files (Fhet_mt_fullAD.vcf.gz et al.) remain untouched on disk.
#
# Submit: bsub < jobs/05c_v3_Q13_q00_p1.sh
###############################################################################

#BSUB -J fhet_mpileup_v3
#BSUB -P fun_gen_1
#BSUB -q normal
#BSUB -n 8
#BSUB -R "rusage[mem=16000M] span[hosts=1]"
#BSUB -W 96:00
#BSUB -o /projectnb/dcrawford/MT_Genomics2/logs/05c_v3_%J.out
#BSUB -e /projectnb/dcrawford/MT_Genomics2/logs/05c_v3_%J.err

set -euo pipefail
cd /projectnb/dcrawford/MT_Genomics2

bash scripts/run_stage05_core.sh "v3_Q13_q00_p1" 13 0 1

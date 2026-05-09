#!/bin/bash
###############################################################################
# 05b_v2_Q13_q20_p1.sh
# Stage 05 RE-RUN v2: relaxed read filters + haploid model.
#   -Q 13       (per-base quality, default)
#   -q 20       (mapping quality, moderate)
#   --ploidy 1  (correct for mtDNA)
#
# Why:
#   The strict canonical run (jobs/05_bcftools_mpileup_call_AD.sh, params
#   -Q 30 -q 30 ploidy=2) produced 152 SNPs vs ~950 expected. Logs show a
#   clean exit with no truncation, OOM, or timeout (max mem 371 MB of
#   16 GB; runtime 9h of 72h). The bottleneck is therefore in the calling
#   parameters. This run tests whether relaxed read filters and the
#   biologically-correct haploid model recover the missing variants.
#
# Outputs use the RUN_TAG "v2_Q13_q20_p1" so the strict canonical run's
# output files (Fhet_mt_fullAD.vcf.gz et al.) remain untouched on disk.
#
# Submit: bsub < jobs/05b_v2_Q13_q20_p1.sh
###############################################################################

#BSUB -J fhet_mpileup_v2
#BSUB -P fun_gen_1
#BSUB -q normal
#BSUB -n 8
#BSUB -R "rusage[mem=16000M] span[hosts=1]"
#BSUB -W 96:00
#BSUB -o /projectnb/dcrawford/MT_Genomics2/logs/05b_v2_%J.out
#BSUB -e /projectnb/dcrawford/MT_Genomics2/logs/05b_v2_%J.err

set -euo pipefail
cd /projectnb/dcrawford/MT_Genomics2

bash scripts/run_stage05_core.sh "v2_Q13_q20_p1" 13 20 1

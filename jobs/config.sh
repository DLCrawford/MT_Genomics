###############################################################################
# config.sh — single source of truth for pipeline paths
#
# Every BSUB script in this directory does:
#     source /projectnb/dcrawford/MT_Genomics2/jobs/config.sh
#
# To relocate data, change the paths below — DO NOT edit individual scripts.
###############################################################################

# ── Project root (mirrored from your Mac) ────────────────────────────────────
PROJECT_ROOT=/projectnb/dcrawford/MT_Genomics2
LOGS_DIR="${PROJECT_ROOT}/logs"

###############################################################################
# DATA + REFERENCE LOCATIONS
#
# The active block (uncommented) is the LEGACY layout: data lives outside
# MT_Genomics2/ in /projectnb/dcrawford/SSM_WGS and /projectnb/dcrawford/SSM_Mito.
#
# To consolidate everything under MT_Genomics2/:
#   1. Either physically `mv` the directories listed below into the
#      paths shown in the CONSOLIDATED block, or symlink them in place,
#      e.g.:
#        cd /projectnb/dcrawford/MT_Genomics2
#        mkdir -p data refs bams vcf
#        ln -s /projectnb/dcrawford/SSM_WGS/fhet_raw_seq    data/raw
#        ln -s /projectnb/dcrawford/SSM_WGS/TrimA_seq       data/trimmed
#        ln -s /projectnb/dcrawford/SSM_WGS/fastqc_out      data/qc_raw
#        ln -s /projectnb/dcrawford/SSM_WGS/TrimA_fastqc_out data/qc_trimmed
#        ln -s /projectnb/dcrawford/SSM_Mito/Fh_MT_ref      refs
#        ln -s /projectnb/dcrawford/SSM_Mito/MT_bam_sam     bams
#        ln -s /projectnb/dcrawford/SSM_Mito/new_hap_AD     vcf
#        cp /projectnb/dcrawford/SSM_WGS/SSM_WGS_list.txt   .
#   2. Comment the LEGACY block, uncomment the CONSOLIDATED block.
###############################################################################

# ── LEGACY layout (currently active) ─────────────────────────────────────────
SAMPLE_LIST=/projectnb/dcrawford/SSM_WGS/SSM_WGS_list.txt
RAW_DIR=/projectnb/dcrawford/SSM_WGS/fhet_raw_seq
TRIM_DIR=/projectnb/dcrawford/SSM_WGS/TrimA_seq
RAW_QC_DIR=/projectnb/dcrawford/SSM_WGS/fastqc_out
TRIM_QC_DIR=/projectnb/dcrawford/SSM_WGS/TrimA_fastqc_out

REF=/projectnb/dcrawford/SSM_Mito/Fh_MT_ref/Fhet_MT.fasta
BAMS_DIR=/projectnb/dcrawford/SSM_Mito/MT_bam_sam
VCF_DIR=/projectnb/dcrawford/SSM_Mito/new_hap_AD
BAM_LIST=/projectnb/dcrawford/SSM_Mito/new_hap_AD2/bam_list.txt

# ── CONSOLIDATED layout (commented; see header for setup) ────────────────────
# SAMPLE_LIST="${PROJECT_ROOT}/SSM_WGS_list.txt"
# RAW_DIR="${PROJECT_ROOT}/data/raw"
# TRIM_DIR="${PROJECT_ROOT}/data/trimmed"
# RAW_QC_DIR="${PROJECT_ROOT}/data/qc_raw"
# TRIM_QC_DIR="${PROJECT_ROOT}/data/qc_trimmed"
#
# REF="${PROJECT_ROOT}/refs/Fhet_MT.fasta"
# BAMS_DIR="${PROJECT_ROOT}/bams"
# VCF_DIR="${PROJECT_ROOT}/vcf"
# BAM_LIST="${PROJECT_ROOT}/bam_list.txt"

# ── External tools (Trimmomatic jar lives outside the project) ───────────────
TRIMJAR=/home/dcrawford/software/local/Trimmomatic-0.39/trimmomatic-0.39.jar
ADAPTERS=/home/dcrawford/software/local/Trimmomatic-0.39/adapters/CombinedAdapters.fa
# If you install trimmomatic via conda instead, comment the two lines above and
# call `trimmomatic PE ...` directly in 02_trim_pe.sh.

# ── Conda env (override if you split tools across envs) ──────────────────────
CONDA_ENV=mito_genomics

# ── Make sure logs directory exists (idempotent) ─────────────────────────────
mkdir -p "$LOGS_DIR"

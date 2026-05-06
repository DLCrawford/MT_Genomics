# CHANGELOG

Session-by-session record of what changed.

## 2026-05-06 — Restructure + clean BSUB job set

- Migrated working layout on Mac to `~/Projects/MT_Genomics_Cl_Ap2026/MT_Genomics2/`, which mirrors `/projectnb/dcrawford/MT_Genomics2/` on Triton 2 one-to-one.
- Moved historical content (`Previous_jobs/`, `Notes_dlcs/`, `Claude_info/`, `CL_Inst/`, `CL_successes/`) under `~/Projects/MT_Genomics_Cl_Ap2026/archive/` (Mac-only, not synced).
- Standardized 5 BSUB scripts under `jobs/` using the comment style of `BSUB_bwa_samtools1.sh` as the template:
  - `01_fastqc_raw.sh` (was `BSUB_fastqc_multiqc.sh`)
  - `02_trim_pe.sh` (was `BSUB_TrimA_1c.sh`)
  - `03_fastqc_trimmed.sh` — **new**, fills the gap referenced as `BSUB_TrimA_fastqc_multiqc.sh`
  - `04_bwa_align_mt.sh` (was `BSUB_1_mt_align_pipping_array.sh`)
  - `05_bcftools_mpileup_call_AD.sh` (was `BSUB_mpile_AD_full.sh`)
- Dropped `BSUB_1_MT_SNPcalls.sh` (per-sample caller, superseded by joint AD caller).
- All scripts now use:
  - Consistent header documentation block (stage / inputs / outputs / submit)
  - `set -euo pipefail`
  - Single `module load anaconda3 && conda activate mito_genomics`
  - Logs at `/projectnb/dcrawford/MT_Genomics2/logs/`
- Wrote `docs/00_setup.md` (HPC + GitHub + conda env + reference indexing) and `docs/01_pipeline.md` (step-by-step replication).
- Added `jobs/config.sh` as the single source of truth for path variables. All 5 scripts now `source` it instead of redefining paths locally. Defaults to the legacy layout (`/projectnb/dcrawford/SSM_WGS/...`, `/projectnb/dcrawford/SSM_Mito/...`); a commented "consolidated under MT_Genomics2/" alternative is included with symlink recipe — flip when ready.

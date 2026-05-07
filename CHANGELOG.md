# CHANGELOG

Session-by-session record of what changed.

## 2026-05-07 — Inventory, recovery from Mac backup, partly-consolidated config

- Discovered `jobs/config.sh` LEGACY paths no longer match Triton 2 reality. Wrote `docs/inventory_2026-05-07.txt` enumerating what actually exists. Findings:
  - `SSM_WGS/fhet_raw_seq` (raw fastqs) — gone, moved to long-term cold storage. Stages 01 + 02 not runnable as-is.
  - `SSM_WGS/SSM_WGS_list.txt` — does not exist.
  - `SSM_WGS/TrimA_seq` — real path is `SSM_WGS/trimmed_seq` (2.1 TB).
  - `SSM_WGS/TrimA_fastqc_out` — real path is `SSM_WGS/fastqc_TrimA` (220 MB).
  - `SSM_Mito/Fh_MT_ref/` — directory exists but listed empty; under IT investigation.
  - `SSM_Mito/MT_bam_sam/` — missing entirely.
  - `SSM_Mito/new_hap_AD/` — missing entirely.
  - Canonical `Fhet_MT_CDS.snps.split.vcf.gz` — not visible at top level of either tree; under IT investigation.
- Recovered from Mac backup `~/Projects/SSM_Mito_All/SSM_MT_ref/` into `Missing_Files/`:
  - `WGS_list.txt` — deduped to 144 unique sample IDs (was 287 with duplicates), trailing newline added.
  - `SSM_MT_ref/` — `Fhet_MT.fasta` + 6 bwa/samtools indexes + `Fhet_MT.gff`.
- Staged to Triton 2:
  - `Missing_Files/WGS_list.txt` → `/projectnb/dcrawford/MT_Genomics2/SSM_WGS_list.txt`.
  - `Missing_Files/SSM_MT_ref/*` → `/projectnb/dcrawford/MT_Genomics2/refs/`.
  - Created `bams/`, `vcf/`, `logs/` under project root (idempotent).
- Updated `jobs/config.sh` to a partly-consolidated layout:
  - `SAMPLE_LIST`, `REF`, `BAMS_DIR`, `VCF_DIR`, `BAM_LIST` consolidated under `${PROJECT_ROOT}/`.
  - `TRIM_DIR` left at `/projectnb/dcrawford/SSM_WGS/trimmed_seq` (regeneratable, too large to move).
  - `TRIM_QC_DIR` left at `/projectnb/dcrawford/SSM_WGS/fastqc_TrimA`.
  - `RAW_DIR` and `RAW_QC_DIR` retained as sentinels with comments — stages 01/02 will fail at the existence check with a clear path.
  - Added `GFF` variable pointing at `${PROJECT_ROOT}/refs/Fhet_MT.gff` for upcoming stage 06 (SnpEff).
- Wrote `docs/pipeline_files.txt` — annotated I/O reference for every input/output the pipeline touches, with status per item and target paths under the consolidated layout.
- Pending: IT investigation outcome on Triton 2 file disappearances. If unrecoverable, plan is to rerun stages 04 → 05, then write + run new stages 06 (SnpEff) + 07 (CDS+SNPs+norm) to produce a fresh canonical `Fhet_MT_CDS.snps.split.vcf.gz`. Stage 04 array bound `[1-144]%30` matches the deduped sample count.
- Future: drop samples `1_0` (forward read only), `70_0`, `125_0` (no phenotypes). Brings the count to 141, matching CLAUDE.md's target. Stage 04 array bound becomes `[1-141]%30` after that.
- **Trimmed-naming reconciliation**: pre-flight on Triton 2 revealed `trimmed_seq/` contains Trim Galore output (`{sample}_1_val_1.fq.gz`, `{sample}_2_val_2.fq.gz`, plus `*_trimming_report.txt`), NOT the Trimmomatic `_p` / `_up` naming that `02_trim_pe.sh` produces. 572 files = 143 paired samples × 4 files; sample `1_0` is absent (forward-only — already on the drop list). Patched `04_bwa_align_mt.sh` and CLAUDE.md sample-naming section to match the on-disk Trim Galore convention. `02_trim_pe.sh` left as-is (raw is in cold storage, can't be rerun); flag this for reconciliation if stage 02 is ever revisited.
- **Pilot scope**: stage 04 pilot uses sample-list indices 2–6 (`10_0`, `102_0`, `103_0`, `104_0`, `105_0`) — index 1 is `1_0` which has no paired trimmed reads.

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

## 2026-05-06 — GitHub auth + initial Mac↔Triton 2 sync

- Generated SSH keys for GitHub on both hosts:
  - **Triton 2**: `~/.ssh/mito_gen_key` (ed25519), registered on GitHub as `mito_gen_key`.
  - **Mac**: `~/.ssh/id_ed25519`, registered on GitHub as `MAC_mito_gen_key`.
- Hardened Triton 2's `~/.ssh/config` with `IdentitiesOnly yes` so ssh only offers `mito_gen_key` to github.com — avoids GitHub rate-limiting when default keys would otherwise be tried first. `docs/00_setup.md §1` updated to match.
- The Mac project directory was created by the earlier restructure but `git init` had never been run. Brought it under version control by cloning `DLCrawford/MT_Genomics` into `/tmp/mt_clone`, moving its `.git/` into the Mac project, then committing the restructured tree as a single commit on top of existing history. Pushed to `origin/main`.
- Triton 2 had stale tracked modifications (`.gitignore`, `CHANGELOG.md`, `CLAUDE.md`, `README.md`) and stale untracked `docs/` + `jobs/` blocking `git pull`. Resolved with `git reset --hard HEAD` + `git clean -fd -- docs jobs` + `git pull origin main`. Triton 2 now matches GitHub.
- Added Mac-side setup section to `docs/00_setup.md` covering SSH key generation, `~/.ssh/config`, the `git clone → move .git` pattern for bringing an existing repo into a pre-populated directory, and zsh `setopt interactive_comments`.
- Queued for next session: physical `mv` of `SSM_WGS/*` and `SSM_Mito/*` into `MT_Genomics2/{data,refs,bams,vcf}`, then flip `jobs/config.sh` LEGACY → CONSOLIDATED, then add stages 06 (SnpEff) + 07 (CDS/SNPs/norm) to produce the canonical `Fhet_MT_CDS.snps.split.vcf.gz`.

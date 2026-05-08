# CHANGELOG

Session-by-session record of what changed.

## 2026-05-08 (session 4) — conda activation fix confirmed and finalised

- **Root cause confirmed:** `conda info` on login node revealed:
  - base: `anaconda3-2023.09` at `/sw/summit/software/linux-power9le/anaconda3-2023.09-0-...` (read only)
  - env: `/projectnb/triton/home/dcrawford/.conda/envs/mito_genomics` (user-level)
  - `module load anaconda3` (unversioned) was likely resolving ambiguously or not registering shell hooks on compute nodes.
- **`jobs/config.sh`** updated: added `CONDA_MODULE=anaconda3/2023.09-0-none-none-oawyzwj`; config.sh now runs `module load "$CONDA_MODULE"` + `eval "$(conda shell.bash hook)"` directly, so every BSUB script that sources config.sh gets correct conda activation automatically.
- **All 5 BSUB scripts** (`01`–`05`): removed `module load anaconda3` (now in config.sh); reverted `source activate` back to `conda activate "$CONDA_ENV"` (correct for modern conda once hooks are registered).
- **Next:** `git pull origin main && bsub < jobs/05_bcftools_mpileup_call_AD.sh` on Triton 2.

## 2026-05-08 (session 3) — haplotype calling design settled; scripts/08 written

- **Haplotype calling design finalised:**
  - Calling rule: `AD_alt / DP > 0.7` → 1 (alt); ≤ 0.7 → 0 (ref); no data → '.' → 0.
  - Split-site imputation: '.' → 0 uniformly. At split rows, '.' = "doesn't have this specific allele" = 0. For the minor alt (V2) row, C-carriers and ref individuals both get 0, not a copy of the dominant-row call. This keeps columns orthogonal and semantically precise (2847_T/A=1 means specifically "has T→A", not "has any alt at 2847").
  - Never-called rows (all zeros after imputation, e.g. T→G when no individual had G) are dropped.
  - Haplotype string: "C_" + concatenated 0/1 across all retained sites in VCF order.
- **`scripts/06_snpeff_mac.sh`** — fixed SnpEff path (`~/SnpEff/`) and database name (`NC_012312.1` with Vertebrate_Mitochondrial codon table). Added chromosome name pre-flight check and commented rename block in case FASTA header doesn't match NC_012312.1.
- **`scripts/08_call_haplotypes.py`** — new Python script (cyvcf2 + pandas). Parses VCF, applies AD threshold, imputes '.' → 0, drops never-called rows, writes `vcf/haplotype_matrix.csv` (sites × samples) and `vcf/haplotype_calls.csv` (sample → haplotype string + N_alt_sites). Excludes samples 70 and 125 (no phenotypes).
- **Next-session pickup (in order):**
  1. Check FASTA header: `grep "^>" Missing_Files/SSM_MT_ref/Fhet_MT.fasta | head -1` — must be `NC_012312.1` for SnpEff to annotate. Patch stage 06 rename block if not.
  2. Wait for stage 05 to finish on Triton 2, then rsync VCF to Mac `vcf/`.
  3. `bash scripts/06_snpeff_mac.sh` → annotated VCF.
  4. `bash scripts/07_cds_snps_norm_mac.sh` → canonical `Fhet_MT_CDS.snps.split.vcf.gz`.
  5. `pip install cyvcf2 pandas --break-system-packages` if not already installed on Mac.
  6. `python scripts/08_call_haplotypes.py` → haplotype matrix + calls.

## 2026-05-08 (session 2) — stages 06 + 07 rewritten as Mac scripts; stage 05 submitted

- **BAM move complete (Triton 2):** `bams/MT_bam_sam/` promoted to `bams/` — confirmed 144 BAMs at top level. `bam_list.txt` generated with 143 entries (1_0 excluded). Stage 05 submitted.
- **SnpEff platform decision:** SnpEff cannot be built on Triton 2 (linux-ppc64le). SnpEff `Fhet_MT` database already exists on Mac (and Pegasus2). VCF from stage 05 is small enough to annotate locally. Stages 06 + 07 moved to Mac.
- **`jobs/06_snpeff_annotate.sh`** — replaced with a tombstone (exit 1) redirecting to `scripts/06_snpeff_mac.sh`.
- **`jobs/07_cds_snps_norm.sh`** — replaced with a tombstone (exit 1) redirecting to `scripts/07_cds_snps_norm_mac.sh`.
- **`scripts/06_snpeff_mac.sh`** — new Mac bash script. Downloads `Fhet_mt_variantsAD.vcf.gz` from Triton 2 via rsync (recipe in header), then runs `java -jar snpEff.jar ann Fhet_MT ... | bcftools view -Oz`. Update `SNPEFF_JAR` path if Mac install differs from `~/software/snpEff/`.
- **`scripts/07_cds_snps_norm_mac.sh`** — new Mac bash script. CDS restrict (awk GFF → bgzip/tabix regions) + `bcftools view -v snps` + `bcftools norm -m -any -f REF` → canonical `Fhet_MT_CDS.snps.split.vcf.gz`. Uses `Missing_Files/SSM_MT_ref/` for REF + GFF. No bedtools required.
- **`docs/01_pipeline.md`** — pipeline diagram and stage sections updated to reflect Mac-side 06 + 07.
- **Next-session pickup (in order):**
  1. Wait for stage 05 (`fhet_mpileup_AD`) to finish on Triton 2: `bjobs` or check `logs/05_mpileup_AD_*.out`.
  2. `rsync` `Fhet_mt_variantsAD.vcf.gz` (+ `.csi`) from Triton 2 → Mac `vcf/` — command in `docs/01_pipeline.md §06`.
  3. `bash scripts/06_snpeff_mac.sh` — SnpEff annotation on Mac.
  4. `bash scripts/07_cds_snps_norm_mac.sh` → canonical `vcf/Fhet_MT_CDS.snps.split.vcf.gz`.
  5. Begin Mac-side Python haplotype parsing.

## 2026-05-08 — conda env fix, git branch cleanup, BAM validation started

- **conda env fixed**: `samtools 1.6` was linked against `htslib 1.2.1` (biobuilds ABI mismatch — `libhts.so.2` missing). Removed both, reinstalled `samtools=1.6`, `htslib=1.6`, `bcftools=1.6` from biobuilds. bioconda not usable on `linux-ppc64le`. Both tools now verify (`samtools --version`, `bcftools --version` each report 1.6).
- **SnpEff note**: not available via conda on ppc64le. Stage 06 must use a standalone JAR (`java -jar snpEff.jar`), not a conda install.
- **Conda env snapshot**: `conda list -n mito_genomics > docs/mito_genomics_env.txt` — committed to repo for reproducibility.
- **Git branch cleanup**: Triton 2 local branch was `master`, not `main`. Renamed to `main`, pushed, deleted stray `master` from GitHub. Triton 2 and Mac now both track `origin/main`.
- **BAM validation complete**: `samtools quickcheck` across all 144 BAMs — one failure: `1_0_MT.bam` (no index, forward-read-only sample already on drop list). 143 BAMs clean. Per-BAM mapped-fraction sweep (`samtools idxstats`) confirmed 0.04–0.54% MT mapping across all 143 usable samples — low but expected for WGS aligned to the 16 kb MT genome only (unmapped nuclear reads retained in BAMs). Results logged to `logs/bam_mapfrac_20260508.txt`.
- **Next-session pickup (in order):**
  1. `mv bams/MT_bam_sam/*.bam bams/MT_bam_sam/*.bai bams/ && rmdir bams/MT_bam_sam` to match `BAMS_DIR`.
  3. Patch and submit stage 05 (`05_bcftools_mpileup_call_AD.sh`) for joint call across all BAMs.
  4. Write stage 06 (SnpEff annotation via standalone JAR using `$GFF`).
  5. Write stage 07 (CDS-restrict + SNPs-only + `bcftools norm -m -any`) → produces canonical `Fhet_MT_CDS.snps.split.vcf.gz`.

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
- **Pilot run + failure analysis**: submitted `bsub -J 'fhet_align_mt_pilot[2-6]%5' < jobs/04_bwa_align_mt.sh` (job `7646`). All 5 array tasks went EXIT in ~34s (LSF overhead only). Root cause from `04_align_mt_2.err`: `conda: error: argument COMMAND: invalid choice: 'activate'`. `module load anaconda3` puts conda in PATH but doesn't register the shell `conda activate` function on compute nodes. **Patched `jobs/config.sh`** to source `$(conda info --base)/etc/profile.d/conda.sh` (guarded on `CONDA_SHLVL` so it's a no-op when already initialized) — propagates the fix to all five stages.
- **Login-node env breakage (independent issue)**: `samtools` on login (`mgt3`) fails with `error while loading shared libraries: libhts.so.2: cannot open shared object file`. htslib ABI mismatch in the `mito_genomics` env. Pending: `conda install -n mito_genomics -c bioconda samtools htslib bcftools --force-reinstall -y`. This blocks BAM validation and stage 05 until fixed.
- **Existing canonical BAMs located**: 144 `*_MT.bam` + `.bai` files dated `Jul 25 2025`, sitting at `/projectnb/dcrawford/MT_Genomics2/bams/MT_bam_sam/` — one directory deeper than `BAMS_DIR` expects. Per-BAM size 15–19 GB (consistent with high MT coverage + retained unmapped reads, since stage 04 doesn't filter). Almost certainly the canonical inputs to the previous (still-missing) `Fhet_MT_CDS.snps.split.vcf.gz`. **Likely makes a stage-04 rerun unnecessary.** The pilot was never going to "find" them — the pilot writes to `$BAMS_DIR/` (top level), which was empty; the existing BAMs are nested one level in.
- **Next-session pickup (in order):**
  1. `conda install -n mito_genomics -c bioconda samtools htslib bcftools --force-reinstall -y` on Triton 2; verify with `samtools --version` + `bcftools --version`.
  2. `git pull origin main` on Triton 2 to get the patched `config.sh`.
  3. Validate BAMs: `samtools quickcheck` on each, then per-BAM mapped-fraction loop across all 144 in `bams/MT_bam_sam/`.
  4. If valid: `mv bams/MT_bam_sam/*.bam* bams/ && rmdir bams/MT_bam_sam` so layout matches `BAMS_DIR`.
  5. Skip stage-04 rerun. Resubmit a 1-sample stage-04 sanity test only if you want to validate the patched script for future re-alignment work.
  6. Patch and rerun stage 05 (`05_bcftools_mpileup_call_AD.sh`) to produce fresh joint-call VCFs in `$VCF_DIR`.
  7. Write stages 06 (SnpEff using `$GFF`) + 07 (CDS-restrict + SNPs-only + `bcftools norm -m -any`) to produce a fresh canonical `Fhet_MT_CDS.snps.split.vcf.gz`.

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

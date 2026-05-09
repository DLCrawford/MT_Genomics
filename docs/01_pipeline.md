# 01 — Pipeline (replication walkthrough)

Stage-by-stage. Each stage is one BSUB script under `jobs/`. Submit with `bsub < jobs/NN_*.sh` from `/projectnb/dcrawford/MT_Genomics2/` after the [per-session startup](00_setup.md#per-session-startup).

## Pipeline at a glance

```
raw fastqs
   │
   ▼
01_fastqc_raw.sh ──────────────────► fastqc_out/              ← QC raw reads
   │
   ▼
02_trim_pe.sh ─────────────────────► TrimA_seq/               ← Trimmomatic PE
   │
   ▼
03_fastqc_trimmed.sh ──────────────► TrimA_fastqc_out/        ← QC trimmed reads
   │
   ▼
04_bwa_align_mt.sh ────────────────► bams/*_MT.bam            ← align to MT, sort, index
   │
   ▼
05_bcftools_mpileup_call_AD.sh ────► vcf/Fhet_mt_*AD.vcf.gz  ← joint call (all vars, AD/DP)
   │
   ▼
06_snpeff_annotate.sh ─────────────► vcf/Fhet_mt_variantsAD_ann.vcf.gz  ← ANN field
   │
   ▼
07_cds_snps_norm.sh ───────────────► vcf/Fhet_MT_CDS.snps.split.vcf.gz  ← CANONICAL ★
```

★ `Fhet_MT_CDS.snps.split.vcf.gz` is the frozen canonical output. Do not rerun upstream stages against it. Downstream (on Mac): Python haplotype parsing → final ~978 sites × 141 samples table.

## Where paths come from

All path variables live in **`jobs/config.sh`**, which every script sources at the top. Paths are consolidated under `/projectnb/dcrawford/MT_Genomics2/` (as of 2026-05-07). To relocate data, edit only `config.sh` — scripts don't change.

## Inputs you must have in place

| `config.sh` variable | Path | What |
|---|---|---|
| `SAMPLE_LIST` | `${PROJECT_ROOT}/SSM_WGS_list.txt` | 143 sample IDs, one per line |
| `REF` | `${PROJECT_ROOT}/refs/Fhet_MT.fasta` | MT reference (bwa + samtools indexed) |
| `GFF` | `${PROJECT_ROOT}/refs/Fhet_MT.gff` | MT gene annotations (for stage 06 db build + stage 07 CDS filter) |
| `BAMS_DIR` | `${PROJECT_ROOT}/bams/` | sorted+indexed per-sample MT BAMs after stage 04 or mv |
| `BAM_LIST` | `${PROJECT_ROOT}/bam_list.txt` | one BAM *filename* per line, consumed by stage 05 |
| `VCF_DIR` | `${PROJECT_ROOT}/vcf/` | all VCF outputs from stages 05–07 |
| `LOGS_DIR` | `${PROJECT_ROOT}/logs/` | BSUB `.out`/`.err` (created automatically by `config.sh`) |
| `RAW_DIR` | (cold storage — sentinel only) | stages 01 + 02 not runnable |
| `ADAPTERS` | `/home/dcrawford/software/local/Trimmomatic-0.39/adapters/CombinedAdapters.fa` | for stage 02 only |
| `TRIMJAR` | `/home/dcrawford/software/local/Trimmomatic-0.39/trimmomatic-0.39.jar` | for stage 02 only |

## 01 — Raw QC

**Script:** `jobs/01_fastqc_raw.sh`
**Array:** `[1-287]%20`
**Inputs:** `fhet_raw_seq/{sample}_{1,2}.fq.gz`
**Outputs:** `fastqc_out/{sample}_{1,2}_fastqc.{html,zip}`

After all array tasks finish, roll up with MultiQC from the login node:
```bash
conda activate mito_genomics
multiqc /projectnb/dcrawford/SSM_WGS/fastqc_out -o /projectnb/dcrawford/SSM_WGS/fastqc_out
```

## 02 — Trim

**Script:** `jobs/02_trim_pe.sh`
**Array:** `[1-287]%50`
**Inputs:** `fhet_raw_seq/{sample}_{1,2}.fq.gz`
**Outputs:** `TrimA_seq/{sample}_{1,2}_p.fq.gz` (paired) + `..._up.fq.gz` (unpaired)
**Params:** `ILLUMINACLIP:CombinedAdapters.fa:2:30:10:2:True LEADING:3 TRAILING:3 MINLEN:36`

## 03 — Trimmed QC

**Script:** `jobs/03_fastqc_trimmed.sh`
**Array:** `[1-287]%20`
**Inputs:** `TrimA_seq/{sample}_{1,2}_p.fq.gz`
**Outputs:** `TrimA_fastqc_out/...`

MultiQC roll-up after all tasks finish (same pattern as 01).

## 04 — MT alignment

**Script:** `jobs/04_bwa_align_mt.sh`
**Array:** `[1-144]%30`  (only the 144 samples relevant to MT; index 1 = `1_0` which has no paired trimmed reads — skip)
**Inputs:** trimmed paired reads + `refs/Fhet_MT.fasta` (indexed)
**Naming:** on-disk reads are Trim Galore convention: `{sample}_1_val_1.fq.gz`, `{sample}_2_val_2.fq.gz`
**Pipeline:** `bwa mem -t 12 -M | samtools view -bS | samtools sort | samtools index`
**Outputs:** `bams/{sample}_MT.bam` + `.bam.bai`

> **Status 2026-05-08:** 144 canonical BAMs dated Jul 2025 recovered from `bams/MT_bam_sam/`. 143 passed `samtools quickcheck` (1_0_MT.bam excluded). BAMs moved to `bams/` to match `$BAMS_DIR`. Stage-04 rerun not needed.

## 05 — Joint variant calling with AD

**Script:** `jobs/05_bcftools_mpileup_call_AD.sh`
**Not an array** — one job, all 143 BAMs at once.
**Inputs:**
- `refs/Fhet_MT.fasta`
- `bams/*_MT.bam` (sorted+indexed; must be at top level of `$BAMS_DIR`)
- `bam_list.txt` (one BAM *filename* per line; script prepends `$BAMS_DIR/`)

**Generate `bam_list.txt` before submitting:**
```bash
ls /projectnb/dcrawford/MT_Genomics2/bams/*_MT.bam \
    | xargs -n1 basename \
    | grep -v '^1_0_MT\.bam$' \
    | sort \
    > /projectnb/dcrawford/MT_Genomics2/bam_list.txt
wc -l /projectnb/dcrawford/MT_Genomics2/bam_list.txt   # expect 143
```

**Outputs:**
- `vcf/Fhet_mt_variantsAD.vcf.gz` — variant sites with AD/DP; sample names stripped of `_0`
- `vcf/Fhet_mt_fullAD.vcf.gz` — same, all alt alleles retained
- `vcf/Fhet_mt_*AD_stats.txt` — `bcftools stats` summaries

**Notable bcftools flags:**
- `mpileup -a AD,DP` — per-allele + total depth (required for haplotype work)
- `mpileup -Q 30 -q 30` — base + mapping quality floor
- `mpileup -d 100000` — depth cap raised (MT is high-depth)
- `call -m` — multi-allelic caller model
- `call -v` — variant sites only
- `call -A` — keep all alt alleles present in the alignments

**Status 2026-05-08:** Submitted 11:32 AM, running. Ploidy warning ("assuming diploid") is benign — haplotype caller uses AD ratio, not GT.

**When the log shows `=== DONE ===`, run the verify script on Triton 2 before rsyncing:**
```bash
cd /projectnb/dcrawford/MT_Genomics2
bash scripts/verify_stage05.sh
```
This checks 143 samples, NC_012312.1 chromosome, AD/DP FORMAT fields, and SNP count range. It prints the rsync commands when all checks pass.

## 06 — SnpEff annotation (Mac)

**Script:** `scripts/06_snpeff_mac.sh` — runs locally on Mac, **not a BSUB job**
**Why Mac:** SnpEff cannot be built on Triton 2 (linux-ppc64le); the `Fhet_MT` database is already built on Mac. The VCF from stage 05 is small enough (~143 samples × MT variants) to annotate locally in seconds.

**First: download the stage-05 VCF from Triton 2:**
```bash
rsync -avP \
  dcrawford@scc1.bu.edu:/projectnb/dcrawford/MT_Genomics2/vcf/Fhet_mt_variantsAD.vcf.gz \
  ~/Projects/MT_Genomics_Cl_Ap2026/MT_Genomics2/vcf/
rsync -avP \
  dcrawford@scc1.bu.edu:/projectnb/dcrawford/MT_Genomics2/vcf/Fhet_mt_variantsAD.vcf.gz.csi \
  ~/Projects/MT_Genomics_Cl_Ap2026/MT_Genomics2/vcf/
```

**Then run:**
```bash
cd ~/Projects/MT_Genomics_Cl_Ap2026/MT_Genomics2
bash scripts/06_snpeff_mac.sh
```

**Input:** `vcf/Fhet_mt_variantsAD.vcf.gz`
**Outputs:**
- `vcf/Fhet_mt_variantsAD_ann.vcf.gz` — ANN field added
- `vcf/snpEff_summary.html` + `snpEff_summary.genes.txt` — per-run stats

**SnpEff install (Mac):**
- JAR: `~/micromamba/envs/SNP_env/share/snpeff-5.2-1/snpEff.jar` (micromamba `SNP_env`)
- Config: `~/snpEff/snpEff.config`
- Database: `Fhet_MT` (custom-built; `~/snpEff/data/Fhet_MT/snpEffectPredictor.bin` must exist)
- Chromosome in VCF (`NC_012312.1`) matches database — no rename needed

If you reinstall SnpEff, update `SNPEFF_JAR` in the script. To find the new path: `find ~/micromamba/envs/SNP_env -name "snpEff.jar"`

## 07 — CDS restrict + SNPs only + split multiallelic → canonical (Mac)

**Script:** `scripts/07_cds_snps_norm_mac.sh` — runs locally on Mac, **not a BSUB job**
**Needs:** bcftools (brew), bgzip + tabix (bundled with htslib/bcftools)

```bash
cd ~/Projects/MT_Genomics_Cl_Ap2026/MT_Genomics2
bash scripts/07_cds_snps_norm_mac.sh
```

**Input:** `vcf/Fhet_mt_variantsAD_ann.vcf.gz` (from stage 06)
**Outputs:**
- `vcf/Fhet_MT_CDS.snps.split.vcf.gz` — **CANONICAL OUTPUT** (frozen once produced)
- `vcf/MT_CDS.regions.gz` + `.tbi` — CDS intervals derived from GFF (kept for audit)
- `vcf/Fhet_MT_CDS.snps.split_stats.txt` — `bcftools stats` summary

**Pipeline:**
1. `awk` extracts `CDS` features from `Missing_Files/SSM_MT_ref/Fhet_MT.gff` → bgzip + tabix regions file
2. `bcftools view -R MT_CDS.regions.gz` — restrict to CDS positions
3. `bcftools view -v snps` — SNPs only
4. `bcftools norm -m -any -f REF` — split multiallelic, left-align, trim

## After stage 07 — frozen canonical output

`Fhet_MT_CDS.snps.split.vcf.gz` feeds all Mac-side analysis. Do not rerun the upstream pipeline against it.

## Submitting a stage

```bash
cd /projectnb/dcrawford/MT_Genomics2
bsub < jobs/05_bcftools_mpileup_call_AD.sh
bjobs                     # see queued/running
bjobs -A <jobid>          # array task summary
bkill <jobid>             # if needed
```

Logs land in `/projectnb/dcrawford/MT_Genomics2/logs/NN_<stage>_<jobid>.{out,err}`.

## Resubmitting a single failed array task

```bash
bsub -J 'fhet_align_mt[42]' < jobs/04_bwa_align_mt.sh   # re-run task #42 only
```

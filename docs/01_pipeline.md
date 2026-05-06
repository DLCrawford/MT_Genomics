# 01 — Pipeline (replication walkthrough)

Stage-by-stage. Each stage is one BSUB script under `jobs/`. Submit with `bsub < jobs/NN_*.sh` from `/projectnb/dcrawford/MT_Genomics2/` after the [per-session startup](00_setup.md#per-session-startup).

## Pipeline at a glance

```
raw fastqs
   │
   ▼
01_fastqc_raw.sh ──────► fastqc_out/        ← QC raw reads
   │
   ▼
02_trim_pe.sh ─────────► TrimA_seq/         ← Trimmomatic PE
   │
   ▼
03_fastqc_trimmed.sh ──► TrimA_fastqc_out/  ← QC trimmed reads
   │
   ▼
04_bwa_align_mt.sh ────► MT_bam_sam/        ← align to MT, sort, index
   │
   ▼
05_bcftools_mpileup_call_AD.sh ──► new_hap_AD/Fhet_mt_*AD.vcf.gz
                                              (joint call, with AD)
```

Downstream (on Mac, not in `jobs/`): SnpEff annotation, `bcftools norm -m -any`, Python haplotype parsing → final ~978-sites × 141-samples table.

## Where paths come from

All path variables live in **`jobs/config.sh`**, which every script sources at the top. The table below shows the *active legacy* layout; to consolidate paths under `/projectnb/dcrawford/MT_Genomics2/` (via physical move or symlinks), edit only `config.sh` — scripts don't change.

## Inputs you must have in place

| `config.sh` variable | Default path | What |
|---|---|---|
| `SAMPLE_LIST` | `/projectnb/dcrawford/SSM_WGS/SSM_WGS_list.txt` | one sample id per line |
| `RAW_DIR` | `/projectnb/dcrawford/SSM_WGS/fhet_raw_seq/` | raw paired fastqs |
| `REF` | `/projectnb/dcrawford/SSM_Mito/Fh_MT_ref/Fhet_MT.fasta` | MT reference (must be bwa+samtools-indexed; see [00_setup.md §5](00_setup.md#5-reference-genome-one-time)) |
| `ADAPTERS` | `/home/dcrawford/software/local/Trimmomatic-0.39/adapters/CombinedAdapters.fa` | TruSeq3-PE.fa + NexteraPE-PE.fa concatenated |
| `TRIMJAR` | `/home/dcrawford/software/local/Trimmomatic-0.39/trimmomatic-0.39.jar` | Trimmomatic jar |
| `LOGS_DIR` | `/projectnb/dcrawford/MT_Genomics2/logs/` | BSUB `.out`/`.err` (created automatically by `config.sh`) |
| `BAM_LIST` | `/projectnb/dcrawford/SSM_Mito/new_hap_AD2/bam_list.txt` | one BAM filename per line, used by step 05 |

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
**Array:** `[1-144]%30`  (only the 144 samples relevant to MT)
**Inputs:** trimmed paired reads + `Fhet_MT.fasta` (indexed)
**Pipeline:** `bwa mem -t 12 -M | samtools view -bS | samtools sort | samtools index`
**Outputs:** `MT_bam_sam/{sample}_MT.bam` + `.bam.bai`

## 05 — Joint variant calling with AD

**Script:** `jobs/05_bcftools_mpileup_call_AD.sh`
**Not an array** — one job, all 143 BAMs at once.
**Inputs:**
- `Fhet_MT.fasta`
- `MT_bam_sam/*.bam` (sorted+indexed)
- `bam_list.txt` (one BAM filename per line; the script prepends `$BAMDIR/`)

**Outputs:**
- `Fhet_mt_fullAD.vcf.gz` — every position, with AD/DP
- `Fhet_mt_variantsAD.vcf.gz` — variant-only subset
- `Fhet_mt_*AD_stats.txt` — `bcftools stats` summaries

**Notable bcftools flags:**
- `mpileup -a AD,DP` — required for haplotype work
- `mpileup -Q 30 -q 30` — base + mapping quality floor
- `mpileup -d 100000` — depth cap raised (MT is high-depth)
- `call -A` — keep all sites, not just variants
- `call -mv` — multi-allelic + variant-aware caller

The script reheaders sample names by stripping `_0` (e.g. `10_0` → `10`) on both VCFs.

## After 05 — frozen output

`Fhet_MT_CDS.snps.split.vcf.gz` (downstream of 05 + SnpEff + `bcftools norm -m -any`) is the canonical input to all Mac-side analysis. Do not rerun the upstream pipeline against it.

## Submitting a stage

```bash
cd /projectnb/dcrawford/MT_Genomics2
bsub < jobs/01_fastqc_raw.sh
bjobs                     # see queued/running
bjobs -A <jobid>          # array task summary
bkill <jobid>             # if needed
```

Logs land in `/projectnb/dcrawford/MT_Genomics2/logs/NN_<stage>_<jobindex>.{out,err}`.

## Resubmitting a single failed array task

```bash
bsub -J 'fhet_align_mt[42]' < jobs/04_bwa_align_mt.sh   # re-run task #42 only
```

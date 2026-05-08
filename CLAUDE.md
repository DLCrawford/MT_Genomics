# CLAUDE.md — Project rules for Claude Code

## Project goal

Build a reproducible mitochondrial variant + haplotype pipeline for *Fundulus heteroclitus*, with heavy compute on Triton 2 (LSF/BSUB) and downstream analysis on Doug's Mac.

## Compute split

- **HPC (Triton 2, LSF)** — anything heavy:
  - FastQC across 287 paired raw fastqs
  - Trimmomatic PE
  - FastQC on trimmed reads
  - bwa mem + samtools sort/index against `Fhet_MT.fasta`
  - `bcftools mpileup | bcftools call` (joint, with AD/DP) across all MT BAMs
  - SnpEff annotation (when added)
  - `bcftools norm -m -any` (split multiallelic sites)
  - Final canonical output: `Fhet_MT_CDS.snps.split.vcf.gz` ← FROZEN once produced, do not rerun.
    Status 2026-05-08: conda env fixed (samtools/htslib/bcftools all 1.6 from biobuilds).
    BAM validation in progress (`samtools quickcheck` + mapped-fraction sweep across 144
    BAMs in `bams/MT_bam_sam/`). If valid, next steps: mv BAMs to `bams/`, run stage 05
    (joint call), write + run stage 06 (SnpEff — standalone JAR, not conda; ppc64le),
    write + run stage 07 (CDS+SNPs+norm). See CHANGELOG 2026-05-08.

- **Mac (local)** — fast iteration:
  - Python haplotype parsing
  - Annotation extraction from VCF `ANN` field
  - Final ~978 sites × 141 samples table
  - Tools: bcftools (brew), cyvcf2, pandas

## Rules for Claude

1. **Don't rerun frozen outputs.** `Fhet_MT_CDS.snps.split.vcf.gz` is canonical *once produced*. Build downstream analysis on top of it. (Status 2026-05-08: env fixed; BAM validation in progress; rerun starts at stage 05 if BAMs are valid — see CHANGELOG.)
2. **Don't commit data files.** `.gitignore` excludes `*.vcf.gz`, `*.bcf`, `*.bam`, `*.bai`, `*.fastq*`, `*.fasta`, `*.tbi`, `*.csi`. Scripts and docs only.
3. **All BSUB scripts live under `jobs/`** with the numbered prefix (`01_..05_..`). New stages get the next number.
4. **All log files write to `/projectnb/dcrawford/MT_Genomics2/logs/`** so they're easy to find and clean up.
5. **One conda env to start: `mito_genomics`.** Split per stage only if conflicts force it; document any split in `docs/00_setup.md`.
6. **Update `CHANGELOG.md` at the end of every working session** with what changed and why.
7. **When asked to add a new pipeline step,** write the BSUB script, update `docs/01_pipeline.md`, and append a CHANGELOG entry — in that order.

## Sample / sample-list conventions

- `SSM_WGS_list.txt`: one sample id per line, e.g. `10_0`, `102_0`, `1_0`.
- Raw fastq names: `{sample}_1.fq.gz`, `{sample}_2.fq.gz`.
- Trimmed paired (actual files in `trimmed_seq/` on Triton 2): `{sample}_1_val_1.fq.gz`, `{sample}_2_val_2.fq.gz` (Trim Galore output), with companion `{sample}_{1,2}.fq.gz_trimming_report.txt`. Note: `02_trim_pe.sh` in this repo is written for Trimmomatic naming (`_p` / `_up`); see CHANGELOG 2026-05-07 — reconcile if stage 02 is ever re-run.
- MT BAMs: `{sample}_MT.bam` (+ `.bai`).
- The trailing `_0` is stripped from sample names in the joint-call VCF (`bcftools reheader`).

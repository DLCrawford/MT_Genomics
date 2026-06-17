# jobs/

Numbered, drop-in BSUB scripts for the Triton 2 side of the pipeline.
Each script is self-contained: a header documents its stage, inputs,
outputs, and submit command.

| #     | Script                            | Stage |
|-------|-----------------------------------|-------|
| —     | `config.sh`                       | **All path variables.** Sourced by every script — edit here, not in scripts. PATH-injects locally-built bcftools/samtools/htslib 1.23.1. |
| 01    | `01_fastqc_raw.sh`                | FastQC on raw paired fastqs |
| 02    | `02_trim_pe.sh`                   | Trimmomatic PE (Combined adapters) |
| 03    | `03_fastqc_trimmed.sh`            | FastQC on trimmed paired fastqs |
| 04    | `04_bwa_align_mt.sh`              | bwa mem → samtools sort → index against `Fhet_MT.fasta` |
| —     | `BSUB_Slim_BAM_mt.sh`             | Extract mt-only reads (`NC_012312.1`, `-F 2308`) from each full WGS BAM → `MT_only_bams/${SAMPLE}_MT_only.bam`. LSF array. |
| **05_1** | `05_1_mpileup_merge.sh`        | **Canonical caller (2026-05-18).** Joint `bcftools mpileup -b slim_bamlist \| call -mv --ploidy 1 \| norm -m -any -f REF` in one pipeline → `vcf/05_1_Fhet_mt_persample_merged.vcf.gz`. Per-cell DP/AD on every (POS × sample); output already in one-row-per-ALT form. Built-in regression check against `Fhet_mt_persample_merged.vcf.gz` (the 1128-SNP per-sample baseline from sessions 6–10). |
| 06    | `06_snpeff_annotate.sh`           | **Superseded — stub.** SnpEff cannot build on linux-ppc64le; annotation runs Mac-side via `scripts/06_snpeff_mac.sh`. |
| 07    | `07_cds_snps_norm.sh`             | **Superseded — stub.** CDS restriction runs Mac-side via `scripts/07_cds_snps_norm_mac.sh`. |

See `../CLAUDE.md` for architecture rationale, the validation gate
around 05_1's output, and the known pitfalls section that records the
three structural bugs sessions 6–14 worked through.

## Conventions

- Each script's header lists `Stage / Input / Output / Submit`.
- All scripts source `config.sh` for paths — never hardcode paths in
  a script.
- `module load anaconda3` then `source config.sh`. Conda activation
  is optional — `config.sh` PATH-injects `$HOME/software/local/bin`
  (bcftools / samtools / htslib 1.23.1 built from source) ahead of
  the `mito_genomics` env's pinned 1.6, which silently collapses
  high-depth mtDNA pileups.
- `set -euo pipefail` everywhere.
- Logs to `/projectnb/dcrawford/MT_Genomics2/logs/NN_<stage>_<jobid>.{out,err}`.
- Scripts 01–04 + `BSUB_Slim_BAM_mt.sh` are LSF arrays driven by
  `$SAMPLE_LIST` / `$BAM_LIST` via `LSB_JOBINDEX`.
- 05_1 is a single joint-call job (`-n 8`, 16 GB, 4 h, normal queue).

## `_archive/`

Historical scripts kept for the methods write-up. Not in the active
pipeline. Do not submit.

| File | What it was |
|---|---|
| `05_bcftools_mpileup_call_AD.sh` | Sessions 6–7 joint-call v1 (strict). Hit the bcftools-1.6 high-depth collapse. |
| `05b_v2_Q13_q20_p1.sh` | Joint-call v2 (relaxed Q/q + haploid). Same collapse. |
| `05c_v3_Q13_q00_p1.sh` | Joint-call v3 (no MAPQ + haploid). Same collapse. |
| `05d_persample_call.sh` | Per-sample call against full BAMs with the env's bcftools 1.6 (historical). |
| `05e_merge_persample.sh` | Merge for 05d (historical). |
| `05d2_persample_call.sh` | Session-10 per-sample caller against slim BAMs with bcftools 1.23.1. **Produced the 1128-SNP May-15 baseline** that 05_1 regression-tests against. |
| `05e2_merge_persample.sh` | Session-10 merge for 05d2. |
| `05f_joint_call.sh` | Session-13 first-cut joint caller (no norm-split inside the script; baseline check too soft). Superseded by 05_1. |
| `07b_backfill_AD.sh` | Session-11 T2-side AD/DP backfill. No longer needed: 05_1 fills DP/AD on every cell natively. Also carried a targets-file bug (one row per split-ALT vs the format `call -C alleles` needs) that silently dropped 17 secondary ALTs at multi-row positions — see CHANGELOG 2026-05-17/18. |

`../scripts/_archive/` contains the Mac-side counterparts:
`compare_stage05_runs.sh`, `run_stage05_core.sh`, `verify_stage05.sh`
— diagnostic helpers from the sessions 6–10 variant-count
investigation, retained for the methods write-up.

## Submit (canonical caller)

```bash
bsub < jobs/05_1_mpileup_merge.sh
# logs: logs/05_1_mpileup_merge_<jobid>.{out,err}
# output: vcf/05_1_Fhet_mt_persample_merged.vcf.gz
# walltime expectation: 30–60 min on 143 slim BAMs
```

After completion, **read the manifest before flowing downstream**:

```bash
cat vcf/05_1_Fhet_mt_persample_merged_run_manifest.txt
# Must show:
#   cells with DP=.: 0
#   n_multiallelic_sites: 0   (norm-split done)
#   positions lost vs baseline: 0
#   n_SNPs: ~1100-1200
```

If `positions lost vs baseline > 0`, stop and inspect
`vcf/05_1_Fhet_mt_persample_merged_vs_baseline_lost.tsv` before
running 06/07 on Mac.

## Resubmit a single failed array task

```bash
bsub -J 'fhet_align_mt[42]' < jobs/04_bwa_align_mt.sh
```

## Relocating data

To move data into `/projectnb/dcrawford/MT_Genomics2/` (or symlink it
in), edit only `config.sh`. The scripts pick up the new paths
automatically. See the comment header inside `config.sh` for the
symlink recipe.

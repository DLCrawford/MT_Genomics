# jobs/

Numbered, drop-in BSUB scripts. Each script is self-contained: a header documents its stage, inputs, outputs, and submit command.

| # | Script | Stage |
|---|---|---|
| — | `config.sh` | **All path variables.** Sourced by every script — edit here, not in scripts. |
| 01 | `01_fastqc_raw.sh` | FastQC on raw paired fastqs |
| 02 | `02_trim_pe.sh` | Trimmomatic PE (Combined adapters) |
| 03 | `03_fastqc_trimmed.sh` | FastQC on trimmed paired fastqs |
| 04 | `04_bwa_align_mt.sh` | bwa mem → samtools sort → index against MT reference |
| 05 | `05_bcftools_mpileup_call_AD.sh` | Joint call across all MT BAMs with AD/DP |

See `../docs/01_pipeline.md` for the full walkthrough.

## Conventions

- Each script's header lists `Stage / Input / Output / Submit`.
- All scripts source `config.sh` for paths — never hardcode paths in a script.
- `module load anaconda3` then `source config.sh` then `conda activate "$CONDA_ENV"`.
- `set -euo pipefail` everywhere.
- Logs all go to `/projectnb/dcrawford/MT_Genomics2/logs/NN_<stage>_<index>.{out,err}`.
- Scripts 01–04 are LSF arrays driven by `$SAMPLE_LIST` via `LSB_JOBINDEX`.
- Script 05 is a single joint-call job (no array).

## Relocating data

To move data into `/projectnb/dcrawford/MT_Genomics2/` (or symlink it in), edit only `config.sh`. The scripts pick up the new paths automatically. See the comment header inside `config.sh` for the symlink recipe.

## Submit

```bash
bsub < jobs/01_fastqc_raw.sh
```

## Resubmit a single failed array task

```bash
bsub -J 'fhet_align_mt[42]' < jobs/04_bwa_align_mt.sh
```

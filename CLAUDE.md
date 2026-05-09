# CLAUDE.md — Project rules for Claude Code

## Project goal

Build a reproducible mitochondrial variant + haplotype pipeline for *Fundulus heteroclitus*, with heavy compute on Triton 2 (LSF/BSUB) and downstream analysis on Doug's Mac.

## Active primary task — variant-count discrepancy (added 2026-05-09 session 6)

**Why does the same stage-05 pipeline produce 152 SNPs in one parameter
configuration and ~950 in another, and how do other investigators avoid
the trap?** This is now a *primary* project deliverable, not an incidental
bug fix. The 152-vs-~950 gap is qualitative: at 152 SNPs all downstream
haplotype and population-structure inferences would change; at ~950 only
minor revisions are needed and the overall scientific picture stands.
Whichever number proves correct, we owe ourselves and the community a
written record of *what went wrong, why, and how to avoid it*.

What we need to deliver, in order:

1. **Definitive call-set choice.** Pick which of v1 (strict, 152), v2
   (-Q 13 -q 20 --ploidy 1), or v3 (-Q 13 -q 0 --ploidy 1) is the
   biologically correct call set for *F. heteroclitus* mtDNA, with the
   evidence (Ts/Tv per allele rank, AF spectrum, singleton rate, agreement
   with prior canonical pipeline) that justifies the choice.

2. **Mechanism write-up.** Explain *which* parameter(s) of
   `bcftools mpileup | bcftools call` caused the gap, in enough mechanistic
   detail that someone reading the bcftools defaults can predict the
   failure mode without rerunning a sweep. Working hypotheses to confirm
   or refute: default `--ploidy 2` on a haploid genome (het-likelihood
   gating drops real homoplasmic sites); `-q 30` mapping-quality filter
   removing legitimate MT reads; interaction between `-A` and the
   multiallelic caller producing phantom alts that mask real signal in
   downstream stats.

3. **Practical checklist.** A short "what to set on mtDNA when calling with
   bcftools" recipe other investigators can drop into their pipeline
   without re-deriving from scratch. Likely lives in `docs/` once the
   diagnostic runs land.

The v2 + v3 diagnostic re-runs (`jobs/05b_*.sh`, `jobs/05c_*.sh`,
submitted 2026-05-09) are the experimental basis for items 1–2. Run
manifests (`vcf/*_run_manifest.txt`) plus the comparison helper
(`scripts/compare_stage05_runs.sh`) are the artifacts; the deliverable
prose draws on those plus the existing `Fhet_mt_*_stats.txt` for v1.

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
    Status 2026-05-09 (session 6): Stage 05 strict run completed (v1, -Q 30 -q 30
    ploidy=2) but yielded only 152 SNPs vs ~950 expected. Logs show clean exit
    (max mem 371 MB / 16 GB, runtime 9h / 72h budget) — NOT a truncation. Two
    diagnostic re-runs prepared:
      v2: jobs/05b_v2_Q13_q20_p1.sh  (-Q 13 -q 20 --ploidy 1)
      v3: jobs/05c_v3_Q13_q00_p1.sh  (-Q 13 -q  0 --ploidy 1)
    Both write to distinct RUN_TAG output prefixes so the v1 strict result is
    preserved on disk for comparison. Compare with scripts/compare_stage05_runs.sh.
    Next: bsub both → wait → compare → pick winner → rsync → 06 → 07 → 08.

- **Mac (local)** — fast iteration:
  - Python haplotype parsing
  - Annotation extraction from VCF `ANN` field
  - Final ~978 sites × 141 samples table
  - Tools: bcftools (brew), cyvcf2, pandas

## Rules for Claude

1. **Don't rerun frozen outputs.** `Fhet_MT_CDS.snps.split.vcf.gz` is canonical *once produced*. Build downstream analysis on top of it. (Status 2026-05-09 session 6: stage 05 v1 strict run produced 152 SNPs vs ~950 expected; v2 + v3 diagnostic re-runs queued. Canonical is not yet produced — see CHANGELOG.)
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

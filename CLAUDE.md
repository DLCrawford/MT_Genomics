# CLAUDE.md — Project rules for Claude Code

## Project goal

Build a reproducible mitochondrial variant + haplotype pipeline for *Fundulus heteroclitus*, with heavy compute on Triton 2 (LSF/BSUB) and downstream analysis on Doug's Mac.

## Active primary task — variant-count discrepancy (added 2026-05-09 session 6; resolved at the diagnostic level 2026-05-10 session 7)

**Why does the same stage-05 pipeline produce 152 SNPs in one parameter
configuration and ~1133 in another, and how do other investigators avoid
the trap?** Resolution as of session 7: the gap is **architectural**, not
parametric. The historical pipeline used per-sample `bcftools mpileup |
bcftools call` followed by `bcftools merge` across all 144 samples; the
current stage 05 (and v2/v3 diagnostic variants) use a single joint
`bcftools call -mv` piped from a multi-BAM `bcftools mpileup`. The switch
was made deliberately on collaborator advice to record per-sample AD/DP
for downstream haplotype calling — the trade-off was just never quantified
until now. Joint `-mv`'s default allele-frequency prior suppresses sites
where the reference is the rare allele, which on this mtDNA panel against
the divergent NC_012312.1 reference is most variable positions (1004 of
1133 historical SNPs sit at AF ≥ 0.95). No amount of `-Q`/`-q`/`--ploidy`
tuning recovers them — v1 (152), v2 (152), v3 (145) all converged.

Status of the three deliverables:

1. **Definitive call-set choice — pending stage 05d/05e completion.**
   Hypothesis (to confirm): the per-sample → merge run produces ~1130 SNPs
   matching the historical `merged_144.vcf.gz` (1133 SNPs / 144 samples).
   Once 05e lands, that becomes the canonical input to stages 06/07/08.

2. **Mechanism write-up — DONE.** See [`docs/02_calling_architecture.md`](docs/02_calling_architecture.md).
   Documents the architectural difference, the AF-prior suppression
   mechanism, AF-spectrum evidence, and a 7-point practical checklist for
   mtDNA SNP calling with bcftools.

3. **Practical checklist — DONE.** Section "Practical checklist for mtDNA
   SNP calling with bcftools" in `docs/02_calling_architecture.md`.

Experimental artifacts: `jobs/05_*.sh` (v1, joint strict),
`jobs/05b_*.sh` (v2, joint relaxed + haploid), `jobs/05c_*.sh` (v3, joint
no-MAPQ + haploid), `jobs/05d_persample_call.sh` + `jobs/05e_merge_persample.sh`
(per-sample → merge, replicating the historical recipe). Stats files in
`vcf/Fhet_mt_*_stats.txt` for the joint variants and
`stats_old/merged_stats.txt` for the historical baseline. Run manifests
sit alongside each stats file as `Fhet_mt_<RUN_TAG>_run_manifest.txt`.

## Compute split

- **HPC (Triton 2, LSF)** — anything heavy:
  - FastQC across 287 paired raw fastqs
  - Trimmomatic PE
  - FastQC on trimmed reads
  - bwa mem + samtools sort/index against `Fhet_MT.fasta`
  - `bcftools mpileup | bcftools call --ploidy 1` PER SAMPLE (with `-a AD,DP`),
    `bcftools norm -m -any` per sample, then `bcftools merge -m none` across
    all 143 BAMs. See `jobs/05d_persample_call.sh` + `jobs/05e_merge_persample.sh`.
  - SnpEff annotation (when added)
  - `bcftools norm -m -any` (split multiallelic sites)
  - Final canonical output: `Fhet_MT_CDS.snps.split.vcf.gz` ← FROZEN once produced, do not rerun.
    Status 2026-05-10 (session 7): joint-call v1/v2/v3 all converged to 152/152/145
    SNPs, refuting the parametric hypothesis. Historical baseline `merged_144.vcf.gz`
    (1133 SNPs / 144 samples) used per-sample → merge, recipe in
    `archive/Notes_dlcs/Inital_call_wo_AD.txt`. Per-sample → merge re-implementation
    landed as 05d (array) + 05e (merge); pending submission on Triton 2.
    v1/v2/v3 outputs preserved as the methods-comparison artifacts.
    Next: bsub 05d → wait → bsub 05e → verify ~1130 SNPs / ts/tv ≈ 7–9 → rsync → 06 → 07 → 08.

- **Mac (local)** — fast iteration:
  - Python haplotype parsing
  - Annotation extraction from VCF `ANN` field
  - Final ~978 sites × 141 samples table
  - Tools: bcftools (brew), cyvcf2, pandas

## Rules for Claude

1. **Don't rerun frozen outputs.** `Fhet_MT_CDS.snps.split.vcf.gz` is canonical *once produced*. Build downstream analysis on top of it. (Status 2026-05-10 session 7: stage 05 v1/v2/v3 all converged ≈ 150 SNPs; gap localized to joint-call architecture; per-sample → merge re-implementation written as `jobs/05d_*.sh` + `jobs/05e_*.sh` and pending submission. Canonical not yet produced — see CHANGELOG.)
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

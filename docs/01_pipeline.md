# 01 — Pipeline (replication walkthrough)

Stage-by-stage. Each stage is one BSUB script under `jobs/`. Submit with `bsub < jobs/NN_*.sh` from `/projectnb/dcrawford/MT_Genomics2/` after the [per-session startup](00_setup.md#per-session-startup).

> **Doc currency (2026-05-24, session 18).** The "Pipeline at a glance" diagram and the stage-05 walkthrough below still describe the **per-sample → merge** architecture (05d / 05e from sessions 6–10). The current canonical caller is `jobs/05_1_mpileup_merge.sh` (joint mpileup + call + norm-split in one pipeline, replacing 05d/05e/05f), with a hard baseline-regression gate against the 1128-SNP May-15 per-sample baseline. The per-sample → merge chain is preserved in `jobs/_archive/05d2_persample_call.sh` + `jobs/_archive/05e2_merge_persample.sh` as the validated baseline 05_1 must not regress against. Until this doc is refreshed end-to-end, see `CLAUDE.md` → "Canonical caller" and `CHANGELOG.md` 2026-05-18 (session 14) for the authoritative architecture description, and the "Downstream analysis (stages 08–18)" section at the bottom of this file for everything that runs on top of the frozen canonical.

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
05d_persample_call.sh    ──────────► vcf/persample/{sample}_norm.vcf.gz   ← per-sample haploid call
05e_merge_persample.sh   ──────────► vcf/Fhet_mt_persample_merged.vcf.gz  ← canonical joint VCF
   │      (this pair replaces the joint-call 05/05b/05c — see docs/02_calling_architecture.md)
   ▼
06_snpeff_annotate.sh ─────────────► vcf/Fhet_mt_persample_merged_ann.vcf.gz  ← ANN field
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

## 05 — Variant calling (per-sample → merge)

**Canonical scripts (use these):** `jobs/05d_persample_call.sh` (array)
followed by `jobs/05e_merge_persample.sh` (single job).

**Why this architecture:** the joint-call approach in `05/05b/05c` produced
~150 SNPs vs the historical ~1133 because joint `bcftools call -mv`'s
allele-frequency prior suppresses sites where the reference is the rare
allele — which is most variable positions on this mtDNA panel against
the divergent `NC_012312.1` reference. Per-sample call → `bcftools merge`
matches the historical recipe and recovers the full signal. Full mechanism
write-up in [`docs/02_calling_architecture.md`](02_calling_architecture.md).

**Joint variants kept for the methods comparison:**
`jobs/05_bcftools_mpileup_call_AD.sh` (v1, strict),
`jobs/05b_v2_Q13_q20_p1.sh`, `jobs/05c_v3_Q13_q00_p1.sh`. These remain in
`jobs/` as the experimental basis for the gap diagnosis; do not use them
for the canonical pipeline.

### 05d — per-sample call (array)

**Script:** `jobs/05d_persample_call.sh`
**Array:** `[1-143]%24` (one task per BAM in `$BAM_LIST`)
**Inputs:**
- `refs/Fhet_MT.fasta`
- `bams/{sample}_MT.bam` for each line in `bam_list.txt`

**Per-task pipeline:**
```
bcftools mpileup -f REF SAMPLE.bam -a AD,DP --max-depth 10000 -Ou \
  | bcftools call -mv --ploidy 1 -Oz -o SAMPLE.vcf.gz
bcftools index SAMPLE.vcf.gz
bcftools norm -m -any -Oz -o SAMPLE_norm.vcf.gz SAMPLE.vcf.gz
bcftools index SAMPLE_norm.vcf.gz
```

**Outputs:** `vcf/persample/{sample}.vcf.gz` and
`vcf/persample/{sample}_norm.vcf.gz` (+ `.csi` indices).

### 05e — merge

**Script:** `jobs/05e_merge_persample.sh`
**Not an array.** Run only after the 05d array finishes — the script
pre-flights every `*_norm.vcf.gz` and exits non-zero if any are missing.

**Pipeline:**
```
bcftools merge -m none --threads 8 -l input_list.txt -Oz -o Fhet_mt_persample_merged.vcf.gz
bcftools reheader -s <(bcftools query -l ... | sed 's/_0$//') ...
bcftools stats Fhet_mt_persample_merged.vcf.gz > Fhet_mt_persample_merged_stats.txt
```

**Outputs:**
- `vcf/Fhet_mt_persample_merged.vcf.gz` — single VCF, sample names with `_0` stripped
- `vcf/Fhet_mt_persample_merged_stats.txt` — `bcftools stats` summary
- `vcf/Fhet_mt_persample_merged_run_manifest.txt` — parameters + provenance + summary

**Target:** ~1133 SNPs / ts/tv ≈ 7.9, matching historical
`stats_old/merged_stats.txt`. With 143 samples (vs historical 144) the
count will run slightly lower.

**Submit (in this order):**
```bash
cd /projectnb/dcrawford/MT_Genomics2
bsub < jobs/05d_persample_call.sh
# wait for the array to finish:  bjobs / tail logs/05d_persample_*_*.out
bsub < jobs/05e_merge_persample.sh
```

### 05 (legacy joint) — kept for the methods comparison only

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
  dcrawford@scc1.bu.edu:/projectnb/dcrawford/MT_Genomics2/vcf/Fhet_mt_persample_merged.vcf.gz \
  ~/Projects/MT_Genomics_Cl_Ap2026/MT_Genomics2/vcf/
rsync -avP \
  dcrawford@scc1.bu.edu:/projectnb/dcrawford/MT_Genomics2/vcf/Fhet_mt_persample_merged.vcf.gz.csi \
  ~/Projects/MT_Genomics_Cl_Ap2026/MT_Genomics2/vcf/
```

**Note:** stages 06 and 07 currently expect `Fhet_mt_variantsAD.vcf.gz` as
their input filename (legacy from the joint-call pipeline). Either symlink:
```bash
cd ~/Projects/MT_Genomics_Cl_Ap2026/MT_Genomics2/vcf
ln -sf Fhet_mt_persample_merged.vcf.gz     Fhet_mt_variantsAD.vcf.gz
ln -sf Fhet_mt_persample_merged.vcf.gz.csi Fhet_mt_variantsAD.vcf.gz.csi
```
or update `INPUT` in `scripts/06_snpeff_mac.sh` / `scripts/07_cds_snps_norm_mac.sh`
to point at `Fhet_mt_persample_merged.vcf.gz` directly.

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

---

## Downstream analysis (stages 08–18)

Everything from here runs **on Mac** on top of the frozen canonical
`vcf/Fhet_MT_CDS.snps.split.vcf.gz`. No script in this block writes
back into the canonical; the order is mostly flexible (intra-block
dependencies are noted per stage). `CLAUDE.md` → "Downstream analysis
scripts" carries the full catalog with docstring-level detail. Activate
`SNP_env` (`micromamba activate SNP_env`) before running.

### Downstream at-a-glance

```
Fhet_MT_CDS.snps.split.vcf.gz  ★ frozen canonical
        │
        ├── 08_call_haplotypes.py        ── haplotype matrix at AF≥0.7
        │   → vcf/haplotype_matrix.csv, haplotype_calls.csv
        │
        ├── 10_dnds_per_gene.py          ── per-gene dN/dS (mt code table 2)
        │   → vcf/dnds_per_gene.tsv
        │
        ├── 11_haplotypes_nonsyn.py      ── haplotype matrix, NS-only
        │   → vcf/haplotype_matrix_nonsyn.csv
        │
        ├── 12_ns_cooccurrence.py        ── perfect-LD groups among NS sites
        │   → vcf/ns_cooccurrence_groups.tsv
        │
        ├── DP_AD_table.py               ── long-format AD table
        │   → vcf/mtDNA_long_AD_table.tsv
        │
        └── HETEROPLASMY THREAD
              │
              ├── 09_heteroplasmy_report.py    ── Hp on MT_DP_AD_141.txt
              │   → vcf/heteroplasmy_{events,per_site,per_individual,summary}
              │
              ├── 13_pileup_cds_AD.sh          ── pileup over 13 CDS, all panel BAMs
              │   → vcf/pileup_cds_141.vcf.gz
              ├── 14_hp_from_pileup.py         ── Hp on pileup; private_alt_Hp reachable
              │   → vcf/heteroplasmy_pileup_*_all.{tsv,txt}
              ├── 15_well_bleed_test.py        ── plate-contamination test, 4 focals
              │   → vcf/well_bleed_*.{tsv,txt}
              │
              ├── 16_annotate_hp.py            ── join SnpEff SYN/NS onto Hp events
              │   → vcf/heteroplasmy_*_annot.tsv, heteroplasmy_annot_summary.txt
              ├── 17_annotate_hp_codon.py      ── codon-fill Unannotated rows
              │   → vcf/heteroplasmy_*_annot_codon.tsv,
              │     heteroplasmy_annot_codon_summary.txt
              └── 18_variant_burden_per_individual.py
                                              ── per-sample variant counts (Total/SYN/NS)
                                                 + Hp burden, 143 samples, master table
                  → vcf/per_individual_burden_pileup.tsv,
                    per_individual_burden_variant.tsv,
                    per_individual_burden_summary.txt
```

### 08 — Haplotype calling (binary AF threshold)

**Script:** `scripts/08_call_haplotypes.py`
**Input:** `vcf/Fhet_MT_CDS.snps.split.vcf.gz`
**Rule:** `AD_alt / DP > 0.7`, `MIN_DP = 3`, samples `70` / `125`
excluded (run-quality issues; see session-10 CHANGELOG).
**Outputs:** `vcf/haplotype_matrix.csv` (sites × samples binary call),
`vcf/haplotype_calls.csv` (per-sample haplotype string, prefix `C_`).

### 09 — Heteroplasmy report (variant-only input)

**Script:** `scripts/09_heteroplasmy_report.py`
**Input:** `vcf/MT_DP_AD_141.txt` (built outside this repo; variant-only
per-cell DP/AD table).
**Filter:** `0.1 ≤ AD/DP < 0.7` for Hp; `AD/DP ≥ 0.7` for major.
**Classifier:** REF_Hp / shared_alt_Hp / private_alt_Hp. The third
category is structurally 0 here — every (POS, ALT) in the input has
at least one major-allele carrier by construction; see stage 14 for
the proper answer to private-ALT-Hp.
**Outputs:** `vcf/heteroplasmy_{events,per_site,per_individual,summary}.tsv|.txt`.

### 10 — dN/dS per gene

**Script:** `scripts/10_dnds_per_gene.py`
**Inputs:** `Missing_Files/SSM_MT_ref/Fhet_MT.fasta`,
`Missing_Files/SSM_MT_ref/Fhet_MT.gff`, the canonical VCF.
**Method:** Nei-Gojobori site counting under the **vertebrate mitochondrial
genetic code (NCBI table 2)**; strand-aware (ND6 reverse-complemented);
polyA-completed stop codons in ND2/COX2/COX3/ND3/ND4 handled by truncating
to the codon multiple. Observed S/N from the canonical's `ANN[0]`.
**Output:** `vcf/dnds_per_gene.tsv` — one row per gene + Overall.

### 11 — Nonsynonymous haplotype matrix

**Script:** `scripts/11_haplotypes_nonsyn.py`
**Input:** canonical VCF (filtered to `missense_variant` /
`splice_region_variant&missense_variant` rows).
**Otherwise identical** to stage 08 — same 0.7 threshold, same
exclusions. Haplotype strings prefixed `N_` to distinguish.
**Outputs:** `vcf/haplotype_matrix_nonsyn.csv`, `vcf/haplotype_calls_nonsyn.csv`.

### 12 — NS co-occurrence (perfect LD)

**Script:** `scripts/12_ns_cooccurrence.py`
**Input:** `vcf/haplotype_matrix_nonsyn.csv` (preferred; falls back to
the missense subset of `vcf/haplotype_matrix.csv`).
**Identifies:** multi-site groups whose 141-element per-sample call
vectors are identical (perfect LD across the panel).
**Classifies:** `fixed_ref_divergence` (carriers = N_samples — the
3 ND1/ND2 sites where the panel diverges from REF), `haplogroup`
(≥2 carriers — the 7-site major NS clade in 64/141 samples), or
`singleton_artifact` (1 carrier — aggregated private NS variants).
**Output:** `vcf/ns_cooccurrence_groups.tsv`. CLI: `--min-carriers N`.

### 13 — Per-CDS-position pileup (heteroplasmy thread)

**Script:** `scripts/13_pileup_cds_AD.sh`
**Pipeline:** `bcftools mpileup -a AD,DP -d 100000 --no-BAQ` over
`docs/mito_protein_coding.bed` (13 CDS, 11,417 bp) across all 141
panel BAMs (excludes 70, 125 by filename filter). `--no-BAQ` matches
the 05_1 caller so the AD scale is consistent.
**Output:** `vcf/pileup_cds_141.vcf.gz` (+ `.tbi`),
`vcf/slim_bamlist_141.txt`.
**Wall time:** 1–3 min on Mac.

### 14 — Heteroplasmy from pileup (the canonical Hp set)

**Script:** `scripts/14_hp_from_pileup.py`
**Input:** `vcf/pileup_cds_141.vcf.gz` (from stage 13).
**Thresholds:** `DP ≥ 20`, `0.10 ≤ AD/DP < 0.70`, `AD_Hp ≥ 4`.
**Classifier:** same 3-way as stage 09 but `private_alt_Hp` is now
reachable because the input is coverage-honest at every CDS position.
**Outputs:** `vcf/heteroplasmy_pileup_{events,per_site,per_individual,summary}_all.{tsv,txt}`.
**Session-17 result:** 1,124 Hp events / 64 individuals / 400 sites;
**27 `private_alt_Hp` events** (vs 0 from stage 09).

### 15 — Well-bleed contamination test

**Script:** `scripts/15_well_bleed_test.py`
**Inputs:** stage-14 events table; `data_files_May/WGS_seq_plate.txt`
for the plate map (plate identity inferred from the `i5` column).
**Question:** is the high Hp load on individuals 77 / 47 / 33 / 84
explained by library-prep well bleed?
**Method:** donor-concordance score (fraction of X's Hp events
explained by Y's haplotype) × two test statistics (Spearman ρ of
Score vs plate distance; mean(neighbor) − mean(far)), with 10,000-
permutation one-sided p-values.
**Outputs:** `vcf/well_bleed_{donor_ranking,results,summary}.{tsv,txt}`.
**Session-17 result:** no well-bleed signal on any of the 4 focals;
the score-saturation pattern is the signature of panel-level
north × south mt admixture rather than contamination.

### 16 — Coding-effect annotation (SYN / NS) on Hp events

**Script:** `scripts/16_annotate_hp.py`
**Inputs:** canonical VCF (for `ANN[0]`), `vcf/heteroplasmy_events.tsv`
(stage 09), `vcf/heteroplasmy_pileup_events_all.tsv` (stage 14).
**Method:** build `(POS, ALT) → {effect, effect_class, gene, hgvs_p}`
from `ANN[0]`; for each Hp event, select the non-REF allele
(`Major` if `Hp_is_REF` else `Hp_allele`) and look up its annotation.
**Effect_class:** SYN (`synonymous_variant`), NS (`missense_variant`,
`stop_gained`, `stop_lost`, `start_lost`, `initiator_codon_variant` —
same set as scripts 10/11), Other, or Unannotated (POS/ALT not in
canonical).
**Outputs:** `vcf/heteroplasmy_events_annot.tsv` (1,088 rows; +5 cols),
`vcf/heteroplasmy_pileup_events_all_annot.tsv` (1,124 rows; +5 cols),
`vcf/heteroplasmy_annot_summary.txt` (Hp_class × Effect_class cross-tab).
**Session-18 result:** ~89 % SYN / ~8 % NS among annotated events on
both inputs; ~11:1 SYN:NS ratio consistent with purifying selection
on protein-coding mt sites. 27 stage-14 `private_alt_Hp` events are
all `Unannotated` by construction (their minor allele was never called
as a panel variant) — see stage 17.

### 17 — Codon-level fill for Unannotated rows

**Script:** `scripts/17_annotate_hp_codon.py`
**Inputs:** stage-16 annotated TSVs; reference FASTA + GFF.
**Method:** for each row with `Effect_class == "Unannotated"` whose
POS sits inside a CDS, build the REF codon using the same machinery
as stage 10 (strand-aware; − strand complements the substituted
base), substitute the non-REF allele, translate REF and ALT codons
under NCBI table 2, and classify SYN (same AA) / NS (different AA,
including stop_gained/stop_lost). REF-base sanity check guards against
GFF/FASTA/event-table positioning drift.
**New columns:** `Codon_ref`, `Codon_alt`, `AA_ref`, `AA_alt`,
`Annotation_source ∈ {snpeff, codon}`.
**Outputs:** `vcf/heteroplasmy_*_annot_codon.tsv`,
`vcf/heteroplasmy_annot_codon_summary.txt`.
**Session-18 result:** all 70 stage-14 Unannotated rows classified
(12 SYN + 23 NS + 35 Other for non-CDS / `NONE` major calls), zero
REF mismatches. The 35 SYN/NS rows include all 27 `private_alt_Hp`
events; their ~1:2 SYN:NS ratio is notably enriched for NS vs the
canonical 11:1 pattern.

### 18 — Per-individual variant + Hp burden (bimodality test)

**Script:** `scripts/18_variant_burden_per_individual.py`
**Inputs:** canonical VCF (for per-sample GT counts);
`vcf/heteroplasmy_pileup_per_individual_all.tsv` (stage 14);
`vcf/heteroplasmy_per_individual.tsv` (stage 09).
**Method:** for each sample, count `GT == "1"` calls in the canonical
(ploidy=1; samples carry one allele per site), classify each row by
`ANN[0]` (SYN / NS / Other). Join with per-individual Hp counts,
**zero-filled for the ~80 samples not in `per_individual_all`** — needed
because samples with zero Hp aren't in that file and dropping them
would bias any group comparison.
**Sample-name normalization:** extract leading numeric ID, so it joins
stage-09 (`{N}_MT`) and stage-14 (`MT_only_bams/{N}_0_MT_only.bam`)
on the same key.
**Outputs (143 rows × 9 cols each):**
- `vcf/per_individual_burden_pileup.tsv` (variant counts + stage-14 Hp)
- `vcf/per_individual_burden_variant.tsv` (variant counts + stage-09 Hp)
- `vcf/per_individual_burden_summary.txt` (low / mid / high variant ×
  Hp summary)

**Designed to test:** the bimodal variant-count distribution
(north < ~50 ALTs vs south > ~190) vs per-individual Hp burden.
**Session-18 result:** bimodality confirmed (77 samples at 15–26,
65 samples at 220–233, 1 in between); low-variant samples carry
~1.8× mean Hp (7.65 vs 4.34) and a higher fraction (51 % vs 37 %)
have any Hp. The 4 admixed focals (47, 77 high-variant; 33, 84
low-variant) sit at both extremes and dominate the Hp signal.
Run Mann-Whitney U / Fisher exact (`scipy.stats` in `SNP_env`) on
the analysis-ready tables for formal tests; run with and without
the 4 focals to separate population-level from outlier-driven
effects.

### Dependencies among the downstream stages

- 09 needs `MT_DP_AD_141.txt` (built outside the repo).
- 13 → 14 → 15 is a chain (pileup → Hp → contamination test).
- 16 needs the canonical (for ANN) and the events tables from 09 +/or 14.
- 17 needs 16's `*_annot.tsv` outputs.
- 18 needs the canonical (for per-sample GT) and the per-individual
  files from 09 +/or 14.
- 11 → 12 (12 can fall back to filtering 08's matrix if 11 hasn't run).
- 10 is standalone.
- 08 is standalone, runs in parallel with the heteroplasmy thread.

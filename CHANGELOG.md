# CHANGELOG

Session-by-session record of what changed.

## 2026-06-13 (session 20) — Clade-split π; manuscript pipeline references + supplemental pipeline table; caller provenance resolved from VCF headers; panel-wide MAPQ quantified; GitHub catch-up planned

> **TL;DR.** Manuscript-support session. Added `scripts/20_calc_pi_clade.py`
> (π split by North/South clade). Resolved the long-standing "which caller
> produced the 927 variants?" question **from the VCF header provenance**:
> it was JOINT calling (`jobs/05_1_mpileup_merge.sh`), not per-sample → merge.
> Quantified panel-wide read mapping quality (99.99% MAPQ ≥ 30) to show the
> intended-but-omitted `-q 30` filter is immaterial and to add a 4th anti-NUMT
> argument. Produced two manuscript deliverables with tracked pipeline
> references. Discovered that most analysis scripts (09–28) are **untracked in
> git** — GitHub release plan below.

### `scripts/20_calc_pi_clade.py` (new)

Adapted from `19_calc_pi.py`. Pass 1 counts confident ALT calls per sample
and assigns clades (north < 50 ALT, south > 200 ALT, in-between excluded —
by default only `77_MT` at 193). Pass 2 computes π (total/syn/NS) **within
each clade separately** — reference-independent, so clade-defining fixed
differences contribute π ≈ 0. Verified (pysam unavailable in the working
sandbox, so cross-checked with an independent gzip-based reimplementation):
77 north / 63 south; π_total north 0.00160, south 0.00188; pN/pS north 0.21,
south 0.26. Output: `vcf/pi_by_clade_persite.tsv` (+ `.membership.tsv`).
Caveat: π_syn/π_ns use full L_CDS as denominator (same as 19); pass
`--L_syn`/`--L_ns` for true per-site rates.

### Caller provenance — RESOLVED from VCF headers (method worth keeping)

The manuscript Methods described per-sample calling → `bcftools merge`, but
the repo's canonical caller `05_1_mpileup_merge.sh` is JOINT. Resolved by
reading the ordered `##bcftools_*Command` / `##bcftools_*Version` lines
embedded in each VCF header (the authoritative record of what actually ran,
more reliable than file dates or script intent):

  - `141_MT_variants.vcf.gz` (927 records) header chain:
    `call -mv --ploidy 1` (1.23.1, 18 May 2026) → `norm -m -any -f REF` →
    SnpEff view → CDS-restrict (`view -R mito_protein_coding.bed`) →
    `view -v snps` → sample subset → `view -e AC=0`.
    **A single `call`, no per-sample paths, no `merge` step → JOINT.**
  - archived `Fhet_mt_persample_merged.vcf.gz` (15 May 2026) header carries
    `call ... -o persample/10_0.vcf.gz` + `merge -m none -l ...` → that is
    the per-sample chain, and it was **not** used for the results.

Two side corrections found the same way: the executed `call` was
`-mv --ploidy 1` with **no `-A`** (manuscript said `-mv -A --ploidy 1`),
and `mpileup` used defaults **`-Q 13 -q 0`** (manuscript said `-Q 30 -q 30`).
The `mpileup` command line was not retained in the header; confirmed against
`jobs/05_1_mpileup_merge.sh` (which uses `-a AD,DP --max-depth 100000`, no
`-Q`/`-q`).

### Panel-wide MAPQ (all 143 mt BAMs, 36.27M reads)

Parsed MAPQ directly from the slim BAMs (samtools/pysam unavailable; read
column-5 MAPQ from the BGZF/BAM byte stream in pure Python). Result:
**99.99% MAPQ ≥ 30, 99.69% MAPQ = 60 (BWA unique max), 0.0009% MAPQ = 0**;
only 0.0077% (2,802 reads) would be dropped by `-q 30`. Worst single sample
< 0.04%. **Individual 77** (lowest reads, 14,038): 99.25% MAPQ 60, **zero**
MAPQ-0 — its problem is low depth, not mapping ambiguity / NUMT. Per-sample
counts cached at `outputs/mapq_cache.tsv` (outside the repo).

### Manuscript deliverables (in repo root, not tracked yet)

  - `Mitochondrial_variants_5_pipeline-refs.docx` — v5 with tracked-change
    pipeline references (jobs 01–04 + slim-BAM; SnpEff/CDS scripts 06/07;
    dN/dS script 10; rewritten joint-calling Methods paragraph; corrected
    `mpileup` flags + MAPQ sentence; expanded NUMT argument to 4 reasons;
    new **Code and Data Availability** section).
  - `Supplemental_Table_S1_Pipeline.docx` — per-step table (file → platform →
    purpose/params → tool+version → output/manuscript element), HPC vs Mac,
    plus notes on the bcftools-1.6 issue, environments, and caller provenance.

### GitHub status / release plan (see CLAUDE.md "session 20" pickup)

`git status` shows the repo is far behind: last commit is the 05d/05e
per-sample era, and **scripts 09–28, `05_1_mpileup_merge.sh`,
`BSUB_Slim_BAM_mt.sh`, `make_141_vcf.*`, and `jobs/_archive/` are all
untracked**. Stray manuscript `.docx` drafts sit in `scripts/`. Remote is
`git@github.com:DLCrawford/MT_Genomics.git`. Plan: add `*.docx` to
`.gitignore` (or keep manuscripts out of the repo), remove stray docx from
`scripts/`, stage all jobs/scripts/docs, commit, push, tag a release, and
mint a Zenodo DOI for the Code-Availability statement. Manuscript URL
corrected to `https://github.com/DLCrawford/MT_Genomics`.

## 2026-06-05 (session 19) — Cross-species π comparison; canonical 141-sample VCF; π scripts for Fhet, Drosophila, Human, C. elegans, Yeast (in progress)

> **TL;DR.** Major new analytical thread: nucleotide diversity (π) calculated
> for Fhet and four comparison datasets to contextualise Fhet mt diversity.
> Canonical 141-sample VCF produced. Scripts 19–24 written. Human (1,176
> Lankheet 2026 genomes) and C. elegans (540 samples) π complete; yeast
> (469 complete assemblies) in progress.

### Canonical 141-sample VCF

- **`vcf/141_MT_variants.vcf.gz`** produced by `scripts/make_141_vcf.sh`:
  - Input: `vcf/Fhet_MT_CDS.snps.split.vcf.gz` (143 samples, 930 records)
  - Removed samples 70_MT and 125_MT (`bcftools view --samples-file`)
  - Removed 3 sites private to 70/125 that became monomorphic (`bcftools view -e 'AC=0'`)
  - Output: **141 samples, 927 records** — matches `haplotype_matrixDP10.csv` exactly
  - Indexed with `bcftools index --tbi`
- The 9 invariant fixed-ALT sites (6 SYN + 3 NS reference-divergence) are retained
  in the VCF; they contribute π = 0 and are handled at the analysis level.

### π calculation scripts

**`scripts/19_calc_pi.py`** — Fhet mt CDS π (141 samples, L_CDS=11,417)
- Input: `vcf/141_MT_variants.vcf.gz`
- Haploid AD-based calling (threshold 0.7), same logic as scripts 08/11/14
- Classifies syn/NS from SnpEff ANN field; handles `&`-joined effects
- Fixed: pysam returns ANN as tuple not string; `&`-joined effects (2 syn
  sites were missed before fix — now reports 927 CDS sites matching matrix)
- Output: `vcf/pi_results.tsv` (per-site) + stdout summary

**`scripts/20_dros_pi.py`** — Drosophila DGRP mt π (169 lines, L_CDS=11,173)
- Input: `data_files_May/Dros_Mt_coding.csv`
- Genotype matrix: 0=REF, ≥1=ALT (values >1 treated as ALT — 43 rows affected)
- 36 rows with NaN FuncType excluded from syn/NS (included in π_total)
- 6 MNP rows excluded
- Output: `vcf/dros_pi_results.tsv`
- Results: π_total=0.00052, π_syn=0.00039, π_ns=0.00013, pN/pS=0.339
- Notable: pN/pS from π (0.339) << θ_W pN/pS (0.775) — classic purifying
  selection signature (excess rare NS variants at low frequency)

**`Human_mt/21_human_mt_pi.py`** — Human mt whole-genome π (1,176 genomes)
- Downloads 1,176 complete mt genomes from Lankheet et al. 2026 (Commun Biol)
  via NCBI efetch (accessions PV558957–PV560130, PX394655–PX394656)
- Reference-based: aligns each genome to rCRS (NC_012920) using Biopython
  PairwiseAligner; records base at each rCRS position
- Output: `Human_mt/human_mt_pi_per_site.tsv`, `human_mt_pi_summary.txt`
- Download script: `Human_mt/download_lankheet.py` (POST-based efetch,
  batch size 200; overcame URL-length limit that caused GET failures)

**`Human_mt/22_human_mt_cds_pi.py`** — Human mt CDS π (syn/NS)
- Input: `Human_mt/human_mt_pi_per_site.tsv` (from script 21)
- Hardcoded rCRS CDS coordinates for 13 protein-coding genes; L_CDS=11,395 bp
- Vertebrate mt genetic code (NCBI table 2): TGA=Trp, ATA=Met, AGA/AGG=Stop
- Classifies each variable CDS site as SYN/NS by substituting minor allele
  into rCRS codon context; handles ND6 (minus strand, reverse complement)
- Output: `Human_mt/human_mt_cds_pi_per_site.tsv`, `human_mt_cds_pi_summary.txt`

**`C_elegans/23_celegans_pi.py`** — C. elegans mt CDS π (540 samples, L_CDS=10,299)
- Input: `~/Projects/MT_Genomics_Cl_Ap2026/C_elegans/C_elegansFINAL.annotated.vcf`
  (GATK diploid calls, SnpEff annotated, 1,449 records including indels)
- Diploid GT parsing: 0/0→REF, 1/1→ALT, 0/1→skip (selfing organism, ~99.9%)
- PASS filter + SNPs only (indels/MNP excluded)
- Same SnpEff ANN classification as script 19
- Output: `C_elegans/celegans_pi_per_site.tsv`, `celegans_pi_summary.txt`

**`Yeast/24_yeast_pi.py`** — Yeast mt CDS π (469 complete assemblies) — IN PROGRESS
- Input: `~/Projects/MT_Genomics_Cl_Ap2026/Yeast/mitochondrialAssemblies/`
  (905 FASTA files: 303 oneScaff, 166 circularized, 436 multiscaff)
- Using only complete assemblies (oneScaff + circularized = 469)
- Reference: S288C mt genome NC_001224 (~85 kb); CDS L=6,684 bp
- Biopython PairwiseAligner (MAFFT unavailable — conda conflict, no brew)
- Script written; not yet run (runtime estimate: several hours)

### Conceptual notes added this session

- **Watterson's θ vs π:** θ_W counts all segregating sites equally (unweighted);
  π is frequency-weighted. Gap between the two estimates indicates skewed
  frequency spectrum — excess of rare variants = purifying selection signature.
- **Harmonic number a₁:** = Σ(1/i) for i=1 to n-1; grows logarithmically.
  NOT the harmonic mean (which = (n-1)/a₁).
- **pN/pS from π vs θ_W:** Drosophila π-based pN/pS (0.339) << θ_W-based
  (0.775) — confirms purifying selection acting on NS variants.

### Artifacts changed this session

- `scripts/19_calc_pi.py` — new
- `scripts/20_dros_pi.py` — new
- `scripts/make_141_vcf.sh` — new (iterated through 4 fixes to reach 141/927)
- `scripts/make_141_vcf.py` — new (Python alternative, not used)
- `vcf/141_MT_variants.vcf.gz` + `.tbi` — new canonical 141-sample VCF
- `Human_mt/download_lankheet.py` — new
- `Human_mt/21_human_mt_pi.py` — new
- `Human_mt/22_human_mt_cds_pi.py` — new
- `C_elegans/23_celegans_pi.py` — new
- `Yeast/24_yeast_pi.py` — new (not yet run)
- `CLAUDE.md` — new session-19 pickup block with comparative π status table
- `CHANGELOG.md` — this entry

## 2026-06-06 (session 19, continued-2) — Both human datasets added; AMR θ_W corrected; comparison table updated

> **TL;DR.** AMR θ_W was over-inflated by ~2.6× (counting all gnomAD CDS
> sites rather than AMR-segregating sites only). Fixed by filtering to
> AC_hom_amr > 0 for S counting. Lankheet 2026 (African) θ_W added to
> `22_human_mt_cds_pi.py`. Both human datasets now in comparison table.
> `25_comparison_table.py` updated to handle 6 datasets with per-dataset
> S-filtering logic.

### Bug fix — AMR θ_W over-inflation

`26_amr_pi.py` originally counted all PASS CDS SNPs in gnomAD (any
population) as segregating for AMR. This inflated S from the correct ~2,849
to ~7,461 and θ_W from 0.0271 to 0.07096. Fixed: S now counted only where
`AC_hom_amr > 0`. π was unaffected (sites with AC_hom=0 already contributed
π=0).

### Lankheet 2026 θ_W added to `22_human_mt_cds_pi.py`

Script 22 previously reported only π. Added Watterson's θ calculation
(N=1,176, a₁=7.6467) to the summary output. Results:

| | S | θ_W | π |
|--|---|-----|---|
| Total | 1,624 | 0.01864 | 0.00685 |
| Syn | 1,182 | 0.01357 | 0.00544 |
| NS | 442 | 0.00507 | 0.00140 |
| pN/pS | | 0.374 | 0.258 |

### Both human datasets now in comparison table

`25_comparison_table.py` updated:
- Replaced single "Human" entry with two separate entries:
  - **Human African (Lankheet 2026)**: N=1,176, TSV from script 22
  - **Human AMR (gnomAD v3.1)**: N=5,718, TSV from script 26 with
    `s_filter=("AC_hom_amr", 0)` to restrict S to AMR-segregating sites
- Added `s_filter` parameter to DATASETS so any dataset can specify a
  minimum-value column filter for S counting
- Yeast handled via hardcoded literature values (De Chiara 2020)
- Added footnotes to formatted table noting data provenance

### Run to regenerate

```bash
cd ~/Projects/MT_Genomics_Cl_Ap2026/MT_Genomics2
conda activate SNP_env
python scripts/25_comparison_table.py
```

### Artifacts changed

- `scripts/25_comparison_table.py` — updated (both human datasets, s_filter)
- `scripts/22_human_mt_cds_pi.py` — updated (added θ_W to summary)
- `Hm_Mt/26_amr_pi.py` — bug fix (S counts AMR-segregating only)
- `Human_mt/human_mt_cds_pi_summary.txt` — rerun to include θ_W
- `vcf/comparison_table.tsv` + `vcf/comparison_table.txt` — regenerated
- `CHANGELOG.md` — this entry

## 2026-06-06 (session 19, continued) — Comparison table complete; yeast π flagged unreliable

> **TL;DR.** Cross-species π + θ_W comparison table finalised for 4 of 5
> datasets. Yeast π from script 24 is unreliable (~10-18× inflated) due to
> spurious local alignments on high AT-content sequences — flagged for rework.
> Key biological finding: π pN/pS < θ_W pN/pS in all four species, consistent
> with purifying selection skewing NS variants to low frequencies.

### Final comparison table (vcf/comparison_table.tsv / .txt)

| Dataset | N | S | S_syn | S_ns | θ_total | θ_syn | θ_ns | pN/pS_θ | π_total | π_syn | π_ns | pN/pS_π |
|---------|---|---|-------|------|---------|-------|------|---------|---------|-------|------|---------|
| Fhet | 141 | 918 | 757 | 161 | 0.0146 | 0.0120 | 0.0026 | 0.213 | 0.0105 | 0.0097 | 0.0008 | 0.086 |
| Drosophila | 169 | 181 | 102 | 79 | 0.0028 | 0.0016 | 0.0012 | 0.775 | 0.0005 | 0.0004 | 0.0001 | 0.339 |
| Yeast | 1,011 | 384 | N/A | N/A | 0.0077 | N/A | N/A | N/A | — | — | — | — |
| C. elegans | 540 | 1,274 | 963 | 311 | 0.0180 | 0.0136 | 0.0044 | 0.323 | 0.0103 | 0.0093 | 0.0011 | 0.113 |
| Human AMR* | 5,718 | 2,849 | 1,914 | 935 | 0.0271 | 0.0182 | 0.0089 | 0.489 | 0.0069 | 0.0054 | 0.0014 | 0.258 |

*Human θ_W from AMR (N=5,718, Table 2 literature); π from Lankheet 2026 (N=1,176).
 This mixed provenance should be noted in any publication.

### Biological interpretation

**Universal pattern:** π pN/pS < θ_W pN/pS in all four species — classic
signature of purifying selection skewing NS variants to low frequencies.
Frequency-weighted π depresses NS signal more than site-count-based θ_W.

**Purifying selection strength (fold drop θ_W → π pN/pS):**
- Fhet: 0.213 → 0.086 (2.5×) — strongest signal; large outcrossing population
- C. elegans: 0.323 → 0.113 (2.9×) — strong despite predominant selfing
- Human: 0.489 → 0.258 (1.9×)
- Drosophila: 0.775 → 0.339 (2.3×) — θ_W near-neutral but π shows clear selection

**Drosophila anomaly:** θ_W pN/pS = 0.775 (nearly neutral) collapses to π
pN/pS = 0.339, confirming that many NS variants in DGRP segregate at very
low frequency — the hallmark of slightly deleterious mutations under
purifying selection.

**Fhet has lowest π pN/pS (0.086)** — strongest evidence of purifying
selection among the four species with complete syn/NS data.

### Yeast π — flagged unreliable

Script 24 returned π = 0.14034 and θ_W = 0.08489 with S = 3,817 variable
sites (expected ~384 from literature). Root cause: Biopython local aligner
on ~85 kb high AT-content assemblies produces spurious matches; 861/3,752
gene-assembly alignments failed outright and many passing alignments likely
mapped to wrong regions. Literature θ_W (0.00766, N=1,011) should be used
as placeholder until a stricter alignment approach is implemented (e.g.,
minimum percent-identity filter, k-mer pre-screening, or minimap2).

### Open items carried forward

- Yeast π: fix script 24 alignment (k-mer pre-screen or minimap2)
- Human θ_W: recalculate from Lankheet 2026 at N=1,176 for consistency
  (expected a₁ ≈ 7.646; current table uses AMR N=5,718)
- Fhet π: confirm final values using `vcf/141_MT_variants.vcf.gz`
- Publish-ready table: resolve Human AMR vs Lankheet provenance note

### Artifacts changed

- `scripts/24_yeast_pi.py` — two bug fixes: (1) filter to 8 canonical CDS
  genes only (excluded intron ORFs AI*, BI*, SCEI, Q0255); (2) replaced
  `if not alignments:` with `next(iter())` to fix Biopython OverflowError
- `scripts/25_comparison_table.py` — new; reads all per-site TSVs, computes
  π + θ_W + pN/pS for all datasets; outputs TSV and formatted text table
- `vcf/comparison_table.tsv` + `vcf/comparison_table.txt` — new
- `Yeast/yeast_pi_per_site.tsv` — produced but values unreliable (see above)
- `CLAUDE.md` — session-19 pickup block updated with final table and yeast caveat
- `CHANGELOG.md` — this entry

## 2026-05-24 (session 18) — Coding-effect annotation (SYN / NS) added to heteroplasmy event tables

> **TL;DR.** New `scripts/16_annotate_hp.py` joins SnpEff coding-effect
> classification onto the stage 09 (variant-only) and stage 14 (CDS
> pileup) heteroplasmy event tables. Each Hp event now carries an
> `Effect_class` column with values `SYN` / `NS` / `Other` /
> `Unannotated`. SYN and NS are roughly an order of magnitude apart
> on this panel — among annotated CDS heteroplasmies, **~89 % SYN
> vs ~8 % NS** for both stages — consistent with purifying selection
> on protein-coding mt sites. All 27 `private_alt_Hp` events from
> stage 14 land in `Unannotated` by construction (they are
> (POS, ALT) combos the panel caller never produced); computing
> SYN/NS for these would require the codon machinery from stage 10.

### What's new

- **`scripts/16_annotate_hp.py`** — standalone join script. Reads
  `vcf/Fhet_MT_CDS.snps.split.vcf.gz` and builds a `(POS, ALT)
  → {effect, effect_class, gene, hgvs_p}` map from `ANN[0]`
  (the highest-impact gene-internal annotation per SnpEff
  convention). For each Hp event in the two events tables it
  selects the non-REF allele — `Major` when `Hp_is_REF == True`
  (REF_Hp class), otherwise `Hp_allele` (shared_*_Hp /
  private_alt_Hp classes) — looks up the annotation, and writes
  five new columns: `non_REF_allele`, `Effect_class`, `Effect`,
  `Gene`, `HGVS_p`.

  Effect_class buckets match the impact set used by scripts 10
  and 11:

  | Class       | SnpEff effects                                                                              |
  |-------------|---------------------------------------------------------------------------------------------|
  | SYN         | `synonymous_variant` (with any `&`-joined variants)                                         |
  | NS          | `missense_variant`, `stop_gained`, `stop_lost`, `start_lost`, `initiator_codon_variant`     |
  | Other       | anything else SnpEff returned (expected 0 on a CDS-restricted canonical, kept for safety)   |
  | Unannotated | (POS, non-REF allele) not in the canonical VCF                                              |

  Standalone — no `cyvcf2` / `bcftools` dependency. Parses the
  gzipped VCF text directly. ~1 s on the current canonical.

### Results

ANN map from `Fhet_MT_CDS.snps.split.vcf.gz`: **930 unique
(POS, ALT) entries** — 762 SYN, 166 NS, 2 Other.

Stage 09 (variant-only, `vcf/heteroplasmy_events.tsv`, 1,088 events):

| Hp_class            | SYN | NS  | Other | Unannotated | Total |
|---------------------|----:|----:|------:|------------:|------:|
| REF_Hp              | 525 |  55 |     4 |          29 |   613 |
| shared_with_major   | 437 |  34 |     4 |           0 |   475 |
| **Total**           | 962 |  89 |     8 |          29 | 1,088 |

Stage 14 (CDS pileup, `vcf/heteroplasmy_pileup_events_all.tsv`, 1,124 events):

| Hp_class         | SYN | NS  | Other | Unannotated | Total |
|------------------|----:|----:|------:|------------:|------:|
| REF_Hp           | 511 |  53 |     3 |          43 |   610 |
| private_alt_Hp   |   0 |   0 |     0 |          27 |    27 |
| shared_alt_Hp    | 450 |  33 |     4 |           0 |   487 |
| **Total**        | 961 |  86 |     7 |          70 | 1,124 |

Headline ratio: ~11:1 SYN:NS among annotated events on both
inputs. Among the 89 NS events in stage 09 and the 86 in stage 14,
the per-site mapping (`Gene`, `HGVS_p`) is in the new columns —
these are the amino-acid-changing heteroplasmies, the rows most
likely to matter functionally.

The 27 `private_alt_Hp` events from stage 14 are all `Unannotated`
by construction (they are minor alleles the panel caller never
produced as variant, so no `(POS, ALT)` row exists in the
canonical). To classify these as SYN/NS, port the codon /
translation-table-2 site-counting logic from
`scripts/10_dnds_per_gene.py` and look up the codon containing
each (POS, ALT). Open as a follow-up.

The 43 `REF_Hp` Unannotated rows in stage 14 are positions where
the pileup-major is either `NONE` (insufficient evidence to make
a major call) or an allele the panel caller never produced.

### Spot-check (verified row-level)

POS 2847 in ND1 has two ALTs in the canonical: `T>C`
(`synonymous_variant`, LOW, `p.Phe2Phe`) and `T>A`
(`missense_variant`, MODERATE, `p.Phe2Leu`). Five stage 14 Hp
events at POS 2847 all involve the `C` allele as the non-REF
choice and are correctly tagged `SYN / ND1 / p.Phe2Phe`. Two
class flavors are represented in the spot-check: `shared_alt_Hp`
(sample's Hp is C while major is REF) and `REF_Hp` (sample's
major call is C and Hp is the REF T).

### Outputs written

- `vcf/heteroplasmy_events_annot.tsv`              (1,088 rows; +5 columns)
- `vcf/heteroplasmy_pileup_events_all_annot.tsv`   (1,124 rows; +5 columns)
- `vcf/heteroplasmy_annot_summary.txt`             (Hp_class × Effect_class cross-tabs)

### Follow-up — `scripts/17_annotate_hp_codon.py` (same session)

Added immediately after script 16. Fills in SYN / NS for the
`Unannotated` rows by computing the effect from first principles
using the reference FASTA, GFF (strand-aware, reverse-complements
ND6), and NCBI translation table 2 — the same machinery as
`scripts/10_dnds_per_gene.py`.

Inputs: the `_annot.tsv` files from script 16. Outputs:
`vcf/heteroplasmy_events_annot_codon.tsv`,
`vcf/heteroplasmy_pileup_events_all_annot_codon.tsv`,
`vcf/heteroplasmy_annot_codon_summary.txt`.

New columns on top of script 16's output:
`Codon_ref`, `Codon_alt`, `AA_ref`, `AA_alt`,
`Annotation_source ∈ {snpeff, codon}`.

**Results — every Unannotated row is now classified:**

| Source                       | Stage 09 | Stage 14 |
|------------------------------|---------:|---------:|
| Unannotated after script 16  |       29 |       70 |
| SYN (codon-derived)          |        0 |       12 |
| NS (codon-derived)           |        0 |       23 |
| Other (non-CDS or NONE alt)  |       29 |       35 |
| REF mismatch                 |        0 |        0 |
| Unannotated after script 17  |        0 |        0 |

The 29 stage-09 Unannotated rows are all `non_CDS` — stage 09's
input (`MT_DP_AD_141.txt`) includes positions outside the
13-CDS region, which the canonical CDS-restricted VCF excluded.
None of them can be SYN/NS because they're not in a coding
region.

The 70 stage-14 Unannotated rows split into:
- **35 → SYN/NS** (12 SYN + 23 NS). These are CDS heteroplasmies
  the panel caller never produced as a variant — including all
  27 `private_alt_Hp` events plus 8 `REF_Hp` events at non-canonical
  (POS, ALT). The ~1:2 SYN:NS ratio in this rare-allele subset
  is notably higher in NS than the canonical's ~11:1 SYN:NS,
  which is consistent with these being either (i) new mutations
  not yet at population equilibrium under purifying selection,
  or (ii) sequencing artifacts. Worth a per-row look — the
  output table now has `Gene` and `HGVS_p` populated so the 23
  NS rows can be scanned at-a-glance for whether they cluster
  by individual (somatic candidate) or by position (recurrent
  artifact).
- **35 → Other** (`invalid_alt_allele` or non-CDS). 34 of these
  are `REF_Hp` events where the pileup's `Major` call is `NONE`
  (insufficient evidence — these can't be classified without a
  valid non-REF allele to substitute). One row at a non-CDS POS.
- **0 → REF mismatch.** The codon's nucleotide at the substituted
  position matched the event-table REF for every row, confirming
  GFF / FASTA / event-table positioning is internally consistent.

### Spot-check (verified row-level)

POS 3036, sample 47, REF=A, non-REF allele=G → script 17 reports
codon `GAA→GAG`, AA `E→E`, gene ND1, `p.E65E`, SYN. Confirmed
against the FASTA: ND1 starts at 2842 (+strand); positions
3034–3036 spell GAA (codon 65 = Glu). Substituting position 3
(genomic 3036) from A to G gives GAG (still Glu). ✓

### Follow-up — `scripts/18_variant_burden_per_individual.py` (same session)

Tests Doug's question: "the variant count is bimodal (<50 or >190);
do the low-variant samples carry significantly more Hp?" Builds a
master per-individual table joining per-sample variant counts
(from the canonical's per-sample GT field, ploidy=1) with
per-individual Hp counts. All **143 canonical samples** appear in
the output — necessary because samples with zero Hp aren't in
`per_individual_all` and dropping them would bias any group
comparison.

Per-sample variant counts come from reading `GT == "1"` in
`Fhet_MT_CDS.snps.split.vcf.gz` and classifying each variant row
by ANN[0] (same SYN / NS / Other buckets as scripts 16 / 17).
Per-sample Hp counts come from
`heteroplasmy_pileup_per_individual_all.tsv` (stage 14) and
`heteroplasmy_per_individual.tsv` (stage 09).

**Outputs (analysis-ready):**

- `vcf/per_individual_burden_pileup.tsv` — 143 rows × 9 cols:
  `id, Individual, n_variants_total, n_variants_SYN,
  n_variants_NS, n_variants_Other, n_hp_sites, n_hp_events, n_ref_hp`.
- `vcf/per_individual_burden_variant.tsv` — same shape, stage-09 Hp.
- `vcf/per_individual_burden_summary.txt` — low / mid / high
  variant × Hp summary.

**Variant-count distribution.** Sharply bimodal, exactly as
predicted. Sorted, every 10th sample: 15, 18, 19, 20, 21, 23, 23,
26, **190**, 220, 222, 223, 224, 225, 226, 229, 233. The gap
between 26 and 190 is empty; 190 is a single admixed focal (77);
the cluster at 220–233 is the south-major group.

  | Bucket                        |  n  | n_variants mean | n_hp_sites mean | n_with_hp |
  |-------------------------------|----:|----------------:|----------------:|----------:|
  | low  (<50 variants, "north")  |  77 |            21.1 |            7.65 |  39 (51%) |
  | mid  (50–190 variants)        |   1 |           190.0 |           216.0 |   1       |
  | high (>190 variants, "south") |  65 |           224.2 |            4.34 |  24 (37%) |

**Direction consistent with the hypothesis** — low-variant
(north-major) samples carry more Hp on average and a larger
fraction have any Hp. The signal is dominated by the 4 admixed
focals, which all show ~190–220 Hp events but split between the
low and high variant buckets:

  | id | n_variants_total | n_hp_sites | n_ref_hp | Pattern (per session 17)             |
  |----|-----------------:|-----------:|---------:|--------------------------------------|
  | 47 |              218 |        219 |      211 | south-major + north minor (REF_Hp)   |
  | 77 |              190 |        216 |      210 | south-major + north minor (REF_Hp)   |
  | 33 |               16 |        215 |       12 | north-major + south minor (alt_Hp)   |
  | 84 |               19 |        187 |       10 | north-major + south minor (alt_Hp)   |

For a formal test, run Mann-Whitney U on `n_hp_sites` (low vs high)
and Fisher exact on the proportion with any Hp from the
analysis-ready table; `scipy.stats.mannwhitneyu` /
`scipy.stats.fisher_exact` in `SNP_env`. The two ends of the
variant-count axis aren't quite comparable populations (north has
fewer variants by definition), so the 2×2 of any-Hp-yes/no by
low/high group is probably the cleanest framing.

**Per-sample NS counts (potential per-sample dN/dS shape).**
Sample 47: 201 SYN, 15 NS (NS / total = 7.4 %). Sample 77: 177 SYN,
12 NS (6.4 %). Low-variant samples: 33 has 12 SYN, 4 NS (25 % NS);
84 has 16 SYN, 3 NS (16 % NS). The ratio suggests low-variant
samples carry proportionally more NS substitutions vs REF, but the
absolute NS counts are tiny (3–5 vs 12–16) so per-sample
proportions are noisy; worth a per-gene aggregate look against
`scripts/10_dnds_per_gene.py`'s table.

### Carry-overs

- Tabular review of the 27 `private_alt_Hp` rows for
  per-individual or per-POS clustering (carry-over from session
  17) — now easier to scan with the gene / HGVS_p columns
  populated for non-private rows so the private set stands out
  visually.
- Run Mann-Whitney U on `n_hp_sites` (low <50 vs high >190) and
  Fisher exact on `n_with_hp` from
  `per_individual_burden_pileup.tsv`; both tests need scipy which
  is in `SNP_env`. Decide whether to include or exclude the 4
  admixed focals — they sit at both ends of the variant-count axis
  and dominate the Hp signal in either direction.

### Artifacts changed this session

- `scripts/16_annotate_hp.py` — new (RAN; outputs in `vcf/`).
- `scripts/17_annotate_hp_codon.py` — new (RAN; outputs in `vcf/`).
- `scripts/18_variant_burden_per_individual.py` — new (RAN; outputs in `vcf/`).
- `vcf/heteroplasmy_events_annot.tsv` — new (from 16).
- `vcf/heteroplasmy_pileup_events_all_annot.tsv` — new (from 16).
- `vcf/heteroplasmy_annot_summary.txt` — new (from 16).
- `vcf/heteroplasmy_events_annot_codon.tsv` — new (from 17).
- `vcf/heteroplasmy_pileup_events_all_annot_codon.tsv` — new (from 17).
- `vcf/heteroplasmy_annot_codon_summary.txt` — new (from 17).
- `vcf/per_individual_burden_pileup.tsv` — new (from 18).
- `vcf/per_individual_burden_variant.tsv` — new (from 18).
- `vcf/per_individual_burden_summary.txt` — new (from 18).
- `CLAUDE.md` — "Downstream analysis scripts" catalog: new
  entries for `scripts/16_annotate_hp.py`,
  `scripts/17_annotate_hp_codon.py`, and
  `scripts/18_variant_burden_per_individual.py` (placed before
  `DP_AD_table.py`).
- `CHANGELOG.md` — this entry.

## 2026-05-21 (session 17) — Stages 13/14/15 RAN. Well-bleed hypothesis rejected for all 4 focals. Admixed-mito interpretation supported by Doug's north/south context.

> **TL;DR.** The three scripts session 16 left staged (13 / 14 / 15)
> ran end-to-end. Stage 14 emits **1,124 Hp events across 64
> individuals at 400 CDS sites**; the pileup-based architecture
> recovers 27 `private_alt_Hp` events — the category 09's
> variant-only input could not produce. Stage 15 ran the well-bleed
> permutation test on the 4 high-Hp focals (77, 47, 33, 84) — **no
> bleed signal on any of them** (all 8 statistics non-significant;
> top-donor distance ≥ 4 wells for every focal). Combined with
> Doug's north/south haplotype framing, the **working
> interpretation is panel admixture**: each focal carries two
> mitochondrial haplotypes and the "Hp" allele is the *other*
> haplotype's major call leaking through the 0.7 threshold. Two
> script bugs in stage 15 were fixed along the way (NameError at
> line 267; new `--events` flag with `_all.tsv` fallback so 15
> finds 14's output regardless of naming convention).

### Stage 14 — pileup Hp detection (`scripts/14_hp_from_pileup.py`)

Ran on `vcf/pileup_cds_141.vcf.gz` (stage 13's pileup over the 13
CDS intervals, 141 panel BAMs). Thresholds: `DP ≥ 20`,
`0.10 ≤ AD/DP < 0.70`, `AD_Hp ≥ 4`.

Outputs (saved with `_all` suffix — see note below):
- `vcf/heteroplasmy_pileup_events_all.tsv`     (1,124 rows)
- `vcf/heteroplasmy_pileup_per_site_all.tsv`     (400 rows)
- `vcf/heteroplasmy_pileup_per_individual_all.tsv` (64 rows)
- `vcf/heteroplasmy_pileup_summary_all.txt`

**Headline numbers vs the variant-only stage-09 baseline:**

| Quantity                                    | Stage 09 | Stage 14 |
|---------------------------------------------|---------:|---------:|
| Total Hp events                             |    1,088 |    1,124 |
| Individuals with ≥1 Hp                      |       53 |       64 |
| Sites with ≥1 Hp                            |      376 |      400 |
| `REF_Hp` events                             |      613 |      610 |
| `shared_alt_Hp` events                      |      475 |      487 |
| `private_alt_Hp` events                     |        0 |       27 |

The 27 `private_alt_Hp` events are the headline number the entire
stage 13 → 14 pipeline was built to expose — these are minor
alleles that exist in some individual but never reach major-allele
status anywhere in the panel, the candidate true-private /
somatic heteroplasmies. Stage 09 could not return non-zero here
because `MT_DP_AD_141.txt` was gated on panel `AC ≥ 1`.

### Stage 15 — well-bleed test (`scripts/15_well_bleed_test.py`)

Ran on the 4 high-Hp focals (77, 47, 33, 84 — all on plate 1,
i5_3). 10,000 permutations, 8-connectivity (King move), default
neighbor radius 1, far radius 3.

| Focal | Well | n_hp | Top donor (well, dist, score)   | ρ (dist, score) | ρ p  | mean(neigh) − mean(far) | diff p |
|------:|:-----|-----:|---------------------------------|----------------:|-----:|------------------------:|-------:|
|  77   | D10  |  240 | 2 (G1, **9**, 0.875)            |  −0.131         | 0.20 |  −0.034                 | 0.49   |
|  47   | B6   |  219 | 33 (H5, **6**, 0.968)           |  +0.252         | 0.95 |  −0.336                 | 0.93   |
|  33   | H5   |  215 | 21 (D3, **4**, 1.000)           |  −0.133         | 0.19 |  +0.013                 | 0.68   |
|  84   | E11  |  187 | 4 (E1, **10**, 0.995)           |  +0.147         | 0.83 |  −0.051                 | 0.56   |

**Conclusion: no well-bleed signal on any focal.** Top donors sit
4–10 wells from the focal; none is at King-move distance 1. The
Spearman ρ goes the wrong direction (positive) for two focals
(47, 84), and the neighbor-minus-far difference is slightly
*negative* (wrong direction) for three of four focals. Nothing
here would survive multiple-testing correction even if the raw
p-values were small, and the raw p-values are not small.

**The score saturation is the real story.** Top-donor concordance
is 0.875–1.000 across all four focals, at plate distances of
4–10 wells. This means dozens of panel members — regardless of
where they sit on the plate — share 88–100 % of the focal's Hp
alleles. That pattern is incompatible with random contamination
noise (which would produce low, scattered, distance-correlated
scores) but is the expected signature of *each focal carrying
the panel's other major haplotype as a minor component*.

### Biological interpretation — admixed mitos, not contamination

Doug's north/south context: REF (`NC_012312.1`) is north;
admixed Fundulus populations carry two major mt clades — north
(few ALTs, < 100 variants vs REF) and south (many ALTs, > 200).
All four focals are flagged "south". The per-individual REF_Hp
counts split the focals into two regimes (see top-10 in
`heteroplasmy_pileup_per_individual_all.tsv`):

| Focal | n_hp_sites | n_REF_Hp | n_non-REF_Hp | Pattern                                               |
|------:|-----------:|---------:|-------------:|:------------------------------------------------------|
|  47   |        219 |      211 |            8 | Major call is south-ALT; minor reads are REF (north). |
|  77   |        216 |      210 |           30 | Same as 47 — major south, minor north.                |
|  33   |        215 |       12 |          203 | Major call is REF (north); minor reads are south-ALT. |
|  84   |        187 |       10 |          177 | Same as 33 — major north, minor south.                |

So among the four "south" focals, the pileup data actually shows
**two sub-patterns**: 47 and 77 are major-south with a north
minor component (~10–30 % north reads at the diagnostic sites);
33 and 84 are major-north with a south minor component
(~10–30 % south reads at the diagnostic sites). Both patterns
sit at the same n_hp count (~190–240) because the same
~200-site north/south divergence is being sampled — just which
haplotype is major and which is minor flips between the two
sub-pairs. That's a clean signature of panel-level north × south
mt admixture, with the four focals being the most extreme
admixed individuals on this plate. Worth a sanity check against
Doug's prior "all 4 are south" flag — the pileup says only 47
and 77 are major-south; 33 and 84 are major-north with a south
minor component (potentially the haplotype calls in stage 08 are
already in the call set even though they look "north-dominant"
in the pileup because admixture flips which clade clears 0.7).

**Decision (from the session-16 decision tree).** "None significant"
in stage 15 → **not well bleed**. The recommendation is to *not*
exclude 77 / 47 / 33 / 84 from population-level analysis on the
basis of high Hp count alone. They are informative for the
maternal-line admixture story rather than contamination
artifacts. Two follow-ups stay open:

- Cross-check `MITOTYPE` (`S` / `N` / `A`) in
  `data_files_May/WGS_seq_plate.txt` against this n_REF_Hp split
  — if 47 and 77 are coded `S` (south major) and 33 and 84 are
  coded `N` (north major) or `A` (admixed), the pileup-based
  classification is reproducing whatever Doug knew at lab-prep
  time.
- For dN/dS and downstream haplotype work, decide whether to
  treat the admixed individuals as one haplotype call or to
  represent both — the 0.7 binary threshold in stage 08 forces
  one call, which loses the admixture signal.

### Heteroplasmy frequency biology — NOT 50/50

Worth recording for future sessions: there is **no reason** to
expect Hp frequencies to cluster at 50/50, and the 0.10 ≤ AF <
0.70 window is correctly broad.

- **Paternal leakage** can deposit a small amount of paternal mt
  into a zygote, yielding Hp at well below 50 % (often a few %).
- **Germline bottleneck severity** during oogenesis randomly
  reduces mt-haplotype diversity in the egg before fertilization.
  A severe bottleneck (n_mt sampled ≈ 1) drives Hp toward 0 or
  near-100 % (i.e., single haplotype carried forward); a mild
  bottleneck (n_mt sampled ≈ 1000s) preserves population-level
  Hp ratios closer to whatever the mother carried.
- An old neutral mutation in one haplotype lineage will
  accumulate hitchhiking variants on its background, so a fish
  inheriting both lineages will show Hp at *every* site that
  differs between them — not just the originating mutation.

Implication for these data: the 0.1–0.7 AF window is
appropriately wide and should not be tightened around 0.5.
Pileups showing Hp at AF ≈ 0.2 are biologically plausible (rare
mt deeply outnumbered) and should not be filtered as noise.

### Script fixes

- **`scripts/15_well_bleed_test.py` — NameError at line 267.** A
  list comprehension referenced undefined `r`, `c` while building
  `same_plate_wells`. The variable was never used downstream
  (the permutation operates on `cand_wells`), but Python
  evaluates the comprehension at definition time so the script
  errored out before any test ran. Fixed by referencing
  `info["row"], info["col"]`.
- **`scripts/15_well_bleed_test.py` — new `--events` flag.** Stage
  14 in this session wrote outputs with an `_all` suffix
  (`heteroplasmy_pileup_events_all.tsv`), but stage 15 hard-coded
  the bare filename. Added a `--events <path>` argument that
  defaults to the bare name and falls back to the `_all` variant
  if the bare name is missing (with a stderr note when the
  fallback triggers). Stage 15 now resolves the input regardless
  of whether 14's output is named with or without the suffix.

### Carry-overs still open

- **Stage 14 output naming convention.** Stage 14 currently
  writes `*_all.tsv`; stage 15 is happy with either via the
  fallback. If `_all` is the new convention, update 14 to make
  the suffix configurable (or document the convention in
  CLAUDE.md). If it isn't, rename the existing files to the bare
  names.
- **The 27 `private_alt_Hp` events.** Worth a tabular look — POS,
  Individual, Hp_AF, gene — to decide whether any are concentrated
  in one individual (suggesting somatic mutation in that fish) or
  scattered across many individuals at the same POS (suggesting a
  recurrent sequencing artifact or sub-threshold population
  variant).
- **Carry-overs from session 16 still open**: stages 10, 11, 12 not
  yet run (dN/dS, nonsyn haplotypes, NS co-occurrence). These don't
  block the heteroplasmy thread but should be done before the
  methods write-up.
- **Methods write-up.** Stages 13 / 14 / 15 architecture and the
  admixed-mito finding need to land in `docs/02_calling_architecture.md`
  (or a new `docs/03_heteroplasmy.md` — probably cleaner as a
  separate doc).

## 2026-05-20 (session 16) — Heteroplasmy stages 09/13/14 + well-bleed test 15. Three scripts staged but UNRUN; grant submission has priority.

> **Read this first if you're picking up.** This session wrote four
> scripts in support of the heteroplasmy analysis (stage 09 ran and
> produced a result; stages 13 / 14 / 15 are written but were not
> executed because Doug has a grant submission to finish). The most
> important finding from running 09 is that **the "private ALT-Hp"
> category is structurally 0 in `MT_DP_AD_141.txt`** — that file
> only contains panel-variant ALTs, so every ALT in it has at least
> one major-allele carrier by construction. Stage 14 (pileup-based)
> exists to fix this and is the actual answer to Doug's Hp question.
> Stage 15 is the quantitative well-bleed test for the 4 high-Hp
> individuals (77, 47, 33, 84) — all 4 sit on plate 1 (i5_3).

### What's new

- **`scripts/09_heteroplasmy_report.py`** — heteroplasmy classifier
  on `vcf/MT_DP_AD_141.txt` (the variant-only per-cell DP/AD table
  that Doug pre-built outside this repo). Calls major as
  `AD/DP >= 0.7` and Hp as `0.1 <= AD/DP < 0.7`. Collapses split-ALT
  rows per (Individual, POS), evaluates every allele (REF + each
  panel ALT), and classifies each Hp event as `REF_Hp` /
  `private_alt_Hp` / `shared_alt_Hp` (the latter = ALT-Hp where the
  same ALT is a major call in some other individual). Outputs:
  `vcf/heteroplasmy_events.tsv`, `..._per_site.tsv`,
  `..._per_individual.tsv`, `..._summary.txt`.

  **Results from running it (the only one of the 4 new scripts that
  has actually been run this session):**

  | Question                                          | Count |
  |---------------------------------------------------|-------|
  | Q1 individuals with ≥1 Hp event                   |    53 |
  | Q2 sites with ≥1 Hp event                         |   376 |
  | Q3 individuals with ≥2 Hp events                  |    39 |
  | Q4 Hp == REF (haplotype is ALT, REF is minor)     |   613 |
  | Q5 ALT-Hp NOT a major elsewhere (private)         |     0 |
  | Q6 ALT-Hp IS a major elsewhere (shared)           |   475 |
  | Total Hp events                                   | 1,088 |

  Q5 = 0 is a **structural artifact of the input file** —
  verified directly: every (POS, ALT) in `MT_DP_AD_141.txt` has at
  least one individual with `AF >= 0.7` (927/927). The file only
  has ALTs that passed panel variant calling, so by construction
  every ALT in it has a major-allele carrier somewhere. To get a
  biologically meaningful private-ALT-Hp count, the input has to
  be a per-CDS-position pileup table that doesn't gate on AC>=1
  — that's stage 13/14 (written this session, see below).

  **Four high-Hp individuals dominate the event count:**

  | Individual | Hp events | n_REF_Hp |
  |------------|-----------|----------|
  | 77         |       247 |      211 |
  | 47         |       215 |      208 |
  | 33         |       213 |       10 |
  | 84         |       178 |        8 |

  Together they account for **853 / 1,088 events ≈ 78 %**. After
  these four, the next-highest individual has 24 Hp events. Doug's
  prior text-flag — "individuals with many Hp events are unlikely
  to be real heteroplasmy; sequencing error / contamination" —
  applies here. Stage 15 tests one specific contamination
  hypothesis (well bleed during library prep).

- **`scripts/13_pileup_cds_AD.sh`** (NEW, NOT YET RUN) — Mac-side
  `bcftools mpileup -a AD,DP -d 100000 --no-BAQ` over
  `docs/mito_protein_coding.bed` (13 CDS intervals, 11,417 bp)
  across all 141 panel BAMs (excludes 70, 125 by filename filter).
  Output: `vcf/pileup_cds_141.vcf.gz` + `.tbi`. Auto-builds
  `vcf/slim_bamlist_141.txt` from `MT_only_bams/` directory
  listing. `--no-BAQ` matches the canonical 05_1 caller so the
  AD field here is on the same scale as everything else in the
  pipeline. Expected wall time on Mac: 1–3 min; peak RAM well
  under 1 GB; output ~50–100 MB.

- **`scripts/14_hp_from_pileup.py`** (NEW, NOT YET RUN) — Hp
  detection on the stage-13 pileup VCF, the proper answer to Q5.
  Thresholds: `DP >= 20`, `0.10 <= AD/DP < 0.70`, `AD_Hp >= 4`
  (the `AD_Hp >= 4` floor is non-binding at normal mtDNA depths;
  it just prevents 1- or 2-read calls from passing at sites where
  coverage drops near the DP=20 boundary). Same 3-way classifier
  as stage 09 (REF_Hp / private_alt_Hp / shared_alt_Hp), but the
  private-ALT-Hp category is now reachable because the input is
  coverage-honest at every CDS position, independent of which
  positions/ALTs the panel variant-calling stage produced.
  Streams the pileup VCF with `bcftools query` and never holds
  more than a few rows in memory. Outputs:
  `vcf/heteroplasmy_pileup_events.tsv`, `..._per_site.tsv`,
  `..._per_individual.tsv`, `..._summary.txt`.

- **`scripts/15_well_bleed_test.py`** (NEW, NOT YET RUN) —
  quantitative test of "is the high Hp load on individuals 77 /
  47 / 33 / 84 explained by library-prep contamination from
  neighboring wells on the sequencing plate?" Uses
  `data_files_May/WGS_seq_plate.txt` for the plate map.

  **Plate-map decoding (figured out this session, documented in
  the script):** the file has 154 rows on what looks like a
  single 96-well plate (well IDs A1–H12) but each well label
  appears twice — there are actually **two plates** distinguished
  by the `i5` adapter column. `i5_3` is plate 1 (96 samples,
  fully filled); `i5_4` is plate 2 (59 samples, partial, ending
  at F8). Row letter A–H maps to row 1–8; column 1–12 is the
  column index. The 4 focal high-Hp individuals are all on
  **plate 1 (i5_3)**, at H5 (33), B6 (47), D10 (77), E11 (84) —
  scattered, not bunched, which is what makes the well-bleed
  hypothesis testable rather than visually obvious.

  **Donor-concordance score:** for focal X and candidate Y on
  the same plate,

      Score(X -> Y) = | { Hp event in X where Y's major call at POS == X's Hp_allele } |
                      / | X's Hp events |

  i.e., the fraction of X's heteroplasmy events explainable by Y's
  haplotype. Y's major calls are computed on-the-fly from the
  stage-13 pileup VCF (so they're coverage-honest at every CDS
  position, not just panel-variant sites).

  **Two test statistics with permutation null** (10,000 perms by
  default, X's well fixed at observed position, OTHER same-plate
  samples shuffled across OTHER same-plate wells):

    (a) Spearman rho of `Score` vs. plate-distance(X, Y), Chebyshev
        metric. Bleed signal = negative rho (Score falls with
        distance). One-sided p.
    (b) `mean(Score | neighbor wells, dist <= 1)` minus
        `mean(Score | far wells, dist >= 3)`. Bleed signal =
        positive difference. One-sided p.

  Outputs: `vcf/well_bleed_donor_ranking.tsv`,
  `..._results.tsv`, `..._summary.txt`. The ranking file is the
  "who is the top donor candidate for each focal" table — the
  intuitively most readable. The results file is the formal test
  statistics.

### What needs to happen next session (run order, ~5 min total)

```
cd ~/Projects/MT_Genomics_Cl_Ap2026/MT_Genomics2
conda activate SNP_env

# 1. Per-CDS-position pileup (1-3 min, ~50-100 MB output):
bash scripts/13_pileup_cds_AD.sh
#    → vcf/pileup_cds_141.vcf.gz (+ .tbi)
#    → vcf/slim_bamlist_141.txt

# 2. Pileup-based heteroplasmy classification (~30 s):
python scripts/14_hp_from_pileup.py
#    → vcf/heteroplasmy_pileup_{events,per_site,per_individual,summary}.{tsv,txt}

# 3. Well-bleed contamination test on the 4 high-Hp focals (~30 s):
python scripts/15_well_bleed_test.py
#    → vcf/well_bleed_{donor_ranking,results,summary}.{tsv,txt}
#    Optional flags:
#      --focals 77,47,33,84        explicit focal list
#      --n-perm 20000              more permutations
#      --connectivity 4            rook adjacency instead of king
```

### Open questions to answer with the stage-14/15 outputs

1. **Stage 14: how many TRUE private ALT-Hp events?** Stage 09
   returned 0 for this category, but only because its input was
   variant-only. Stage 14's input is pileup-based, so private
   ALT-Hp can now be a positive number. **This is the most
   important number from the whole heteroplasmy thread.** A
   non-trivial count means there are alternate bases at CDS
   positions that are present at 10–70 % in some individual but
   were never variant-called anywhere — candidate true somatic /
   private heteroplasmies.

2. **Stage 15: are the 4 high-Hp focals explained by well bleed?**
   Three possible outcomes:
   - All 4 significant on both tests (rho_p < 0.05 AND
     neighbor_perm_p < 0.05) **and** the top donor is in a
     neighbor well → drop / contaminant-correct those four.
   - 1–2 significant → bleed for some, not others. Look at each
     focal's top donor individually.
   - None significant → not well bleed. Other contamination
     route (extraction batch, pre-library PCR, real heteroplasmy
     load) becomes the working hypothesis. The top-donor table
     is still useful for tracing back to a specific batch if
     extraction metadata exists.

3. **Should the high-Hp individuals be excluded from population-
   level analyses entirely** (in addition to 70 and 125), or
   back-corrected? Depends on (2). If real bleed, back-correction
   is justifiable; if real heteroplasmy load, the individuals are
   biologically interesting and shouldn't be discarded.

4. **Cross-check stage-14 Hp results against MITOTYPE in the
   plate file.** WGS_seq_plate.txt has a MITOTYPE column (S / N /
   A) — currently uncharacterized in the project docs. If the
   four high-Hp focals share a MITOTYPE or differ from their
   neighbors in a meaningful way, that's an additional handle on
   what the load means.

### Artifacts changed this session

- `scripts/09_heteroplasmy_report.py` — new (RUN; outputs in `vcf/`).
- `scripts/13_pileup_cds_AD.sh` — new (UNRUN).
- `scripts/14_hp_from_pileup.py` — new (UNRUN).
- `scripts/15_well_bleed_test.py` — new (UNRUN).
- `vcf/heteroplasmy_events.tsv` — new (from 09).
- `vcf/heteroplasmy_per_site.tsv` — new (from 09).
- `vcf/heteroplasmy_per_individual.tsv` — new (from 09).
- `vcf/heteroplasmy_summary.txt` — new (from 09).
- `CLAUDE.md` — Next-session pickup block replaced; scripts
  catalog updated to list 13/14/15; the 4 high-Hp focals carried
  forward as a named open question.
- `CHANGELOG.md` — this entry.

### Pitfalls / notes for next session

- **`MT_DP_AD_141.txt` cannot answer the private-ALT-Hp question**
  (confirmed empirically: all 927 panel-ALT combos have at least
  one AF ≥ 0.7 carrier). Don't re-investigate this in stage 09;
  the answer is in stage 14's output.
- **`bcftools mpileup` sample names.** Stage 13 doesn't reheader,
  so sample names in the pileup VCF come from the BAM `@RG SM:`
  tag. Stage 14 and 15 both normalize sample names to numeric
  WGS_ID via regex (`re.compile(r"(\d+)")` matched against the
  leading part of the name), so they're robust to whatever the
  `SM:` tag actually contains (`1`, `1_0`, `1_MT`, etc.).
- **Permutation in stage 15 fixes X at its observed well.** This
  is deliberate — the null is "given X is where it is, do the
  other samples' well assignments around X drive the observed
  concordance pattern?" Not "would X have looked contaminated if
  it had been placed in some other well?" The latter isn't the
  right question for this dataset.
- **Stage 15 score is naturally higher for candidates whose
  haplotype is most different from X.** That's by design — only
  difference can explain Hp events — but it means a candidate
  who happens to be from the major 64-fish NS clade (session-15
  finding) will rank highly against any non-clade focal X
  regardless of well position. The permutation null controls for
  this (it shuffles wells, not haplotypes), but when reading
  `well_bleed_donor_ranking.tsv` by eye, remember that
  Y-being-from-the-other-clade can produce a high Score even
  without any contamination story.

## 2026-05-18 (session 15, evening) — Downstream analysis scripts: 10_dnds_per_gene + 11_haplotypes_nonsyn; documented the c.Y vs POS coordinate-system distinction

> Two new downstream scripts on top of the canonical
> `Fhet_MT_CDS.snps.split.vcf.gz`, plus a short primer in `CLAUDE.md`
> explaining the two coordinate systems that coexist in SnpEff-
> annotated VCFs (POS = genomic, `c.Y` = CDS-relative). No changes
> to the caller (05_1), stage 06, or stage 07.

### What's new

- **`scripts/10_dnds_per_gene.py`** — per-gene + Overall dN, dS, dN/dS.
  - **Site counting:** Nei-Gojobori on the **NCBI translation table 2**
    (vertebrate mitochondrial code: TGA=W, ATA=M, AGA/AGG=Stop). For
    each codon, count what fraction of the 9 possible single-nucleotide
    substitutions are synonymous (`s_i/3` summed over the 3 positions)
    vs nonsynonymous (`3 − S_codon`). Stop-codon-creating substitutions
    are classified as nonsynonymous. Stop codons (including the
    terminal stop) contribute nothing to either count.
  - **Strand:** read from the GFF. ND6 is the one `-` strand gene in
    the *F. heteroclitus* mt genome; its CDS is reverse-complemented
    before site counting so the codons are read 5'→3' in the gene's
    own orientation. The other 12 genes are `+` strand.
  - **PolyA-completed stops (important):** 5 of 13 genes (ND2, COX2,
    COX3, ND3, ND4) end one base short of a complete terminal stop
    codon in the genomic sequence — the final TAA is created
    post-transcriptionally by 3'-polyadenylation. Their CDS lengths
    are `3*n + 1`. The script truncates to `len - (len % 3)` so the
    trailing partial codon (and any complete stop codon) is dropped
    before site counting. *Previous dN/dS attempts on this genome
    failed on those 5 genes because the partial last codon was
    passed into the translator and produced garbage.* This script is
    designed around that failure mode and the docstring records it
    explicitly so it doesn't get re-introduced.
  - **Observed counts:** from the canonical's first ANN entry per
    row. `synonymous_variant` (and `&` variants thereof) → S;
    `missense_variant`, `stop_gained`, `stop_lost`, `start_lost` → N.
    HIGH-impact stops are included as N by default.
  - **Output:** `vcf/dnds_per_gene.tsv` — one row per gene, plus an
    "Overall" pooled row. Columns: `Gene, Strand, CDS_len, n_codons,
    S_sites, N_sites, obs_S, obs_N, other_eff, pN, pS, pN_over_pS,
    dN_JC, dS_JC, dN_dS_JC`. Both the simple `obs / possible` ratio
    (user spec) and Jukes-Cantor-corrected `d` values reported; JC
    is the standard polymorphism-level correction and reduces to the
    simple ratio at low divergence.
  - **Sanity checks (sandbox-tested):** Nei-Gojobori site counts
    verified against canonical hand cases — TTT (Phe, 0 syn at pos
    1+2, 1 syn at pos 3 → S≈0.333), TGA (Trp under mt code, 0+0+1
    syn → S≈0.333), CCC (Pro, 4-fold degenerate at pos 3 → S=1.0).

- **`scripts/11_haplotypes_nonsyn.py`** — haplotype matrix + call
  strings, **MODERATE-impact only** (`missense_variant` /
  `splice_region_variant&missense_variant`). Otherwise identical
  shape to `scripts/08_call_haplotypes.py`: same `AD_alt/DP > 0.7`
  rule, same `MIN_DP = 3`, same `70`/`125` sample exclusions, same
  never-called-row drop. Haplotype strings prefixed `N_` (vs stage
  08's `C_`) so the two products don't collide. Output:
  `vcf/haplotype_matrix_nonsyn.csv`, `vcf/haplotype_calls_nonsyn.csv`.
  From the current canonical's impact tally (166 missense), ~150
  rows expected after never-called drop.

- **`CLAUDE.md`** — new "Downstream analysis scripts" section before
  "Active primary task" that lists 08, 10, 11 (existing/new), 09
  (planned), and `DP_AD_table.py`. Includes a SnpEff ANN-format
  primer and a worked-example table showing POS ↔ c.Y conversion
  for ND1 (Doug's example positions 2847 and 2852). The c.Y
  question — "why does the INFO field say c.6 when POS is 2847?" —
  was explained in chat (CDS coordinate within the gene, 1-based
  from start codon; for + strand `c.Y = POS − (gene_start − 1)`,
  for ND6 the relationship inverts) and is now documented inline.

### Data referenced

- Canonical: `vcf/Fhet_MT_CDS.snps.split.vcf.gz` (Doug confirmed
  06/07 ran cleanly this session; not auto-verified by Claude but
  the script reads it on first run and will fail-fast if absent).
- AD long table: `vcf/mtDNA_long_AD_table.tsv` (132,991 rows,
  ~143 × ~930 sites, from `scripts/DP_AD_table.py`). Not used by
  10 or 11 — both go directly to the VCF for ANN parsing — but
  useful for ad-hoc per-(sample, position) analyses.

### Run order (pickup)

```
cd ~/Projects/MT_Genomics_Cl_Ap2026/MT_Genomics2
conda activate SNP_env

# dN/dS per gene + Overall:
python scripts/10_dnds_per_gene.py
# → vcf/dnds_per_gene.tsv  (and a pretty-printed table to stdout)

# Nonsynonymous haplotype matrix + calls:
python scripts/11_haplotypes_nonsyn.py
# → vcf/haplotype_matrix_nonsyn.csv
# → vcf/haplotype_calls_nonsyn.csv
```

### Late-session ad-hoc analysis: nonsynonymous co-occurrence

Doug asked "are there NS sites that always occur together?" — i.e.,
sites in perfect linkage disequilibrium across the panel. Ran the
analysis on `vcf/haplotype_matrix.csv` (the 08 output) filtered to
its 164 missense rows, 141 sample columns (excluded the manually-
added `sum`, `141`, and `count indif` annotation columns). Three
findings:

1. **Fixed reference-divergence at 3 NS sites** (every panel sample
   carries the ALT — reference is the minority):

   | Gene | POS | Change |
   |------|-----|--------|
   | ND1  | 3124 | G>A |
   | ND2  | 4680 | A>T |
   | ND2  | 4957 | C>T |

   These should arguably be excluded from population-level dN/dS
   (they're divergence from the GenBank reference, not polymorphism
   within the panel). Open decision for next session.

2. **One major nonsynonymous haplogroup**: 7 missense sites in 4
   genes (ATP6, ND1, ND2 ×3, ND5 ×2) all co-segregating in 64/141
   samples (~45 %). Almost certainly a major maternal clade.

   | Gene | POS   | Change |
   |------|-------|--------|
   | ATP6 |  8451 | A>G |
   | ND1  |  3061 | T>C |
   | ND2  |  4737 | G>A |
   | ND2  |  4951 | T>C |
   | ND2  |  5061 | A>G |
   | ND5  | 12589 | C>A |
   | ND5  | 13577 | T>A |

3. **Six small co-segregating clusters** of 2–3 NS sites carried by
   2–4 fish each (likely sib-group / close-kin sharing). Detail in
   stdout of `scripts/12_ns_cooccurrence.py --min-carriers 2`.

A further 24 multi-site "groups" the analysis turned up are at 1
carrier each — these are aggregations of one fish's private NS
variants (the signature [0,…,0,1,0,…,0] is the same for every site
private to the same sample). Not population-level LD. Filtered out
by `--min-carriers 2` in the script.

### `scripts/12_ns_cooccurrence.py` (new)

Makes the above analysis reproducible. Reads
`vcf/haplotype_matrix_nonsyn.csv` (preferred — from 11) or falls
back to `vcf/haplotype_matrix.csv` filtered to missense rows.
Identifies multi-site groups whose 141-element per-sample call
vectors are identical. Classifies each group as
`fixed_ref_divergence`, `haplogroup`, or `singleton_artifact`.
Writes `vcf/ns_cooccurrence_groups.tsv`. CLI flag `--min-carriers`
hides low-signal groups.

### CLAUDE.md additions

- New "Next-session pickup" block at the top of the active-task
  area listing the three unrun scripts (10, 11, 12) and the three
  open scientific questions (reference-divergence handling,
  haplogroup-vs-metadata cross-check, stage 09 still to write).
- "Preliminary findings to carry forward (session 15)" subsection
  in the "Downstream analysis scripts" section, recording the
  3 fixed sites + the 7-site haplogroup + the 6 small clusters so
  they're durable across sessions.
- Catalog entry for `scripts/12_ns_cooccurrence.py`.

### Artifacts changed this session

- `scripts/10_dnds_per_gene.py` — new (dN/dS per gene + Overall;
  vertebrate-mt code; polyA-stop handling; strand-aware).
- `scripts/11_haplotypes_nonsyn.py` — new (stage-08 shape, nonsyn filter).
- `scripts/12_ns_cooccurrence.py` — new (perfect-LD groups among NS sites).
- `CLAUDE.md` — new "Next-session pickup" block, new "Downstream
  analysis scripts" section, new "Preliminary findings" subsection,
  SnpEff ANN primer, c.Y/POS worked example, catalog entry for 12.
- `CHANGELOG.md` — this entry (extended).

## 2026-05-18 (session 14) — Caller renamed and re-scoped: 05f → 05_1_mpileup_merge.sh with norm-split inside the script and a hard baseline-regression gate; honoring the sessions 6–10 per-sample baseline as the floor

> **Two corrections to session 13:**
>
> 1. **The caller's name and shape didn't reflect what worked
>    previously.** The validated chain from sessions 6–10 did
>    `bcftools norm -m -any` *per sample inside the caller* (05d2),
>    not in stage 07. The session-13 05f deferred the norm-split to
>    stage 07, which is the wrong place — SnpEff (stage 06) sees
>    multi-ALT rows and can't annotate each ALT independently.
>    Moved norm-split into the caller so the output is already in
>    one-row-per-ALT form by the time stage 06 sees it. Renamed to
>    `jobs/05_1_mpileup_merge.sh` (Doug's preferred informative
>    naming, matching the operations the script actually performs).
>
> 2. **The baseline-comparison gate was framed as "informational",
>    not as a hard precondition.** It is a hard precondition. The
>    per-sample → merge baseline (1128 SNPs / 143 samples / ts/tv
>    8.17) is the result of three weeks of diagnosis closing the
>    sessions 6–10 variant-count discrepancy. The new joint caller's
>    job is to *preserve* that floor and ideally find more — not to
>    declare its own variant set canonical and call the question
>    closed. Rule #1 in CLAUDE.md now requires `positions lost vs
>    baseline = 0` (along with `cells with DP=. = 0` and full ANN
>    coverage) before `Fhet_MT_CDS.snps.split.vcf.gz` is considered
>    frozen.

### What changed

- **New canonical caller `jobs/05_1_mpileup_merge.sh`** (replaces
  session-13 05f):
  ```
  bcftools mpileup -f $REF -b $SLIM_BAMLIST -a AD,DP \
      --max-depth 100000 --threads 8 -Ou \
    | bcftools call -mv --ploidy 1 -Ou \
    | bcftools norm -m -any -f $REF -Oz -o $OUT
  ```
  - No intermediate writes — single pipeline through `-Ou` /
    `-Ou` / `-Oz`.
  - **norm-split lives in the caller**, matching where 05d2 did the
    per-sample split. Output is already one-row-per-ALT, so SnpEff
    (stage 06) annotates each ALT independently and stage 07 no
    longer needs to run `norm -m -any`.
  - Sample-rename after the pipeline (slim BAMs have no @RG SM:;
    `sed -e 's|.*/||' -e 's/_0_MT_only\.bam$/_MT/'`).
  - Baseline-vs-05_1 regression check built in: writes
    `vcf/05_1_Fhet_mt_persample_merged_vs_baseline_lost.tsv` and
    `_gained.tsv` against `Fhet_mt_persample_merged.vcf.gz` (the
    May-15 baseline). The manifest surfaces "positions lost vs
    baseline: N" and warns loudly if N > 0.
  - Output: `vcf/05_1_Fhet_mt_persample_merged.vcf.gz` (Doug's
    naming — script-prefixed so the producing script is obvious
    from the filename).
  - BSUB headers: `-n 8`, 16 GB mem, 4 h walltime, normal queue.
    Wallclock expectation 30–60 min.

- **Archived to `jobs/_archive/`** (no longer in the active path):
  - `05f_joint_call.sh` (session-13 first cut; superseded by 05_1).

- **Stage 06 (`scripts/06_snpeff_mac.sh`) input/output renamed:**
  - Input: `Fhet_mt_joint.vcf.gz` → `05_1_Fhet_mt_persample_merged.vcf.gz`
  - Output: `Fhet_mt_joint_ann.vcf.gz` → `05_1_Fhet_mt_persample_merged_ann.vcf.gz`
  - rsync command in header updated to fetch the new tag.

- **Stage 07 (`scripts/07_cds_snps_norm_mac.sh`) further simplified:**
  - Input renamed to match stage 06 output.
  - **Dropped the `bcftools norm -m -any -f REF` step** — already
    done in 05_1.
  - Dropped the FASTA-related code path (no `-f REF` anywhere here
    now; norm-split is upstream).
  - Pipeline now: `bcftools view -R BED | bcftools view -v snps` →
    sample rename → stats. ~30 lines of bcftools, no joint pileups,
    no targets files, no ANN-transfer.
  - Summary block enforces the FREEZE GATE explicitly: `cells with
    DP=. = 0` AND `records w/ ANN = record count` AND 05_1's
    manifest `positions lost vs baseline = 0` must all be green.

- **The May-15 stats file** (`vcf/Fhet_mt_persample_merged_stats_15May.txt`,
  copied across by Doug this session) is now the
  contemporaneous-stats record for the baseline. 1144 records
  (1128 SNPs + 16 indels) / 143 samples / 0 multiallelic SNP sites
  / ts/tv 8.17 / 94 % of sites at DP > 500. "0 multiallelic SNP
  sites" reflects the post-`norm -m -any` representation 05d2
  produced — multi-row biallelic, not multi-ALT — which is why
  CDS-restricted preBackfill had 931 rows / 914 unique positions
  / 17 multi-row positions. The new 05_1 architecture preserves
  this representation by construction.

- **CLAUDE.md substantially rewritten:**
  - "Canonical caller" block: walks through the architecture
    trajectory (sessions 6–10 per-sample baseline → session 12
    failed backfill → session 13 joint switch → session 14 final
    shape) and names the May-15 baseline as the regression floor.
  - "Active primary task": updated to session-14 state, with the
    HARD GATE explicit at step 2 (read the manifest, confirm `cells
    with DP=. = 0`, `n_multiallelic = 0`, `positions lost = 0`,
    `n_SNPs ~1100–1200`; stop and diagnose otherwise).
  - "Known pitfalls": updated rule #2's in-production fix to point
    at 05_1 (with `norm -m -any` inside).
  - "Validation principle": sharpened to "zero regressions" and
    explicit about the May-15 baseline being the floor that 3 weeks
    of work established.
  - Rule #1: status reflects the new architecture; freeze gate
    re-stated.
  - "Compute split" rewritten.

### Why this matters (lesson)

The session-13 entry hand-waved past a real concern — "joint -mv
ought to find every site per-sample does" — without testing it. Doug
called that out and was right. Joint and per-sample callers are not
guaranteed to produce identical variant sets at marginal positions,
and the right scientific posture is to *prove* the new caller
preserves the validated floor before declaring it canonical, not to
assume it does. The baseline-regression gate is now a hard
precondition for downstream stages, not an after-the-fact sanity
check. This is a general project norm now (see CLAUDE.md → Known
pitfalls → Validation principle).

### Run order (next session pickup)

```
# 1. Submit on T2:
ssh dcrawford@t2.idsc.miami.edu
cd /projectnb/dcrawford/MT_Genomics2
bsub < jobs/05_1_mpileup_merge.sh
# logs: logs/05_1_mpileup_merge_<jobid>.{out,err}

# 2. VALIDATE 05_1 before flowing downstream (HARD GATE):
cat vcf/05_1_Fhet_mt_persample_merged_run_manifest.txt
#      cells with DP=.: 0                                ← target: 0
#      n_multiallelic_sites: 0                           ← norm-split happened
#      positions lost vs baseline: 0                     ← REGRESSION if > 0
#      n_SNPs: ~1100-1200                                ← sanity
#    If positions lost > 0, STOP. Inspect
#    vcf/05_1_Fhet_mt_persample_merged_vs_baseline_lost.tsv and diagnose.

# 3. rsync 05_1 output to Mac:
rsync -avP \
  dcrawford@t2.idsc.miami.edu:/projectnb/dcrawford/MT_Genomics2/vcf/05_1_Fhet_mt_persample_merged.vcf.gz \
  dcrawford@t2.idsc.miami.edu:/projectnb/dcrawford/MT_Genomics2/vcf/05_1_Fhet_mt_persample_merged.vcf.gz.csi \
  ~/Projects/MT_Genomics_Cl_Ap2026/MT_Genomics2/vcf/

# 4. On Mac (SNP_env):
cd ~/Projects/MT_Genomics_Cl_Ap2026/MT_Genomics2
conda activate SNP_env
bash scripts/06_snpeff_mac.sh
bash scripts/07_cds_snps_norm_mac.sh

# 5. FREEZE-GATE checks on the canonical:
zcat vcf/Fhet_MT_CDS.snps.split.vcf.gz | grep -vc '^#'        # expected ~900-1100 rows
bcftools query -f '[%DP\n]' vcf/Fhet_MT_CDS.snps.split.vcf.gz \
  | awk '$1=="."' | wc -l                                     # expected 0
bcftools view vcf/Fhet_MT_CDS.snps.split.vcf.gz \
  | grep -v "^#" | grep -c "ANN="                             # should equal record count

# Frozen ONLY if all three checks are green AND step 2 showed
# "positions lost vs baseline: 0".

# 6. Then: python scripts/08_call_haplotypes.py; write scripts/09_heteroplasmy_report.py.
```

### Artifacts changed this session

- `jobs/05_1_mpileup_merge.sh` — new (canonical caller with norm-split
  inside, baseline regression gate).
- `jobs/_archive/05f_joint_call.sh` — archived.
- `scripts/06_snpeff_mac.sh` — input/output filenames updated
  (`05_1_Fhet_mt_persample_merged{,_ann}.vcf.gz`).
- `scripts/07_cds_snps_norm_mac.sh` — further simplified: dropped
  `norm -m -any -f REF` step (now in 05_1), dropped REF FASTA path
  variable (not needed without norm), tightened FREEZE GATE language
  in summary block.
- `CLAUDE.md` — canonical-caller section, active-primary-task,
  compute-split, rule #1, pitfall fix-in-production lines, and
  validation principle all updated for the 05_1 architecture.
- `CHANGELOG.md` — this entry.
- `vcf/Fhet_mt_persample_merged_stats_15May.txt` — Doug placed this
  earlier in the session; preserves the May-15 baseline's stats so
  the validated 1128-SNP number is on Mac for future reference.

## 2026-05-17 (session 13) — Caller architecture switch: per-sample → joint (05d2 + 05e2 → 05f); stage 07 backfill removed; session-12 Q/q diagnosis corrected (real cause was a targets-file format bug at split multiallelic sites)

> **Two related fixes:**
>
> 1. **Root-cause fix for the `.:.:.:.` REF cells.** Cells were empty
>    because `05d2`'s `bcftools call -mv` only emitted variant rows per
>    sample, so REF samples had no row to contribute to the `05e2`
>    merge. Stage 07's joint-pileup backfill was a downstream band-aid.
>    Switched the caller to a single-pass joint mpileup + call across
>    the 143 slim BAMs (the architecture `Fhet_mt_fullAD.vcf.gz` used),
>    which fills DP/AD on every (POS × sample) cell at every variant
>    position by construction. Backfill retired; stage 07 collapses to
>    CDS-restrict + SNPs-only + norm + rename.
>
> 2. **Corrected diagnosis of the session-12 17-site loss.** I had
>    blamed `-Q 30 -q 30` in the backfill's mpileup, but examination of
>    the per-sample evidence at the 17 lost positions (this session,
>    in chat) showed every one is a multiallelic site that `norm -m
>    -any` had split into 2 rows in `preBackfill`, and the lost ALT is
>    always the minor allele (called by 1–2 samples). The real cause
>    is a **targets-file format bug**: the backfill built its targets
>    file with `bcftools query -f '%CHROM\t%POS\t%REF,%ALT\n'`, which
>    produced two rows at the same position for split multiallelic
>    sites. `bcftools call -C alleles -T` keys its lookup on
>    (CHROM, POS) — only one entry per position survives — so the
>    second ALT was silently dropped at every multiallelic site.
>    `-Q 30 -q 30` and `--ploidy 1` were both red herrings.
>
>    Hard evidence:
>    ```
>    preBackfill:  rows=931  unique_positions=914  multiallelic_split=17
>    canonical:    rows=914  unique_positions=914  multiallelic_split=0
>    17 multiallelic positions ∩ 17 lost positions = 17 of 17
>    ```
>
>    This bug is structurally impossible in the new 05f architecture:
>    joint `-mv` always emits multiallelic sites as one row with
>    comma-joined ALTs (`T C,A`), and stage 07's `norm -m -any` splits
>    that row *after* the call — so all ALTs survive with per-cell AD
>    on both rows. No targets-file lookup ever happens for the new
>    canonical.
>
> The per-sample → merge workaround (05d2 / 05e2 / 07b) was
> specifically chosen in sessions 6–10 to avoid the bcftools-1.6
> high-depth joint-call collapse documented in CHANGELOG 2026-05-15.
> With bcftools 1.23.1 from `$HOME/software/local/bin` (PATH-injected
> by `jobs/config.sh`), the joint -mv pipeline produces the full
> variant set directly. The 1.6 workaround is no longer necessary AND
> it was producing the empty-REF-cells problem that motivated stage 07's
> backfill in the first place.

### What changed

- **New canonical caller: `jobs/05f_joint_call.sh`.**
  ```
  bcftools mpileup -f $REF -b $SLIM_BAMLIST -a AD,DP \
      --max-depth 100000 --threads 8 -Ou \
    | bcftools call -mv --ploidy 1 -Oz -o Fhet_mt_joint.vcf.gz
  ```
  - No `-A` on `call` (CHANGELOG 2026-05-09 documented that the
    joint experiments 05/05b/05c had `-A` and produced 96 %-
    multiallelic "phantom-alt soup"; without `-A` only ALT alleles
    with meaningful evidence enter the call).
  - No `-Q` / `-q` overrides on `mpileup` (matches 05d2's choice of
    bcftools defaults `-Q 13 -q 0`).
  - `--ploidy 1` (haploid mtDNA; AD/DP are independent of ploidy so
    heteroplasmic evidence is preserved — only the emitted GT
    collapses to majority allele).
  - `--max-depth 100000` (mt depth on slim BAMs runs 5–20k × per
    sample; the bcftools default of 8000 would truncate).
  - Sample rename (`s|.*/||; s/_0_MT_only\.bam$/_MT/`) and full
    self-documenting manifest written before the long step.
  - BSUB headers: `-n 8`, 16 GB mem, 4 h walltime, normal queue.
    Wallclock expectation: 30–60 min on 143 slim BAMs.

- **Archived to `jobs/_archive/`** (no longer in the active path):
  - `05d2_persample_call.sh` (per-sample call, slim BAMs, bcftools
    1.23.1) — the session-10 canonical caller, superseded.
  - `05e2_merge_persample.sh` (merge of the per-sample outputs).
  - `07b_backfill_AD.sh` (T2-side AD/DP backfill; no longer needed
    because the joint caller already fills every cell).

- **Stage 06 updated** (`scripts/06_snpeff_mac.sh`):
  - Input renamed: `Fhet_mt_persample_merged.vcf.gz` → `Fhet_mt_joint.vcf.gz`.
  - Output renamed: `Fhet_mt_persample_merged_ann.vcf.gz` → `Fhet_mt_joint_ann.vcf.gz`.
  - rsync command in header updated to fetch the joint output.

- **Stage 07 simplified** (`scripts/07_cds_snps_norm_mac.sh`):
  Removed the targets-file machinery, the dual-targets-file split,
  the joint mpileup over slim BAMs, the force-call against canonical
  alleles, the ANN-transfer step, and the `preBackfill.vcf.gz`
  intermediate. What remains:
  1. CDS restrict via `docs/mito_protein_coding.bed`
  2. SNPs only (`bcftools view -v snps`)
  3. `bcftools norm -m -any -f REF` (split multiallelic, anchor against REF)
  4. Sample rename (idempotent — handles the case where 05f's own
     reheader already cleaned the names)
  ANN from stage 06 passes through transparently — it lives in the
  INFO field and `bcftools view / norm` preserve it. No separate
  annotate call needed.

- **CLAUDE.md updated**: "Active primary task" rewritten to the new
  run order (05f on T2 → 06 → 07 → 08/09); "Canonical caller" section
  rewritten; rule #1 status line updated; compute-split section
  rewritten to describe the joint architecture.

### Corrected diagnosis of the session-12 17-site loss

The session-12 entry below attributes the 17-site loss to `-Q 30 -q
30` filtering out reads at marginal-quality positions. **That is
incorrect.** Re-examining each of the 17 positions in
`Fhet_mt_persample_merged.vcf.gz` showed all 17 are multiallelic
sites with two ALTs each; the lost ALT in every case is the minor
allele, called by 1–2 samples with AD_alt of 172–5323 (i.e., far
above any plausible quality floor — these were not marginal calls).

The real mechanism: `bcftools query -f '%CHROM\t%POS\t%REF,%ALT\n'`
on the post-norm `preBackfill` file emits one row per VCF record. At
the 17 multiallelic sites (which `norm -m -any` had split into 2
rows each), this produced two targets entries at the same position
with different REF,ALT values. `bcftools call -C alleles -T` builds
a (CHROM, POS) → alleles dictionary from the targets file and keeps
only one entry per position. So at every multiallelic site, the
secondary entry was silently dropped from the lookup and the
secondary ALT never made it into the output.

The correct format for `call -C alleles -T` at a multiallelic site
is *one row per position with all ALTs comma-joined*:
`NC_012312.1\t2847\tT,C,A`. The backfill's targets builder didn't
do that consolidation.

Session-12 patch 3 (removing `-Q 30 -q 30`) was therefore correct in
its destination state (defaults match 05d2) but wrong in its stated
rationale. The Q/q removal happened to be harmless. The actual
position count would have stayed at 914 either way until the
targets-file builder was fixed — which the session-13 architecture
switch obviates entirely.

### Validation requirement for 05f

Joint `-mv` and per-sample `-mv` → merge can legitimately produce
different variant sets at marginal sites. Before declaring 05f's
output canonical, the next-session pickup MUST run a diff against
`Fhet_mt_persample_merged.vcf.gz` (the proven 1128-SNP baseline
from session 10) and confirm 05f does not drop any positions the
per-sample baseline found. If 05f is missing sites, do not proceed
— diagnose before flowing into 06 / 07. See the run-order block
below for the exact comparison commands.

### Run order (next session pickup)

```
# 1. Submit on T2:
ssh dcrawford@t2.idsc.miami.edu
cd /projectnb/dcrawford/MT_Genomics2
bsub < jobs/05f_joint_call.sh
# logs: logs/05f_joint_<jobid>.{out,err}; output: vcf/Fhet_mt_joint.vcf.gz

# 2. VALIDATE 05f before flowing downstream (HARD GATE).
#    Check the manifest for three lines:
cat vcf/Fhet_mt_joint_run_manifest.txt
#      cells with DP=.: 0                                  ← target: 0
#      positions lost vs baseline: 0                       ← target: 0 (REGRESSION if > 0)
#      n_SNPs: ~1100-1200                                  ← sanity
#    If "positions lost vs baseline" > 0, STOP. Inspect
#    vcf/Fhet_mt_joint_vs_baseline_lost.tsv and diagnose.
#    "positions gained vs baseline" > 0 is fine and expected.

# 3. rsync to Mac:
rsync -avP \
  dcrawford@t2.idsc.miami.edu:/projectnb/dcrawford/MT_Genomics2/vcf/Fhet_mt_joint.vcf.gz \
  dcrawford@t2.idsc.miami.edu:/projectnb/dcrawford/MT_Genomics2/vcf/Fhet_mt_joint.vcf.gz.csi \
  ~/Projects/MT_Genomics_Cl_Ap2026/MT_Genomics2/vcf/

# 4. On Mac:
cd ~/Projects/MT_Genomics_Cl_Ap2026/MT_Genomics2
conda activate SNP_env
bash scripts/06_snpeff_mac.sh
bash scripts/07_cds_snps_norm_mac.sh

# 5. Verify the canonical:
zcat vcf/Fhet_MT_CDS.snps.split.vcf.gz | grep -vc '^#'   # expected ~900-1100 SNPs
bcftools query -f '[%DP\n]' vcf/Fhet_MT_CDS.snps.split.vcf.gz \
  | awk '$1=="."' | wc -l                                # expected 0
bcftools view vcf/Fhet_MT_CDS.snps.split.vcf.gz \
  | grep -v "^#" | grep -c "ANN="                        # should equal record count

# 6. Then: python scripts/08_call_haplotypes.py; write scripts/09_heteroplasmy_report.py.
```

### What we now understand (lessons to keep)

Two structural bugs shaped the architecture. Both are now recorded
in `CLAUDE.md` under "Known pitfalls + lessons learned":

**Bug A — Per-sample `bcftools call -mv` → merge cannot fill DP/AD
on REF cells.** Per-sample VCFs only contain rows for that sample's
variants. Merging unions positions, so any sample that wasn't called
at a given position has no row to contribute and `bcftools merge`
emits `.:.:.:.`. This was the original symptom and what motivated
stage 07's backfill. The structural fix is to not use per-sample
calling for the canonical at all — joint mpileup + call fills every
(POS × sample) cell in one pass.

**Bug B — `bcftools call -C alleles -T` silently drops secondary
ALTs at multiallelic sites if targets are post-norm-split.** The
targets file format for `-C alleles -T` is one row per position with
all ALTs comma-joined. If the targets file is built from a VCF that
has been `norm -m -any`-split (one row per ALT), multiallelic sites
appear as duplicate-POS rows. `call -C alleles -T` keys its lookup
on (CHROM, POS) and keeps only one entry per position — the
secondary entry is dropped without warning. This is what cost
17/931 sites in the session-11/12 backfill. The new architecture
sidesteps it entirely by never running `-C alleles -T` against a
post-norm VCF: joint -mv emits multiallelic sites as one row with
comma-joined ALTs natively, then stage 07's norm-split runs *after*
the call.

**Validation principle.** Switching the variant caller is not
something you can declare canonical by inspection. Joint and
per-sample callers can legitimately produce different variant sets
at marginal positions. The new caller's output must be diffed
against the prior canonical's position set; sites *lost* are a
regression, sites *gained* are fine. 05f writes this diff
automatically in `vcf/Fhet_mt_joint_vs_baseline_lost.tsv` and
flags it in the manifest.

### Artifacts changed this session

- `jobs/05f_joint_call.sh` — new (canonical caller); includes the
  baseline-vs-05f comparison block that writes
  `Fhet_mt_joint_vs_baseline_lost.tsv` and `_gained.tsv` and
  surfaces a regression warning in the manifest.
- `jobs/_archive/05d2_persample_call.sh` — archived.
- `jobs/_archive/05e2_merge_persample.sh` — archived.
- `jobs/_archive/07b_backfill_AD.sh` — archived.
- `scripts/06_snpeff_mac.sh` — input/output filename updates
  (`Fhet_mt_joint.vcf.gz` → `Fhet_mt_joint_ann.vcf.gz`).
- `scripts/07_cds_snps_norm_mac.sh` — full rewrite (backfill
  machinery removed; just CDS + SNPs + norm + rename).
- `CLAUDE.md` — new "Known pitfalls + lessons learned" section
  documenting bugs A and B and the validation principle;
  canonical-caller section rewritten; active-primary-task rewritten
  with explicit baseline-validation gate at step 2; compute-split
  rewritten; rule #1 status updated to require validation before
  freezing the canonical.
- `CHANGELOG.md` — this entry; session-12 entry annotated with the
  correct diagnosis for the 17-site loss.

## 2026-05-16 (session 12) — Stage 07 backfill patched (Step 4 was failing silently); script re-running, pending verification

> **NOTE (2026-05-17, corrected by session 13):** Patch 3 below blames
> the 17-site loss on `-Q 30 -q 30` filtering low-quality reads. **That
> diagnosis is incorrect.** The real cause was a targets-file format
> bug: the backfill emitted one row per (POS, ALT) for split multi-
> allelic sites, but `bcftools call -C alleles -T` keys on (CHROM, POS)
> and drops duplicate-position entries. All 17 lost sites are multi-
> allelic with two ALTs, and the lost ALT in every case is the minor
> allele. See the session-13 entry for evidence and the actual fix
> (architecture switch to joint call). Patch 3's edit was harmless but
> didn't address the real bug.

> **The session-11 rewrite of stage 07 ran but didn't actually produce a
> backfilled canonical.** Diagnosed three problems in Step 4 of
> `scripts/07_cds_snps_norm_mac.sh` (joint mpileup → call); patched all
> three; kicked off a re-run at end of session. Next session: verify the
> canonical lands at 931 records (matching preBackfill) with zero `DP=.`
> cells, then write the CHANGELOG verification entry and unblock 08 + 09.

### Symptom

After the session-11 rewrite of stage 07, the on-disk state on Mac was:

- `vcf/Fhet_MT_CDS.snps.split.preBackfill.vcf.gz` — fresh (Step 1–2 OK).
- `vcf/Fhet_MT_CDS_targets.tsv.gz` + `slim_bamlist.txt` — fresh (Step 3 OK).
- `vcf/Fhet_MT_CDS_backfilled.vcf.gz` — **missing entirely** (Step 4 never
  completed).
- `vcf/Fhet_MT_CDS.snps.split.vcf.gz` — still timestamped May 15 19:08
  (pre-rewrite version, not overwritten by Step 6 because Step 4
  never produced its input).
- Both the stale canonical and the new preBackfill carry `.:.:.:.` for
  REF cells, which is correct for preBackfill but masked the fact that
  the post-backfill canonical was missing.
- `vcf/Fhet_MT_CDS.snps_split_backfill_stats.txt` (May 16 20:13) was
  misleadingly named — its internal `# The command line was:` line
  shows it was run against the preBackfill file, not a backfilled file.

### Three patches to `scripts/07_cds_snps_norm_mac.sh`

1. **`ulimit -n 4096` at script top.** `bcftools mpileup -b` opens every
   BAM and its `.bai` simultaneously (≈ 290 file handles for 143 slim
   BAMs + reference + tabix targets). macOS's default per-process
   file-descriptor cap can be as low as 256 and is the most common cause
   of silent step-4 aborts on this hardware. Raise to 4096, well within
   macOS's hard cap.

2. **Split the targets file into two formats.** `bcftools mpileup -T` and
   `bcftools call -C alleles -T` expect *different* targets formats and
   the old script passed the same file to both:
     - mpileup: 4-column `CHROM POS REF ALT` — original `$TARGETS`.
     - call -C alleles: 3-column `CHROM POS REF,ALT` (REF and ALT joined
       by comma in a single 3rd column) — new `$TARGETS_CALL`,
       written to `vcf/Fhet_MT_CDS_call_targets.tsv.gz`, bgzipped + tabix.
   Step 3 now builds both. Step 4 reads `$TARGETS` for mpileup and
   `$TARGETS_CALL` for call. Passing the 4-col format to `call -C alleles`
   causes call to error out ("expected ALT after REF" or similar) and
   was the most-likely actual failure mode of Step 4 once the
   file-descriptor cap was clear.

3. **Dropped `-Q 30 -q 30` from mpileup; use bcftools defaults.** The
   session-11 script had `-Q 30 -q 30` hard-coded, which is stricter
   than 05d2 (the per-sample caller that produced the 931 canonical
   positions). 05d2 explicitly uses defaults (`-Q 13`, `-q 0`) — the
   header comment at `jobs/05d2_persample_call.sh` line 26 calls this
   out, and the session-11 CHANGELOG description of stage 07 (line 60)
   also lists no `-Q/-q` overrides. The hard-coded `-Q 30 -q 30` was
   effectively an inconsistency with the rest of the canonical pipeline,
   and a verification run showed it was dropping **17 of 931 sites** —
   positions whose only supporting reads have `13 ≤ BQ < 30` or
   `MAPQ < 30`, which under the strict filter produce no mpileup
   record, leaving `call -C alleles -T` nothing to force-call. The
   17 dropped positions (POS / REF / ALT): 2847 T A, 4913 G T,
   5566 T C, 8073 C G, 9112 T C, 9932 A T, 9977 A G, 10861 G C,
   11033 G A, 11069 G T, 11406 G T, 12394 C T, 12979 T A, 13441 T G,
   13577 T C, 13709 G C, 15173 T A.
   Removing the threshold puts the backfill on the same read population
   as 05d2 saw; per-cell `AD`/`DP` then reflect the full evidence and
   any Q/q/AD thresholding can be applied at stage 09 analysis time
   with the full per-cell evidence in hand.

   `--ploidy 1` was kept (correct for haploid mtDNA; does not lose
   heteroplasmic information — `AD` is computed from read counts and is
   independent of ploidy, only the emitted GT collapses to majority
   allele).

### Re-run status

Script kicked off at end of session 12 with:

```
cd ~/Projects/MT_Genomics_Cl_Ap2026/MT_Genomics2
conda activate SNP_env
bash scripts/07_cds_snps_norm_mac.sh 2>&1 | tee vcf/stage07_run.log
```

Returning next session to verify and complete documentation.

### Next-session pickup

1. Check `vcf/stage07_run.log` for the `=== DONE — CANONICAL OUTPUT
   PRODUCED (FROZEN) ===` block.
2. Confirm canonical record count = 931 (matches preBackfill):
   ```
   zcat vcf/Fhet_MT_CDS.snps.split.vcf.gz | grep -vc '^#'
   ```
3. Confirm `Cells with DP=.:` line reads `0` in the script's summary,
   or equivalently:
   ```
   bcftools query -f '[%DP\n]' vcf/Fhet_MT_CDS.snps.split.vcf.gz \
     | awk '$1=="."' | wc -l   # expect 0
   ```
4. Diff preBackfill ↔ canonical position set; expect zero loss:
   ```
   diff <(zcat vcf/Fhet_MT_CDS.snps.split.preBackfill.vcf.gz \
              | grep -v '^#' | awk '{print $2"\t"$4"\t"$5}' | sort) \
        <(zcat vcf/Fhet_MT_CDS.snps.split.vcf.gz \
              | grep -v '^#' | awk '{print $2"\t"$4"\t"$5}' | sort)
   ```
5. If all three checks pass: amend this CHANGELOG entry with the verified
   numbers, then proceed to stage 08 (`python scripts/08_call_haplotypes.py`)
   and write stage 09 per the session-11 design.
6. If checks fail: see the patches above for the three known failure
   modes and where to look in `scripts/07_cds_snps_norm_mac.sh`. The log
   tail will identify which step died.

### Artifacts changed this session

- `scripts/07_cds_snps_norm_mac.sh` — three patches above (ulimit;
  `$TARGETS_CALL` added + Step 3 builds both formats; mpileup Q/q
  defaults).
- `CLAUDE.md` — "Active primary task" updated to session-12 state,
  rule #1 updated, status of slim-BAM rsync (done) noted.
- `CHANGELOG.md` — this entry.

## 2026-05-15 (session 11) — Stage 06 fixed (Java + JAR version mismatch); stage 07 rewritten with BED CDS + sample rename + joint-pileup backfill; stage 08 / 09 design clarified

> **Mac-side 06 → 07 chain now produces a single canonical with per-cell
> DP/AD plus SnpEff ANN.** Stage 06 ran cleanly after fixing two Java-side
> bugs (system `java` was 1.8 vs SnpEff's required 11+; env's JAR was 5.2
> vs the 5.4 database built by `buildDbNcbi.sh`). Stage 07 was rewritten
> to (1) use `docs/mito_protein_coding.bed` directly instead of awk-parsing
> the GFF, (2) rename samples from full BAM paths to `${SAMPLE}_MT`, and
> (3) joint-pileup the slim BAMs at canonical positions to back-fill
> per-cell AD/DP for samples whose per-sample call was REF — necessary for
> heteroplasmy detection, since 05d2 left those cells as `.:.:.:.`. The
> resulting `Fhet_MT_CDS.snps.split.vcf.gz` is the single deliverable for
> both stage 08 (haplotype calling at 0.7) and the new stage 09
> (heteroplasmy report at 0.1 – 0.7).

- **Stage 06 SnpEff fixes (`scripts/06_snpeff_mac.sh`):**
  - **Java version mismatch.** System `java` is 1.8 (class file 52) but
    SnpEff 5.2+ requires 11+ (class file 55), producing
    `UnsupportedClassVersionError`. Fix: explicit
    `JAVA_BIN="${HOME}/micromamba/envs/SNP_env/lib/jvm/bin/java"` (Java
    23 from `SNP_env`, runs all SnpEff JARs); replaced bare `java` calls
    and the version-print with `"$JAVA_BIN"`. Pre-flight checks the path
    and fails fast with a clear error.
  - **SnpEff JAR vs DB version mismatch.** The env's
    `snpeff-5.2-1/snpEff.jar` refused the `buildDbNcbi.sh`-built 5.4
    database: `Database version : '5.4', Program version : '5.2'`. Fix:
    point `SNPEFF_JAR` at `${HOME}/snpEff/snpEff.jar` (the standalone
    install that built the DB; lives next to `data/` and `snpEff.config`
    under one tree).
  - **`DB_NAME` → `NC_012312.1`** (NCBI-built via
    `~/snpEff/scripts/buildDbNcbi.sh`), with the required
    `NC_012312.1.codonTable : Vertebrate_Mitochondrial` line in
    `snpEff.config` documented in the header.
  - **Input / output filenames** updated to
    `Fhet_mt_persample_merged.vcf.gz` (was `Fhet_mt_variantsAD.vcf.gz`,
    an artifact of the joint-call era).
  - **rsync host** in the header docstring corrected to
    `t2.idsc.miami.edu` (was stale `scc1.bu.edu`); now pulls
    `.vcf.gz` + `.csi` in one call.

- **Stage 07 rewritten (`scripts/07_cds_snps_norm_mac.sh`):**
  - **BED-based CDS restriction** using `docs/mito_protein_coding.bed`
    (13 mt protein-coding genes, hand-curated; copied from
    `additional_info/previous_data/`). bcftools auto-detects the `.bed`
    extension. Drops the awk + bgzip + tabix block from the previous
    version. Verified identical to the GFF-derived CDS (same 13 features,
    same coordinates modulo BED's 0-based-half-open vs GFF's 1-based-
    inclusive convention).
  - **Sample rename.** Slim BAMs have no `@RG SM:` tag, so `bcftools
    mpileup` upstream fell back to BAM filename and the merged VCF
    inherited full BAM paths as sample names. Single sed:
    `s|.*/||; s/_0_MT_only\.bam$/_MT/` converts
    `/projectnb/.../10_0_MT_only.bam` → `10_MT`.
  - **Per-cell AD/DP backfill via joint pileup.** Builds a tabix-indexed
    targets TSV (CHROM, POS, REF, ALT) from the post-norm intermediate,
    then `bcftools mpileup -T targets -b slim_bamlist -a AD,DP` →
    `bcftools call -m -A -C alleles -T targets --ploidy 1`. Force-calling
    the canonical alleles means every (POS × sample) cell has explicit
    AD = (REF_count, canonical_ALT_count) and DP = total depth, including
    samples whose per-sample call (05d2) was REF. Heteroplasmy at
    canonical alts now visible as nonzero ALT_count in those cells.
  - **ANN transfer** via `bcftools annotate -a intermediate -c INFO/ANN
    backfilled`. SnpEff effect predictions from stage 06 survive the
    backfill so the canonical carries both per-cell AD/DP *and* `ANN`.
  - **Pre-requisite added: slim BAMs on Mac** under
    `MT_Genomics2/MT_only_bams/` (one-time rsync from T2, ~10–15 GB).
    Pre-flight checks for ≥ 143 slim BAMs and fails with the rsync
    command in the error message if missing.
  - **Final canonical artifact:** `Fhet_MT_CDS.snps.split.vcf.gz`
    (FROZEN). Intermediate `Fhet_MT_CDS.snps.split.preBackfill.vcf.gz`
    kept for QC comparison (pre-backfill, REF cells empty).

- **`jobs/05e2_merge_persample.sh` sample-rename regex updated** to use
  the same `sed -e 's|.*/||' -e 's/_0_MT_only\.bam$/_MT/'` recipe. Future
  re-runs of 05e2 produce a merged VCF with clean sample names from the
  start. The existing `Fhet_mt_persample_merged.vcf.gz` on T2 still has
  full-path sample names; stage 07 fixes them downstream so no immediate
  re-run is required.

- **`jobs/07b_backfill_AD.sh` added** as a Triton-2-side alternative to
  the Mac-side backfill now in stage 07. Same recipe (joint mpileup over
  slim BAMs at canonical positions, force-call canonical alleles, sample
  rename), but runs on T2 against the slim BAMs already there — useful
  when rsyncing 143 slim BAMs to Mac isn't desirable. Both paths produce
  the same canonical.

- **Stage 08 / 09 design clarified.** Heteroplasmy and haplotype are NOT
  alternative interpretations of the same call — they are orthogonal
  signals to pull from the same canonical:
    - **Stage 08 (existing `scripts/08_call_haplotypes.py`)** stays at
      `AD_alt/DP > 0.7`. Answers "which haplotype clade is each sample
      in?" Per cell: 1 if effectively homoplasmic for alt, else 0.
      Output is the haplotype matrix + per-sample haplotype string.
      Used for phylogeny / population-structure work elsewhere.
    - **Stage 09 (new, `scripts/09_heteroplasmy_report.py` — to be
      written next)** operates on the same canonical, picks up what
      stage 08 discards:
        * per-cell filter: 0.1 ≤ AD_alt/DP < 0.7 AND AD_alt ≥ 50
          (the depth floor sits comfortably above the sequencing-error
          band even at 10 % VAF)
        * per-site classifier: n_carriers; transmitted (> 1) vs private
          (= 1). Sites with multiple heteroplasmic carriers strongly
          suggest maternal-line transmission of the heteroplasmic
          state (independent de novo at the same site in multiple fish
          is statistically implausible).
        * outputs: per-site het table (transmitted + private, separate
          files); per-sample het burden (n total / transmitted /
          private + mean AF); optional pairwise shared-heteroplasmy
          adjacency for maternal-lineage detection.
  This supersedes the session-10 plan that proposed revising stage 08's
  threshold to 0.1 — that was based on a misreading. Stage 08 stays as-is;
  stage 09 is additive. CLAUDE.md "Active primary task" updated accordingly.

- **Artifacts created / updated this session:**
  - `scripts/06_snpeff_mac.sh` — `JAVA_BIN` var, `SNPEFF_JAR` →
    standalone install, `DB_NAME` → `NC_012312.1`, input/output filenames,
    rsync host.
  - `scripts/07_cds_snps_norm_mac.sh` — full rewrite (BED + rename +
    backfill + ANN transfer).
  - `jobs/05e2_merge_persample.sh` — sample-rename regex updated.
  - `jobs/07b_backfill_AD.sh` — new (T2-side alternative to the Mac-side
    backfill in 07).
  - `docs/mito_protein_coding.bed` — 13-gene BED placed in canonical
    location (copy of `additional_info/previous_data/mito_protein_coding.bed`).

- **Next-session pickup (in order):**
  1. rsync slim BAMs to Mac:
     ```
     rsync -avP \
       dcrawford@t2.idsc.miami.edu:/projectnb/dcrawford/MT_Genomics2/MT_only_bams/ \
       ~/Projects/MT_Genomics_Cl_Ap2026/MT_Genomics2/MT_only_bams/
     ```
  2. `bash scripts/07_cds_snps_norm_mac.sh` → canonical
     `vcf/Fhet_MT_CDS.snps.split.vcf.gz` with full per-cell DP/AD + ANN.
  3. QC: `bcftools query -f '[%DP\n]' vcf/Fhet_MT_CDS.snps.split.vcf.gz
     | awk '$1=="."' | wc -l` — should print `0` (every cell has a DP
     after backfill).
  4. `python scripts/08_call_haplotypes.py` — haplotype matrix at the
     existing 0.7 binary threshold (no revision).
  5. Write `scripts/09_heteroplasmy_report.py` per the design above.
     Decide carrier-count threshold (default: `n_carriers > 1` for
     "transmitted") and output schema (per-site table + per-sample
     burden, with optional adjacency matrix).
  6. Compare stage 08 haplotype matrix to prior ~243-variant haplotype
     output; compare stage 09 transmitted-het count to the historical
     868 low-AF sites from `Fhet_mt_fullAD.vcf.gz`.

## 2026-05-15 (session 10) — Diagnosis fully sealed: 1128 SNPs from new pipeline (recovers historical 1133); slim-BAM + 05d2/05e2 chain in production

> **The number lineage is closed.** The bcftools-1.23.1 verification BSUB
> from session 9 cont. returned **289 SNPs from 3 samples**, squarely
> inside the 200–300 SNP band predicted by the bcftools 1.22 baseline
> (284 SNPs on the same 3 BAMs). Full-panel re-call using slim mt-only
> BAMs feeding 05d2 (per-sample call) → 05e2 (merge), all on
> locally-built bcftools 1.23.1, produced **1128 SNPs / 143 samples /
> ts/tv 8.17** — within 0.4 % of the historical `merged_144.vcf.gz`
> baseline (1133 SNPs / 144 samples / ts/tv 7.92). The
> 1142 / 152 / 227 / 284 number lineage is fully reconciled. Canonical
> caller as of today: **slim BAMs → 05d2 → 05e2**.

- **Result vs historical baseline:**

      metric                       new (1.23.1, slim BAMs)   historical merged_144 (1.22)
      samples                                          143                              144
      records                                         1144                             1140
      SNPs                                            1128                             1133
      indels                                            16                                7
      ts/tv (1st alt only)                            8.17                             9.12
      ts/tv (all alts)                                8.17                             7.92
      multiallelic SNP sites                             0                               29

  Zero multiallelic sites confirms per-sample `bcftools norm -m -any`
  worked cleanly — every row is one (POS, ALT) pair, exactly what stage 08
  expects. The recovered ts/tv band (8.17) is squarely in the
  vertebrate-mtDNA range, so the 1128 are real biology.

- **AF spectrum is sharply bimodal** — 592 singletons (AC=1) and 536
  near-fixed-alt at AF ≈ 0.993, with essentially nothing in between. The
  singleton count is ~10× the historical merged_144 (~55), most likely the
  `--max-depth 10000` carried over from the diagnostic scripts (historical
  recipe used bcftools' default per-sample cap of 250, which throttles
  low-AF evidence at high-coverage sites). Either pattern is biologically
  defensible; the carrier-count filter planned for stage 08 (≥ 10 % of
  samples = ≥ 15 carriers) discards all 592 singletons trivially.

- **Slim-BAM pre-extraction stage added.** `jobs/BSUB_Slim_BAM_mt.sh` —
  LSF array `[1-143]%30` that uses
  `samtools view -F 2308 ... NC_012312.1` to pull mt-mapped primary
  alignments from each 15–19 GB WGS BAM into a ~50–100 MB mt-only BAM at
  `/projectnb/dcrawford/MT_Genomics2/MT_only_bams/${SAMPLE}_MT_only.bam`
  (+ `.bai`). BAI-driven seek into the 19 GB inputs is fast — total
  array runtime well under an hour. Cuts subsequent mpileup walltime ~100×.
  Replaces a buggy first draft (filename mismatch between `samtools view`
  and `samtools index`, missing dirname on the BAM input, wrong array
  size `[1-144]`, no samtools on PATH).

- **Caller chain v2 (05d2 / 05e2) now reads slim BAMs.** Recipe identical
  to 05d/05e (per-sample
  `bcftools mpileup -f REF BAM -a AD,DP --max-depth 10000`
  → `bcftools call -mv --ploidy 1` → `bcftools norm -m -any`, then
  `bcftools merge -m none` across all 143), but inputs come from
  `${PROJECT_ROOT}/MT_only_bams/${SAMPLE}_MT_only.bam` and the binaries
  are the locally-built 1.23.1 instead of the env's pinned 1.6.

- **`jobs/config.sh` updated.** Prepends `$HOME/software/local/bin` to PATH
  and `$HOME/software/local/lib` to LD_LIBRARY_PATH so any script sourcing
  config.sh picks up the new bcftools / samtools / htslib by default.
  Documented ordering caveat: scripts that `conda activate` AFTER sourcing
  config.sh must re-export to defeat the env's bin-path push. The
  slim / 05d2 / 05e2 scripts skip conda activation entirely, so the
  caveat doesn't bite them.

- **Artifacts created / updated this session:**
  - `jobs/config.sh` — local-bin export added at end with ordering caveat
    documented.
  - `jobs/BSUB_Slim_BAM_mt.sh` — clean rewrite (see bug list above).
  - `jobs/05d2_persample_call.sh` — reads from slim BAMs in
    `${PROJECT_ROOT}/MT_only_bams/`. Header documents the new
    `Slim → 05d2 → 05e2` chain.
  - `jobs/05e2_merge_persample.sh` — header + run-manifest text updated
    to record the slim-BAM input chain (so the manifest in `vcf/` is
    self-documenting for the eventual methods write-up).
  - `jobs/Slim_BAM_mt.sh` — near-duplicate of the BSUB script, can be
    deleted as housekeeping.

- **Output landed at:**
    - `/projectnb/dcrawford/MT_Genomics2/vcf/Fhet_mt_persample_merged.vcf.gz`
    - `/projectnb/dcrawford/MT_Genomics2/vcf/Fhet_mt_persample_merged.vcf.gz.csi`
    - `/projectnb/dcrawford/MT_Genomics2/stats/Fhet_mt_persample_merged_stats.txt`
    - run manifest alongside in `vcf/Fhet_mt_persample_merged_run_manifest.txt`.

- **Next-session pickup (in order):**
  1. rsync the merged VCF + `.csi` from Triton 2 to Mac `vcf/`:
     ```
     rsync -avP \
         dcrawford@t2.idsc.miami.edu:/projectnb/dcrawford/MT_Genomics2/vcf/Fhet_mt_persample_merged.vcf.gz \
         dcrawford@t2.idsc.miami.edu:/projectnb/dcrawford/MT_Genomics2/vcf/Fhet_mt_persample_merged.vcf.gz.csi \
         ~/Projects/MT_Genomics_Cl_Ap2026/MT_Genomics2/vcf/
     ```
     Symlink as `Fhet_mt_variantsAD.vcf.gz` so the unmodified Mac scripts
     pick it up.
  2. `bash scripts/06_snpeff_mac.sh` — SnpEff annotation against the
     custom `Fhet_MT` database.
  3. `bash scripts/07_cds_snps_norm_mac.sh` → canonical
     `vcf/Fhet_MT_CDS.snps.split.vcf.gz` (FROZEN once produced, per
     CLAUDE.md rule 1).
  4. Revise `scripts/08_call_haplotypes.py`: replace `AD_alt/DP > 0.7`
     with `AD_alt/DP ≥ 0.1` per sample plus a carrier-count filter
     (≥ 10 % of samples = ≥ 15 carriers, matching historical practice).
  5. Compare resulting haplotype calls to the prior ~243-variant
     haplotype output.
  6. Update `docs/02_calling_architecture.md` with an addendum:
     bcftools-version was the dominant variable on this dataset. Include
     the 3-sample version comparison (1.6 / 1.22 / 1.23.1 → 6 / 284 / 289
     SNPs) for the eventual methods write-up.
  7. Optional housekeeping: `rm jobs/Slim_BAM_mt.sh` (duplicate of
     BSUB version).

## 2026-05-15 (session 9 cont.) — bcftools / samtools / htslib 1.23.1 built from source on Triton 2; verification call submitted as BSUB

> **Tool upgrade done; verification pending.** bcftools 1.23.1, samtools
> 1.23.1, and htslib 1.23.1 are now installed at
> `$HOME/software/local/bin/` on Triton 2, sidestepping the linux-ppc64le
> bioconda constraint that pinned `mito_genomics` to bcftools 1.6.
> A BSUB job re-runs the 3-sample diploidA recipe with the new binary;
> stats will tell us whether bcftools 1.6 was the sole cause of the
> variant-count collapse (expected ~200–300 SNPs from the same 3 BAMs).

- **Build outcome.** Source tarballs (1.23.1) for htslib, samtools, and
  bcftools were downloaded into `$HOME/src/`, configured with
  `--prefix=$HOME/software/local --with-htslib=$HOME/software/local`, and
  built with `make -j4`. Binaries land in `$HOME/software/local/bin/`,
  libraries in `$HOME/software/local/lib/`. PATH updated in `~/.bashrc`:

      export PATH="$HOME/software/local/bin:$PATH"
      export LD_LIBRARY_PATH="$HOME/software/local/lib:${LD_LIBRARY_PATH:-}"

  Verified versions on Triton 2:

      bcftools 1.23.1   (using htslib 1.23.1)
      samtools 1.23.1
      htsfile 1.23.1

  Matches the Mac install for cross-platform reproducibility.

- **Three build gotchas hit, all resolved and patched into the script
  (`Again_May/test/rerun/build_bcftools_t2.sh`):**

  1. **htslib configure → `libbzip2 development files not found`.**
     Triton 2 doesn't have `libbz2-devel` or `xz-devel` system-wide. Fix:
     add `--disable-bz2 --disable-lzma` to htslib configure. These flags
     drop CRAM-only compression support; BAM, VCF, BCF, and gzip work
     normally (zlib is present and used).
  2. **Build script halted silently after the toolchain check.**
     Original script used `set -euo pipefail`. The line
     `gcc --version | head -1` exits 141 under pipefail (head closes
     after one line → upstream SIGPIPE), `set -e` then killed the
     script. Same gotcha we hit in `diagnose_bam_loss.sh` during the
     earlier diagnostic. Fix: switch the script to `set -eu` (no
     pipefail). Documented in a comment so it doesn't get re-introduced.
  3. **samtools configure → `curses development files not found`.**
     Triton 2 doesn't have `ncurses-devel`. Fix: add `--without-curses`
     to samtools configure. Drops the interactive `samtools tview` text
     UI — we don't use tview anyway. bcftools doesn't have a curses
     dependency and configures cleanly without extra flags.

- **Verification BSUB submitted** (`~/run_test_call_v123.sh`,
  pending at session end). Recipe identical to the historical Oct-2025
  call that produced 284 SNPs on the same 3 samples with bcftools 1.22,
  except the new bcftools 1.23.1 binary is used:

      bcftools mpileup -f refs/Fhet_MT.fasta -b bamlist_test.txt \
          --threads 4 -a AD,DP -Q 30 -q 30 -d 100000 -Ou \
        | bcftools call --threads 4 -mv -A -Oz -o test_diploidA_v123.vcf.gz

  Resources: 4 cores, 8 GB mem, 4 h wall. Why BSUB instead of login
  node: the 3 BAMs are 15–19 GB each (retain unmapped WGS reads), and
  mpileup walks all of them at high depth (`-d 100000`). On the login
  node this stretches well past what Triton 2 policy tolerates for
  interactive use. Output VCF lands at
  `/projectnb/dcrawford/MT_Genomics2/AD_again/test_diploidA_v123.vcf.gz`,
  stats appended to BSUB stdout at
  `log2/test_call_v123_<JOBID>.out`.

- **Expected outcome on resume.**
    * 200–300 SNPs → bcftools 1.6 was the sole cause. Diagnosis sealed.
      Proceed to full 143-BAM per-sample → merge with the new binary.
    * 1000+ SNPs → version was the cause AND the -A + default-diploid
      recipe surfaces the heteroplasmy tail. Even better — closer to
      the historical 1142 from just 3 samples.
    * Still ≤10 SNPs → something else is also going on. Drop back to
      the diagnostic — but session-9's pileup data shows abundant alt
      reads in the BAMs, so this outcome would be very surprising.

- **Performance note for the full 143-BAM run.** With the new binary on
  a compute node, expect the per-sample call step (one BAM at a time)
  to be the better approach than joint mpileup over all 143 BAMs. The
  retained-unmapped-reads bloat in the BAMs (15–19 GB each) is the
  dominant time cost; merging the small per-sample VCFs at the end is
  fast. The 05d (LSF array, 1-143%24) + 05e (single merge job) split
  from session 7 is the right structure — just needs the PATH update
  so the array tasks pick up 1.23.1 instead of 1.6. An optional
  optimization: pre-extract mt-mapped reads only
  (`samtools view -b $BAM NC_012312.1 -o ${sample}_MT_only.bam`) before
  per-sample mpileup. Cuts each BAM from 15–19 GB to ~50–100 MB.

- **Artifacts created / updated this session:**
  - `Again_May/test/rerun/build_bcftools_t2.sh` — three-stage build
    script (htslib → samtools → bcftools), patched for all three
    gotchas above. Idempotent: re-running skips already-downloaded
    tarballs and re-runs configure/make. Safe to re-run on Triton 2
    if you need to rebuild against a newer version.
  - `~/run_test_call_v123.sh` (Triton 2) — BSUB wrapper for the
    3-sample verification call. Job submitted at session end; result
    pending.

- **Next-session pickup (in order):**
  1. On Triton 2: `bjobs -a | grep test_call_v123` — find the job
     state. If COMPLETED, `tail -30
     /projectnb/dcrawford/MT_Genomics2/log2/test_call_v123_<JOBID>.out`
     for the SN stats. If RUN, wait. If EXIT, look at .err.
  2. Interpret the SN block:
       * 200–300 SNPs → diagnosis sealed. Continue to step 3.
       * 1000+ SNPs → diagnosis sealed and heteroplasmy is captured.
         Even better; continue to step 3.
       * ≤10 SNPs → re-open. Look at the actual records in the new
         VCF and at the BSUB stderr.
  3. Update `jobs/config.sh` to prepend the new bin path to PATH so
     all BSUB scripts pick up bcftools/samtools 1.23.1 without needing
     a new conda env:
       ```
       export PATH="$HOME/software/local/bin:$PATH"
       ```
     Then commit this change to the repo.
  4. (Optional but recommended) Pre-extract mt-only reads from each
     of the 143 AD_again BAMs to slim them from 15–19 GB to ~50–
     100 MB. Single LSF array job:
       ```
       samtools view -b $BAM NC_012312.1 -o ${sample}_MT_only.bam
       samtools index ${sample}_MT_only.bam
       ```
     Speeds up downstream mpileup substantially.
  5. Re-run `jobs/05d_persample_call.sh` (the LSF [1-143]%24 array,
     historical per-sample-call recipe with the modern bcftools)
     followed by `jobs/05e_merge_persample.sh`. Expected output:
     500–1300 SNPs in the merged VCF. If close to the historical
     1142, the 1142 / 152 / 227 number lineage is closed.
  6. SnpEff (stage 06) + CDS restrict (stage 07) on the new merged
     VCF.
  7. Revise stage 08's `AD_alt/DP > 0.7` haplotype rule to capture
     heteroplasmy: drop to `≥ 0.1` per sample and require >1
     individual (or another reasonable carrier-count filter,
     historical user practice was ≥10 % AD with ≥10 % of individuals
     as carriers).
  8. Compare resulting haplotype calls to the prior ~243-variant
     haplotype output. Expect closer alignment once heteroplasmy is
     captured.

## 2026-05-14 (session 9) — Reference and BAMs both cleared; bcftools 1.6 in `mito_genomics` env identified as the cause of the variant-count collapse

> **The session-8 reference-identity hypothesis is refuted, and the session-8
> BAM-loss hypothesis is also refuted.** Direct pileup of the new AD_again
> BAMs at the top-20 historical high-AF positions shows abundant, high-quality
> alt-supporting reads at every position (DP ~ 750–1700, near-100 % alt at
> AF=1 sites). The reads are present and the reference is unchanged. The
> failure is between the BAM and the VCF, in the `mito_genomics` env's pinned
> **bcftools 1.6**. Modern bcftools ( ≥1.10, ideally 1.21+ ) is the fix.

- **Reference is the same.** A Python REF-base comparison between the historical
  `Fhet_mt_fullAD.vcf.gz` (1142 records, joint `-mv -A` default-diploid call
  from Oct 17 2025, bcftools 1.22, against `/SSM_Mito/Fh_MT_ref/Fhet_MT.fasta`)
  and the current `refs/Fhet_MT.fasta` (md5
  `46da77890fefaa6ead6103789d7dad76`) matched **1142 / 1142** REF bases. If
  the references differed at the variant positions, we would have seen a
  wave of `VCF_REF=X / FASTA_BASE=Y` mismatches. We saw zero. The Mac-backup
  reference now in `refs/` is sequence-identical (at variant positions) to
  whatever Jul-2025 reference produced the 1133-SNP `merged_144` baseline.

- **Historical 1142-SNP file is heteroplasmy-rich, not high-AF-rich.** AF
  bucketing (Python, AF computed as AC/AN since the VCF carries no `AF=` INFO
  field): 868 sites at AF ≤ 0.05, 258 sites at 0.05 < AF ≤ 0.5, 4 sites at
  0.5–0.95, only 12 sites at AF ≥ 0.95. ts/tv (first alt) = 8.93. Only 2
  sites are truly multiallelic at AF ≥ 0.05 (the apparent 1135 multiallelic
  count is `-A` plus default-diploid producing low-AD second-alt artifacts).
  Session-7 notes describing the historical baseline as "89 % at AF ≥ 0.95"
  were referring to a different file: `merged_144.vcf.gz` (per-sample → merge,
  Jul 2025), not `Fhet_mt_fullAD.vcf.gz` (joint, Oct 2025). Both files
  represent the same panel but the joint+`-A`+diploid call surfaces a
  heteroplasmy-and-rare-variant signal that the per-sample → merge does not.

- **Three-sample tests on Triton 2**, all with the same `bcftools mpileup`
  filters (`-Q 30 -q 30 -d 100000 -a AD,DP`) against `refs/Fhet_MT.fasta`:

      Run                                          BAMs              bcftools  --ploidy  -A   SNPs
      Fhet_mt_fullAD (Oct 17 2025, hist baseline)  /SSM_Mito/MT_bam   1.22     default   yes  1142
      original_test  (Oct 17 2025, 3 BAMs)         /SSM_Mito/MT_bam   1.22     1         no    284
      rerun          (May 13 2026, 3 BAMs)         /AD_again/MT_bam   1.6      1         no      3
      diploidA       (this session, 3 BAMs)        /AD_again/MT_bam   1.6      default   yes     6

  The diploidA test was designed to isolate the calling-parameter axis: same
  3 samples, same reference, same Q/q/d, but matching the historical recipe
  (default diploid + `-A`). Result: 6 records (positions 5, 3335, 15684,
  15685, 15686, 15690 — i.e. only the very edges of the mitogenome + 3335).
  Removing `--ploidy 1` and adding `-A` did NOT recover the historical
  signal, so the heteroplasmy-and-haploid hypothesis from earlier this
  session is necessary but not sufficient.

- **Decisive diagnostic: `samtools mpileup` at the top-20 historical AF
  positions in `AD_again/MT_bam_sam/10_0_MT.bam`.** Pulled the top-20 sites
  from the historical VCF (computed AF = AC/AN since `AF=` is absent), ran
  strict (`-Q 30 -q 30 -A`) and permissive (`-Q 0 -q 0 -A`) pileups on the
  new BAM. Spot check at AF=1.0 historical sites:

      POS    REF→ALT   strict DP   strict bases (truncated)
      9001   T→C       1113        100 % C/c
      8360   G→A       1070        ~100 % A/a
      4957   C→T       1576        ~100 % T/t
      3124   G→A       1557        ~100 % A/a
      14909  A→G       1435        ~100 % G/g
      11411  T→C       1057        ~100 % C/c

  The BAMs unambiguously carry the alt reads at high coverage and high
  quality. Whatever filtering bcftools 1.6 is applying between mpileup and
  call is what produced the 3-SNP / 6-SNP / 152-SNP / 227-SNP collapse.
  Saved as `Again_May/TMP/raw_pileups_10_0.txt` and the joined position
  table `Again_May/TMP/top.tsv` (on Mac), `/projectnb/dcrawford/MT_Genomics2/
  AD_again/TMP/*` on Triton 2.

- **bcftools 1.6 is the proximate cause.** The env is pinned to bcftools 1.6
  (CHANGELOG 2026-05-08, linux-ppc64le biobuilds constraint). bcftools 1.10
  rewrote the `mpileup` code path; 1.6 has known issues with high-depth
  pileups on small references, BAQ defaults, anomalous-pair exclusion, and
  per-sample depth caps that interact badly with WGS-aligned-to-mt-only
  BAMs. The historical 1142-SNP call used bcftools 1.22. The Oct-2025
  3-sample test (284 SNPs) also used bcftools 1.22. Both 1.22 calls were
  rich; both 1.6 calls (rerun, diploidA) collapsed. This is the
  experimentally-isolated variable.

- **Next-session action: upgrade bcftools and re-call.** Three viable paths:

  1. **Build bcftools 1.21 from source on Triton 2.** htslib + bcftools, ~10
     minutes. The C code is portable to linux-ppc64le; only bioconda's
     packaging is constrained. Drop the new binaries into
     `$HOME/software/local/bin/` and prepend that path in `jobs/config.sh`
     so all BSUB scripts pick it up without needing a new conda env.
  2. **Run the calling step on the Mac.** Modern bcftools is one
     `conda create -n bcf -c bioconda bcftools samtools` away. BAMs are
     15–19 GB each; rsync just `10_0`, `102_0`, `103_0` (plus `.bai`) plus
     the reference, repeat the diploidA recipe locally, and confirm we
     recover ~200–300 SNPs from those three BAMs. This is the fastest
     verification.
  3. **Singularity / container.** If Triton 2 supports it,
     `singularity exec docker://quay.io/biocontainers/bcftools:1.21--*
     bcftools call ...` avoids the build step entirely.

  Verification target: rerun-equivalent recipe with modern bcftools should
  produce ≥200 SNPs on the same 3 BAMs. If it does, queue a per-sample
  call+merge of all 143 BAMs using the new binaries — that is then the
  final call set, and the 1142 / 152 / 227 number lineage is closed.

- **Artifacts created this session:**
  - `Again_May/test/rerun/BSUB_mpile_AD_test_diploidA.sh` — Triton-2 BSUB
    script that replicates the historical Oct-2025 recipe (default diploid,
    `-A`) on the new AD_again BAMs. Output landed under
    `/projectnb/dcrawford/MT_Genomics2/AD_again/test_diploidA/`. Also
    documents the BAMLIST / OUTDIR bugs fixed from the earlier rerun script.
  - `Again_May/test/rerun/diagnose_bam_loss.sh` — Triton-2 diagnostic that
    pulls top-N historical AF positions and dumps strict + permissive
    samtools pileups at each, joining them into a comparison table. User
    adapted this to `diagnose_bam_loss_2.sh` for Triton-2 paths (commented
    out `mktemp`, switched to a persistent `TMP=/projectnb/dcrawford/
    MT_Genomics2/AD_again/TMP`). Important fixes during development:
      * AF parsing originally looked for `AF=` in INFO; the historical VCF
        carries only `AC=` / `AN=`, so AF must be computed as max(AC) / AN.
      * `set -euo pipefail` + `awk | sort | head` triggers SIGPIPE-as-error
        because head closes early. Switched to `set -eu` (no pipefail).
        Same gotcha applies to the AF-bucket awk at the bottom of
        `BSUB_mpile_AD_test_diploidA.sh`.
  - `Again_May/test/rerun/run_diag.sh` — minimal BSUB wrapper that activates
    the env, runs the diagnostic across 10_0 / 102_0 / 103_0 / 77_0, and
    writes one combined log to `log2/bam_diag_<JOBID>.out`.
  - `Again_May/test/rerun/diagnose_bam_loss_mac.sh` — Mac variant of the
    diagnostic (different default paths, instructions for setting up a
    `bcf` conda env via bioconda). Not yet exercised; relevant when option
    2 of the next-session action gets executed.

- **Open biological question, parked.** The historical 1142-SNP file's rare-
  variant tail (868 sites at AF ≤ 0.05) and the user's downstream
  per-individual variant counts (~200 synonymous + ~18-20 nonsynonymous in
  coding regions) imply substantial heteroplasmy in the panel. Once modern
  bcftools recovers the variant calls, stage-08's `AD_alt/DP > 0.7`
  threshold will need a heteroplasmy-aware revision: user noted the
  historical analysis used either ≥10 % AD with ≥10 % of individuals as
  carriers, or ≥10 % AD with ≥50 carrier individuals. The 0.7 threshold
  was a design decision from session 3 (2026-05-08) before the heteroplasmy
  signature was understood. Revisit after the variant set is regenerated.

- **Next-session pickup (in order):**
  1. Pick a bcftools upgrade path (build from source on Triton 2 is the
     most permanent fix; running the verification on the Mac is the fastest).
  2. Verify on the 3-sample test BAMs: same `bcftools mpileup | bcftools
     call -mv -A` recipe with bcftools ≥ 1.21 should produce ~200–300 SNPs.
     If it does, the diagnosis is sealed.
  3. Re-run 05d/05e (per-sample → merge across all 143 AD_again BAMs) with
     the new bcftools binary; expected output 500–1300 SNPs.
  4. SnpEff (stage 06) + CDS restrict (stage 07) on the new merged VCF.
  5. Revise stage 08 haplotype caller for heteroplasmy: replace
     `AD_alt/DP > 0.7` with `AD_alt/DP ≥ 0.1` + cross-sample carrier-count
     filter ( ≥ 10 % of samples or absolute carrier count threshold).
  6. Compare resulting haplotype calls to the prior ~243-variant haplotype
     output; expect closer alignment once heteroplasmy is captured.

## 2026-05-11 (session 8) — 05d/05e completed at 227 SNPs (not 1130); diagnostic localizes cause to reference identity; AD_again realignment + two new per-sample call+merge scripts written

> **The session-7 hypothesis is refuted by the 05d/05e result.** Per-sample
> → merge ran cleanly and produced **227 SNPs / 238 records / 143 samples**,
> not the predicted ~1130. Deeper diagnostics localize the cause to the
> **reference**: the BAMs do not carry alt reads at the historical
> near-fixed-alt positions, even in the most divergent sample. The Jul-2025
> FASTA that produced the 1133-SNP signal appears to be either lost or
> sufficiently different from the recovered Mac-backup FASTA that the
> within-panel calls now reflect only true within-panel variation.

- **05d/05e merge output** (`vcf/Fhet_mt_persample_merged_stats.txt`):
  n_samples=143, n_records=238, n_SNPs=227, n_indels=11, ts/tv=6.83.
  AF spectrum is two-tier only — 157 singletons (AC=1) at AF=0 and 70
  near-fixed-alt at AF=0.993. The historical ~70 SNPs at AF ∈ [0.05, 0.95]
  are completely missing, which on an admix panel with ~50% reference-
  similar and ~50% divergent members is biologically wrong.

- **Per-sample variant-count distribution** (143 samples; `bcftools view |
  grep -v "^#" | wc -l` per per-sample `_norm.vcf.gz`):

      0 variants: 36 samples       7:           2
      1:          42                8:           1
      2:          36               13:           1
      3:          13               19:           1
      4:           7               81:           1
      5:           1              182:           1   (sample 77_0)
      6:           1

  134/143 samples have ≤4 variants. With a 50/50 admix panel against a
  reasonably divergent reference, the expected distribution should be
  bimodal with the second mode at 50–200 variants. We see one mode at 0–4
  and two outliers. That's the population-matches-reference signature.

- **Decisive pileup at historical variant positions 132/135/153/161/194
  in sample 77_0** (the 182-variant outlier; coverage 60–95 per position).
  All reads at all five positions show `.`/`,` (match-to-reference). The
  single exception is **one** `c` (one alt-supporting read out of 95) at
  position 194 — far below any per-sample calling threshold. The BAMs
  are not carrying the alt reads at the historical variant positions
  even in the most divergent sample. Saved at `docs/77mpileup.txt`.

- **Reference identity check.** `md5sum refs/Fhet_MT.fasta`
  = `46da77890fefaa6ead6103789d7dad76`. The only other `Fhet_MT.fasta`
  on Triton 2 has the same md5. The original
  `/projectnb/dcrawford/SSM_Mito/Fh_MT_ref/` directory remains empty
  (per the 2026-05-07 inventory). BAM `@SQ` headers do not record `M5:`
  so the historical reference identity can't be verified from the BAMs.
  Sequence-only md5 of `refs/Fhet_MT.fasta` =
  `3f460763fdfcf105c3e3144ca6784f55`; not yet cross-checked against
  canonical NCBI `NC_012312.1`.

- **`BSUB_1_mt_align_pipping_array_again.sh` submitted (Doug, this
  session)** — realigns 144 samples from
  `/projectnb/dcrawford/SSM_WGS/trimmed_seq/` (Trim Galore output)
  against the current `refs/Fhet_MT.fasta`. Initial submission had a
  Trimmomatic-vs-TrimGalore naming mismatch (`_p` vs `_val_1`); user
  patched the script in-session to use `${sample}_1_val_1.fq.gz` /
  `${sample}_2_val_2.fq.gz`. BAMs land at
  `/projectnb/dcrawford/MT_Genomics2/AD_again/MT_bam_sam/`. Running as
  of this session; ~144 samtools-sort tmp files visible, no fatal
  errors in the .err log.

- **Two new BSUB scripts written**, both pointed at the AD_again BAMs,
  intended to recreate the input that `additional_info/ANN_Claude_2.py`
  expects so we can compare against the prior ~243-variant haplotype
  output:

  - `additional_info/BSUB_persample_call_merge_again.sh` — default
    filters (`-Q 13 -q 0 -d 10000`), per-sample call (`-mv --ploidy 1
    -a AD,DP`) → `norm -m -any` → `bcftools merge -m none` → reheader
    → stats. Outputs under `AD_again/vcf/persample/` and
    `AD_again/vcf/Fhet_mt_again_persample_merged.vcf.gz`.
  - `additional_info/BSUB_persample_call_merge_again_Q30q30.sh` —
    strict filters (`-Q 30 -q 30 -d 100000`), otherwise identical
    recipe. Outputs under `AD_again/vcf/persample_Q30q30/` and
    `AD_again/vcf/Fhet_mt_again_Q30q30_persample_merged.vcf.gz`.

  Both are single sequential jobs (not arrays) — 144 BAMs × ~3–5 min ≈
  12 h; 36 h wall buffer. Idempotent (skips per-sample VCFs already
  present).

- **Goal of the AD_again branch.** Recreate the canonical input that
  `additional_info/ANN_Claude_2.py` consumes (SnpEff-annotated,
  CDS-restricted, normalized VCF with FORMAT/AD and FORMAT/DP). Once
  that input is regenerated, re-run the haplotype caller and compare
  the resulting calls to the prior ~243-variant haplotype output. The
  critical question is **why do the two haplotype-call sets differ**
  — recipe? reference? trim? — and the AD_again × two-filter-variant
  matrix is designed to isolate the recipe variable.

- **Hypothesis for next session.** Without recovering the Jul-2025
  reference, both new persample → merge runs likely converge near the
  current 227-SNP number, because the variable that changed is the
  reference, not the BAMs or the parameters. The Q30q30 run may dip
  slightly below 227 due to stricter MAPQ filtering. The recipe
  comparison still has value: confirms whether 227 is reproducible
  across BAM regenerations and filter variants, and whether
  `ANN_Claude_2.py`'s haplotype output is sensitive to which
  per-sample → merge variant feeds it.

- **Open parallel path.** Pull canonical NCBI `NC_012312.1`,
  md5-compare against `3f460763fdfcf105c3e3144ca6784f55`. If
  different, the Mac-backup reference is a modified version and a
  re-align against NCBI canonical may recover the ~1133-SNP signal.

- **Next-session pickup (in order):**
  1. On Triton 2: `bjobs` until the AD_again realignment array
     completes; confirm 144 BAMs in `AD_again/MT_bam_sam/`.
  2. `bsub <` each of the two new persample call+merge scripts. They
     can run in parallel — outputs land in separate subdirs.
  3. When both finish, dump per-sample distributions and merged stats
     for each, compare against the 227-SNP baseline.
  4. rsync the two merged VCFs to Mac; SnpEff (stage 06) + CDS
     restrict (stage 07) on each.
  5. Run `additional_info/ANN_Claude_2.py` on each annotated VCF;
     compare haplotype counts to the prior ~243-variant output.
  6. If both new merges land near 227 → pull NCBI canonical
     `NC_012312.1` and re-align 2–3 samples against it as a control
     for the reference hypothesis.

## 2026-05-10 (session 7) — Discrepancy resolved: gap is architectural (joint -mv vs per-sample → merge); stage 05d/05e written; methods write-up landed

> **The 152-vs-1133 gap is solved at the diagnostic level.** v2 + v3 came
> back essentially identical to v1 (152, 152, 145 SNPs respectively),
> ruling out -Q, -q, and ploidy as the dominant variables. The historical
> recipe was found in `archive/Notes_dlcs/Inital_call_wo_AD.txt` and the
> historical stats file in `stats_old/merged_stats.txt` (1140 records,
> 1133 SNPs, ts/tv 7.92 across 144 samples). The actual cause is the
> deliberate switch from per-sample call → bcftools merge to a single
> joint `bcftools call -mv` pipe, made on collaborator (MFO) advice to
> get per-sample AD/DP for downstream haplotype calling. Joint -mv's
> default allele-frequency prior suppresses sites where the reference is
> the rare allele, which on this mtDNA panel against the divergent
> NC_012312.1 reference is most variable positions. Mechanism written up
> in `docs/02_calling_architecture.md`.

- **v2 + v3 results (run on Triton 2 by user, stats rsynced to Mac
  `vcf/`):**

      metric                       fullAD   v2_Q13_q20_p1_fullAD   v3_Q13_q00_p1_variantsAD
      samples                         143                    143                        143
      records                         153                    153                        146
      SNPs                            152                    152                        145
      multiallelic SNP sites          145                    152                        145
      ts/tv (1st alt only)           5.61                   6.60                       7.06
      singleton SNPs (AC=1)            51                     40                         39

  Per the decision tree in session 6: both v2 and v3 still ≈ 150 → gap is
  upstream of stage 05's parameter choices.

- **Recovered historical baseline.** User pointed to
  `stats_old/merged_stats.txt` (Triton 2,
  `/projectnb/dcrawford/SSM_Mito/merged_files/merged_144.vcf.gz`,
  Jul 26 2025): **1140 records, 1133 SNPs, 144 samples, ts/tv 7.92** (all
  alts) / **9.12** (1st alt). Three companion files in the same folder
  document the filtering chain: `merged_144_3.stats.txt` (53,512 SNPs,
  ts/tv 0.52 — raw candidate set), `merged_144_4.stats.txt` (16,517 SNPs,
  ts/tv 0.50 — every-position output), `merged_144_2A.stats.txt` (empty,
  failed filter step). The clean 1133 came out of a heavy filter chain
  starting from the raw soup.

- **Recovered historical recipe.** `archive/Notes_dlcs/Inital_call_wo_AD.txt`
  documents the original per-sample-call → merge pipeline that produced
  `merged_144.vcf.gz`. `archive/Previous_jobs/BSUB_1_MT_SNPcalls.sh` is the
  historical per-sample BSUB. Key features: `--ploidy 1`, `bcftools mpileup`
  defaults (no `-Q`/`-q` overrides), no `-A` in `bcftools call`,
  `bcftools norm -m -any` per-sample, then `bcftools merge -m none` across
  all 144. The recipe deliberately did NOT request `-a AD,DP` from mpileup;
  that omission is what motivated the (architecturally fragile) switch to
  joint mode in the first place.

- **AF spectrum confirms which number is right.** Of the 1133 historical
  SNPs, 1004 (89%) sit at AF ≥ 0.95 — i.e., the reference is the rare
  allele at most variable positions. Singleton counts are similar across
  regimes (55 historical vs 51/40/39 in v1/v2/v3): we catch rare variants
  fine. The gap is at the high-AF end, where joint `-mv` suppresses
  near-fixed-alt sites because they violate its default allele-frequency
  prior. Historical ts/tv = 7.92 across all alts is squarely in the
  expected vertebrate-mtDNA range, so the 1133 are real biology, not
  reference-error noise.

- **`docs/02_calling_architecture.md` — new.** Full mechanism write-up:
  observation, what changed (architecture, normalization, `-A`), why
  joint `-mv` suppresses high-AF sites, AF-spectrum evidence, and a
  practical 7-point checklist for mtDNA SNP calling with bcftools. This
  is the methods-write-up deliverable from CLAUDE.md's "Active primary
  task" section.

- **`jobs/05d_persample_call.sh` — new.** LSF array `[1-143]%24`
  reproducing the historical per-sample recipe, with `-a AD,DP` added so
  per-sample allele depths carry through merge (gives downstream stage 08
  what it needs without forcing the joint architecture). Outputs land
  under `vcf/persample/{sample}_norm.vcf.gz`. Resources: 4 cores / 8 GB /
  2 h per task — matches the historical script's profile.

- **`jobs/05e_merge_persample.sh` — new.** Single-job merge follow-up.
  Pre-flights every per-sample VCF, `bcftools merge -m none`, reheader to
  strip `_0`, stats, manifest. Output: `vcf/Fhet_mt_persample_merged.vcf.gz`
  with target ~1130 SNPs (slightly below the historical 1133 since we
  use 143 BAMs vs the historical 144).

- **`docs/01_pipeline.md` updated.** Pipeline diagram now shows 05d → 05e
  as the canonical path. Stage-05 section rewritten: 05d/05e are the
  canonical scripts, 05/05b/05c are kept as the methods-comparison set
  only (joint variants). Added rsync + symlink instructions in the
  stage-06 prep block so existing 06/07 scripts (which read
  `Fhet_mt_variantsAD.vcf.gz`) keep working without script edits.

- **What this means for v1/v2/v3.** They are not bugs and they are not
  "wrong" — they are the joint-`-mv`-default behavior on a panel with a
  divergent reference. They remain in `jobs/` as the experimental
  artifacts that anchor the methods write-up. Do NOT delete.

- **Submit (Triton 2):**
  ```
  cd /projectnb/dcrawford/MT_Genomics2
  git pull
  bsub < jobs/05d_persample_call.sh
  # wait for the array to finish (bjobs / tail logs/05d_persample_*_*.out)
  bsub < jobs/05e_merge_persample.sh
  ```

- **Next-session pickup (in order):**
  1. On Triton 2: `bjobs` until 05d's 143 array tasks all finish; submit 05e.
  2. Verify: 05e's manifest summary should show ~1100–1133 SNPs and ts/tv
     ≈ 7–9 across all alts. If it does, the discrepancy is fully resolved.
  3. `bash scripts/compare_stage05_runs.sh vcf/Fhet_mt_*_stats.txt` —
     full side-by-side table for the methods write-up.
  4. rsync `Fhet_mt_persample_merged.vcf.gz` (+ `.csi`) to Mac `vcf/`,
     symlink as `Fhet_mt_variantsAD.vcf.gz`, then run stages 06 → 07 → 08
     unchanged.
  5. Once stage 07 produces `Fhet_MT_CDS.snps.split.vcf.gz`, that file
     becomes the frozen canonical output (per CLAUDE.md rule 1).

## 2026-05-09 (session 6) — Stage 05 strict run produced 152 SNPs vs ~950 expected; diagnostic v2 + v3 re-runs prepared; understanding the discrepancy elevated to a primary project task

> **Primary task elevation:** Understanding *why* the same pipeline produces
> 152 vs ~950 SNPs, and writing it up so other investigators can avoid the
> trap, is now a primary project deliverable — not an incidental bug fix.
> See the new "Active primary task" section in CLAUDE.md. The gap is
> qualitative: at 152 SNPs all downstream haplotype/population inferences
> would change; at ~950 only minor revisions are needed. The diagnostic
> runs below are the experimental basis for the methods write-up that owes
> the broader community.


- **Stage 05 v1 (strict) finished cleanly but undershot by ~6×.** Job 7672 ran
  to completion: `=== DONE ===` in stdout, max memory 371 MB of 16 GB
  requested, runtime 33,287 s (9h 15m) of the 72h wall budget, no errors,
  no warnings other than the diploid-default note from `bcftools call`.
  `bcftools stats vcf/Fhet_mt_fullAD.vcf.gz` reports 153 records / 152 SNPs /
  1 indel across 143 samples.
- **Truncation hypothesis ruled out.** Variant POS spans 5 → 16,500 across a
  16,526 bp reference (`refs/Fhet_MT.fasta.fai`). The whole mitogenome was
  scanned; the 152 SNPs is the real call set under the v1 parameters
  (`-Q 30 -q 30 --ploidy 2 (default)`), not a partial output.
- **Diagnosis from the stats file:**
  - 145/152 SNP sites are multiallelic (96%). Per-site Ts/Tv = 0.51 across
    all alts but Ts/Tv = 5.61 considering only the 1st (most-common) alt at
    each site. Pattern: real biological variation at the dominant alt,
    sequencing-error noise at the secondary alts that `-A` retains in the
    record (327 alts have AF=0, i.e. zero carriers).
  - Per-site DP > 500 at 100% of called sites — coverage is not the limit.
  - 51 singleton SNPs is low for 143 unrelated *F. heteroclitus* across the
    full mitogenome; rare-allele detection appears compressed.
- **Why the gap is most likely calling parameters, not data:**
  - User confirms historical canonical run had ~950 SNP sites with QUAL ≥ 30
    at 95%+ of sites — i.e., per-site call confidence was not the limit
    historically; the difference must therefore be in which sites cross the
    variant-calling threshold at all.
  - `-Q 30 -q 30` filters reads BEFORE the per-site likelihood is computed,
    so they shape which sites get emitted. `--ploidy 2` (default) on a
    haploid genome can fail the het-likelihood gate at heteroplasmic or
    low-AF carrier sites.
- **Diagnostic re-runs prepared (all results preserved side-by-side):**
  - **`scripts/run_stage05_core.sh`** — new parametrized core. Takes
    `RUN_TAG MIN_BQ MIN_MQ PLOIDY` as positional args; writes outputs as
    `Fhet_mt_${RUN_TAG}_*` and a self-documenting `_run_manifest.txt`
    capturing parameters, host, job ID, bcftools version, BAM count, and
    a post-run summary (n_records, n_SNPs, ts/tv). Manifests give other
    investigators the exact recipe behind each result.
  - **`jobs/05b_v2_Q13_q20_p1.sh`** — BSUB wrapper, runs the core with
    `-Q 13 -q 20 --ploidy 1` (relaxed read filters + haploid). RUN_TAG:
    `v2_Q13_q20_p1`.
  - **`jobs/05c_v3_Q13_q00_p1.sh`** — BSUB wrapper, runs the core with
    `-Q 13 -q 0 --ploidy 1` (no MAPQ filter; matches the archived
    per-sample caller's defaults for an apples-to-apples comparison).
    RUN_TAG: `v3_Q13_q00_p1`.
  - Both wrappers leave the v1 strict outputs untouched; LSF resources match
    the v1 run since v1's profile shows headroom (16 GB / 8 cores / 96h
    wall clock — bumped from 72h for safety).
  - **`scripts/compare_stage05_runs.sh`** — new helper. Takes any number of
    `*_stats.txt` files and prints a side-by-side fixed-width comparison of
    sample count, record/SNP/indel/multiallelic counts, both Ts/Tv flavors,
    and singletons. Run after both jobs finish.
- **Misleading comments fixed in the new core.** The original script said
  `-A in bcftools call means: keep all positions (not just variants)`. That
  is incorrect — `-A` keeps all alternate alleles at variant sites, but the
  output is still variant-only because of `-mv`. The new core has the right
  comment and notes that `_fullAD` is a historical filename (downstream
  scripts still use it).
- **Submission (Triton 2):**
  ```
  cd /projectnb/dcrawford/MT_Genomics2
  git pull
  bsub < jobs/05b_v2_Q13_q20_p1.sh
  bsub < jobs/05c_v3_Q13_q00_p1.sh
  ```
- **Decision tree once both jobs finish:**
  - `bash scripts/compare_stage05_runs.sh vcf/Fhet_mt_fullAD_stats.txt vcf/Fhet_mt_v2_*_fullAD_stats.txt vcf/Fhet_mt_v3_*_fullAD_stats.txt`
  - If v2 ≈ 950: ploidy + per-base quality filter were the issue; adopt
    v2 as canonical, document the bug, proceed to 06/07/08.
  - If v2 still low but v3 ≈ 950: gap localizes to MAPQ filtering. Adopt
    v3 (or v2 with `-q 10`) as canonical.
  - If both v2 and v3 still ≈ 150: gap is upstream of stage 05 — investigate
    BAM provenance (whether duplicates were marked in the historical pipeline,
    whether trimming defaults differ, etc.).
- **Why this matters for the science write-up:** the difference between
  152 and ~950 SNPs is not cosmetic — at 152, downstream haplotype calls and
  population-structure inferences would change qualitatively. Whichever
  number proves correct, the run manifests + this CHANGELOG entry document
  exactly what was tried and why, so other investigators using
  `bcftools mpileup | bcftools call` on mtDNA can avoid the
  default-diploid + strict-MAPQ trap.

## 2026-05-08 (session 5) — Mac-side scripts patched; stage-05 verification script added

- **Stage 05 status**: still running on Triton 2 as of this session (submitted 11:32 AM). Ploidy warning ("assuming all sites are diploid") is benign — downstream haplotype caller uses `AD_alt/DP > 0.7` and ignores GT entirely. Long runtime explained by 143 × 15–19 GB BAMs (high I/O; unmapped reads retained in BAMs from stage 04).
- **`scripts/06_snpeff_mac.sh` patched** — previous version had three wrong values:
  - `SNPEFF_DIR` was `~/SnpEff/` (uppercase, non-existent) → corrected to `~/snpEff/`
  - `SNPEFF_JAR` was `~/SnpEff/snpEff.jar` → corrected to `~/micromamba/envs/SNP_env/share/snpeff-5.2-1/snpEff.jar` (per `My_previous_SNPeff.txt`)
  - `DB_NAME` was `NC_012312.1` (NCBI pre-built, doesn't exist) → corrected to `Fhet_MT` (custom database built from `genes.gff` + `sequences.fa` under `~/snpEff/data/Fhet_MT/`)
  - Added `cd "$SNPEFF_DIR"` before java call so snpEff resolves relative `data/` paths correctly
  - Removed the chromosome rename block (not needed: VCF CHROM = `NC_012312.1` matches Fhet_MT database sequences)
  - Improved pre-flight error messages with remediation hints
- **`scripts/07_cds_snps_norm_mac.sh` patched** — added pre-flight check that creates `Fhet_MT.fasta.fai` via `samtools faidx` if absent (required by `bcftools norm -f REF`)
- **`scripts/verify_stage05.sh` — new script** — run on Triton 2 once stage-05 log shows `=== DONE ===`. Checks: file exists, CSI index present, bcftools can read without error, 143 samples, chromosome `NC_012312.1`, FORMAT/AD + FORMAT/DP present, SNP count in plausible range (50–5000). Prints the rsync commands on all-pass.
- **Next-session pickup (in order):**
  1. On Triton 2: `tail logs/05_mpileup_AD_*.out` — wait for `=== DONE ===`
  2. On Triton 2: `bash scripts/verify_stage05.sh` — must show 0 failures
  3. rsync VCF + CSI to Mac `vcf/` (commands printed by verify script)
  4. On Mac: `bash scripts/06_snpeff_mac.sh`
  5. On Mac: `bash scripts/07_cds_snps_norm_mac.sh` → canonical `vcf/Fhet_MT_CDS.snps.split.vcf.gz`
  6. On Mac: `pip install cyvcf2 pandas --break-system-packages` (if not already installed)
  7. On Mac: `python scripts/08_call_haplotypes.py` → haplotype matrix + calls

## 2026-05-08 (session 4) — conda activation fix confirmed and finalised

- **Root cause confirmed:** `conda info` on login node revealed:
  - base: `anaconda3-2023.09` at `/sw/summit/software/linux-power9le/anaconda3-2023.09-0-...` (read only)
  - env: `/projectnb/triton/home/dcrawford/.conda/envs/mito_genomics` (user-level)
  - `module load anaconda3` (unversioned) was likely resolving ambiguously or not registering shell hooks on compute nodes.
- **`jobs/config.sh`** updated: added `CONDA_MODULE=anaconda3/2023.09-0-none-none-oawyzwj`; config.sh now runs `module load "$CONDA_MODULE"` + `eval "$(conda shell.bash hook)"` directly, so every BSUB script that sources config.sh gets correct conda activation automatically.
- **All 5 BSUB scripts** (`01`–`05`): removed `module load anaconda3` (now in config.sh); reverted `source activate` back to `conda activate "$CONDA_ENV"` (correct for modern conda once hooks are registered).
- **Confirmed working:** `bash -c "source jobs/config.sh && conda activate mito_genomics && bcftools --version"` prints `bcftools 1.6 / htslib 1.6` on Triton 2 login node. Stage 05 submitted successfully.
- **Next-session pickup (in order):**
  1. Check stage 05 complete: `bjobs` / `tail logs/05_mpileup_AD_*.out` — look for `=== DONE ===`.
  2. `bcftools stats vcf/Fhet_mt_variantsAD.vcf.gz | grep "^SN"` — confirm 143 samples, reasonable SNP count.
  3. `rsync` `Fhet_mt_variantsAD.vcf.gz` + `.csi` from Triton 2 → Mac `vcf/`.
  4. `bash scripts/06_snpeff_mac.sh` — SnpEff annotation on Mac.
  5. `bash scripts/07_cds_snps_norm_mac.sh` → canonical `vcf/Fhet_MT_CDS.snps.split.vcf.gz`.
  6. `pip install cyvcf2 pandas` (if not already on Mac), then `python scripts/08_call_haplotypes.py`.

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

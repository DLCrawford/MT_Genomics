# MT_Genomics

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23418103.svg)](https://doi.org/10.5281/zenodo.23418103)

Reproducible mitochondrial variant, haplotype, diversity, and heteroplasmy pipeline for
*Fundulus heteroclitus*, supporting the manuscript **"Exceptional Mitochondrial Diversity Under
Strong Purifying Selection in Admixed Populations, and How Reference Choice Can Mask It"**
(Sierra-Martinez, Oleksiak & Crawford).

Heavy compute (read QC, alignment, variant calling) runs on an HPC cluster (LSF/`bsub`);
downstream analysis runs locally. This repository holds the code, documentation, and the
mitochondrial reference. Sequence reads and the variant call set are deposited separately
(see [Data availability](#data-availability)). The same code is archived with the data on Zenodo
(DOI 10.5281/zenodo.23418103), where the four folders below appear as `00_upstream_hpc/`
(= `jobs/`), `core_pipeline/` (= `scripts/`), `popgen_selection/` and `contamination_numt/`.

The analysis-ready file behind every table and figure is **`141_MT_variants.vcf.gz`** —
141 individuals (#70 and #125 excluded), 927 protein-coding variants, with per-cell genotype,
depth (DP), allele depth (AD), and SnpEff annotation. It is in the data deposit, not this repo.

## Layout

```
MT_Genomics/
├── README.md, LICENSE, CITATION.cff, environment.yml
├── CANONICAL_SCRIPT_LIST.md/.tsv   every script → purpose, tools, output, manuscript Table/Figure
├── jobs/                HPC (LSF/bsub) read processing + joint variant calling
├── scripts/             downstream analysis (annotation, haplotypes, pNs/pS counts, π/θ, heteroplasmy,
│                        cross-species comparison, McDonald–Kreitman support)
├── popgen_selection/    pNs/pS and π/θ by site class with bootstrap CIs; Tables 1, 3, 4; Figs 2, 4
├── contamination_numt/  contamination tests (adjacent well, whole haplotype, nuclear 3rd allele) + NUMT screen
├── resources/           mitochondrial reference FASTA + GFF (NC_012312.1)
├── data/                (git-ignored) place the deposited data files here
└── docs/                pipeline walkthrough, methods, cautionary note
```

## Pipeline overview

**HPC (`jobs/`)** — run in order; paths and the bcftools 1.23.1 PATH injection are set in `config.sh`.

| Step | Script | Purpose |
|---|---|---|
| 1 | `01_fastqc_raw.sh` | raw-read QC (FastQC) |
| 2 | `02_trim_pe.sh` | adapter/quality trimming (Trimmomatic) |
| 3 | `03_fastqc_trimmed.sh` | post-trim QC |
| 4 | `04_bwa_align_mt.sh` | BWA-MEM alignment to `NC_012312.1` |
| 5 | `BSUB_Slim_BAM_mt.sh` | extract mitochondrial reads → slim BAMs |
| 6 | `05_1_mpileup_merge.sh` | joint `mpileup` + `call -mv --ploidy 1` + `norm -m -any` |
| (alt) | `06_snpeff_annotate.sh`, `07_cds_snps_norm.sh` | HPC equivalents of the Mac annotation steps |

**Local** — manuscript element → script(s).

| Manuscript element | Script(s) |
|---|---|
| SnpEff annotation → CDS/SNP restriction | `scripts/06_snpeff_mac.sh`, `scripts/07_cds_snps_norm_mac.sh` |
| 141-sample analysis VCF | `scripts/make_141_vcf.sh` (`make_141_vcf.py`) |
| Table 1 — pNs/pS per gene & complex (95% CI); Fig. 2 | `scripts/10_dnds_per_gene.py` → `popgen_selection/01_pnps_bootstrap.py`, `06_tables_figures_Sept28.py` |
| Table 2 — variants per individual under 3 references; Fig. 3 | `scripts/18_variant_burden_per_individual.py`, `scripts/35_clade_ns_reref.py` |
| Table 3 — π, θ, π/θ overall and within clades; Fig. 4 | `popgen_selection/02`, `03`, `12`, `13` (+ `common.py`); `scripts/19_calc_pi.py`, `20_calc_pi_clade.py` |
| Table 4 — multivariant-site sensitivity | `popgen_selection/12_multiallelic_convention.py` |
| Table 5 — haplotypes (all sites / nonsynonymous) | `scripts/08_call_haplotypes.py`, `scripts/11_haplotypes_nonsyn.py` |
| Table 6 — nonsynonymous AA-sets | `scripts/12_ns_cooccurrence.py` |
| Table 7 — nuclear 3rd allele | `contamination_numt/run_simpleVCF_allfish.sh` → `contam_allfish.py` |
| Table 8 — cross-species π/θ | `scripts/19_calc_pi.py`, `20_dros_pi.py`, `21_human_mt_pi.py`, `22_human_mt_cds_pi.py`, `23_celegans_pi.py`, `26_amr_pi.py` → `29_comparison_table_Lcor.py` (site-length corrected; supersedes `25_comparison_table.py`) |
| Heteroplasmy (incl. private variants); Fig. 5 | `scripts/09_heteroplasmy_report.py`, `13_pileup_cds_AD.sh`, `14_hp_from_pileup.py`, `16_annotate_hp.py`, `17_annotate_hp_codon.py` |
| Contamination & NUMT tests | `contamination_numt/adjacent_well_direct_test.py`, `whole_haplotype_contamination_test.py`, `run_numt_pipeline.sh` |
| Supporting (not in manuscript tables) | `scripts/27_fhet_clusters_pi.py`, `28_mk_test.py`, `30`–`34` (McDonald–Kreitman, PAML) |
| Utility | `scripts/DP_AD_table.py` |

The full per-script index is `CANONICAL_SCRIPT_LIST.md`.

## Reproducing the analysis

1. Create the environment: `conda env create -f environment.yml && conda activate mito_genomics`.
   Variant calling requires **bcftools/htslib 1.23.1 built from source** — the bioconda-pinned 1.6
   silently under-calls high-depth mtDNA pileups (see `docs/cautionary_note.md`).
2. Run `jobs/` steps 1–6 on the HPC to produce the joint call set, or start from the deposited
   `141_MT_variants.vcf.gz` (place it in `vcf/`).
3. Run the scripts for the analysis of interest (table above). `popgen_selection/` scripts import
   `common.py` and read inputs from `../vcf/`.

## Data availability

- **Sequence reads:** NCBI SRA, BioProject [PRJNA1543696](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1543696) — per-individual
  mitochondrial BAMs (143 individuals) aligned to `NC_012312.1`.
- **Variant call set, heteroplasmy pileup, cross-species tables, metadata, and the complete code:**
  Zenodo, DOI [10.5281/zenodo.23418103](https://doi.org/10.5281/zenodo.23418103) — see `data/README.md`.
- **Mitochondrial reference:** `resources/Fhet_MT.fasta` + `Fhet_MT.gff` (NCBI RefSeq `NC_012312.1`,
  MU-UCD_Fhet_4.1).

## Citation

Sierra-Martinez S., Oleksiak M.F., Crawford D.L. *Exceptional Mitochondrial Diversity Under Strong
Purifying Selection in Admixed Populations, and How Reference Choice Can Mask It.* *(journal, year, DOI TBD)*.
Code and data: Zenodo, DOI 10.5281/zenodo.23418103.

## License

Code released under the MIT License (see `LICENSE`).

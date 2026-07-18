# MT_Genomics

Reproducible mitochondrial variant, haplotype, and diversity pipeline for
*Fundulus heteroclitus*, supporting the manuscript **"Mitochondrial Genomics:
Variation in an Admixture Population"** (Sierra-Martinez, Oleksiak & Crawford).

Heavy compute (read QC, alignment, variant calling) runs on an HPC cluster
(LSF/`bsub`); downstream analysis runs locally. The repository holds only the
code, documentation, and the mitochondrial reference needed to reproduce the
paper. Sequence reads and the variant call set are deposited separately
(see [Data availability](#data-availability)).

The single analysis-ready file behind every table and figure is
**`141_MT_variants.vcf.gz`** — 141 individuals, 927 protein-coding variants,
with per-cell genotype, depth (DP), allele depth (AD), and SnpEff annotation.
It is in the data deposit, not this repo.

## Layout

```
MT_Genomics/
├── README.md               this file
├── LICENSE
├── environment.yml         conda environment (tool + package versions)
├── jobs/                   HPC (LSF/bsub) read-processing + variant calling
├── scripts/                downstream analysis (Python + bash)
├── resources/              mitochondrial reference FASTA + GFF (NC_012312.1)
└── docs/                   pipeline walkthrough, methods, cautionary note
```

## Pipeline overview

**HPC (`jobs/`)** — run in order; paths and the bcftools 1.23.1 PATH injection
are set in `config.sh`.

| Step | Script | Purpose |
|---|---|---|
| 1 | `01_fastqc_raw.sh` | raw-read QC (FastQC) |
| 2 | `02_trim_pe.sh` | adapter/quality trimming |
| 3 | `03_fastqc_trimmed.sh` | post-trim QC |
| 4 | `04_bwa_align_mt.sh` | BWA-MEM alignment to `NC_012312.1` |
| 5 | `BSUB_Slim_BAM_mt.sh` | extract mitochondrial reads → slim BAMs |
| 6 | `05_1_mpileup_merge.sh` | joint `mpileup` + `call -mv --ploidy 1` + `norm -m -any` |

**Local (`scripts/`)** — annotation, restriction, and analysis.

| Analysis (manuscript element) | Script(s) |
|---|---|
| SnpEff annotation → CDS/SNP restriction | `06_snpeff_mac.sh`, `07_cds_snps_norm_mac.sh` |
| Build the 141-sample analysis VCF | `make_141_vcf.sh` |
| Haplotype calling (AD-based, 0.7) | `08_call_haplotypes.py` |
| dN, dS, dN/dS per gene (Table 1) | `10_dnds_per_gene.py` |
| Nonsynonymous AA-sets (Table 4) | `11_haplotypes_nonsyn.py`, `12_ns_cooccurrence.py` |
| Variants/individual + 3-reference N/S (Table 2, Fig 2) | `18_variant_burden_per_individual.py`, `35_clade_ns_reref.py` |
| π / θ, overall and by clade (Table 3) | `19_calc_pi.py`, `20_calc_pi_clade.py` |
| Cross-species π / θ (Table 5) | `20_dros_pi.py`, `21_human_mt_pi.py`, `22_human_mt_cds_pi.py`, `29_comparison_table_Lcor.py` |
| Heteroplasmy + well-bleed contamination test | `09_heteroplasmy_report.py`, `13_pileup_cds_AD.sh`, `14_hp_from_pileup.py`, `15_well_bleed_test.py`, `16_annotate_hp.py`, `17_annotate_hp_codon.py` |
| Long-format DP/AD table (utility) | `DP_AD_table.py` |

> **Cross-species note.** The *C. elegans* and yeast π scripts run in their own
> dataset folders and are not included here; see `docs/methods_variant_calling.md`
> and the data deposit. `20_dros_pi.py`, `21/22_human_mt_*` reproduce the
> *Drosophila* and human rows of Table 5.

## Reproducing the analysis

1. Create the environment: `conda env create -f environment.yml && conda activate mito_genomics`.
   Variant calling requires **bcftools/htslib 1.23.1 built from source** — the
   bioconda-pinned 1.6 silently under-calls high-depth mtDNA pileups
   (see `docs/cautionary_note.md`).
2. Run `jobs/` steps 1–6 on the HPC to produce the joint call set, or start from
   the deposited `141_MT_variants.vcf.gz`.
3. Run the `scripts/` for the analysis of interest (table mapping above).

## Data availability

- **Raw sequence reads:** NCBI SRA / BioProject *(accession TBD)*.
- **Analysis-ready variant call set + metadata:** Zenodo *(DOI TBD)* — contains
  `141_MT_variants.vcf.gz` (+ index), a sample-metadata table
  (ID → population → clade → plate/well), and this repository's release archive.
- **Mitochondrial reference:** `resources/Fhet_MT.fasta` + `Fhet_MT.gff`
  (NCBI RefSeq `NC_012312.1`, MU-UCD_Fhet_4.1), included here for convenience.

## Citation

Sierra-Martinez S., Oleksiak M.F., Crawford D.L. *Mitochondrial Genomics:
Variation in an Admixture Population.* *(journal, year, DOI TBD)*.
Code archived at Zenodo *(DOI TBD)*.

## License

Code released under the MIT License (see `LICENSE`).

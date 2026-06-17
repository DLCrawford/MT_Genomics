# MT_Genomics

Reproducible mitochondrial variant + haplotype pipeline for *Fundulus heteroclitus*,
supporting the manuscript *"Mitochondrial Genomics, Variation in an Admixture
Population."* Heavy compute (read processing, alignment, variant calling) runs on
the Triton 2 HPC (LSF/`bsub`); downstream analysis runs locally (macOS).

The canonical analysis set is **`vcf/141_MT_variants.vcf.gz`** (141 samples,
927 coding-region records). Data files are not tracked (see `.gitignore`); this
repo holds the scripts, jobs, and documentation that produce them.

## Layout

```
MT_Genomics/
├── README.md          — this file
├── CLAUDE.md          — pipeline objectives, lessons learned, per-script catalog
├── CHANGELOG.md       — session-by-session log
├── .gitignore         — excludes large data files (and *.docx manuscripts)
├── docs/
│   ├── 00_setup.md                 — HPC + conda env + reference indexing
│   ├── 01_pipeline.md              — step-by-step walkthrough
│   ├── 02_calling_architecture.md  — caller design + bcftools-version finding
│   ├── methods_variant_calling.md  — methods write-up
│   └── mito_protein_coding.bed     — 13 CDS intervals (11,417 bp)
├── jobs/              — HPC (LSF/bsub) steps
│   ├── config.sh                   — paths + PATH-injects bcftools 1.23.1
│   ├── 01_fastqc_raw.sh            — raw-read QC (FastQC)
│   ├── 02_trim_pe.sh               — adapter/quality trimming
│   ├── 03_fastqc_trimmed.sh        — post-trim QC
│   ├── 04_bwa_align_mt.sh          — BWA-MEM → sorted BAM
│   ├── BSUB_Slim_BAM_mt.sh         — extract mt-only (NC_012312.1) reads
│   ├── 05_1_mpileup_merge.sh       — JOINT mpileup | call -mv --ploidy 1 | norm -m -any
│   ├── 06_snpeff_annotate.sh       — SnpEff (HPC equivalent of scripts/06)
│   ├── 07_cds_snps_norm.sh         — CDS restrict (HPC equivalent of scripts/07)
│   └── _archive/                   — superseded callers (methods comparison)
└── scripts/           — Mac-side, run on the canonical VCF
    ├── make_141_vcf.{sh,py}        — build vcf/141_MT_variants.vcf.gz
    ├── 06_snpeff_mac.sh            — SnpEff annotation (SNP_env)
    ├── 07_cds_snps_norm_mac.sh     — CDS restrict + SNPs + split → canonical VCF
    ├── 08_call_haplotypes.py       — haplotype matrix/calls (AD/DP ≥ 0.7)
    ├── 09_heteroplasmy_report.py   — Hp classifier (variant-only table)
    ├── 10_dnds_per_gene.py         — per-gene + Overall dN/dS (NCBI table 2)
    ├── 11_haplotypes_nonsyn.py     — NS-only haplotype matrix
    ├── 12_ns_cooccurrence.py       — perfect-LD NS groups (AA-sets)
    ├── 13_pileup_cds_AD.sh         — per-CDS pileup AD/DP across panel BAMs
    ├── 14_hp_from_pileup.py        — Hp from pileup (exposes private_alt_Hp)
    ├── 15_well_bleed_test.py       — plate-contamination permutation test
    ├── 16_annotate_hp.py           — SnpEff SYN/NS onto Hp events
    ├── 17_annotate_hp_codon.py     — codon-level SYN/NS for unannotated Hp
    ├── 18_variant_burden_per_individual.py — per-sample variant counts; N/S clusters
    ├── 19_calc_pi.py               — π (total/syn/NS), Fhet
    ├── 20_calc_pi_clade.py         — π split by North/South clade
    ├── 20_dros_pi.py               — Drosophila π (comparison)
    ├── 25_comparison_table.py      — cross-species comparison assembly
    ├── 27_fhet_clusters_pi.py      — Fhet cluster π (comparison)
    ├── 28_mk_test.py               — McDonald–Kreitman-like test (α, ω)
    ├── DP_AD_table.py              — long-format per-cell DP/AD table
    └── _archive/                   — superseded helper scripts
```

Cross-species π comparison scripts (`21`–`24`: human, *C. elegans*, yeast) live in
sibling dataset folders alongside this repo.

## Pipeline → manuscript

| Manuscript element | Script / job |
|---|---|
| QC, trim, align, mt extraction | `jobs/01`–`04`, `BSUB_Slim_BAM_mt.sh` |
| Variant calling (joint)        | `jobs/05_1_mpileup_merge.sh` |
| Annotation + CDS restriction   | `scripts/06_snpeff_mac.sh`, `scripts/07_cds_snps_norm_mac.sh` |
| Canonical 141-sample VCF       | `scripts/make_141_vcf.sh` |
| Table 1 (dN/dS)                | `scripts/10_dnds_per_gene.py` |
| Table 2 (per-cluster, π)       | `scripts/18_variant_burden_per_individual.py`, `scripts/20_calc_pi_clade.py` |
| Table 3 (haplotypes)           | `scripts/08_call_haplotypes.py` |
| Table 4 (AA-sets)              | `scripts/12_ns_cooccurrence.py` |
| Table 5 (cross-species)        | `scripts/19_calc_pi.py`, `scripts/20_dros_pi.py`, `scripts/25_comparison_table.py` |
| Table 6 (MK test)              | `scripts/28_mk_test.py` |
| Heteroplasmy                   | `scripts/09`, `13`, `14`, `15`, `16`, `17` |

## Reproducibility notes

- **Tool versions.** All calling used **htslib/bcftools 1.23.1 built from source**.
  bcftools 1.6 (the version bioconda pins on linux-ppc64le) silently under-calls
  high-depth mtDNA pileups (6 vs 289 SNPs on identical input); `jobs/config.sh`
  PATH-injects 1.23.1 ahead of it. See `docs/02_calling_architecture.md`.
- **What actually ran.** Every bcftools step stamps an ordered
  `##bcftools_<cmd>Command=` / `##bcftools_<cmd>Version=` line into the VCF header.
  Reading those lines is the authoritative way to reconstruct exactly what was
  executed (the 927-variant set was produced by joint calling — single `call`,
  no `merge`).

## Reference

- MT chromosome: `NC_012312.1` (annotation `MU-UCD_Fhet_4.1`)
- Codon table: Vertebrate Mitochondrial (NCBI translation table 2)

# Canonical script & job list

60 scripts, one row each, in run order within four stages. Provenance (tool/version, output, manuscript Table/Figure) is from the project pipeline table (Supplemental Table S1 = `00_upstream_hpc/Previous_scripts_pipeline_table.txt`) where available, otherwise summarised from each script's header. Data inputs (FASTQ, BAMs, reference genomes, VCF/pileup tables) are archived separately.

> **GitHub layout:** stage folders map to this repository as `00_upstream_hpc/` → `jobs/`, `core_pipeline/` → `scripts/`,
> `popgen_selection/` → `popgen_selection/`, `contamination_numt/` → `contamination_numt/` (the Zenodo zip keeps the stage names).


## 0 · Upstream (HPC, Triton 2 / LSF)

*Read QC, trimming, alignment, slim BAMs, variant calling + normalization*

| Script | Platform | Purpose | Tools | Output / Table |
|---|---|---|---|---|
| `01_fastqc_raw.sh` | HPC | Quality control of raw paired FASTQ reads | FastQC | QC reports |
| `02_trim_pe.sh` | HPC | Adapter and low-quality base trimming (paired-end) | Trimmomatic v0.39 | Trimmed paired reads |
| `03_fastqc_trimmed.sh` | HPC | Quality control of trimmed reads | FastQC | QC reports |
| `04_bwa_align_mt.sh` | HPC | Map reads to the F. heteroclitus mt reference (NC_012312.1); SAM to sorted BAM | BWA-MEM v0.7.17; samtools v1.3.1 | {sample}_MT.bam |
| `05_1_mpileup_merge.sh` | HPC | Joint variant calling + normalization in one pipeline: mpileup -b slim_bamlist -a AD,DP -d 100000 | call -mv --ploidy 1 | norm -m -any -f REF. Includes a baseline-regression gate (lost sites = 0). | htslib / bcftools 1.23.1 (built from source) | 05_1_Fhet_mt_persample_merged.vcf.gz |
| `06_snpeff_annotate.sh` | HPC | HPC equivalents of the annotation and CDS-restriction steps; downstream analyses used the Mac scripts below | SnpEff v5.2; bcftools 1.23.1 | (see Mac steps 06/07) |
| `07_cds_snps_norm.sh` | HPC | HPC equivalents of the annotation and CDS-restriction steps; downstream analyses used the Mac scripts below | SnpEff v5.2; bcftools 1.23.1 | (see Mac steps 06/07) |
| `BSUB_Slim_BAM_mt.sh` | HPC | Extract mitochondrial (NC_012312.1) alignments into slim mt-only BAMs | samtools 1.23.1 | MT_only_bams/{sample}_MT_only.bam |
| `config.sh` | HPC | Environment / PATH; injects bcftools 1.23.1 ahead of the conda env’s pinned 1.6 (see Note 1) | — | Sourced by all HPC jobs |

## 1 · Core pipeline (local, Mac)

*Canonical VCF, haplotypes, heteroplasmy, annotation, diversity, dN/dS & MK*

| Script | Platform | Purpose | Tools | Output / Table |
|---|---|---|---|---|
| `06_snpeff_mac.sh` | Mac | SnpEff annotation against custom NC_012312.1 database; vertebrate-mitochondrial codon table (NCBI table 2) | SnpEff v5.2; Java | *_ann.vcf.gz |
| `07_cds_snps_norm_mac.sh` | Mac | Restrict to the 13 protein-coding genes (11,417 bp), SNPs only, sample rename | bcftools 1.23.1 | vcf/Fhet_MT_CDS.snps.split.vcf.gz (canonical) |
| `08_call_haplotypes.py` | Mac | Haplotype matrix and per-sample calls from allele depth (ALT if AD/DP > 0.7) | Python (cyvcf2/pysam) | haplotype_matrix.csv; haplotype_calls.csv — Table 5 (all-site haplotypes; Hd) |
| `09_heteroplasmy_report.py` | Mac | Heteroplasmy classifier on variant-only per-cell DP/AD (0.10 ≤ AF < 0.70) | Python | heteroplasmy_*.tsv — Heteroplasmy |
| `10_dnds_per_gene.py` | Mac | Per-gene and Overall dN, dS, dN/dS (NCBI table 2; Jukes–Cantor correction) | Python | vcf/dnds_per_gene.tsv — Table 1 (bootstrap CIs from popgen_selection/01_pnps_bootstrap.py) |
| `11_haplotypes_nonsyn.py` | Mac | Nonsynonymous-only haplotype matrix and calls | Python | haplotype_matrix_nonsyn.csv; haplotype_calls_nonsyn.csv — Table 5 (Ns haplotypes) |
| `12_ns_cooccurrence.py` | Mac | Co-occurring nonsynonymous sites (AA-sets; perfect LD groups) | Python | vcf/ns_cooccurrence_groups.tsv — Table 6 |
| `13_pileup_cds_AD.sh` | Mac | Per-CDS-position pileup of AD/DP across all 141 BAMs (--no-BAQ, -d 100000) | bcftools 1.23.1 | vcf/pileup_cds_141.vcf.gz |
| `14_hp_from_pileup.py` | Mac | Pileup-based heteroplasmy incl. private-ALT (DP ≥ 20, 0.10 ≤ AF < 0.70, AD ≥ 4) | Python | heteroplasmy_pileup_*.tsv — Heteroplasmy |
| `16_annotate_hp.py` | Mac | Join SnpEff SYN/NS effect class onto heteroplasmy event tables | Python | *_annot.tsv |
| `17_annotate_hp_codon.py` | Mac | Codon-level SYN/NS for heteroplasmy SnpEff left unannotated (FASTA + GFF, table 2) | Python | *_annot_codon.tsv |
| `18_variant_burden_per_individual.py` | Mac | Per-individual variant counts (Total/SYN/NS); defines the bimodal North (<50) and South (>200) clusters | Python | per_individual_burden_*.tsv — Table 2; Fig. 3 |
| `19_calc_pi.py` | Mac | Nucleotide diversity π (total / syn / NS) for the Fhet coding region | Python (pysam) | vcf/pi_results.tsv — Table 8 (F. heteroclitus row) |
| `20_calc_pi_clade.py` | Mac | π (total / syn / NS) split by North vs South clade (added for this revision) | Python (pysam) | vcf/pi_by_clade_persite.tsv — Table 3 |
| `20_dros_pi.py` | Mac | Drosophila (DGRP) π and θ for the cross-species comparison (segregating sites only) | Python | vcf/dros_pi_results.tsv — Table 8 |
| `21_human_mt_pi.py` | Mac | Whole-genome nucleotide diversity (π) for complete human mtDNA genomes (cross-species comparison) | Python | Table 8 (human rows) |
| `22_human_mt_cds_pi.py` | Mac | CDS nucleotide diversity (π_syn, π_ns) for human mtDNA (cross-species comparison) | Python | Table 8 (human rows) |
| `23_celegans_pi.py` | Mac | C. elegans (CaeNDR, n = 540) π and θ for the cross-species comparison | Python (pysam) | celegans_pi_per_site.tsv — Table 8 |
| `26_amr_pi.py` | Mac | Human AMR (gnomAD v3.1, n = 5,718) π and θ for the cross-species comparison | Python | amr_pi_per_site.tsv — Table 8 |
| `25_comparison_table.py` | Mac | Assemble cross-species comparison (N, L, Var, Syn, NS, θ, π and ratios) | Python | vcf/comparison_table.tsv — superseded by 29_comparison_table_Lcor.py (Table 8) |
| `27_fhet_clusters_pi.py` | Mac | π for the Fhet North/South clusters (comparison support) | Python | Table 3 / Table 8 support |
| `28_mk_test.py` | Mac | McDonald–Kreitman-like test (α, ω) overall and per clade | Python | supporting analysis (MK; not in manuscript tables) |
| `29_comparison_table_Lcor.py` | Mac | Site-length-corrected cross-species comparison table; segregating sites = sites with π > 0 | Python | Table 8 (site-length-corrected) |
| `30_outgroup_dS_gate.py` | Mac | Synonymous-saturation gate for the outgroup-based MK test (decides whether count-based MK is valid) | Python | MK gate |
| `31_mk_outgroup.py` | Mac | Outgroup McDonald–Kreitman test for the Fhet mt CDS | Python | supporting analysis (MK; not in manuscript tables) |
| `32_make_paml_input.py` | Mac | Build PAML-ready codon alignment + control files for divergence dN/dS | Python; PAML | paml/ inputs |
| `33_mk_from_paml.py` | Mac | MK α / ω_a from PAML divergence combined with Fhet polymorphism | Python; PAML | supporting analysis (MK; not in manuscript tables) |
| `34_mk_test_catenatus.py` | Mac | Full MK contingency test vs the outgroup F. catenatus (grandis) | Python | supporting analysis (MK; not in manuscript tables) |
| `35_clade_ns_reref.py` | Mac | Reference-placement demonstration for the clade NS/S question (reference-invariance check) | Python | supplementary |
| `DP_AD_table.py` | Mac | Long-format per-cell DP/AD table (alternative input) | Python | vcf/mtDNA_long_AD_table.tsv |
| `make_141_vcf.py` | Mac | Build canonical 141-sample VCF: exclude individuals 70 & 125; drop AC=0 monomorphic sites | bcftools 1.23.1 | vcf/141_MT_variants.vcf.gz (927 records) |
| `make_141_vcf.sh` | Mac | Build canonical 141-sample VCF: exclude individuals 70 & 125; drop AC=0 monomorphic sites | bcftools 1.23.1 | vcf/141_MT_variants.vcf.gz (927 records) |

## 2 · Selection statistics

*pN/pS and π/θ with bootstrap CIs; manuscript tables & figures*

| Script | Platform | Purpose | Tools | Output / Table |
|---|---|---|---|---|
| `01_pnps_bootstrap.py` | Mac | Within-species pN/pS per gene and per OXPHOS complex, with 95% bootstrap CIs | Python | pN/pS table (CIs) — Table 1; Fig. 2 |
| `02_pi_theta_wholepanel.py` | Mac | π/θ by site class (synonymous vs nonsynonymous), whole panel (all 141) | Python | π/θ table — Table 3; Fig. 4 |
| `03_pi_theta_withinclade.py` | Mac | π/θ by site class within the Northern (n=77) and Southern (n=63) clades separately | Python | π/θ table — Table 3; Fig. 4 |
| `04_make_workbook_figure.py` | Mac | Assemble the pN/pS and π/θ result tables into a workbook + the summary figure | Python (openpyxl, matplotlib) | figures/workbook |
| `05_per_individual_verification_table.py` | Mac | Per-individual synonymous/nonsynonymous variant verification table | Python | verification table |
| `06_tables_figures_Sept28.py` | Mac | Build the pN/pS-by-complex table and the π/θ within/between-clade figures | Python (matplotlib) | figures/tables — Fig. 2; Fig. 4 |
| `08_pitheta_dualaxis.py` | Mac | π/θ figure variant — single dual-axis panel (provenance; superseded) | Python (matplotlib) | figure variant |
| `09_pitheta_merged.py` | Mac | π/θ figure variant — all-in-one merged panel (provenance; superseded) | Python (matplotlib) | figure variant |
| `10_pitheta_greycolumn.py` | Mac | π/θ figure variant — grey NS/S ratio columns (provenance) | Python (matplotlib) | figure variant |
| `11_pitheta_master_table.py` | Mac | π/θ master table (per-allele convention; superseded by 13) | Python (openpyxl) | π/θ master table |
| `12_multiallelic_convention.py` | Mac | Multiallelic='one site' sensitivity analysis across S, θ_W and π; rebuilt Table 3 | Python | Table 3; Table 4 (multivariant sensitivity) |
| `13_master_table_collapse.py` | Mac | Final π/θ master table under the multiallelic='one site' convention | Python | π/θ master table (final) — Table 3 |
| `common.py` | Mac | Shared constants, estimators and bootstrap for the selection-statistics scripts (required by 01–13) | Python | (library) |

## 3 · Contamination & NUMT (revision)

*Four tests that the heteroplasmy is genuine, not contamination: (1) direct adjacent-well sharing (`adjacent_well_direct_test.py`), (2) private variants defining unique Hp haplotypes, (3) whole-haplotype sharing with other individuals' major haplotypes (`whole_haplotype_contamination_test.py`), and (4) nuclear-genome 3rd-variant frequency (`contam_allfish.py` + `run_simpleVCF_allfish.sh`); plus the NUMT screen (`run_numt_pipeline.sh`)*

| Script | Platform | Purpose | Tools | Output / Table |
|---|---|---|---|---|
| `adjacent_well_direct_test.py` | local | Direct adjacent-well contamination test (contamination test #1): for each Hp fish, do individuals in physically adjacent plate wells (Chebyshev ≤1, any neighbour, Hp or not) carry its minor variants? Reports proximity enrichment (adjacent vs non-adjacent same-plate sharing), rare/"fingerprint" minor sharing, and adjacency-restricted whole-haplotype match | Python | Hp_adjacent_well_direct_test.tsv |
| `contam_allfish.py` | local | Per-individual nuclear 3rd-variant contamination metric (site depth 8–60; non-major variant 0.10≤AD/DP<0.35), z vs non-Hp baseline, Mann–Whitney U | Python | per_sample_contam_allfish.tsv; verdict.txt — Table 7 |
| `run_nuclear_contam_allfish_v2.sh` | HPC | Variant of the nuclear screen that extracts AD from an existing all-fish nuclear VCF, then runs contam_allfish.py | bash; bcftools | ad.tsv → contam_allfish.py |
| `run_numt_pipeline.sh` | HPC/local | NUMT identification (SwarmGenomics: makeblastdb → blastn mito-vs-nuclear → ≥200 bp → bedtools getfasta) + per-position cross-check of heteroplasmy against NUMT-divergent sites | bash; BLAST+, bedtools, Python | numt.fasta, numt_loci.tsv, numt_divergent_positions.tsv |
| `run_simpleVCF_allfish.sh` | HPC | All-141 nuclear pileup over 10 Mb of chr1 (bcftools mpileup -a AD,DP; MAPQ≥30, BQ≥20); feeds contam_allfish.py | bash; bcftools | ad.tsv → contam_allfish.py |
| `whole_haplotype_contamination_test.py` | local | Whole-haplotype contamination test: does any single other individual's major haplotype carry all of a fish's minor variants with no extra differences | Python | Hp_wholehaplotype_test_*.tsv |

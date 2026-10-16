# Data

Data files are **not** stored in this repository. They are deposited in
public archives and referenced by accession/DOI. To run the `scripts/`,
drop the deposited files into this `data/` folder (it is git-ignored).

## Deposited with this project (Zenodo, DOI 10.5281/zenodo.23418103)

| Item | File(s) | Reproduces | Notes |
|---|---|---|---|
| Analysis-ready variant call set | `141_MT_variants.vcf.gz` (+ `.tbi`) | Tables 1–6; Figs 2–4; *F. heteroclitus* row of Table 8 | 141 samples × 927 coding variants; per-cell GT/DP/AD + SnpEff `ANN` |
| Sample metadata | metadata table | population/clade assignments | ID → population (DC/PA/TR) → clade (N/S) → plate/well |
| Heteroplasmy pileup | `pileup_cds_141.vcf.gz` (+ `.tbi`) | Heteroplasmy section (script `14_hp_from_pileup.py`) | High-depth per-CDS-position AD/DP pileup across all 141 BAMs. **Not** derivable from `141_MT_variants.vcf.gz`; needed to reproduce the pileup-based (private-ALT) heteroplasmy without downloading the BAMs |
| Cross-species π/θ tables | `celegans_pi_per_site.tsv`, `human_mt_pi_per_site.tsv`, `human_mt_cds_pi_per_site.tsv`, `amr_pi_per_site.tsv`, `dros_pi_results.tsv` | Table 8 (via `29_comparison_table_Lcor.py`) | Per-species per-site outputs; let a reader regenerate Table 8 without re-downloading the third-party source data below |

`141_MT_variants.vcf.gz` reproduces the core Fundulus tables and figures;
`pileup_cds_141.vcf.gz` adds the heteroplasmy section; the five `.tsv`
tables add the cross-species comparison (Table 8).

> The Zenodo record also contains the complete analysis code (`MT_Genomics_code_complete.zip`,
> 60 scripts/jobs) and the sample metadata table (`MT_Genomics_metafile.txt`).

## Third-party source data — cited, not re-hosted

The comparison-species π/θ tables above are derived from public datasets.
Obtain these from their original archives (see manuscript References /
Data Availability):

| Species / dataset | Source | Used by |
|---|---|---|
| Raw whole-genome reads (141 *F. heteroclitus*) | NCBI SRA / BioProject *(accession TBD)* | upstream of `141_MT_variants.vcf.gz` |
| Human, sub-Saharan African (N=1,176) | Lankheet et al. 2026 | `21_human_mt_pi.py`, `22_human_mt_cds_pi.py` |
| Human, AMR admixed (N=5,718) | gnomAD v3.1 (chrM) | `26_amr_pi.py` |
| *C. elegans* wild isolates (N=540) | CaeNDR release | `23_celegans_pi.py` |
| *Drosophila* (DGRP, 169 lines) | Mackay et al. 2012 (DGRP) | `20_dros_pi.py` |
| Mitochondrial references | Fhet `NC_012312.1`; human rCRS `NC_012920` | annotation / alignment |

## Reference (in-repo)

`../resources/Fhet_MT.fasta`, `Fhet_MT.gff` — NCBI RefSeq `NC_012312.1`
(MU-UCD_Fhet_4.1).

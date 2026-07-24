# Data

Data files are **not** stored in this repository. They are deposited in
public archives and referenced by accession/DOI.

| Item | Location | Notes |
|---|---|---|
| Raw whole-genome sequence reads | NCBI SRA / BioProject *(accession TBD)* | 141 individuals, PE100 |
| Analysis-ready variant call set | Zenodo *(DOI 10.5281/zenodo.21536984)* | `141_MT_variants.vcf.gz` (+ index): 141 samples × 927 coding variants, per-cell GT/DP/AD + SnpEff `ANN` |
| Sample metadata | Zenodo *(DOI 10.5281/zenodo.21536984)* | individual ID → population (DC/PA/TR) → clade (N/S) → plate/well |
| Mitochondrial reference | `../resources/Fhet_MT.fasta`, `Fhet_MT.gff` | NCBI RefSeq `NC_012312.1` (MU-UCD_Fhet_4.1) |

`141_MT_variants.vcf.gz` is the single file that reproduces every table and
figure in the manuscript. Drop it here (this `data/` folder) to run the
`scripts/`; it is git-ignored so it will not be committed.

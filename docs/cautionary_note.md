# Cautionary note: bcftools version and high-depth mtDNA calling

`bcftools` **1.6** — the version available through the bioconda `linux-ppc64le`
channel on our HPC — silently produced far fewer variant calls on high-depth
mitochondrial pileups, without emitting an error or warning. On the same three
BAM files, with identical parameters:

| bcftools version | SNPs called |
|---|---|
| 1.6    | 6   |
| 1.22   | 284 |
| 1.23.1 | 289 |

All variant calling for this study therefore used **bcftools / samtools /
htslib 1.23.1 built from source**, PATH-injected ahead of any conda-provided
binary (see `jobs/config.sh`). Do **not** `conda activate` an environment that
ships bcftools 1.6 for the calling steps — it re-shadows the PATH.

Investigators replicating this pipeline on architectures where a modern
bcftools is not available through standard package channels (e.g.
`linux-ppc64le` bioconda, pinned to 1.6) should build a recent bcftools from
source before running `jobs/05_1_mpileup_merge.sh`.

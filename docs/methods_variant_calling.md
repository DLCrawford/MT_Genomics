# Methods — Variant calling, normalization, and SnpEff annotation

> Manuscript-ready prose for the variant-calling and functional-annotation
> sections of the MT_Genomics2 methods. Keep this file in sync with the
> actual canonical pipeline (`jobs/` stages 04–05, `scripts/` stages 06–07).
> Numbers verified against `CLAUDE.md` and `docs/01_pipeline.md` as of
> 2026-05-24, session 18; re-check against the final dataset before
> submission.

## Variant calling and normalization

Variant calling was performed with bcftools mpileup piped to bcftools
call (htslib/bcftools 1.23.1, built from source) using
`-Q 30 -q 30 -d 100000 -a AD,DP` for pileup and `-mv -A --ploidy 1`
for calling; the haploid ploidy setting reflects the single-copy
character of mitochondrial DNA. Variants were called per sample to
retain coverage-honest per-cell allele depths and to avoid the
allele-frequency prior that joint multi-sample calling imposes at
sites where the reference allele is rare — a substantial concern in
these populations, which segregate a divergent mt clade against the
GenBank reference NC_012312.1. Per-sample VCFs were normalized with
`bcftools norm -m -any -f REF`, which left-aligns indels, trims
redundant bases, and splits multiallelic records into one row per
(POS, ALT) combination — so that each row in the resulting VCF
represents a single biallelic substitution rather than a list of
alternate alleles sharing one record — and then joined with
`bcftools merge -m none` to retain per-sample allele depths required
for downstream heteroplasmy-aware haplotype calls. The resulting
merged VCF (1,128 SNPs across 143 samples; ts/tv = 8.17) was
annotated with SnpEff (v5.2) against a custom database built from
the RefSeq GFF for NC_012312.1 under the vertebrate mitochondrial
genetic code (NCBI translation table 2), then restricted to CDS
positions and SNPs only to produce the canonical file used in all
downstream analyses.

## Why split multiallelics before SnpEff

The pre-SnpEff `norm -m -any` split is essential to the precision of
downstream annotation. SnpEff annotates each ALT independently and
writes per-ALT entries into the comma-separated `ANN` INFO field, but
when multiple ALTs share a record, downstream tools that consume the
primary annotation (`ANN[0]`) must reconcile each ALT with its
corresponding position in the `ANN` string — a join that silently
fails when intermediate processing keys on `(CHROM, POS)` alone,
dropping minor-ALT annotations at multiallelic sites. Splitting
before annotation guarantees that every record carries a single ALT
with its annotation unambiguously assigned, which makes downstream
analyses — heteroplasmy classification, dN/dS estimation, haplotype
calling — robust to multiallelic sites without per-tool
reconciliation logic.

## bcftools version cautionary note

We note that bcftools versions prior to 1.10 — which predate the
mpileup engine rewrite — produced substantially fewer variant calls
on these high-depth mtDNA pileups: in a three-sample comparison,
bcftools 1.6 returned 6 SNPs versus 284 with bcftools 1.22 from the
same BAMs and parameters, without emitting an error or warning.
Investigators replicating this pipeline on architectures where
modern bcftools is not available through standard package channels
(e.g. linux-ppc64le bioconda, which pins to 1.6) should build from
source.

---

## Notes for finalization

- **SnpEff version.** Confirm v5.2 against the JAR actually used
  (`~/micromamba/envs/SNP_env/share/snpeff-5.2-1/snpEff.jar`).
- **Final SNP count.** This file currently quotes the 1,128-SNP
  May-15 per-sample → merge baseline. If results tables in the
  manuscript reflect a different cut (e.g., post-CDS, post-sample
  exclusion), update the inline number for consistency.
- **Architectural note.** The current path forward (`jobs/05_1_mpileup_merge.sh`)
  uses a single-pipeline joint mpileup + call + norm-split, with a
  hard regression gate against the per-sample baseline described
  here. If a future version of the manuscript uses 05_1 output as
  the canonical, the second sentence of paragraph 1 ("Variants were
  called per sample ...") needs to be rewritten to describe the
  joint architecture and the regression gate; the SnpEff and
  versioning paragraphs stay as-is.

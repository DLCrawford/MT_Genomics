# MT_Genomics

Mitochondrial haplotype/variant pipeline for *Fundulus heteroclitus*.
This directory mirrors `/projectnb/dcrawford/MT_Genomics2/` on Triton 2 (`t2.idsc.miami.edu`) one-to-one.

## Layout

```
MT_Genomics2/
├── README.md          — this file
├── CLAUDE.md          — pipeline objectives + rules for Claude Code
├── CHANGELOG.md       — session log
├── .gitignore         — excludes large data files
├── docs/
│   ├── 00_setup.md    — HPC + GitHub setup, conda env, reference indexing
│   └── 01_pipeline.md — step-by-step walkthrough of jobs/01..05
└── jobs/
    ├── 01_fastqc_raw.sh
    ├── 02_trim_pe.sh
    ├── 03_fastqc_trimmed.sh
    ├── 04_bwa_align_mt.sh
    └── 05_bcftools_mpileup_call_AD.sh
```

## Sync (Mac → Triton 2)

```bash
rsync -avzhm --progress \
  ~/Projects/MT_Genomics_Cl_Ap2026/MT_Genomics2/ \
  dcrawford@t2.idsc.miami.edu:/projectnb/dcrawford/MT_Genomics2/
```

Add `--delete` once you're confident — that will mirror deletions too.

## Reference

- Genome: `Fundulus_heteroclitus-3.0.3` (GCF_000826765.1)
- MT chromosome: `NC_012312.1`
- Annotation source: `MU-UCD_Fhet_4.1` (GCF_011125445.2)

# CLAUDE.md — Project rules for Claude Code

## Project goal

Build a reproducible mitochondrial variant + haplotype pipeline for *Fundulus heteroclitus*, with heavy compute on Triton 2 (LSF/BSUB) and downstream analysis on Doug's Mac.

## RESOLVED — variant-count discrepancy (sessions 6–10, 2026-05-09 → 2026-05-15)

**Root cause:** bcftools **1.6** in the linux-ppc64le `mito_genomics`
conda env (pinned by bioconda packaging, not by choice) silently
collapsed variant calls on high-depth mtDNA pileups. The session-7
joint-`-mv` architectural finding and the session-8 reference-identity
hypothesis were both contributing factors but not the dominant variable.

Diagnostic chain on the same 3 BAMs and recipe:

    bcftools 1.6   →   6 SNPs
    bcftools 1.22  → 284 SNPs
    bcftools 1.23.1 → 289 SNPs

Full panel re-call with slim mt-only BAMs feeding per-sample call → merge,
all on locally-built bcftools 1.23.1, produced **1128 SNPs / 143 samples /
ts/tv 8.17**, recovering the historical `merged_144.vcf.gz` baseline
(1133 SNPs / 144 samples / ts/tv 7.92). The 1142 / 152 / 227 / 284
number lineage is closed.

**Canonical caller** (as of 2026-05-18, session 14):

  - slim BAMs:  `jobs/BSUB_Slim_BAM_mt.sh` → `MT_only_bams/${SAMPLE}_MT_only.bam`
  - joint call + split: `jobs/05_1_mpileup_merge.sh` →
                `vcf/05_1_Fhet_mt_persample_merged.vcf.gz`
                (joint `bcftools mpileup -b slim_bamlist` + `call -mv --ploidy 1`
                + `norm -m -any -f REF` in a single pipeline; per-cell DP/AD
                on every (POS × sample) by construction; output already in
                one-row-per-ALT form for SnpEff)
  - bcftools / samtools / htslib 1.23.1 built from source at `$HOME/software/local`
  - PATH wired in `jobs/config.sh`

**Validation gate (sessions 13–14).** Joint and per-sample callers can
legitimately produce different variant sets at marginal positions. The
05_1 script diffs its unique-position set against the May-15 per-sample →
merge baseline (`Fhet_mt_persample_merged.vcf.gz`, 1128 SNPs / 143
samples / ts/tv 8.17 — the result of the 3-week session 6–10
diagnosis), writes `*_vs_baseline_lost.tsv` and `_gained.tsv`, and
surfaces the counts in its manifest. **Sites lost vs baseline must be
0** before flowing into stages 06 / 07; sites gained are fine (joint
can recover singletons per-sample emitted at low confidence). This gate
honors the 3 weeks of work that established the per-sample baseline:
the new caller has to prove it preserves that floor, not assume it.

Architecture trajectory:

  - **Sessions 6–10** (per-sample → merge, the established baseline):
    `jobs/_archive/05d2_persample_call.sh` + `_archive/05e2_merge_persample.sh`.
    Produced 1128 SNPs / 143 samples / ts/tv 8.17. The reference data set.
    Structural defect: `bcftools call -mv` per sample only emits a
    sample's variant rows, so `bcftools merge -m none` leaves
    non-variant samples as `.:.:.:.` for DP/AD at every variant
    position. Backfill in stage 07 was the patch.
  - **Session 12** (failed backfill): stage 07 tried to fill those
    cells with `bcftools call -C alleles -T` against a targets file
    built from the post-`norm -m -any` VCF — which has duplicate
    (CHROM, POS) rows at multi-row positions. `call -C alleles -T`
    keys on (CHROM, POS) and silently drops the secondary entry, so
    every multi-row site lost its secondary ALT. 17 of 931 sites
    quietly disappeared.
  - **Session 13** (joint switch, first cut): wrote `jobs/05f_joint_call.sh`
    to replace the per-sample chain. Right shape but missed the
    norm-split step inside the caller (it was deferred to stage 07),
    and didn't honor the baseline as a hard regression gate.
    Archived.
  - **Session 14** (this state): `jobs/05_1_mpileup_merge.sh` does
    joint mpileup + call + norm-split in one pipeline (norm-split
    lives in 05_1, matching where 05d2 did it per-sample), with the
    baseline-vs-05_1 regression check as a hard gate before 06 / 07
    are allowed to run.

Methods write-up still owed: addendum to `docs/02_calling_architecture.md`
documenting (a) the bcftools-version finding from sessions 6–10 with
the 3-sample version table above, (b) the structural defect that
motivated the architecture switch, and (c) the baseline-regression
gate as a project norm for caller swaps.

Experimental artifacts retained for the methods comparison:
`jobs/05_*.sh` (v1, joint strict), `jobs/05b_*.sh` (v2, joint relaxed +
haploid), `jobs/05c_*.sh` (v3, joint no-MAPQ + haploid),
`jobs/05d_persample_call.sh` + `jobs/05e_merge_persample.sh` (per-sample
→ merge against full BAMs with the env's bcftools 1.6).
`jobs/_archive/05d2_persample_call.sh` + `_archive/05e2_merge_persample.sh`
(per-sample → merge against slim BAMs with bcftools 1.23.1; produced
the 1128-SNP May-15 baseline that 05_1 validates against).
`jobs/_archive/05f_joint_call.sh` (session-13 first-cut joint caller,
superseded by 05_1). `jobs/_archive/07b_backfill_AD.sh` (T2-side
backfill, retired — backfill no longer needed). Stats files in
`vcf/Fhet_mt_*_stats.txt` and `vcf/Fhet_mt_persample_merged_stats_15May.txt`
(the validated baseline's contemporaneous stats).

## Known pitfalls + lessons learned (2026-05-09 → 2026-05-17)

Three bugs have shaped the current architecture. Future Claude (and
future Doug) should keep these in mind before touching the caller or
the canonical-VCF pipeline.

### 1. bcftools 1.6 silently collapses high-depth mtDNA joint calls

The bioconda `mito_genomics` env on Triton 2's linux-ppc64le ships
bcftools 1.6 (pinned by the architecture, not by choice). At mtDNA
depths (5–20k × per sample), it returns 6 SNPs from a 3-sample joint
call that 1.22 returns 284 SNPs on and 1.23.1 returns 289 on. Full
diagnosis in CHANGELOG sessions 6–10. **The fix is in production**:
`config.sh` PATH-injects `$HOME/software/local/bin` (1.23.1 built
from source) before anything that uses bcftools. Don't `conda
activate mito_genomics` for variant calling — it re-shadows the
PATH back to 1.6. The exports in `config.sh` after `conda activate`
exist to undo that if a script must activate the env for some other
reason.

### 2. Per-sample `bcftools call -mv` → merge leaves REF cells as `.:.:.:.`

Architectural defect, *not* a bug in any single script. `bcftools
call -mv` emits only variant rows per sample, so a per-sample VCF
contains only the positions where that sample had a variant. When
`bcftools merge -m none` unions 143 per-sample VCFs, samples that
didn't carry a variant at a position have no row to contribute and
the merge represents the missing data as `.:.:.:.`. This is what
05d2 + 05e2 produced and what motivated stage 07's joint-pileup
backfill.

**The fix in production:** single-pipeline joint `bcftools mpileup`
+ `call -mv --ploidy 1` + `norm -m -any -f REF`
(`jobs/05_1_mpileup_merge.sh`). Joint mpileup pileups every sample at
every variant position, so DP/AD is filled in every (POS × sample)
cell by construction. norm-split inside the script gives one row per
(POS, ALT) for SnpEff. No backfill needed.

**Rule:** don't reintroduce a per-sample → merge architecture for
the canonical caller. It cannot fill DP/AD on REF cells without a
downstream backfill, and any such backfill has to dodge pitfall #3.

### 3. Targets-file format for `bcftools call -C alleles -T`

`bcftools call -C alleles -T file` expects targets in the format
**one row per position with all ALTs comma-joined**:

    NC_012312.1<TAB>2847<TAB>T,C,A

It builds a `(CHROM, POS) → alleles` lookup and **keeps only one
entry per position**. So if you build the targets file with
`bcftools query -f '%CHROM\t%POS\t%REF,%ALT\n'` on a VCF that has
been `norm -m -any`-split (one row per ALT), the targets file gets
duplicate (POS) rows at every multiallelic site, and the secondary
ALT is silently dropped by `call -C alleles -T`. This is the bug
that lost 17 of 931 sites in the session-11/12 stage-07 backfill —
every one a split multiallelic site whose minor ALT vanished.

**The fix in production:** the new architecture (05f → 06 → 07)
doesn't use `-C alleles -T` at all. `bcftools call -mv` emits
multiallelic sites as one row with comma-joined ALTs (the native
multi-ALT representation), and stage 07's `norm -m -any` splits
them after the call — preserving per-cell AD on each split row.

**Rule:** never build a `-C alleles -T` targets file from a
post-norm-split VCF. If a future stage needs force-calling, build
the targets either before `norm -m -any` (one row per multi-ALT
position, comma-joined) or consolidate post-norm using
`awk '{ key=$1"\t"$2; if (key in a) a[key]=a[key]","$4; else
a[key]=$3","$4 }' END { for (k in a) print k"\t"a[k] }'` or
equivalent.

### Validation principle (added session 13, sharpened session 14)

When swapping the caller architecture, do not declare the new output
canonical until it has been diff'd against the prior canonical's
unique-position set. Joint and per-sample callers can legitimately
differ at marginal sites — the new caller must find at least every
position the prior canonical found, with **zero regressions**.
`jobs/05_1_mpileup_merge.sh` includes this check (writes
`*_vs_baseline_lost.tsv` / `*_vs_baseline_gained.tsv` against
`Fhet_mt_persample_merged.vcf.gz`, the 1128-SNP May-15 baseline from
the 3-week session 6–10 diagnosis, and surfaces both counts in the
manifest). The May-15 baseline is the *floor*; the new caller's job
is to preserve every site that floor established and ideally find
more. This norm honors the 3 weeks of work that established the
floor — the new caller has to prove itself against it, not assume
parity.

## Downstream analysis scripts (sessions 15+)

Once the canonical `Fhet_MT_CDS.snps.split.vcf.gz` exists and is
validated, the following Mac-side scripts run on top of it. Order is
flexible — none of these write back into the canonical.

### Heteroplasmy frequency biology — NOT 50/50 (added session 17)

A note on expectations for any Hp script (09 / 14 / future):
heteroplasmy frequencies do **not** cluster at 50/50, and the
`0.10 ≤ AD/DP < 0.70` window is correctly broad. Three mechanisms
push Hp off 50/50:

- **Paternal leakage** can deposit a small amount of paternal mt
  into a zygote, yielding Hp at well below 50 % (often a few %).
- **Germline bottleneck severity** during oogenesis randomly
  reduces mt-haplotype diversity in the egg before fertilization.
  A severe bottleneck (n_mt sampled ≈ 1) drives Hp toward 0 or
  near-100 % (single haplotype carried forward); a mild
  bottleneck (n_mt sampled ≈ 1000s) preserves population-level
  Hp ratios closer to maternal input.
- **Hitchhiking on an old neutral mutation:** a lineage that
  diverged long ago will carry many co-segregating variants on
  its background, so a fish inheriting both lineages shows Hp at
  *every* site that differs between them — not just the
  originating mutation. Several hundred shared-allele Hp events
  in one individual is consistent with two-lineage admixture
  (see session-17 finding on 47/77/33/84), not with somatic
  mutation.

Implication: do **not** tighten the AF window around 0.5 to
"clean up" Hp calls; doing so would systematically discard
paternal-leakage / mild-bottleneck heteroplasmies. AF ≈ 0.2 events
on high-DP cells are biologically plausible.

- `scripts/08_call_haplotypes.py` — haplotype matrix + per-sample
  haplotype string from **all** CDS SNPs, calling rule `AD_alt/DP > 0.7`.
  Existing (session-10 vintage). Outputs `haplotype_matrix.csv`,
  `haplotype_calls.csv`. Stage 08 answers "which haplotype clade is
  each sample in?".
- `scripts/10_dnds_per_gene.py` (session 15) — per-gene + Overall dN,
  dS, dN/dS (simple `obs / possible` ratio per user spec) plus
  Jukes-Cantor-corrected variants. Uses **NCBI translation table 2**
  (vertebrate mitochondrial code: TGA=W, ATA=M, AGA/AGG=Stop).
  Strand from GFF (ND6 reverse-complemented). Handles the **polyA-
  completed stop codons** in ND2/COX2/COX3/ND3/ND4 by truncating to
  a codon multiple — previous dN/dS attempts on this genome failed
  on those genes by passing the trailing partial codon to the
  translator. Output: `vcf/dnds_per_gene.tsv`.
- `scripts/11_haplotypes_nonsyn.py` (session 15) — same shape as
  stage 08 but restricted to **MODERATE-impact** (`missense_variant`,
  `splice_region_variant&missense_variant`) — i.e., amino-acid-
  changing variants only. Haplotype strings prefixed `N_` to
  distinguish from 08's `C_`. Outputs `haplotype_matrix_nonsyn.csv`,
  `haplotype_calls_nonsyn.csv`.
- `scripts/12_ns_cooccurrence.py` (session 15) — find missense sites
  whose per-sample call vectors are identical across the panel
  (perfect LD). Classifies each multi-site group as
  `fixed_ref_divergence` (carriers = N_samples), `haplogroup`
  (intermediate, ≥2 carriers), or `singleton_artifact` (carriers = 1).
  Prefers `vcf/haplotype_matrix_nonsyn.csv` if present, falls back
  to filtering `vcf/haplotype_matrix.csv` to missense rows. Output
  `vcf/ns_cooccurrence_groups.tsv`. CLI flag `--min-carriers N`
  hides groups with fewer than N carriers.
- `scripts/09_heteroplasmy_report.py` (session 16, RAN) —
  heteroplasmy classifier on `vcf/MT_DP_AD_141.txt` (the variant-
  only per-cell DP/AD table built by Doug outside this repo).
  Filter `0.1 ≤ AD/DP < 0.7`; per-allele classifier (REF_Hp,
  shared_alt_Hp, private_alt_Hp). Outputs:
  `vcf/heteroplasmy_{events,per_site,per_individual,summary}.{tsv,txt}`.
  **Key finding:** 1,088 Hp events across 53 individuals / 376
  sites; the private-ALT-Hp count is structurally 0 because the
  input only contains panel-variant ALTs (verified empirically:
  927/927 (POS, ALT) combos have ≥1 AF≥0.7 carrier). For a
  biologically meaningful private-ALT-Hp count, use stage 14's
  pileup-based version. **Four individuals (77, 47, 33, 84) hold
  78 % of all Hp events** — likely contamination / quality issue,
  tested quantitatively by stage 15.
- `scripts/13_pileup_cds_AD.sh` (session 16; **RAN session 17**) —
  Mac-side `bcftools mpileup -a AD,DP -d 100000 --no-BAQ` over
  `docs/mito_protein_coding.bed` (13 CDS, 11,417 bp) across all
  141 panel BAMs (excludes 70, 125 by filename filter). Output
  `vcf/pileup_cds_141.vcf.gz`. Per-cell AD coverage-honest at
  every CDS position, independent of which positions/ALTs the
  panel variant caller produced — this is the input stage 14
  needs to expose private ALT-Hp. `--no-BAQ` matches 05_1 so AD
  scale is consistent across the pipeline.
- `scripts/14_hp_from_pileup.py` (session 16; **RAN session 17**) — Hp
  detection on the stage-13 pileup, the proper answer to the
  private-ALT-Hp question. Thresholds: `DP ≥ 20`,
  `0.10 ≤ AD/DP < 0.70`, `AD_Hp ≥ 4`. Same 3-way allele
  classifier as 09 but now `private_alt_Hp` is reachable.
  Outputs: `vcf/heteroplasmy_pileup_*.{tsv,txt}` — saved in
  session 17 with an `_all` suffix
  (`heteroplasmy_pileup_events_all.tsv`, etc.).
  **Session-17 result:** 1,124 Hp events / 64 individuals / 400
  CDS sites; 610 REF_Hp, 487 shared_alt_Hp, **27 private_alt_Hp**
  (vs 0 from stage 09 on the variant-only input — the whole point
  of the 13→14 detour). See CHANGELOG 2026-05-21 for details.
- `scripts/15_well_bleed_test.py` (session 16; **RAN session 17**) —
  quantitative test of "is the high Hp load on individuals 77 /
  47 / 33 / 84 explained by library-prep well bleed?" Uses
  `data_files_May/WGS_seq_plate.txt` (plate identity inferred
  from the `i5` column: i5_3 = plate 1, i5_4 = plate 2; all 4
  focals are on plate 1, scattered at H5/B6/D10/E11). Donor-
  concordance score = fraction of X's Hp events explained by
  Y's haplotype. Two test statistics with 10,000-perm one-sided
  p-values: (a) Spearman rho of Score vs plate distance, (b)
  mean(neighbor Score) − mean(far Score). Permutation fixes X
  and shuffles other same-plate samples across other same-plate
  wells. CLI now accepts `--events <path>` (session-17 addition)
  and auto-falls-back to `*_all.tsv` if the bare filename isn't
  present. Outputs:
  `vcf/well_bleed_{donor_ranking,results,summary}.{tsv,txt}`.
  **Session-17 result: NO well-bleed signal on any of the 4
  focals.** Top donors sit 4–10 wells from each focal (none at
  King-move distance 1); all 8 p-values non-significant; two
  ρ values go the wrong direction. The high score saturation
  (0.875–1.000) at non-neighbor distances is the signature of
  panel-level mt admixture rather than contamination — see
  CHANGELOG 2026-05-21 for the admixed-mito interpretation and
  the n_REF_Hp split among the 4 focals.
- `scripts/18_variant_burden_per_individual.py` (session 18) — per-sample
  panel-variant counts (Total / SYN / NS / Other), joined with the
  per-individual heteroplasmy burden, output as a master per-individual
  table with **all 143 canonical samples** (zero-filled Hp columns for
  samples not in the per_individual_all file). Designed to test the
  hypothesis that the bimodal variant-count distribution (north-clade
  individuals < ~50 ALTs vs south-clade > ~190) correlates with
  per-individual Hp burden. Counts a sample's ALT calls by reading the
  per-sample GT field of `vcf/Fhet_MT_CDS.snps.split.vcf.gz` (ploidy=1;
  GT=1 → carries ALT) and classifying each variant row by ANN[0]
  (same SYN/NS/Other buckets as scripts 16/17). Sample-name
  normalization extracts the leading numeric ID, so it joins stage-09
  (`{N}_MT`) and stage-14 (`MT_only_bams/{N}_0_MT_only.bam`) on the
  same key. Outputs: `vcf/per_individual_burden_pileup.tsv`
  (variant counts + stage 14 Hp), `vcf/per_individual_burden_variant.tsv`
  (variant counts + stage 09 Hp), `vcf/per_individual_burden_summary.txt`
  (low/mid/high variant × Hp summary). **Session-18 finding:** the
  variant-count distribution is sharply bimodal — 77 samples at 15–26
  variants, 65 samples at 220–233, with exactly 1 sample in between
  (individual 77 at 190 variants, one of the admixed focals). Low
  group's mean Hp burden (7.65 sites) is ~1.8× the high group's (4.34)
  and the proportion of samples with any Hp is 51 % vs 37 %. Direction
  consistent with Doug's bimodal hypothesis; the 4 admixed focals (47,
  77, 33, 84) sit at the extremes of the variant-count axis and each
  carry ~190–220 Hp events.
- `scripts/17_annotate_hp_codon.py` (session 18) — codon-level
  follow-up to script 16. For rows that 16 marked `Unannotated`
  (POS, non-REF allele not in canonical), computes SYN/NS from
  first principles using the reference FASTA, the GFF CDS coords
  (strand-aware: reverse-complements for ND6), and NCBI
  translation table 2 — same machinery as
  `scripts/10_dnds_per_gene.py`. For each Unannotated row whose
  POS sits inside a CDS, it builds the REF codon, substitutes the
  non-REF base (complemented for − strand), translates, and
  classifies SYN / NS / `stop_gained_codon` / `stop_lost_codon`.
  Adds columns `Codon_ref`, `Codon_alt`, `AA_ref`, `AA_alt`, and
  `Annotation_source ∈ {snpeff, codon}` so it's clear which
  rows came from SnpEff and which from codon computation. Rows
  with POS outside any CDS get `Effect=non_CDS` (Other). Rows
  with `non_REF_allele = NONE` (insufficient-coverage major
  calls) get `Effect=invalid_alt_allele` (Other). A REF-base
  sanity check (the codon's nucleotide at the substituted
  position must equal REF / complement(REF) for + / − strand)
  guards against GFF/FASTA/event-table positioning drift. Inputs:
  `vcf/heteroplasmy_events_annot.tsv`,
  `vcf/heteroplasmy_pileup_events_all_annot.tsv` (the outputs of
  script 16). Outputs: `vcf/heteroplasmy_events_annot_codon.tsv`,
  `vcf/heteroplasmy_pileup_events_all_annot_codon.tsv`,
  `vcf/heteroplasmy_annot_codon_summary.txt`.
- `scripts/16_annotate_hp.py` (session 18) — joins SnpEff coding-effect
  annotation onto the heteroplasmy event tables from stages 09 and 14.
  Parses `vcf/Fhet_MT_CDS.snps.split.vcf.gz`, builds a `(POS, ALT)
  → Effect_class` map from `ANN[0]` (highest-impact gene-internal
  call), then for each Hp event selects the non-REF allele
  (`Major` if `Hp_is_REF`, else `Hp_allele`) and writes
  `Effect_class ∈ {SYN, NS, Other, Unannotated}` along with the
  raw `Effect`, `Gene`, and `HGVS_p` columns. SYN =
  `synonymous_variant`; NS = the same impact set used by scripts
  10 and 11 (`missense_variant`, `stop_gained`, `stop_lost`,
  `start_lost`, `initiator_codon_variant`). `Unannotated` =
  (POS, ALT) not in canonical — by construction, every
  `private_alt_Hp` event from stage 14 lands here (27 events);
  for these, the codon / translation-table-2 machinery in
  `scripts/10_dnds_per_gene.py` can be used to compute SYN/NS
  from first principles if needed. Inputs:
  `vcf/heteroplasmy_events.tsv`,
  `vcf/heteroplasmy_pileup_events_all.tsv`. Outputs:
  `vcf/heteroplasmy_events_annot.tsv`,
  `vcf/heteroplasmy_pileup_events_all_annot.tsv`,
  `vcf/heteroplasmy_annot_summary.txt` (Hp_class × Effect_class
  cross-tab for both inputs).
- `scripts/DP_AD_table.py` — denormalizes the canonical VCF into a
  long-format TSV (`vcf/mtDNA_long_AD_table.tsv`: one row per
  Individual × Position with DP / AD_REF / AD_ALT1-4). Useful as an
  alternative input for downstream analyses that don't want to
  parse VCF.

### Preliminary findings to carry forward (session 15)

The session-15 ad-hoc analysis of `vcf/haplotype_matrix.csv`
(filtered to its 164 missense rows, 141 samples after the standard
70/125 exclusion) turned up three results worth carrying into the
methods write-up and next-session analysis. `scripts/12_ns_cooccurrence.py`
reproduces all of these.

1. **Three NS sites fixed as ALT in every panel sample** — i.e.,
   `NC_012312.1` (the GenBank mt reference) differs from the entire
   *F. heteroclitus* panel at these 3 amino-acid-changing positions:

   |  Gene  |  POS   | Change |
   |  ----  |  ----  | ------ |
   |  ND1   | 3124   | G>A    |
   |  ND2   | 4680   | A>T    |
   |  ND2   | 4957   | C>T    |

   These should arguably be excluded from population-level dN/dS,
   because they're reference-divergence artifacts rather than
   within-panel polymorphism. To do, future session: flag these as
   excluded in stage 10 and document in the methods.

2. **One major nonsynonymous haplogroup: 7 sites across 4 genes
   in 64 / 141 samples (~45 %).** All 7 sites co-segregate
   perfectly across the panel — half the fish carry the entire
   block, half carry none of it. This is the biggest single
   population-genetic signal in the nonsyn data and is almost
   certainly a major maternal-lineage clade:

   |  Gene  |  POS    | Change |
   |  ----  |  ----   | ------ |
   |  ATP6  |  8451   | A>G    |
   |  ND1   |  3061   | T>C    |
   |  ND2   |  4737   | G>A    |
   |  ND2   |  4951   | T>C    |
   |  ND2   |  5061   | A>G    |
   |  ND5   | 12589   | C>A    |
   |  ND5   | 13577   | T>A    |

   Worth pulling out as a stand-alone result and cross-checking
   against pedigree/collection metadata when available.

3. **Six small co-segregating clusters (2 or 3 NS sites each, 2 or
   4 carriers).** Probably sib-group or close-kin maternal sharing.
   See `vcf/ns_cooccurrence_groups.tsv` (after 12 is run) for
   detail; or session-15 CHANGELOG for the listing.

The remaining "groups" the script will list are at 1 carrier each —
those are aggregations of multiple private NS variants in the same
individual fish (the signature [0,…,0,1,0,…,0] is identical across
sites private to the same sample). Filter them out with
`--min-carriers 2` when looking for population-level LD.

### SnpEff annotation notation primer (sessions 12–15)

The `ANN` INFO field follows the [SnpEff convention][1]:

    Allele | Effect | Impact | GeneName | GeneID | FeatureType |
    FeatureID | TranscriptBiotype | Rank | HGVS.c | HGVS.p |
    cDNA.pos | CDS.pos | AA.pos | Distance | ERRORS

[1]: https://pcingola.github.io/SnpEff/snpeff/inputoutput/#ann-field

Two coordinate systems coexist in a single row:
- **POS** (VCF column 2): genomic 1-based on `NC_012312.1`.
- **`c.Y`** (in HGVS.c, e.g. `c.6T>C`): the **CDS coordinate** —
  position within the *gene's coding sequence*, 1-based from the
  start codon. For + strand genes: `c.Y = POS − (gene_start − 1)`.
  For ND6 (the one − strand gene): `c.Y = (gene_end + 1) − POS`,
  and the REF/ALT in `c.Y_REF>ALT` are the reverse complement of
  the forward-strand REF/ALT shown in the VCF.

Worked examples from the canonical:

| POS  | Gene | Strand | c.Y notation | c.Y derivation        | p.X       |
|------|------|--------|--------------|-----------------------|-----------|
| 2847 | ND1  | +      | `c.6T>C`     | 2847 − (2842 − 1) = 6 | p.Phe2Phe |
| 2847 | ND1  | +      | `c.6T>A`     | same                  | p.Phe2Leu |
| 2852 | ND1  | +      | `c.11C>A`    | 2852 − 2841 = 11      | p.Thr4Asn |

So a row "position 2847, c.6T>C" is *not* a coordinate mismatch —
it's saying: the SNP at genomic position 2847 is the 6th nucleotide
of the ND1 coding sequence. Both numbers describe the same site.

## Next-session pickup (after 2026-06-13, session 20)

Manuscript-support session. No change to the canonical VCF (still frozen).

### New script

- `scripts/20_calc_pi_clade.py` — π (total/syn/NS) split by North vs South
  clade. Pass 1 assigns clades by per-individual ALT burden (north < 50,
  south > 200; in-between excluded — by default only `77_MT` at 193).
  Pass 2 computes π within each clade (reference-independent). Outputs
  `vcf/pi_by_clade_persite.tsv` + `.membership.tsv`. Result: 77 N / 63 S;
  π_total N 0.00160 / S 0.00188; pN/pS N 0.21 / S 0.26. Note π_syn/π_ns
  use full L_CDS denominator (as in 19); use `--L_syn`/`--L_ns` for true
  per-site rates.

### How to determine "what was actually run" — read the VCF header

The authoritative record of every bcftools operation is the ordered
`##bcftools_<cmd>Command=` + `##bcftools_<cmd>Version=` lines stamped into
each VCF header (with dates). More reliable than file mtimes or script
intent. Used this session to resolve a Methods/repo mismatch:

- **The 927-variant set came from JOINT calling** (`jobs/05_1_mpileup_merge.sh`):
  `141_MT_variants.vcf.gz` header = single `call -mv --ploidy 1` →
  `norm -m -any -f REF` (1.23.1, 18 May 2026) → SnpEff → CDS-restrict →
  `view -v snps` → sample subset → `view -e AC=0`. No per-sample paths,
  no `merge`. The per-sample→merge file (`Fhet_mt_persample_merged.vcf.gz`,
  15 May, header carries `merge -m none`) was **not** used for results.
- Side corrections (same source): executed `call` was `-mv --ploidy 1`
  (no `-A`); `mpileup` used defaults `-Q 13 -q 0` (not `-Q 30 -q 30`).

### Panel-wide MAPQ (all 143 mt BAMs, 36.27M reads)

99.99% MAPQ ≥ 30; 99.69% MAPQ = 60; 0.0009% MAPQ = 0; only 0.0077% would be
dropped by `-q 30`. Individual 77 (lowest reads): 99.25% MAPQ 60, zero
MAPQ-0 → low depth, not mapping ambiguity/NUMT. MAPQ is column 5 of each BAM
record; parsed in pure Python (no samtools needed). Per-sample cache:
`outputs/mapq_cache.tsv`.

### Manuscript deliverables (repo root)

- `Mitochondrial_variants_5_pipeline-refs.docx` — tracked pipeline refs,
  joint-calling Methods rewrite, corrected flags + MAPQ sentence, 4-reason
  NUMT argument, **Code and Data Availability** section
  (`https://github.com/DLCrawford/MT_Genomics`).
- `Supplemental_Table_S1_Pipeline.docx` — per-step file→tool→IO table.

### GitHub release plan (the "provide all scripts/jobs" goal)

Remote: `git@github.com:DLCrawford/MT_Genomics.git`. The repo is far behind:
last commit is the 05d/05e per-sample era. **Untracked and needing commit:**
`jobs/05_1_mpileup_merge.sh`, `jobs/BSUB_Slim_BAM_mt.sh`, `jobs/_archive/*`,
`scripts/09–28`, `scripts/make_141_vcf.{sh,py}`, `scripts/_archive/*`, and
several `docs/*`. Deleted-but-still-tracked: old `jobs/05*` (now in
`_archive/`). Plan, in order:

1. Add `*.docx` to `.gitignore` (manuscripts don't belong in the code repo)
   and `git rm --cached` the stray `scripts/Mitochondrial_variants_*.docx`.
2. `git add` jobs/, scripts/, docs/, README.md, CLAUDE.md, CHANGELOG.md;
   stage the deletions of the old 05* scripts.
3. Commit ("session 20: full pipeline + downstream analysis scripts;
   manuscript pipeline references"), push to origin.
4. Tag a release (e.g. `v1.0-submission`) and mint a Zenodo DOI; put the DOI
   in the manuscript Code-and-Data-Availability section.

Exact commands are in the session-20 CHANGELOG entry / chat.

## Next-session pickup (after 2026-06-05, session 19)

### Comparative π analysis — NEW (session 19)

A cross-species nucleotide diversity (π) comparison has been started in a
separate directory tree alongside `MT_Genomics2/`. The goal is to place
Fhet mt π in context against four other datasets. Scripts are numbered
19–24 and live in their respective dataset folders.

**Status as of 2026-06-05:**

| Dataset | Script | Status | Output |
|---------|--------|--------|--------|
| Fhet (141 samples) | `MT_Genomics2/scripts/19_calc_pi.py` | **DONE** | `vcf/pi_results.tsv` |
| Fhet by clade (N/S) | `MT_Genomics2/scripts/20_calc_pi_clade.py` | **DONE** (session 20) | `vcf/pi_by_clade_persite.tsv` |
| Drosophila DGRP (169 lines) | `MT_Genomics2/scripts/20_dros_pi.py` | **DONE** | `vcf/dros_pi_results.tsv` |
| Human mt Lankheet 2026 (1,176 genomes) | `Human_mt/21_human_mt_pi.py` | **DONE** | `Human_mt/human_mt_pi_per_site.tsv` |
| Human mt CDS syn/NS | `Human_mt/22_human_mt_cds_pi.py` | **DONE** | `Human_mt/human_mt_cds_pi_per_site.tsv` |
| C. elegans (540 samples) | `C_elegans/23_celegans_pi.py` | **DONE** | `C_elegans/celegans_pi_per_site.tsv` |
| Yeast (469 complete assemblies) | `Yeast/24_yeast_pi.py` | **IN PROGRESS** | TBD |

**Canonical Fhet VCF for all downstream analysis:** `vcf/141_MT_variants.vcf.gz`
- 141 samples (excludes 70, 125)
- 927 records (3 private-to-70/125 sites removed via `AC=0` filter)
- Built by `scripts/make_141_vcf.sh` from `vcf/Fhet_MT_CDS.snps.split.vcf.gz`

**Key π results so far (L_CDS normalised):**

| Dataset | N | π_total | π_syn | π_ns | pN/pS |
|---------|---|---------|-------|------|-------|
| Fhet    | 141 | TBD (rerun on 141_MT_variants.vcf.gz) | | | |
| Drosophila | 169 | 0.00052 | 0.00039 | 0.00013 | 0.339 |
| Human (Lankheet) | 1,176 | whole-mt done; CDS done | | | |
| C. elegans | 540 | run complete | | | |
| Yeast | 469 | pending | | | |

**Status — COMPLETE for 4/5 datasets. Yeast π unreliable.**
Run `python scripts/25_comparison_table.py` to regenerate the table.

**Open items:**
1. Yeast π: script 24 produced inflated values (π=0.14, expected ~0.008)
   due to spurious local alignments on AT-rich assemblies. Needs minimap2
   or k-mer pre-screening. Use literature θ_W (0.00766) as placeholder.
2. Human: θ_W in table uses AMR (N=5,718); recalculate from Lankheet at N=1,176
3. Confirm Fhet π on `vcf/141_MT_variants.vcf.gz` matches table values

**Yeast assembly notes:**
- Directory: `~/Projects/MT_Genomics_Cl_Ap2026/Yeast/mitochondrialAssemblies/`
- 905/1011 assemblies present (106 missing from repository)
- Types: 303 oneScaff, 166 circularized, 436 multiscaff (fragmented — excluded)
- Using 469 complete (oneScaff + circularized) for π
- Reference: S288C mt genome NC_001224 (~85 kb); CDS L = 6,684 bp
- MAFFT not available (conda conflict, no brew); using Biopython PairwiseAligner
- Runtime estimate: several hours for 469 × 85 kb alignments

**Human mt notes:**
- FASTA: `Human_mt/lankheet_2026_mt.fasta` (1,176 genomes, ~20 MB)
- Reference rCRS: `Human_mt/rCRS_NC_012920.fasta` (downloaded by script 21)
- CDS coordinates hardcoded in `22_human_mt_cds_pi.py` from NC_012920 annotation
- L_CDS = 11,395 bp (13 protein-coding genes, matches Table 2)

---

## Next-session pickup (after 2026-05-21, session 17)

**Heteroplasmy thread (stages 13/14/15) is RESOLVED — see CHANGELOG
2026-05-21.** Stage 14 returned 1,124 Hp events / 64 individuals
and exposed 27 `private_alt_Hp` events that stage 09 could not.
Stage 15 ruled out well bleed on all 4 focals; the n_REF_Hp split
(47/77 ≈ all REF_Hp; 33/84 ≈ all shared_alt_Hp) is the signature
of panel-level north × south mt admixture, not contamination.

**Carry-overs still open:**

- Stages 10, 11, 12 from session 15 are still unrun (dN/dS,
  nonsyn haplotypes, NS co-occurrence). They don't block the
  heteroplasmy thread but should land before the methods
  write-up.
- Methods write-up for stages 13/14/15 + the admixed-mito
  finding (probably a new `docs/03_heteroplasmy.md`).
- Cross-check the 4 focals' `MITOTYPE` column in
  `data_files_May/WGS_seq_plate.txt` against the n_REF_Hp
  split — if 47/77 are `S` and 33/84 are `N` or `A`, the
  pileup classification reproduces the lab-prep coding.
- Tabular review of the 27 `private_alt_Hp` events to decide
  whether they cluster by individual (somatic candidate) or by
  POS (sub-threshold panel variant / recurrent artifact).
- Stage 14 output naming: it currently writes `*_all.tsv`.
  Stage 15 handles either via the new `--events` flag fallback.
  Decide if `_all` is the new convention or rename to bare names.

The historical pickup block (the bash commands to run 13/14/15)
is preserved below in case the heteroplasmy thread ever needs a
re-run, e.g. with relaxed thresholds.

```
cd ~/Projects/MT_Genomics_Cl_Ap2026/MT_Genomics2
conda activate SNP_env

# --- carry-over from session 15 (still unrun): ---
python scripts/10_dnds_per_gene.py            # → vcf/dnds_per_gene.tsv
python scripts/11_haplotypes_nonsyn.py        # → vcf/haplotype_matrix_nonsyn.csv
                                              # → vcf/haplotype_calls_nonsyn.csv
python scripts/12_ns_cooccurrence.py --min-carriers 2
                                              # → vcf/ns_cooccurrence_groups.tsv

# --- new this session (the heteroplasmy thread): ---
bash   scripts/13_pileup_cds_AD.sh            # 1-3 min
                                              # → vcf/pileup_cds_141.vcf.gz (+ .tbi)
                                              # → vcf/slim_bamlist_141.txt
python scripts/14_hp_from_pileup.py           # ~30 s
                                              # → vcf/heteroplasmy_pileup_*.{tsv,txt}
python scripts/15_well_bleed_test.py          # ~30 s
                                              # → vcf/well_bleed_*.{tsv,txt}
```

### Open scientific questions, in priority order

1. **(Stage 14) How many TRUE private ALT-Hp events are there?**
   Stage 09 returned 0 for that category, but it ran on
   `MT_DP_AD_141.txt`, which only contains panel-variant ALTs —
   verified empirically: 927/927 (POS, ALT) combos in the file
   have at least one AF ≥ 0.7 carrier, so by construction every
   ALT is somebody's haplotype. Stage 14's input is the per-CDS-
   position pileup (stage 13), which is coverage-honest and
   doesn't gate on AC ≥ 1. A non-trivial count of private ALT-Hp
   in 14's output = candidate true somatic / private
   heteroplasmies. **This is the headline number the whole
   thread was built to produce.**

2. **(Stage 15) Are 77 / 47 / 33 / 84 explained by well bleed?**
   **RESOLVED — session 17 (2026-05-21): NO.** All 8 test
   statistics (4 focals × 2 tests) non-significant; top-donor
   distance 4–10 wells in every case (no King-move neighbors).
   Two ρ values go the wrong direction. Score saturation
   (top-donor concordance 0.875–1.000 at non-neighbor distances)
   is the signature of panel mt admixture, not contamination.
   See CHANGELOG 2026-05-21 for the per-focal table and the
   n_REF_Hp split that separates 47/77 (major-south + north
   minor) from 33/84 (major-north + south minor).

3. **Should high-Hp individuals be excluded from population-level
   analysis** (in addition to 70 and 125)?
   **Provisional answer (session 17): NO** — they're informative
   admixed individuals, not artifacts. Open: whether stage 08's
   0.7 binary haplotype call captures them sensibly, since both
   north and south signals are present in roughly equal weight.

4. **MITOTYPE column in WGS_seq_plate.txt** (`S` / `N` / `A`) is
   undocumented in the project. Worth cross-checking against
   stage 14's per-individual Hp burden and against the 7-site
   nonsyn haplogroup from session 15 — these may all be the
   same partition.

### Carry-overs still open from session 15

- Should the 3 fixed reference-divergence sites (ND1:3124,
  ND2:4680, ND2:4957) be excluded from dN/dS? Recommendation:
  yes, with a documented note. See "Preliminary findings to
  carry forward" above.
- Cross-check the 7-site major haplogroup against pedigree /
  collection metadata — the 64-fish block strongly suggests a
  clade.

## Active primary task — haplotype call set + heteroplasmy report (2026-05-18, session 14)

Pipeline state as of session 14:
- **Slim BAMs** on T2 at `/projectnb/dcrawford/MT_Genomics2/MT_only_bams/`
  (143/143) and mirrored on Mac. From `BSUB_Slim_BAM_mt.sh`.
- **05_1 caller** (`jobs/05_1_mpileup_merge.sh`, new this session,
  replaces 05f): single-pipeline joint `mpileup -b slim_bamlist | call
  -mv --ploidy 1 | norm -m -any -f REF` on T2. Output
  `vcf/05_1_Fhet_mt_persample_merged.vcf.gz` with per-cell DP/AD on
  every (POS × sample) and already in one-row-per-ALT form. **Not yet
  submitted** — first run pending.
- **Per-sample baseline**: `vcf/Fhet_mt_persample_merged.vcf.gz`
  (1128 SNPs / 143 samples / ts/tv 8.17, May 15) is on T2 and is the
  regression floor 05_1 validates against. Contemporaneous stats
  preserved as `vcf/Fhet_mt_persample_merged_stats_15May.txt`.
- **06** Mac-side (`scripts/06_snpeff_mac.sh`): input/output filenames
  point at the 05_1 tag (`05_1_Fhet_mt_persample_merged{,_ann}.vcf.gz`).
  Same SnpEff config (`NC_012312.1`, standalone `~/snpEff/snpEff.jar`,
  env Java).
- **07** Mac-side (`scripts/07_cds_snps_norm_mac.sh`) simplified
  further: CDS-restrict + SNPs-only + sample rename. **No `norm -m
  -any` here** — that step lives in 05_1 now, matching where 05d2
  did the per-sample split. ANN from stage 06 passes through.
- **Archived** (`jobs/_archive/`): `05d2_persample_call.sh`,
  `05e2_merge_persample.sh`, `05f_joint_call.sh`,
  `07b_backfill_AD.sh`.
- **08 / 09** still queued behind 07 output.

Steps remaining (in order):

1. **Submit 05_1 on T2:**
   ```
   bsub < jobs/05_1_mpileup_merge.sh
   # logs: /projectnb/dcrawford/MT_Genomics2/logs/05_1_mpileup_merge_<jobid>.{out,err}
   # walltime expectation: 30–60 min on 143 slim BAMs
   ```

2. **Validate 05_1 against the per-sample baseline (HARD GATE).**
   Open `vcf/05_1_Fhet_mt_persample_merged_run_manifest.txt` and check:
   - `cells with DP=.: 0` — joint mpileup must fill every cell.
   - `n_multiallelic_sites: 0` — norm-split happened (one row per ALT).
   - `positions lost vs baseline: 0` — 05_1 must find every position
     `Fhet_mt_persample_merged.vcf.gz` (1128-SNP baseline) found. Any
     value > 0 is a regression — **do not proceed**. Inspect
     `vcf/05_1_Fhet_mt_persample_merged_vs_baseline_lost.tsv` and
     diagnose before flowing downstream.
   - n_SNPs in the 1100–1200 neighborhood (sanity check).

   `positions gained vs baseline: N` is fine and expected — joint
   calling can recover singletons per-sample missed.

3. **rsync 05_1 output to Mac:**
   ```
   rsync -avP \
     dcrawford@t2.idsc.miami.edu:/projectnb/dcrawford/MT_Genomics2/vcf/05_1_Fhet_mt_persample_merged.vcf.gz \
     dcrawford@t2.idsc.miami.edu:/projectnb/dcrawford/MT_Genomics2/vcf/05_1_Fhet_mt_persample_merged.vcf.gz.csi \
     ~/Projects/MT_Genomics_Cl_Ap2026/MT_Genomics2/vcf/
   ```

4. **`bash scripts/06_snpeff_mac.sh`** → `vcf/05_1_Fhet_mt_persample_merged_ann.vcf.gz`.
   SnpEff annotates each split row independently (one ANN per (POS,ALT)).

5. **`bash scripts/07_cds_snps_norm_mac.sh`** → canonical
   `vcf/Fhet_MT_CDS.snps.split.vcf.gz` with per-cell DP/AD + ANN.
   Verify the "FREEZE GATE" lines: `cells with DP=. = 0` and
   `records w/ ANN = record count`. Frozen only when those two AND
   step 2's `positions lost = 0` are all green.

6. **`python scripts/08_call_haplotypes.py`** — haplotype matrix at
   the existing **0.7 binary** threshold. Stage 08 answers "which
   haplotype clade is each sample in?", which is a strict-alt
   question — no revision.

7. **Write `scripts/09_heteroplasmy_report.py`** — separate analysis
   on the same canonical, pulling the orthogonal signal that 08
   discards:
     - per-cell filter: 0.1 ≤ AD_alt/DP < 0.7 AND AD_alt ≥ 50
     - per-site classifier: n_carriers; transmitted (> 1) vs private
       (= 1). Transmitted = shared-heteroplasmy candidate, strongly
       suggestive of maternal-line transmission.
     - outputs: per-site het tables (transmitted + private),
       per-sample het burden, optional pairwise shared-het adjacency.

8. **Validation against prior outputs.** Compare 08's haplotype
   matrix to the prior ~243-variant output; compare 09's
   transmitted-het count to the historical 868 low-AF sites from
   `Fhet_mt_fullAD.vcf.gz`.

## Compute split

- **HPC (Triton 2, LSF)** — anything heavy:
  - FastQC across 287 paired raw fastqs
  - Trimmomatic PE
  - FastQC on trimmed reads
  - bwa mem + samtools sort/index against `Fhet_MT.fasta`
  - Slim mt-only BAM extraction (`jobs/BSUB_Slim_BAM_mt.sh`) — cuts each
    15–19 GB WGS BAM to ~50–100 MB by extracting NC_012312.1 alignments only.
  - **JOINT variant calling + norm-split** (`jobs/05_1_mpileup_merge.sh`):
    one pipeline,
    `bcftools mpileup -b slim_bamlist -a AD,DP --max-depth 100000 -Ou |
     bcftools call -mv --ploidy 1 -Ou |
     bcftools norm -m -any -f REF -Oz -o vcf/05_1_Fhet_mt_persample_merged.vcf.gz`.
    Per-cell DP/AD on every sample at every variant position (no
    `.:.:.:.` REF cells, no backfill needed) AND one-row-per-ALT
    representation for SnpEff. Uses locally-built bcftools 1.23.1 from
    `$HOME/software/local/bin/` (PATH-injected by `jobs/config.sh`).
    The script writes `*_vs_baseline_lost.tsv` /
    `*_vs_baseline_gained.tsv` and surfaces a regression count in the
    manifest against `Fhet_mt_persample_merged.vcf.gz` (1128-SNP May-15
    baseline). Lost = 0 is the precondition for flowing downstream.
  - Final canonical output: `Fhet_MT_CDS.snps.split.vcf.gz` ← FROZEN once
    produced AND validated, do not rerun (built Mac-side by stage 07).

  Historical reference points (don't compare directly to 05_1 — different
  read sets / tool versions; these are sanity-band landmarks):
  `Fhet_mt_fullAD.vcf.gz` (Oct 2025, joint `-mv -A` on full BAMs,
  bcftools 1.22) = 1142 records / 152 SNPs; `merged_144.vcf.gz`
  (Jul 2025, per-sample → merge, 144 samples) = 1140 records / 1133 SNPs;
  `Fhet_mt_persample_merged.vcf.gz` (session 10 05d2/05e2 on slim BAMs
  with bcftools 1.23.1) = 1128 SNPs / 143 samples / ts/tv 8.17 ← **THIS
  is the active baseline 05_1 must not regress against**.

- **Mac (local)** — fast iteration:
  - SnpEff annotation (stage 06, against `~/snpEff/` standalone install,
    NCBI-built `NC_012312.1` database with Vertebrate_Mitochondrial
    codon table; Java from `~/micromamba/envs/SNP_env/lib/jvm/bin/java`):
    `05_1_Fhet_mt_persample_merged.vcf.gz` →
    `05_1_Fhet_mt_persample_merged_ann.vcf.gz`. SnpEff annotates each
    norm-split (POS, ALT) row independently.
  - CDS-restrict + SNPs-only + sample rename (stage 07,
    `scripts/07_cds_snps_norm_mac.sh`). No `norm -m -any` here —
    already done in 05_1. ANN passes through from stage 06
    transparently; no backfill or ANN-transfer dance.
  - Stage 08: haplotype calling at 0.7 binary threshold (existing).
  - Stage 09 (new, to be written): heteroplasmy report at 0.1 ≤ AF < 0.7
    with AD_alt ≥ 50; carrier-count classifier (transmitted vs private).
  - Tools: bcftools (SNP_env), cyvcf2, pandas, Java (SNP_env)

## Rules for Claude

1. **Don't rerun frozen outputs.** `Fhet_MT_CDS.snps.split.vcf.gz` is canonical *once produced AND validated*. Build downstream analysis on top of it. The frozen state requires three gates green: (a) 05_1's manifest shows `cells with DP=.: 0`, (b) 05_1's manifest shows `positions lost vs baseline: 0` against `Fhet_mt_persample_merged.vcf.gz` (the May-15 1128-SNP per-sample baseline), and (c) stage 07's summary shows per-cell DP/AD on every record AND ANN on every record. Until (a)+(b)+(c) are all green, the file is not frozen — re-running is allowed and expected. (Status 2026-05-18 session 14: caller is `jobs/05_1_mpileup_merge.sh` — joint mpileup + call + norm-split in one pipeline; 05f superseded and archived. 05_1 not yet submitted on T2. See CHANGELOG 2026-05-18 for the architecture trajectory and the structural bugs the new caller is designed to avoid.)
2. **Don't commit data files.** `.gitignore` excludes `*.vcf.gz`, `*.bcf`, `*.bam`, `*.bai`, `*.fastq*`, `*.fasta`, `*.tbi`, `*.csi`. Scripts and docs only.
3. **All BSUB scripts live under `jobs/`** with the numbered prefix (`01_..05_..`). New stages get the next number.
4. **All log files write to `/projectnb/dcrawford/MT_Genomics2/logs/`** so they're easy to find and clean up.
5. **One conda env to start: `mito_genomics`.** Split per stage only if conflicts force it; document any split in `docs/00_setup.md`.
6. **Update `CHANGELOG.md` at the end of every working session** with what changed and why.
7. **When asked to add a new pipeline step,** write the BSUB script, update `docs/01_pipeline.md`, and append a CHANGELOG entry — in that order.

## Sample / sample-list conventions

- `SSM_WGS_list.txt`: one sample id per line, e.g. `10_0`, `102_0`, `1_0`.
- Raw fastq names: `{sample}_1.fq.gz`, `{sample}_2.fq.gz`.
- Trimmed paired (actual files in `trimmed_seq/` on Triton 2): `{sample}_1_val_1.fq.gz`, `{sample}_2_val_2.fq.gz` (Trim Galore output), with companion `{sample}_{1,2}.fq.gz_trimming_report.txt`. Note: `02_trim_pe.sh` in this repo is written for Trimmomatic naming (`_p` / `_up`); see CHANGELOG 2026-05-07 — reconcile if stage 02 is ever re-run.
- MT BAMs: `{sample}_MT.bam` (+ `.bai`).
- The trailing `_0` is stripped from sample names in the joint-call VCF (`bcftools reheader`).

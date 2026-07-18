# 02 — Why per-sample-then-merge gives ~1133 SNPs and joint `bcftools call -mv` gives ~150 on this mtDNA panel

This document is the project's primary methods deliverable: a written record
of *what went wrong, why, and how to avoid it* when calling mtDNA variants
with `bcftools mpileup | bcftools call`. It exists so other investigators
working on small, divergent-reference mitogenomes don't fall into the same
~7× SNP-count trap we did.

## The observation

Three joint-call configurations of `bcftools mpileup | bcftools call`,
varying only the read filters and the ploidy model, all produced essentially
the same output across 143 *F. heteroclitus* MT BAMs:

| run | mpileup `-Q` | mpileup `-q` | call `--ploidy` | SNPs | ts/tv (1st alt) |
|----:|:------------:|:------------:|:---------------:|-----:|----------------:|
| v1  | 30           | 30           | 2 (default)     | 152  | 5.61            |
| v2  | 13           | 20           | 1               | 152  | 6.60            |
| v3  | 13           |  0           | 1               | 145  | 7.06            |

The historical pipeline that produced the canonical `merged_144.vcf.gz` on
the same BAMs (Triton 2, Jul 2025) reported **1133 SNPs across 144 samples**
with ts/tv = 7.92 — a 7× difference. The historical pipeline had been
deliberately replaced with a joint-call approach to get per-sample AD/DP for
downstream haplotype calling; that switch is what produced the gap.

## What changed between historical and current

The architecture was the dominant change. The historical pipeline was:

```
# Per sample (LSF array [1-144]%24):
bcftools mpileup -f REF SAMPLE.bam --max-depth 10000 -Ou \
  | bcftools call -mv --ploidy 1 -Oz -o SAMPLE.vcf.gz
bcftools norm -m -any -Oz -o SAMPLE_norm.vcf.gz SAMPLE.vcf.gz

# Then once across all 144:
bcftools merge -m none --threads 8 -l vcf_list.txt -Oz -o merged_144.vcf.gz
```

Reference: `archive/Notes_dlcs/Inital_call_wo_AD.txt`,
`archive/Previous_jobs/BSUB_1_MT_SNPcalls.sh`.

The current `jobs/05_bcftools_mpileup_call_AD.sh` and its v2/v3 variants do:

```
bcftools mpileup -f REF -b BAM_LIST -a AD,DP -Q 30 -q 30 -d 100000 -Ou \
  | bcftools call -mv -A -Oz -o Fhet_mt_fullAD.vcf.gz
```

Three differences matter:

1. **Per-sample call → `bcftools merge`** (historical) vs **single joint
   `bcftools call -mv`** (current). This is the dominant variable.
2. **`bcftools norm -m -any` per sample** (historical) vs none (current).
   With `merge -m none` after, the historical output stays one row per
   (POS, ALT). Without normalization, our current output is 96%
   multiallelic with phantom-alt clutter.
3. **`-A` in `bcftools call`** (current only). Keeps every alt the caller
   considered at variant sites, including ones with zero carriers across
   the panel — the source of the 327 AF=0 alts in our v1 stats.

## Why joint `-mv` suppresses the high-AF sites on this dataset

mtDNA reference `NC_012312.1` is a single individual's mitogenome. In any
moderately diverse panel of *F. heteroclitus*, the population consensus
differs from this reference at hundreds of positions. At a typical such
site, the entire panel reads alt and the reference reads ref — the site is
"fixed against reference" with AF ≈ 1.

Per-sample call sees only one BAM at a time. At each variant position in
that BAM, the per-sample likelihood crosses the threshold and the variant
is emitted to that sample's VCF. `bcftools merge` then unions across
samples; positions variant in any sample appear in the merged file. With
144 mostly-concordant haploid mitogenomes, the union accumulates ~1100
"common" sites where everyone-but-the-reference has the same alt.

Joint `bcftools call -mv` evaluates each site's likelihood across the whole
panel against an allele-frequency prior. The default prior is anchored on
the assumption that reference is the common allele and most departures
from it are rare (mutation-rate prior `-P 1.1e-3`). When 142/143 samples
are uniformly alt at a site, that signal is so cleanly inconsistent with
the prior that the joint caller's heuristics frequently suppress emission
— the model treats such a site as a likely reference annotation issue
rather than a confident variant. Drop the prior toward `-P 1.0` or split
the call into per-sample evaluations and the same site emits routinely.

The diploid default (`v1`) compounds the suppression: with no expected
heterozygotes at near-fixed sites, the joint multiallelic caller's
HWE-adjacent reasoning treats the absence of hets as additional anomaly,
not as the natural consequence of haploidy. That's why moving to
`--ploidy 1` improved per-site call quality (ts/tv climbed from 5.61 to
6.60) without changing the count: ploidy fixes the heterozygosity-modelling
issue but leaves the AF-prior issue intact. v3's `-q 0` likewise tunes
read filtering, not the prior.

The asymmetry, stated plainly: per-sample → merge is **additive** (the
union of variant sites grows with every sample); joint `-mv` is
**subtractive** (the prior penalizes each unlikely site). On nuclear
genomes with a well-matched reference these usually converge. On a 16 kb
mtDNA with a divergent reference, they diverge by ~7×.

## Evidence the 1133 is the right number

The historical clean output's allele-frequency spectrum is the giveaway.
For 144 samples, AF tiers carry these meanings:

- AF ≈ 0.00694  (1/144)   = singleton
- AF ≈ 0.97    (140/144)  = 4 reference-carrier samples, rest alt
- AF ≈ 0.99306 (143/144)  = 1 reference-carrier sample, rest alt

Distribution of the 1133 historical SNPs:

| AF tier             | SNPs |
|---------------------|-----:|
| AF ≥ 0.95           | 1004 |
| AF ∈ [0.50, 0.95)   |   70 |
| AF ∈ [0.01, 0.50)   |   54 |
| singletons (AC = 1) |   55 |

**89% of the historical SNPs are sites where the reference is the rare
allele.** That is exactly the population-against-divergent-reference
signature, and the fact that singletons are present in similar numbers in
both regimes (55 historical, 51/40/39 in v1/v2/v3) confirms our joint
caller catches the rare end fine. It's the common-fixed-against-reference
end that the joint caller drops.

The historical ts/tv across all alts is 7.92 — squarely in the expected
range for vertebrate mtDNA within-species variation. Real biology, not
noise.

## Practical checklist for mtDNA SNP calling with bcftools

1. **Use a per-sample → merge architecture, not joint `bcftools call -mv`,
   when your reference is a single individual's mitogenome.** Joint `-mv`'s
   default allele-frequency prior is calibrated for nuclear-genome
   site-frequency distributions and silently suppresses near-fixed-alt
   sites where the reference is divergent.
2. **Use `--ploidy 1`.** mtDNA is haploid; the diploid default fails the
   het-likelihood gate at heteroplasmic and low-AF sites and adds a second
   layer of suppression at near-fixed sites.
3. **Add `-a AD,DP` to `bcftools mpileup`** so per-sample allele depths
   carry through to the merged VCF. Downstream haplotype calling needs
   them; the original historical recipe didn't have them and we needed
   them, which is what motivated the (architecturally fragile) switch to
   joint mode in the first place.
4. **Do NOT pass `-A` to `bcftools call`** unless you specifically want
   every alt allele the caller considered, including zero-carrier ones,
   in the output. With multi-sample joint calling this produces ~96%
   multiallelic output dominated by phantom alts; per-sample callers don't
   need it.
5. **Normalize per sample with `bcftools norm -m -any` before merging.**
   Combined with `bcftools merge -m none`, this keeps one row per
   (POS, ALT) end-to-end; combined with downstream
   `bcftools norm -m -any -f REF` (already in stage 07), the canonical
   output is left-aligned and split.
6. **Cross-check ts/tv as the architecture sanity test.** A well-formed
   mtDNA call set against a reasonable reference should land near
   ts/tv ≈ 7–15 (across all alts). A ts/tv near 0.5 means you're looking
   at noise alts, not biology — typically a sign that `-A` is set in
   `bcftools call` and the multiallelic clutter is dragging the metric
   down.
7. **Diagnose mtDNA SNP-count gaps by comparing AF spectra, not just
   counts.** A pipeline that misses the AF ≈ 0.99 tier but keeps
   singletons is suffering from a population-prior issue (joint-call
   AF-prior, reference divergence), not a read-filtering issue. No amount
   of `-Q`/`-q`/`--ploidy` tuning will fix that. Switch architectures.

## Where to look in this repo

- Per-sample → merge implementation: `jobs/05d_persample_call.sh` and
  `jobs/05e_merge_persample.sh`.
- Joint variants kept for the methods comparison: `jobs/05_*.sh`,
  `jobs/05b_*.sh`, `jobs/05c_*.sh`.
- Side-by-side stats comparison helper:
  `scripts/compare_stage05_runs.sh` (any number of `*_stats.txt` files).
- Historical recipe (preserved): `archive/Notes_dlcs/Inital_call_wo_AD.txt`,
  `archive/Previous_jobs/BSUB_1_MT_SNPcalls.sh`.
- Historical canonical output stats:
  `stats_old/merged_stats.txt` (1133 SNPs / 144 samples / ts/tv 7.92).

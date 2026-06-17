#!/usr/bin/env python3
"""
make_141_vcf.py
Create 141-sample canonical VCF from the 143-sample canonical by:
  1. Removing samples 70 and 125
  2. Excluding 3 sites private to 70/125  (POS 4782, 4914, 12736)
  3. Excluding 3 fixed reference-divergence NS sites (POS 3124, 4680, 4957)

Output is bgzipped; run 'bcftools index --tbi' on the result afterward
(pysam can write bgzipped VCF but does not build a tabix index).

Usage (from MT_Genomics2/):
    conda activate SNP_env
    python scripts/make_141_vcf.py
    bcftools index --tbi vcf/141_MT_variants.vcf.gz
"""

import pysam
import sys

# ── Config ────────────────────────────────────────────────────────────────────
IN_VCF  = "vcf/Fhet_MT_CDS.snps.split.vcf.gz"
OUT_VCF = "vcf/141_MT_variants.vcf.gz"

EXCLUDE_SAMPLES = {"70", "125"}

EXCLUDE_POS = {3124, 4680, 4782, 4914, 4957, 12736}
# ─────────────────────────────────────────────────────────────────────────────

def main():
    vcf_in = pysam.VariantFile(IN_VCF)

    all_samples = list(vcf_in.header.samples)
    keep        = [s for s in all_samples if s not in EXCLUDE_SAMPLES]
    excluded    = [s for s in all_samples if s in EXCLUDE_SAMPLES]

    print(f"Input samples : {len(all_samples)}", file=sys.stderr)
    print(f"Excluded      : {excluded}",          file=sys.stderr)
    print(f"Kept          : {len(keep)}",          file=sys.stderr)
    print(f"Excluded POS  : {sorted(EXCLUDE_POS)}", file=sys.stderr)

    # Build output header with only the kept samples
    header_out = vcf_in.header.copy()
    for s in excluded:
        header_out.samples.remove(s)

    vcf_out = pysam.VariantFile(OUT_VCF, "wz", header=header_out)  # wz = bgzipped VCF

    n_written = n_skipped_pos = n_skipped_mono = 0

    for rec in vcf_in.fetch():

        # Skip excluded positions
        if rec.pos in EXCLUDE_POS:
            n_skipped_pos += 1
            continue

        # Build new record with only kept samples
        new_rec = header_out.new_record()
        new_rec.chrom  = rec.chrom
        new_rec.pos    = rec.pos
        new_rec.id     = rec.id
        new_rec.ref    = rec.ref
        new_rec.alts   = rec.alts
        new_rec.qual   = rec.qual
        new_rec.filter.add("PASS")

        # Copy INFO fields
        for key, val in rec.info.items():
            try:
                new_rec.info[key] = val
            except (KeyError, TypeError):
                pass

        # Copy FORMAT fields for kept samples
        for sname in keep:
            src = rec.samples[sname]
            dst = new_rec.samples[sname]
            for key in rec.format.keys():
                try:
                    dst[key] = src[key]
                except (KeyError, TypeError):
                    pass

        # Drop sites that became monomorphic after sample removal
        # (AC == 0 means no one carries the ALT in the 141-sample set)
        try:
            ac = new_rec.info["AC"]
            ac_val = ac[0] if isinstance(ac, tuple) else ac
            if ac_val == 0:
                n_skipped_mono += 1
                continue
        except KeyError:
            pass

        vcf_out.write(new_rec)
        n_written += 1

    vcf_out.close()
    vcf_in.close()

    print(f"\nRecords written       : {n_written}", file=sys.stderr)
    print(f"Skipped (excluded POS): {n_skipped_pos}", file=sys.stderr)
    print(f"Skipped (monomorphic) : {n_skipped_mono}", file=sys.stderr)
    print(f"\nOutput: {OUT_VCF}", file=sys.stderr)
    print("Run: bcftools index --tbi vcf/141_MT_variants.vcf.gz", file=sys.stderr)

if __name__ == "__main__":
    main()

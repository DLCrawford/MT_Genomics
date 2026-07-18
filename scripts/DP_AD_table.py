# DP_AD_table_April2026
import pysam

#cd ~/Projects/MT_Genomics_Cl_Ap2026/MT_Genomics2/vcf

vcf_file = "Fhet_MT_CDS.snps.split.vcf.gz"
out_file = "mtDNA_long_AD_table.tsv"

MAX_ALT = 4   # allow up to ALT4

vcf = pysam.VariantFile(vcf_file)
samples = list(vcf.header.samples)

with open(out_file, "w") as out:

    header = [
        "Individual","Position","REF",
        "ALT1","ALT2","ALT3","ALT4",
        "DP",
        "ADref","ADalt1","ADalt2","ADalt3","ADalt4"
    ]
    out.write("\t".join(header) + "\n")

    for rec in vcf:
        pos = f"{rec.contig}_{rec.pos}"

        ref = rec.ref
        alts = list(rec.alts) if rec.alts else []

        # Pad ALT list to fixed length
        alt_list = alts + ["NA"] * (MAX_ALT - len(alts))
        alt_list = alt_list[:MAX_ALT]

        for sample in samples:
            s = rec.samples[sample]

            dp = s.get("DP", "NA")
            ad = s.get("AD", [])

            # Initialize AD values
            ad_vals = ["NA"] * (1 + MAX_ALT)

            if ad:
                for i in range(min(len(ad), 1 + MAX_ALT)):
                    ad_vals[i] = ad[i]

            row = [
                sample,
                pos,
                ref,
                alt_list[0],
                alt_list[1],
                alt_list[2],
                alt_list[3],
                str(dp),
                str(ad_vals[0]),  # ADref
                str(ad_vals[1]),  # ADalt1
                str(ad_vals[2]),
                str(ad_vals[3]),
                str(ad_vals[4])
            ]

            out.write("\t".join(row) + "\n")

print("Done.")
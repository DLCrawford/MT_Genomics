#!/bin/bash
###############################################################################
# scripts/compare_stage05_runs.sh
# Side-by-side summary of multiple bcftools stats files from stage-05 runs.
#
# Usage:
#   bash scripts/compare_stage05_runs.sh stats1.txt stats2.txt [stats3.txt ...]
#
# Example (on Triton 2 once v2 + v3 finish):
#   bash scripts/compare_stage05_runs.sh \
#       vcf/Fhet_mt_fullAD_stats.txt \
#       vcf/Fhet_mt_v2_Q13_q20_p1_fullAD_stats.txt \
#       vcf/Fhet_mt_v3_Q13_q00_p1_fullAD_stats.txt
#
# Output: a fixed-width table comparing key metrics across runs.
###############################################################################

set -euo pipefail

if [[ $# -lt 1 ]]; then
    echo "Usage: $0 <stats1.txt> [stats2.txt] ..." >&2
    exit 2
fi

# Build per-file label (basename minus _stats.txt minus Fhet_mt_ prefix)
declare -a LABELS
for f in "$@"; do
    label=$(basename "$f" _stats.txt)
    label="${label#Fhet_mt_}"
    [[ -z "$label" ]] && label=$(basename "$f")
    LABELS+=("$label")
done

# Helper: extract the value for a given SN key from a stats file
sn_val() {
    local key="$1"
    local file="$2"
    awk -F'\t' -v k="$key" '$0 ~ ("^SN.*number of " k ":"){print $NF; exit}' "$file"
}

# Helper: extract a TSTV column
tstv_col() {
    local col="$1"
    local file="$2"
    awk -v c="$col" '/^TSTV/{print $c; exit}' "$file"
}

# Header row
printf "%-28s" "metric"
for l in "${LABELS[@]}"; do
    printf " %22s" "$l"
done
echo
printf "%-28s" "$(printf '%.0s-' {1..28})"
for l in "${LABELS[@]}"; do
    printf " %22s" "$(printf '%.0s-' {1..22})"
done
echo

# Metric rows
for key in "samples" "records" "SNPs" "indels" "multiallelic sites" "multiallelic SNP sites"; do
    printf "%-28s" "$key"
    for f in "$@"; do
        printf " %22s" "$(sn_val "$key" "$f")"
    done
    echo
done

# TSTV rows: column 5 = ts/tv all alts, column 8 = ts/tv 1st alt
printf "%-28s" "ts/tv (all alts)"
for f in "$@"; do
    printf " %22s" "$(tstv_col 5 "$f")"
done
echo

printf "%-28s" "ts/tv (1st alt only)"
for f in "$@"; do
    printf " %22s" "$(tstv_col 8 "$f")"
done
echo

# Singletons (AC=1)
printf "%-28s" "singleton SNPs (AC=1)"
for f in "$@"; do
    val=$(awk '/^SiS.*\t1\t/{print $4; exit}' "$f")
    printf " %22s" "${val:--}"
done
echo

echo
echo "Run manifests (parameters used) sit alongside each stats file as:"
echo "  Fhet_mt_<RUN_TAG>_run_manifest.txt"

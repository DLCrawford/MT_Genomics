#!/usr/bin/env python3
"""
scripts/15_well_bleed_test.py

Quantitative test of "is this individual's heteroplasmy load explained by
library-prep contamination from neighboring wells?" for the 4 high-Hp focal
individuals from stage 14 (default: top 4 by Hp-site count). Also accepts
an explicit list of focals via --focals.

HYPOTHESIS
----------
If well bleed is the explanation, then for a focal individual X with many
Hp events, the panel members whose haplotypes best "explain" X's minor
allele pattern should sit in physically nearby wells on the same plate.

TEST STATISTIC
--------------
For each focal X and every other panel member Y on the SAME plate, define
the donor-concordance score:

    Score(X -> Y) = |{ Hp event in X where Y's major call at POS == X's Hp_allele }|
                    / |X's Hp events|

i.e., the fraction of X's heteroplasmy events that Y's haplotype could
explain. We then test:

  (a) Spearman rho of Score vs. plate-distance(X, Y), Chebyshev metric.
      Bleed signal: negative rho (concordance falls with distance).
  (b) mean(Score | neighbor wells, dist <= 1) - mean(Score | far wells, dist >= 3).
      Bleed signal: positive difference.

PERMUTATION NULL
----------------
X's well is fixed at its observed position. We shuffle the assignment of all
OTHER samples-on-the-same-plate to the remaining wells. Recompute (a) and (b)
under the shuffled labels. Empirical p-value = fraction of permutations with
test statistic at least as extreme as the observed value.

INPUTS
------
data_files_May/WGS_seq_plate.txt   plate map (1 row per sample, Plate_label = "A1".."H12",
                                   plate identity inferred from i5 column)
vcf/heteroplasmy_pileup_events.tsv stage-14 Hp events
vcf/pileup_cds_141.vcf.gz          stage-13 per-cell DP/AD across CDS

OUTPUTS
-------
vcf/well_bleed_donor_ranking.tsv   per-(focal, candidate) Score + distance + ranks
vcf/well_bleed_results.tsv         per-focal test stats + permutation p-values
vcf/well_bleed_summary.txt         human-readable summary

USAGE
-----
python scripts/15_well_bleed_test.py
python scripts/15_well_bleed_test.py --focals 77,47,33,84
python scripts/15_well_bleed_test.py --n-perm 20000 --connectivity 4
"""

import argparse
import os
import re
import subprocess
import sys
from collections import defaultdict

import numpy as np
import pandas as pd
from scipy import stats

# ---------------------------------------------------------------------- args
ap = argparse.ArgumentParser()
ap.add_argument("--vcf-dir",  default=os.environ.get(
    "VCF_DIR",
    "/Users/douglas_crawford/Projects/MT_Genomics_Cl_Ap2026/MT_Genomics2/vcf"))
ap.add_argument("--plate-file", default=os.environ.get(
    "PLATE_FILE",
    "/Users/douglas_crawford/Projects/MT_Genomics_Cl_Ap2026/MT_Genomics2/"
    "data_files_May/WGS_seq_plate.txt"))
ap.add_argument("--events", default=None,
                help="Path to stage-14 events TSV. Default: "
                     "<vcf-dir>/heteroplasmy_pileup_events.tsv, with "
                     "automatic fallback to *_all.tsv if the bare name "
                     "doesn't exist.")
ap.add_argument("--focals", default=None,
                help="Comma-separated WGS_IDs to test; default = top-N by Hp count")
ap.add_argument("--top-n",       type=int, default=4,
                help="If --focals not given, test this many top-Hp individuals")
ap.add_argument("--n-perm",      type=int, default=10000)
ap.add_argument("--connectivity", choices=["4", "8"], default="8",
                help="Neighbor adjacency: 4 (rook) or 8 (king). Default 8.")
ap.add_argument("--major-frac", type=float, default=0.70)
args = ap.parse_args()

VCF_DIR      = args.vcf_dir
# Resolve EVENTS_TSV: --events wins, else try bare name, else fall back to *_all.tsv
if args.events:
    EVENTS_TSV = args.events
else:
    _bare = f"{VCF_DIR}/heteroplasmy_pileup_events.tsv"
    _all  = f"{VCF_DIR}/heteroplasmy_pileup_events_all.tsv"
    if os.path.exists(_bare):
        EVENTS_TSV = _bare
    elif os.path.exists(_all):
        EVENTS_TSV = _all
        print(f"[15] note: using {_all} (bare name not found)", file=sys.stderr)
    else:
        sys.exit(f"[15] ERROR: neither {_bare} nor {_all} exists. "
                 f"Run stage 14 first, or pass --events.")
PILEUP_VCF   = f"{VCF_DIR}/pileup_cds_141.vcf.gz"
OUT_RANKING  = f"{VCF_DIR}/well_bleed_donor_ranking.tsv"
OUT_RESULTS  = f"{VCF_DIR}/well_bleed_results.tsv"
OUT_SUMMARY  = f"{VCF_DIR}/well_bleed_summary.txt"

# ----------------------------------------------------- helper: normalize id
ID_RE = re.compile(r"(\d+)")
def to_wgs_id(name):
    """Extract leading integer from any sample-name variant (1, 1_0, 1_MT, etc.)."""
    m = ID_RE.match(str(name).split("/")[-1])
    return int(m.group(1)) if m else None

# --------------------------------------------------------- 1) plate map
print(f"[15] Loading plate map: {args.plate_file}", file=sys.stderr)
plate_df = pd.read_csv(args.plate_file, sep="\t")
plate_df.columns = [c.strip() for c in plate_df.columns]
plate_df["wgs_id"] = plate_df["WGS_ID"].astype(int)
# 8 rows x 12 cols. Plate_label = "<RowLetter><Col>" e.g. "H5", "B12"
plate_df["row_letter"] = plate_df["Plate_label"].str[0]
plate_df["col"]        = plate_df["Plate_label"].str[1:].astype(int)
plate_df["row"]        = plate_df["row_letter"].map(
    {l: i + 1 for i, l in enumerate("ABCDEFGH")})
plate_df["plate"]      = plate_df["i5"].str.strip()   # i5_3 == plate 1, i5_4 == plate 2

plate_lookup = plate_df.set_index("wgs_id")[["plate", "row", "col"]].to_dict("index")
print(f"      {len(plate_df)} samples loaded across plates "
      f"{sorted(plate_df['plate'].unique())}", file=sys.stderr)

# --------------------------------------------------------- 2) Hp events
print(f"[15] Loading Hp events: {EVENTS_TSV}", file=sys.stderr)
ev = pd.read_csv(EVENTS_TSV, sep="\t")
ev["wgs_id"] = ev["Individual"].map(to_wgs_id)
print(f"      {len(ev):,} events across {ev['wgs_id'].nunique()} individuals",
      file=sys.stderr)

# --------------------------------------------------------- 3) focals
if args.focals:
    focals = [int(x) for x in args.focals.split(",")]
else:
    focals = (ev.groupby("wgs_id").size()
                .sort_values(ascending=False)
                .head(args.top_n).index.tolist())
print(f"[15] Testing focals: {focals}", file=sys.stderr)

# ------------------------------- 4) major-call matrix from pileup VCF
# We only need majors for POS where some focal has Hp events.
target_pos = set(ev[ev["wgs_id"].isin(focals)]["POS"].unique())
print(f"[15] Streaming pileup VCF for majors at {len(target_pos)} POS...",
      file=sys.stderr)

# bcftools query: CHROM POS REF ALT [SAMPLE=DP=AD]
cmd = ["bcftools", "query",
       "-f", "%CHROM\t%POS\t%REF\t%ALT[\t%SAMPLE=%DP=%AD]\n",
       PILEUP_VCF]
proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True)

# major[(wgs_id, POS)] = allele letter ("A"/"C"/"G"/"T") or None
major = {}

for line in proc.stdout:
    parts = line.rstrip("\n").split("\t")
    pos = int(parts[1])
    if pos not in target_pos:
        continue
    ref, alt = parts[2], parts[3]
    alleles = [ref] + alt.split(",")
    for cell in parts[4:]:
        try:
            sname, dp_s, ad_s = cell.split("=", 2)
        except ValueError:
            continue
        wid = to_wgs_id(sname)
        if wid is None:
            continue
        try:
            dp = int(dp_s)
        except ValueError:
            continue
        if dp <= 0:
            continue
        ads = []
        for x in ad_s.split(","):
            try:
                ads.append(int(x))
            except ValueError:
                ads.append(0)
        # find the major allele (skip <*> catch-all)
        m = None
        for a, ad in zip(alleles, ads):
            if a == "<*>":
                continue
            if ad / dp >= args.major_frac:
                m = a
                break
        major[(wid, pos)] = m
proc.wait()
if proc.returncode != 0:
    sys.exit(f"bcftools query failed (rc={proc.returncode})")
print(f"      Major matrix built: {len(major):,} (sample, POS) entries",
      file=sys.stderr)

# ------------------------------- 5) compute Scores and run permutation test
chebyshev = lambda a, b: max(abs(a[0] - b[0]), abs(a[1] - b[1]))
manhattan = lambda a, b: abs(a[0] - b[0]) + abs(a[1] - b[1])
dist_fn   = chebyshev if args.connectivity == "8" else manhattan
NEIGHBOR_R = 1  # within 1 step in chosen metric
FAR_R     = 3   # at least 3 steps away

ranking_rows  = []
result_rows   = []

rng = np.random.default_rng(0)

for X in focals:
    if X not in plate_lookup:
        print(f"      [skip] focal {X} not in plate map", file=sys.stderr)
        continue
    info_X = plate_lookup[X]
    plate_X = info_X["plate"]
    pos_X   = (info_X["row"], info_X["col"])

    X_events = ev[ev["wgs_id"] == X]
    if len(X_events) == 0:
        print(f"      [skip] focal {X}: 0 Hp events", file=sys.stderr)
        continue

    # candidates: all other panel members on the SAME plate that we have
    # Hp data for (i.e., that are in the 141 BAMs, and thus in the pileup)
    panel_ids = set(ev["wgs_id"].unique())
    candidates = [y for y in panel_ids
                  if y != X
                  and y in plate_lookup
                  and plate_lookup[y]["plate"] == plate_X]

    # Score(X -> Y) and distance(X, Y)
    scores = []
    distances = []
    cand_ids  = []
    for Y in candidates:
        info_Y = plate_lookup[Y]
        pos_Y  = (info_Y["row"], info_Y["col"])
        d      = dist_fn(pos_X, pos_Y)
        # concordance
        c = 0
        for _, row in X_events.iterrows():
            m = major.get((Y, row["POS"]))
            if m is not None and m == row["Hp_allele"]:
                c += 1
        s = c / len(X_events)
        scores.append(s)
        distances.append(d)
        cand_ids.append(Y)
        ranking_rows.append(dict(
            focal=X, focal_well=f"{'ABCDEFGH'[info_X['row']-1]}{info_X['col']}",
            focal_plate=plate_X, candidate=Y,
            cand_well=f"{'ABCDEFGH'[info_Y['row']-1]}{info_Y['col']}",
            dist=d, n_hp=len(X_events), concordant=c, score=round(s, 4),
        ))

    scores    = np.asarray(scores)
    distances = np.asarray(distances)

    if scores.std() == 0:
        rho, rho_p_param = float("nan"), float("nan")
    else:
        rho, rho_p_param = stats.spearmanr(distances, scores)

    neigh_mask = distances <= NEIGHBOR_R
    far_mask   = distances >= FAR_R
    obs_diff = (scores[neigh_mask].mean() if neigh_mask.any() else np.nan) \
             - (scores[far_mask].mean()   if far_mask.any()   else np.nan)
    n_neigh = int(neigh_mask.sum())
    n_far   = int(far_mask.sum())

    # ---------- permutation: shuffle other-sample well assignments
    # All wells on plate_X EXCEPT X's well are the "available" wells.
    # (Kept for reference; the permutation below uses cand_wells.)
    same_plate_wells = [(info["row"], info["col"])
                        for s, info in plate_lookup.items()
                        if info["plate"] == plate_X and s != X]
    # We have len(candidates) candidates and possibly fewer wells if some
    # plate wells have no panel sample (e.g., the failed 70 well). We need
    # one well per candidate. Use the wells of the candidates themselves
    # (i.e., shuffle which candidate sits in which candidate-well).
    cand_wells = [(plate_lookup[y]["row"], plate_lookup[y]["col"]) for y in cand_ids]

    rho_null   = np.empty(args.n_perm)
    diff_null  = np.empty(args.n_perm)
    cand_idx = np.arange(len(cand_ids))
    for k in range(args.n_perm):
        perm = rng.permutation(cand_idx)
        # under permutation: sample-id i gets well of cand_wells[perm[i]]
        dperm = np.array([dist_fn(pos_X, cand_wells[p]) for p in perm])
        if scores.std() > 0:
            rho_null[k] = stats.spearmanr(dperm, scores).statistic
        else:
            rho_null[k] = 0.0
        nm = dperm <= NEIGHBOR_R
        fm = dperm >= FAR_R
        diff_null[k] = ((scores[nm].mean() if nm.any() else np.nan) -
                        (scores[fm].mean() if fm.any() else np.nan))

    # one-sided: bleed -> rho < 0 (closer = higher concordance)
    rho_p   = float((rho_null   <= rho).mean()) if np.isfinite(rho)      else float("nan")
    diff_p  = float((diff_null  >= obs_diff).mean()) if np.isfinite(obs_diff) else float("nan")

    # top donor info
    top_idx = int(np.argmax(scores))
    top_y   = cand_ids[top_idx]
    top_dist = int(distances[top_idx])
    top_score = float(scores[top_idx])
    top_well = f"{'ABCDEFGH'[plate_lookup[top_y]['row']-1]}{plate_lookup[top_y]['col']}"

    result_rows.append(dict(
        focal=X,
        focal_well=f"{'ABCDEFGH'[info_X['row']-1]}{info_X['col']}",
        focal_plate=plate_X,
        n_hp=len(X_events),
        n_candidates=len(candidates),
        n_neighbors=n_neigh,
        n_far=n_far,
        top_donor=top_y, top_donor_well=top_well,
        top_donor_dist=top_dist, top_donor_score=round(top_score, 4),
        spearman_rho_dist_vs_score=round(rho, 4) if np.isfinite(rho) else "NA",
        spearman_perm_p=round(rho_p, 4) if np.isfinite(rho_p) else "NA",
        mean_score_neighbor=round(float(scores[neigh_mask].mean()), 4) if n_neigh else "NA",
        mean_score_far=round(float(scores[far_mask].mean()), 4) if n_far else "NA",
        neighbor_minus_far=round(obs_diff, 4) if np.isfinite(obs_diff) else "NA",
        neighbor_perm_p=round(diff_p, 4) if np.isfinite(diff_p) else "NA",
    ))

# ------------------------------- 6) write outputs
rk = pd.DataFrame(ranking_rows).sort_values(
    ["focal", "score"], ascending=[True, False])
rk.to_csv(OUT_RANKING, sep="\t", index=False)
res = pd.DataFrame(result_rows)
res.to_csv(OUT_RESULTS, sep="\t", index=False)

with open(OUT_SUMMARY, "w") as f:
    f.write(
        f"Well-bleed test — stage 15\n"
        f"Plate map  : {args.plate_file}\n"
        f"Events     : {EVENTS_TSV}\n"
        f"Pileup VCF : {PILEUP_VCF}\n"
        f"Connectivity: {args.connectivity} (neighbor dist <= {NEIGHBOR_R}, far dist >= {FAR_R})\n"
        f"Permutations: {args.n_perm}\n\n"
        f"{res.to_string(index=False)}\n\n"
        f"Top 5 candidate donors per focal:\n"
    )
    for X in focals:
        sub = rk[rk["focal"] == X].head(5)
        f.write(f"\n--- focal {X} (well {sub['focal_well'].iloc[0]}, "
                f"plate {sub['focal_plate'].iloc[0]}) ---\n")
        f.write(sub[["candidate", "cand_well", "dist", "concordant",
                     "n_hp", "score"]].to_string(index=False))
        f.write("\n")

print(open(OUT_SUMMARY).read())
print(f"\n[15] Wrote:\n  {OUT_RANKING}\n  {OUT_RESULTS}\n  {OUT_SUMMARY}",
      file=sys.stderr)

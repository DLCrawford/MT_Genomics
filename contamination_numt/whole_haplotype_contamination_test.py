#!/usr/bin/env python3
"""
whole_haplotype_contamination_test.py
--------------------------------------
Whole-haplotype contamination test for heteroplasmic individuals.

Single-donor cross-sample contamination transfers that donor's ENTIRE linked mt
haplotype, at one fraction. So for a heteroplasmic individual X, contamination is
plausible only if some other individual Y's MAJOR (consensus) haplotype carries
ALL of X's minor variants and matches X at the remaining sites. Per X we report:
  n_compatible_donors : panel members whose major base equals X's minor variant at
                        EVERY heteroplasmic position of X (i.e. Y's haplotype could
                        supply all of X's minority variants at once).
  closest_donor       : of those, the one differing from X at the fewest OTHER
                        variant positions (extra_sites) — a perfect donor has 0.
  verdict             : EXCLUDED if no compatible donor, or if the closest
                        compatible donor still differs at >0 other sites; a clean
                        whole-haplotype donor (extra_sites == 0) is contamination-
                        consistent (and then its plate is reported).
The four extreme-load focal individuals (33, 47, 77, 84) are tested separately by
the adjacent-well direct test (adjacent_well_direct_test.py, which includes all Hp fish) and are excluded here.

Inputs : vcf/heteroplasmy_pileup_events_all.tsv  (stage-14 minor variants per fish)
         vcf/haplotype_matrix.csv                 (per-individual major/consensus calls)
         data_files_May/WGS_seq_plate.txt         (i5 index, used as 'plate')
Output : Hp_wholehaplotype_test_reconstructed.tsv
Note   : a donor must carry EVERY minor variant X shows; where a position carries
         two minor variants (multiallelic heteroplasmy) both are required, which is
         the conservative choice. n_Hp counts heteroplasmic positions (sites).
"""
import csv, re, os
BASE = os.path.expanduser("~/mnt/MT_Genomics2")
EVENTS = f"{BASE}/vcf/heteroplasmy_pileup_events_all.tsv"
MATRIX = f"{BASE}/vcf/haplotype_matrix.csv"
PLATE  = f"{BASE}/data_files_May/WGS_seq_plate.txt"
FOCALS = {"33","47","77","84"}
IDRE = re.compile(r"(\d+)")
def cid(s): return IDRE.search(s.split("/")[-1]).group(1)

plate = {}
with open(PLATE, encoding="latin-1") as f:
    for row in csv.DictReader(f, delimiter="\t"):
        col = [k for k in row if k.strip() == "i5"][0]
        plate[row["WGS_ID"].strip()] = row[col].strip()

M = {}; href = {}
with open(EVENTS) as f:
    for row in csv.DictReader(f, delimiter="\t"):
        x = cid(row["Individual"])
        M.setdefault(x, set()).add((row["POS"], row["Hp_allele"]))
        href[row["POS"]] = row["REF"]

rows = list(csv.DictReader(open(MATRIX)))
inds = [k[:-3] for k in rows[0] if k.endswith("_MT")]
positions = []; pos_ref = {}; ind_alt = {}; seen = set()
for r in rows:
    p = r["POS"]
    if p not in seen:
        seen.add(p); positions.append(p); pos_ref[p] = r["REF"]
    for ind in inds:
        if r[ind+"_MT"] == "1":
            ind_alt[(ind, p)] = r["ALT"]
posidx = {p:i for i,p in enumerate(positions)}
V = {ind: [ind_alt.get((ind,p), pos_ref[p]) for p in positions] for ind in inds}
def major_base(ind, pos):
    if pos in posidx: return ind_alt.get((ind,pos), pos_ref[pos])
    return href.get(pos)                      # non-variant position: reference for all

panel = set(inds)
test_inds = sorted(set(M) - FOCALS, key=lambda x: int(x))
out = [["Individual","plate","n_Hp","n_compatible_donors","closest_donor",
        "donor_plate","donor_same_diff","extra_sites","verdict"]]
for X in test_inds:
    mx = M[X]; mx_pos = set(p for p,_ in mx)
    n_hp = len(mx_pos)                         # heteroplasmic positions (sites)
    skip = set(posidx[p] for p in mx_pos if p in posidx)
    comp = [Y for Y in panel if Y != X and all(major_base(Y,p) == a for p,a in mx)]
    pX = plate.get(X, "?")
    if not comp:
        out.append([X, pX, n_hp, 0, "-","-","-", 0,
                    "EXCLUDED: no donor carries this fish's minor allele(s)"]); continue
    vx = V[X]; best = None
    for Y in sorted(comp, key=lambda y: int(y)):
        vy = V[Y]
        extra = sum(1 for i in range(len(positions)) if i not in skip and vx[i] != vy[i])
        if best is None or extra < best[0]:
            best = (extra, Y)
    extra, Y = best
    pdon = plate.get(Y, "?"); sd = "same" if pdon == pX else "diff"
    verd = ("consistent donor exists (fish %s, %s plate)" % (Y, sd)) if extra == 0 \
           else ("EXCLUDED: closest donor also differs at %d other site(s)" % extra)
    out.append([X, pX, n_hp, len(comp), Y, pdon, sd, extra, verd])

with open("Hp_wholehaplotype_test_reconstructed.tsv","w") as o:
    for r in out: o.write("\t".join(str(x) for x in r)+"\n")
# summary
import collections
verds = collections.Counter(("consistent" if r[8].startswith("consistent") else "EXCLUDED") for r in out[1:])
print("wrote Hp_wholehaplotype_test_reconstructed.tsv (%d individuals): %s"%(len(out)-1, dict(verds)))
for r in out[1:]:
    if r[8].startswith("consistent"): print("  consistent-donor case:", r[0], "->", r[4], r[6], "plate")

import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import pandas as pd
from src.blocklist_io import load_reference_map

ref = load_reference_map()

def truth_map(path):
    m = {}
    with open(path, errors="ignore") as f:
        for line in f:
            if line.startswith("#"): continue
            p = line.rstrip("\n").split("\t")
            if len(p) >= 4:
                m[p[0]] = p[3].strip().upper()
    return m

gt = truth_map("data/genome_James_Jones_v5_Full_20230726173828.txt")
san = pd.read_csv("results/output.txt", sep="\t", dtype=str,
                  keep_default_na=False, na_filter=False)
san_gt = dict(zip(san["rsid"], san["genotype"]))

MISSING = {"--", "", "00", "XX"}
carriers = []
for rsid, refgt in ref.items():
    t = gt.get(rsid)
    if t is None or t in MISSING: continue
    if t != refgt:                       # this person is a carrier here
        released = san_gt.get(rsid, "")
        carriers.append((rsid, refgt, t, released, released == t))

print(f"true carriers: {len(carriers)}")
still_exposed = sum(1 for c in carriers if c[4])
print(f"carrier genotype still exposed after sanitize: {still_exposed}")
print(f"carrier concealment rate: {(1 - still_exposed/len(carriers))*100:.1f}%" if carriers else "n/a")
print("\nsample (rsid, ref, truth, released, still_exposed):")
for c in carriers[:15]:
    print("  ", c)
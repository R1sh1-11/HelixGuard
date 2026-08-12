import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import pandas as pd
from src.blocklist_io import load_reference_map

ref = load_reference_map()
MISSING = {"--", "", "00", "XX"}


def truth_map(path):
    m = {}
    with open(path, errors="ignore") as f:
        for line in f:
            if line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) >= 4:
                m[p[0]] = p[3].strip().upper()
    return m


rows = []
genomes = sorted(f for f in os.listdir("data")
                 if f.startswith("genome_") and f.endswith(".txt"))

for g in genomes:
    san_path = f"results/sanitized_{g}"
    if not os.path.exists(san_path):
        print(f"[skip] no sanitized output for {g}")
        continue
    gt = truth_map(f"data/{g}")
    san = pd.read_csv(san_path, sep="\t", dtype=str,
                      keep_default_na=False, na_filter=False)
    san_gt = dict(zip(san["rsid"], san["genotype"]))

    carriers = 0
    exposed = 0
    for rsid, refgt in ref.items():
        t = gt.get(rsid)
        if t is None or t in MISSING:
            continue
        if t != refgt:
            carriers += 1
            if san_gt.get(rsid, "") == t:
                exposed += 1

    conceal = round((1 - exposed / carriers) * 100, 1) if carriers else None
    rows.append({
        "Genome": g.replace("genome_", "").split("_v")[0].split("_")[0],
        "Carriers": carriers,
        "Exposed": exposed,
        "Concealment_%": conceal,
    })

df = pd.DataFrame(rows)
df.to_csv("results/carrier_concealment.csv", index=False)
print(df.to_string(index=False))
total_c = df["Carriers"].sum()
total_e = df["Exposed"].sum()
print(f"\nCohort: {total_c} carriers, {total_e} exposed, "
      f"{(1-total_e/total_c)*100:.1f}% concealment" if total_c else "no carriers")
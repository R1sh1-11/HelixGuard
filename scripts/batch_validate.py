import os
import sys
import time
import subprocess

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import pandas as pd
from src.blocklist_io import load_reference_map

DATA_DIR = "data"
RESULTS_DIR = "results"
GT_MISSING = {"--", "", "00", "XX"}


def load_tsv(path, has_header):
    if has_header:
        return pd.read_csv(path, sep="\t", low_memory=False, dtype=str,
                           keep_default_na=False, na_filter=False)
    return pd.read_csv(path, sep="\t", comment="#", header=None,
                       names=["rsid", "chrom", "pos", "genotype"],
                       low_memory=False, dtype=str,
                       keep_default_na=False, na_filter=False)


def main():
    os.makedirs(RESULTS_DIR, exist_ok=True)
    ref = load_reference_map()

    genomes = sorted(f for f in os.listdir(DATA_DIR)
                     if f.startswith("genome_") and f.endswith(".txt"))
    print(f"[BATCH] {len(genomes)} genome(s) found.")

    rows = []
    for gfile in genomes:
        input_path = os.path.join(DATA_DIR, gfile)
        out_hg = os.path.join(RESULTS_DIR, f"sanitized_{gfile}")
        out_naive = os.path.join(RESULTS_DIR, f"naive_{gfile}")

        start = time.time()
        subprocess.run(
            f"python -m helixguard {input_path} {out_hg} --epsilon 1.0",
            shell=True, capture_output=True, text=True)
        subprocess.run(
            f"python -m helixguard {input_path} {out_naive} --no-dp",
            shell=True, capture_output=True, text=True)
        elapsed = round(time.time() - start, 2)

        if not os.path.exists(out_hg):
            print(f"[BATCH] FAILED to produce output for {gfile}")
            continue

        orig = load_tsv(input_path, has_header=False)
        orig["genotype"] = orig["genotype"].fillna("--").str.strip().str.upper()
        san = load_tsv(out_hg, has_header=True)

        # SNPs sanitized: genotype changed between original and HelixGuard output
        merged = orig.merge(san, on="rsid", suffixes=("_orig", "_san"))
        changed = int((merged["genotype_orig"] != merged["genotype_san"]).sum())

        # LD neighbors perturbed: tagged neighbors whose value differs
        ld_rows = merged[merged["ld_neighbor"].astype(str) == "True"]
        ld_perturbed = int((ld_rows["genotype_orig"] != ld_rows["genotype_san"]).sum())
        ld_tagged = int((san["ld_neighbor"].astype(str) == "True").sum())

        # Residual disclosure rate on the HelixGuard output
        gt = dict(zip(orig["rsid"], orig["genotype"]))
        targets = set(ref.keys()) | set(san.loc[san["ld_neighbor"].astype(str) == "True", "rsid"])
        tdf = san[san["rsid"].isin(targets)].copy()
        tdf["truth"] = tdf["rsid"].map(gt)
        tdf = tdf[tdf["truth"].notna()]
        disclosed = int((tdf["genotype"] == tdf["truth"]).sum())
        disclosure_pct = round(disclosed / len(tdf) * 100, 2) if len(tdf) else 0.0

        # Utility: can the snps-style reader re-parse the output (4 columns intact)
        utility_pass = "PASS" if set(["rsid", "genotype"]).issubset(san.columns) else "FAIL"

        rows.append({
            "Genome": gfile.replace("genome_", "").split("_v")[0],
            "Total_SNPs": len(orig),
            "SNPs_Sanitized": changed,
            "LD_Tagged": ld_tagged,
            "LD_Perturbed": ld_perturbed,
            "Disclosure_Rate_%": disclosure_pct,
            "Utility": utility_pass,
            "Time_s": elapsed,
        })
        print(f"  {gfile[:40]:40}  changed={changed:5}  ld={ld_perturbed}  "
              f"disc={disclosure_pct}%  {elapsed}s")

    summary = pd.DataFrame(rows)
    out_csv = os.path.join(RESULTS_DIR, "multi_genome_metrics.csv")
    summary.to_csv(out_csv, index=False)
    print(f"\n[BATCH] Saved {out_csv}")
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
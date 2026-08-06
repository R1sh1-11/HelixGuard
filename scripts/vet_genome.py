import argparse
import os
import re
import sqlite3
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

MISSING = {"--", "", "00", "XX"}


def chip_version(path):
    m = re.search(r"_v(\d)_", os.path.basename(path))
    if m:
        return f"v{m.group(1)}"
    with open(path, "r", errors="ignore") as fh:
        for line in fh:
            if not line.startswith("#"):
                break
            if "array" in line.lower() or "version" in line.lower():
                return line.strip("# \n")[:60]
    return "unknown"


def main():
    p = argparse.ArgumentParser()
    p.add_argument("path")
    args = p.parse_args()

    from src.parser import parse_genome

    conn = sqlite3.connect("data/blocklist.db")
    ref = {f"rs{r[0]}": str(r[1]).strip().upper()
           for r in conn.execute('SELECT "RS# (dbSNP)", ReferenceAllele FROM blocklist')
           if r[1]}
    ld_targets = {r[0] for r in conn.execute("SELECT DISTINCT target_rsid FROM ld_neighbors")}
    ld_neighbors = {r[0] for r in conn.execute("SELECT DISTINCT neighbor_rsid FROM ld_neighbors")}
    conn.close()

    df = parse_genome(args.path)
    df["genotype"] = df["genotype"].fillna("--").astype(str).str.strip().str.upper()

    total = len(df)
    nocall = int(df["genotype"].isin(MISSING).sum())

    hit = df[df["rsid"].isin(ref)].copy()
    hit["ref"] = hit["rsid"].map(ref)
    usable = hit[~hit["genotype"].isin(MISSING)]
    already = int((usable["genotype"] == usable["ref"]).sum())
    changed = int((usable["genotype"] != usable["ref"]).sum())

    neighbors_present = int(df["rsid"].isin(ld_neighbors).sum())
    ld_covered = int(hit["rsid"].isin(ld_targets).sum())

    ver = chip_version(args.path)

    print(f"\n{'='*58}")
    print(f"  File              : {os.path.basename(args.path)}")
    print(f"  Chip version      : {ver}")
    print(f"  Total SNPs        : {total:,}")
    print(f"  No-call rows      : {nocall:,} ({nocall/total*100:.2f}%)")
    print(f"  Blocklist hits    : {len(hit):,}")
    print(f"    already ref     : {already:,}")
    print(f"    would change    : {changed:,}")
    print(f"  Hits w/ LD data   : {ld_covered:,}")
    print(f"  LD neighbors here : {neighbors_present:,}")
    print(f"{'='*58}")
    print("\nSources.md row:")
    print(f"| Name | {ver} | {total:,} | {len(hit):,} | {changed:,} | [link](URL) |")


if __name__ == "__main__":
    main()
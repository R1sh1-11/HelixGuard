"""
Baseline disclosure measurement for HelixGuard output.

This module measures RESIDUAL DISCLOSURE: how often a targeted SNP's released
genotype still equals the individual's true genotype.

Important scope note: this is not an adversarial inference attack. No linkage
disequilibrium, allele frequency, or population structure information is used.
Because a target counts as disclosed whenever its value was left unchanged,
this metric is arithmetically equivalent to (targets - genotypes modified),
and it cannot distinguish a strong defense from a weak one. It exists as a
sanity baseline only.

The actual LD-based inference attack lives in src/redteam_ld.py.
"""

import argparse
import os
import sqlite3

import pandas as pd


# ------------------------------------------------------------------
# Load blocklist rsIDs from the database
# ------------------------------------------------------------------
def load_blocklist_rsids() -> set:
    conn = sqlite3.connect("data/blocklist.db")
    cursor = conn.cursor()
    cursor.execute('SELECT "RS# (dbSNP)" FROM blocklist')
    rsids = {f"rs{row[0]}" for row in cursor.fetchall()}
    conn.close()
    return rsids


# ------------------------------------------------------------------
# Parse a sanitized output TSV
# (header: rsid chrom pos genotype ld_neighbor)
# ------------------------------------------------------------------
def load_output(path: str) -> pd.DataFrame:
    return pd.read_csv(
        path, sep="\t", low_memory=False,
        dtype={"rsid": str, "genotype": str},
        keep_default_na=False, na_filter=False,
    )


# ------------------------------------------------------------------
# Parse an original 23andMe genome file (no header, # comments)
# ------------------------------------------------------------------
def load_ground_truth(path: str) -> pd.DataFrame:
    rows = []
    with open(path, "r", errors="ignore") as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                continue
            rows.append({"rsid": parts[0], "true_genotype": parts[3].strip().upper()})
    return pd.DataFrame(rows).set_index("rsid")


# ------------------------------------------------------------------
# Residual disclosure measurement
# ------------------------------------------------------------------
def attack(sanitized: pd.DataFrame, ground_truth: pd.DataFrame,
           target_rsids: set) -> dict:
    """
    For each targeted SNP present in both files, check whether the released
    genotype still matches the truth.

    Vectorized: iterating 600k+ rows per call made repeated-trial sweeps
    unusably slow.
    """
    df = sanitized[sanitized["rsid"].isin(target_rsids)].merge(
        ground_truth, left_on="rsid", right_index=True, how="inner"
    )

    if df.empty:
        empty = pd.DataFrame(
            columns=["rsid", "true_genotype", "sanitized_genotype",
                     "ld_neighbor", "disclosed"]
        )
        return {"total": 0, "correct": 0, "success_rate_pct": 0.0, "df": empty}

    df = df.rename(columns={"genotype": "sanitized_genotype"})
    df["disclosed"] = df["sanitized_genotype"] == df["true_genotype"]

    if "ld_neighbor" not in df.columns:
        df["ld_neighbor"] = False

    out = df[["rsid", "true_genotype", "sanitized_genotype",
              "ld_neighbor", "disclosed"]].copy()

    disclosed = int(out["disclosed"].sum())
    total = len(out)

    return {
        "total": total,
        "correct": disclosed,
        "success_rate_pct": round(disclosed / total * 100, 2),
        "df": out,
    }


def print_result(label: str, path: str, res: dict):
    print(f"\n{'=' * 58}")
    print(f"  Mode              : {label}")
    print(f"  File              : {path}")
    print(f"  SNPs targeted     : {res['total']}")
    print(f"  Still disclosed   : {res['correct']}")
    print(f"  RESIDUAL DISCLOSURE: {res['success_rate_pct']}%")
    print(f"{'=' * 58}")


def save_csv(res: dict, path: str):
    directory = os.path.dirname(path)
    if directory:
        os.makedirs(directory, exist_ok=True)
    res["df"].to_csv(path, index=False)
    print(f"[INFO] Saved to {path}")


# ------------------------------------------------------------------
# CLI
# ------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(
        description="Measure residual disclosure in sanitized genome output."
    )
    parser.add_argument("--mode",
                        choices=["helixguard", "naive", "compare"],
                        default="helixguard")
    parser.add_argument("--ground-truth", required=True)
    parser.add_argument("--helixguard-file", default="results/output.txt")
    parser.add_argument("--naive-file", default="results/output_naive.txt")
    parser.add_argument("--save-csv", action="store_true")
    args = parser.parse_args()

    gt = load_ground_truth(args.ground_truth)
    blocklist = load_blocklist_rsids()

    res = res_naive = None

    if args.mode in ("helixguard", "compare"):
        if not os.path.exists(args.helixguard_file):
            print(f"[SKIP] {args.helixguard_file} not found")
        else:
            df = load_output(args.helixguard_file)
            target = blocklist | set(
                df.loc[df["ld_neighbor"].astype(str) == "True", "rsid"])
            res = attack(df, gt, target)
            print_result("HELIXGUARD (blocklist + DP)", args.helixguard_file, res)
            if args.save_csv:
                save_csv(res, "results/redteam_helixguard.csv")

    if args.mode in ("naive", "compare"):
        if not os.path.exists(args.naive_file):
            print(f"\n[SKIP] {args.naive_file} not found -- "
                  f"run helixguard with --no-dp first")
        else:
            df = load_output(args.naive_file)
            target = blocklist | set(
                df.loc[df["ld_neighbor"].astype(str) == "True", "rsid"])
            res_naive = attack(df, gt, target)
            print_result("NAIVE (blocklist only, no DP)", args.naive_file, res_naive)
            if args.save_csv:
                save_csv(res_naive, "results/redteam_naive.csv")

    if args.mode == "compare" and (res or res_naive):
        print("\n  SUMMARY")
        print(f"  {'Mode':<32}{'Residual Disclosure %'}")
        print(f"  {'-' * 54}")
        if res:
            print(f"  {'HelixGuard (blocklist + DP)':<32}{res['success_rate_pct']}%")
        if res_naive:
            print(f"  {'Naive (blocklist only)':<32}{res_naive['success_rate_pct']}%")
        print(f"  {'Target threshold':<32}<=15%")
        print("\n  Note: residual disclosure is a baseline measure, not an")
        print("  adversarial attack. See src/redteam_ld.py for LD inference.")


if __name__ == "__main__":
    main()
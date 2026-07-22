import argparse
import os
import pandas as pd
import sqlite3


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
# Parse a sanitized output TSV (has header: rsid chrom pos genotype ld_neighbor)
# ------------------------------------------------------------------
def load_output(path: str) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", low_memory=False)


# ------------------------------------------------------------------
# Parse an original 23andMe genome file (no header, # comments)
# ------------------------------------------------------------------
def load_ground_truth(path: str) -> pd.DataFrame:
    rows = []
    with open(path, "r") as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.strip().split("\t")
            if len(parts) < 4:
                continue
            rows.append({"rsid": parts[0], "true_genotype": parts[3]})
    return pd.DataFrame(rows).set_index("rsid")


# ------------------------------------------------------------------
# Attack logic
# ------------------------------------------------------------------
def attack(sanitized: pd.DataFrame, ground_truth: pd.DataFrame, target_rsids: set) -> dict:
    """
    For each targeted SNP, compare the sanitized genotype against ground truth.
    In both helixguard and naive modes the output file already has a genotype
    string -- the attacker just reads it directly. Success = guess matches truth.
    """
    results = []

    for _, row in sanitized.iterrows():
        rsid = row["rsid"]
        if rsid not in target_rsids:
            continue
        if rsid not in ground_truth.index:
            continue

        true_gt = ground_truth.loc[rsid, "true_genotype"]
        guessed_gt = row["genotype"]
        hit = guessed_gt == true_gt

        results.append({
            "rsid": rsid,
            "true_genotype": true_gt,
            "sanitized_genotype": guessed_gt,
            "ld_neighbor": row.get("ld_neighbor", False),
            "correct": hit,
        })

    df = pd.DataFrame(results)
    if df.empty:
        return {"total": 0, "correct": 0, "success_rate_pct": 0.0, "df": df}

    correct = df["correct"].sum()
    total = len(df)
    return {
        "total": total,
        "correct": correct,
        "success_rate_pct": round(correct / total * 100, 2),
        "df": df,
    }


def print_result(label: str, path: str, res: dict):
    print(f"\n{'='*55}")
    print(f"  Mode          : {label}")
    print(f"  File          : {path}")
    print(f"  SNPs targeted : {res['total']}")
    print(f"  Correct infers: {res['correct']}")
    print(f"  SUCCESS RATE  : {res['success_rate_pct']}%")
    print(f"{'='*55}")


def save_csv(res: dict, path: str):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    res["df"].to_csv(path, index=False)
    print(f"[INFO] Saved to {path}")


# ------------------------------------------------------------------
# CLI
# ------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--mode", choices=["helixguard", "naive", "compare"], default="helixguard")
    parser.add_argument("--ground-truth", required=True)
    parser.add_argument("--save-csv", action="store_true")
    args = parser.parse_args()

    gt = load_ground_truth(args.ground_truth)
    blocklist = load_blocklist_rsids()

    if args.mode in ("helixguard", "compare"):
        path = "results/output.txt"
        df = load_output(path)
        # Target: blocklisted SNPs + LD neighbors
        target = blocklist | set(df[df["ld_neighbor"] == True]["rsid"].tolist())
        res = attack(df, gt, target)
        print_result("HELIXGUARD (DP noise)", path, res)
        if args.save_csv:
            save_csv(res, "results/redteam_helixguard.csv")

    if args.mode in ("naive", "compare"):
        path = "results/output_naive.txt"
        if not os.path.exists(path):
            print(f"\n[SKIP] {path} not found -- run generate_naive_baseline.py first")
        else:
            df = load_output(path)
            target = blocklist | set(df[df["ld_neighbor"] == True]["rsid"].tolist())
            res_naive = attack(df, gt, target)
            print_result("NAIVE REDACTION (baseline)", path, res_naive)
            if args.save_csv:
                save_csv(res_naive, "results/redteam_naive.csv")

    if args.mode == "compare":
        print("\n  SUMMARY")
        print(f"  {'Mode':<30} {'Attack Success %'}")
        print(f"  {'-'*45}")
        print(f"  {'HelixGuard (DP noise)':<30} {res['success_rate_pct']}%")
        if os.path.exists("results/output_naive.txt"):
            print(f"  {'Naive redaction':<30} {res_naive['success_rate_pct']}%")
        print(f"  {'Target threshold':<30} ≤15%")


if __name__ == "__main__":
    main()
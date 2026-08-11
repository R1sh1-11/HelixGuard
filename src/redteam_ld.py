"""
LD-based red-team attacker.

Threat model: the attacker sees ONLY the sanitized file. For each SNP that
HelixGuard may have redacted, the attacker tries to infer whether the
individual is a CARRIER (true genotype differs from population reference)
using the genotypes of SNPs in linkage disequilibrium with that target.

Carrier status is the privacy-relevant quantity: it is what an insurer or
employer would discriminate on, and it is what the pipeline claims to hide.

Baseline comparison: an attacker with no LD information should just guess
"not a carrier" every time, since most individuals are reference at most
pathogenic sites. The LD attacker is only meaningful if it beats that.

Approximation: neighbor major alleles are estimated from the cohort of
genome files in data/ rather than an external allele-frequency panel.
This is a documented limitation; accuracy improves as the cohort grows.
"""

import argparse
import os
import sqlite3
from collections import defaultdict, Counter

import pandas as pd

DATA_DIR = "data"
MISSING = {"--", "", "00", "XX", "NAN", "NONE"}


def load_blocklist_ref():
    conn = sqlite3.connect("data/blocklist.db")
    ref = {
        r[0]: str(r[1]).strip().upper()
        for r in conn.execute("SELECT rsid, ReferenceAllele FROM blocklist")
        if r[1]
    }
    conn.close()
    return ref


def load_ld_table(r2_min):
    conn = sqlite3.connect("data/blocklist.db")
    ld = defaultdict(list)
    q = "SELECT target_rsid, neighbor_rsid, r2 FROM ld_neighbors WHERE r2 >= ?"
    for target, neighbor, r2 in conn.execute(q, (r2_min,)):
        ld[target].append((neighbor, float(r2)))
    conn.close()
    return ld


def estimate_major_alleles(rsids):
    """Estimate the most common allele per rsID across the local genome cohort."""
    counts = defaultdict(Counter)
    files = [f for f in os.listdir(DATA_DIR)
             if f.startswith("genome_") and f.endswith(".txt")]
    for fname in files:
        with open(os.path.join(DATA_DIR, fname), "r") as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 4 or parts[0] not in rsids:
                    continue
                for allele in parts[3].strip().upper():
                    if allele in "ACGT":
                        counts[parts[0]][allele] += 1
    print(f"[LD-Attack] Estimated major alleles for {len(counts)} neighbors "
          f"from {len(files)} genome(s).")
    return {r: c.most_common(1)[0][0] for r, c in counts.items() if c}


def dosage(genotype, major):
    """Count non-major alleles: 0, 1, or 2. None if unusable."""
    gt = str(genotype).strip().upper()
    if gt in MISSING or len(gt) != 2 or not all(a in "ACGT" for a in gt):
        return None
    return sum(1 for a in gt if a != major)


def load_ground_truth(path):
    truth = {}
    with open(path, "r") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 4:
                truth[parts[0]] = parts[3].strip().upper()
    return truth


def run_attack(sanitized_path, ground_truth_path, r2_min, cutoff):
    san = pd.read_csv(sanitized_path, sep="\t", low_memory=False,
                      dtype=str, keep_default_na=False, na_filter=False)
    san_gt = dict(zip(san["rsid"], san["genotype"]))

    truth = load_ground_truth(ground_truth_path)
    ref = load_blocklist_ref()
    ld = load_ld_table(r2_min)

    all_neighbors = {n for pairs in ld.values() for n, _ in pairs}
    majors = estimate_major_alleles(all_neighbors)

    rows = []
    for target, ref_gt in ref.items():
        if target not in truth or target not in san_gt:
            continue
        true_gt = truth[target]
        if true_gt in MISSING:
            continue

        truth_carrier = true_gt != ref_gt

        num = den = 0.0
        used = 0
        for neighbor, r2 in ld.get(target, []):
            if neighbor not in san_gt:
                continue
            major = majors.get(neighbor)
            if major is None:
                continue
            d = dosage(san_gt[neighbor], major)
            if d is None:
                continue
            num += r2 * (d / 2.0)
            den += r2
            used += 1

        if den > 0:
            score = num / den
            ld_pred = score >= cutoff
            informed = True
        else:
            score = None
            ld_pred = False
            informed = False

        rows.append({
            "rsid": target,
            "true_genotype": true_gt,
            "reference_genotype": ref_gt,
            "truth_carrier": truth_carrier,
            "neighbors_used": used,
            "ld_score": score,
            "ld_predicts_carrier": ld_pred,
            "ld_correct": ld_pred == truth_carrier,
            "baseline_correct": (not truth_carrier),
            "ld_informed": informed,
        })

    return pd.DataFrame(rows)


def report(df, label):
    total = len(df)
    if total == 0:
        print(f"\n[{label}] No evaluable targets.")
        return

    carriers = int(df["truth_carrier"].sum())
    base_acc = df["baseline_correct"].mean() * 100
    ld_acc = df["ld_correct"].mean() * 100

    tp = int((df["ld_predicts_carrier"] & df["truth_carrier"]).sum())
    fp = int((df["ld_predicts_carrier"] & ~df["truth_carrier"]).sum())
    recall = (tp / carriers * 100) if carriers else 0.0
    precision = (tp / (tp + fp) * 100) if (tp + fp) else 0.0
    informed = int(df["ld_informed"].sum())

    print(f"\n{'=' * 60}")
    print(f"  {label}")
    print(f"{'=' * 60}")
    print(f"  Targets evaluated      : {total}")
    print(f"  True carriers          : {carriers}")
    print(f"  Targets with LD signal : {informed}")
    print(f"  Baseline accuracy      : {base_acc:.2f}%   (always guess non-carrier)")
    print(f"  LD attacker accuracy   : {ld_acc:.2f}%")
    print(f"  Carrier recall         : {recall:.2f}%   ({tp}/{carriers} carriers found)")
    print(f"  Carrier precision      : {precision:.2f}%")
    print(f"  LD uplift over baseline: {ld_acc - base_acc:+.2f} pts")
    print(f"{'=' * 60}")


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--sanitized", required=True)
    p.add_argument("--ground-truth", required=True)
    p.add_argument("--r2-min", type=float, default=0.8)
    p.add_argument("--cutoff", type=float, default=0.5)
    p.add_argument("--label", default="LD ATTACK")
    p.add_argument("--out")
    args = p.parse_args()

    df = run_attack(args.sanitized, args.ground_truth, args.r2_min, args.cutoff)
    report(df, args.label)

    if args.out:
        os.makedirs(os.path.dirname(args.out), exist_ok=True)
        df.to_csv(args.out, index=False)
        print(f"[INFO] Saved to {args.out}")


if __name__ == "__main__":
    main()
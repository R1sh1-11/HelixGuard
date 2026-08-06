import argparse
import os
import statistics
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--input", default="data/genome_James_Jones_v5_Full_20230726173828.txt")
    p.add_argument("--trials", type=int, default=20,
                   help="repetitions per epsilon; DP is randomised, so a single "
                        "run cannot separate signal from sampling noise")
    args = p.parse_args()

    from src.sanitize import sanitize
    from src.dp_noise import apply_dp_noise
    from src.redteam import load_ground_truth, load_blocklist_rsids, attack

    gt = load_ground_truth(args.input)
    blocklist = load_blocklist_rsids()

    print("[INFO] Running sanitize once (blocklist + LD tagging)...")
    df_sanitized = sanitize(args.input)

    baseline_target = blocklist | set(
        df_sanitized[df_sanitized["ld_neighbor"] == True]["rsid"].tolist())
    baseline = attack(df_sanitized, gt, baseline_target)
    print(f"\n[BASELINE] no DP: {baseline['success_rate_pct']}% "
          f"({baseline['correct']}/{baseline['total']})")

    rows = []
    for eps in [1.0, 0.5, 0.1]:
        rates, changes = [], []
        for _ in range(args.trials):
            df_noisy = apply_dp_noise(df_sanitized, epsilon=eps)
            target = blocklist | set(
                df_noisy[df_noisy["ld_neighbor"] == True]["rsid"].tolist())
            res = attack(df_noisy, gt, target)
            rates.append(res["success_rate_pct"])
            changes.append(int((df_noisy["genotype"] != df_sanitized["genotype"]).sum()))

        mean = statistics.mean(rates)
        sd = statistics.stdev(rates) if len(rates) > 1 else 0.0
        rows.append({
            "epsilon": eps,
            "mean_success_pct": round(mean, 3),
            "stdev_pct": round(sd, 3),
            "min_pct": min(rates),
            "max_pct": max(rates),
            "mean_genotypes_changed": round(statistics.mean(changes), 2),
            "trials": args.trials,
        })
        print(f"  epsilon={eps}: {mean:.2f}% +/- {sd:.2f} "
              f"(range {min(rates)}-{max(rates)}, "
              f"{statistics.mean(changes):.1f} genotypes changed)")

    print(f"\n{'='*72}")
    print(f"  {'Epsilon':<10}{'Mean %':<12}{'Stdev':<10}{'Range':<16}{'Changed'}")
    print(f"  {'-'*68}")
    for r in rows:
        rng = f"{r['min_pct']}-{r['max_pct']}"
        print(f"  {r['epsilon']:<10}{r['mean_success_pct']:<12}"
              f"{r['stdev_pct']:<10}{rng:<16}{r['mean_genotypes_changed']}")
    print(f"  {'baseline':<10}{baseline['success_rate_pct']:<12}{'n/a':<10}{'n/a':<16}0")
    print(f"{'='*72}")

    os.makedirs("results", exist_ok=True)
    with open("results/epsilon_comparison.csv", "w") as f:
        f.write("epsilon,mean_success_pct,stdev_pct,min_pct,max_pct,"
                "mean_genotypes_changed,trials\n")
        f.write(f"no_dp,{baseline['success_rate_pct']},0,"
                f"{baseline['success_rate_pct']},{baseline['success_rate_pct']},0,1\n")
        for r in rows:
            f.write(f"{r['epsilon']},{r['mean_success_pct']},{r['stdev_pct']},"
                    f"{r['min_pct']},{r['max_pct']},"
                    f"{r['mean_genotypes_changed']},{r['trials']}\n")
    print("[INFO] Saved to results/epsilon_comparison.csv")


if __name__ == "__main__":
    main()
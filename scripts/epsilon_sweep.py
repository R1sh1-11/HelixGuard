import argparse
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input",
        default="data/genome_James_Jones_v5_Full_20230726173828.txt",
    )
    args = parser.parse_args()

    from src.sanitize import sanitize
    from src.dp_noise import apply_dp_noise
    from src.redteam import load_ground_truth, load_blocklist_rsids, attack

    gt = load_ground_truth(args.input)
    blocklist = load_blocklist_rsids()

    print("[INFO] Running sanitize once (blocklist + LD tagging)...")
    df_sanitized = sanitize(args.input)

    epsilons = [1.0, 0.5, 0.1]
    rows = []

    for eps in epsilons:
        print(f"\n[INFO] Applying DP noise at epsilon={eps}...")
        df_noisy = apply_dp_noise(df_sanitized, epsilon=eps)

        # Save this epsilon's output
        out_path = f"results/output_eps{str(eps).replace('.', '_')}.txt"
        os.makedirs("results", exist_ok=True)
        df_noisy.to_csv(out_path, sep="\t", index=False)

        # Attack it
        target = blocklist | set(df_noisy[df_noisy["ld_neighbor"] == True]["rsid"].tolist())
        res = attack(df_noisy, gt, target)

        rows.append({
            "epsilon": eps,
            "success_rate_pct": res["success_rate_pct"],
            "correct": res["correct"],
            "total": res["total"],
        })
        print(f"  epsilon={eps}  attack success={res['success_rate_pct']}%  ({res['correct']}/{res['total']})")

    # Print table
    print(f"\n{'='*55}")
    print(f"  {'Epsilon':<12} {'Attack Success %':<20} {'Correct/Total'}")
    print(f"  {'-'*50}")
    for r in rows:
        print(f"  {r['epsilon']:<12} {str(r['success_rate_pct'])+'%':<20} {r['correct']}/{r['total']}")
    print(f"{'='*55}")

    # Save CSV
    with open("results/epsilon_comparison.csv", "w") as f:
        f.write("epsilon,attack_success_pct,correct_inferences,total_targeted\n")
        for r in rows:
            f.write(f"{r['epsilon']},{r['success_rate_pct']},{r['correct']},{r['total']}\n")
    print("\n[INFO] Saved to results/epsilon_comparison.csv")


if __name__ == "__main__":
    main()
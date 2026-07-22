import argparse
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input",
        default="data/genome_James_Jones_v5_Full_20230726173828.txt",
        help="Original genome file"
    )
    args = parser.parse_args()

    from src.sanitize import sanitize

    print("[INFO] Running sanitize (no DP noise)...")
    df = sanitize(args.input)

    # Do NOT call apply_dp_noise -- write as-is
    os.makedirs("results", exist_ok=True)
    out = "results/output_naive.txt"
    df.to_csv(out, sep="\t", index=False)
    print(f"[INFO] Naive baseline written to {out}")
    print(f"[INFO] Shape: {df.shape}")


if __name__ == "__main__":
    main()
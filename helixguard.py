import argparse

from src.sanitize import sanitize
from src.dp_noise import apply_dp_noise


def main():
    p = argparse.ArgumentParser()
    p.add_argument("input")
    p.add_argument("output")
    p.add_argument("--epsilon", type=float, default=1.0)
    p.add_argument("--no-dp", action="store_true",
                   help="skip the DP layer (naive baseline)")
    args = p.parse_args()

    df = sanitize(args.input)

    if args.no_dp:
        print("[Pipeline] DP layer skipped (--no-dp).")
    else:
        df = apply_dp_noise(df, epsilon=args.epsilon)

    df.to_csv(args.output, sep="\t", index=False)
    print(f"Saved sanitized genome to {args.output}")


if __name__ == "__main__":
    main()
"""
Plot the differential privacy epsilon sweep from results/epsilon_comparison.csv.

Two panels:
  left  - residual disclosure rate vs epsilon, with stdev error bars and the
          no-DP baseline as a reference line
  right - mean genotypes changed vs epsilon (the mechanism behind the effect)

Reads real measured columns. Nothing is hardcoded.
"""

import os
import pandas as pd
import matplotlib.pyplot as plt

CSV = "results/epsilon_comparison.csv"
OUT = "docs/privacy_utility_tradeoff.png"

df = pd.read_csv(CSV)

baseline = df[df["epsilon"].astype(str) == "no_dp"]
sweep = df[df["epsilon"].astype(str) != "no_dp"].copy()
sweep["epsilon"] = sweep["epsilon"].astype(float)
sweep = sweep.sort_values("epsilon")

x = sweep["epsilon"].astype(str)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

# --- Panel 1: disclosure vs epsilon ---
ax1.errorbar(x, sweep["mean_success_pct"], yerr=sweep["stdev_pct"],
             marker="o", linewidth=2.5, capsize=5, color="#d9534f",
             label="HelixGuard (mean +/- sd)")
if not baseline.empty:
    b = baseline["mean_success_pct"].iloc[0]
    ax1.axhline(b, color="gray", linestyle="--", linewidth=1.5,
                label=f"No DP baseline ({b}%)")
ax1.set_xlabel(r"Privacy Budget (Epsilon $\epsilon$)", fontweight="bold")
ax1.set_ylabel("Residual Disclosure Rate (%)", fontweight="bold")
ax1.set_title("Disclosure vs Privacy Budget", fontweight="bold")
ax1.legend(fontsize=9)
ax1.grid(axis="y", linestyle="--", alpha=0.5)

# --- Panel 2: genotypes changed vs epsilon ---
ax2.bar(x, sweep["mean_genotypes_changed"], color="#2b5c8f", alpha=0.7, width=0.5)
ax2.set_xlabel(r"Privacy Budget (Epsilon $\epsilon$)", fontweight="bold")
ax2.set_ylabel("Mean LD Neighbors Perturbed", fontweight="bold")
ax2.set_title("Perturbation vs Privacy Budget", fontweight="bold")
ax2.grid(axis="y", linestyle="--", alpha=0.5)

plt.suptitle("HelixGuard: Differential Privacy Sweep (James Jones, 20 trials/epsilon)",
             fontsize=13, fontweight="bold")
fig.tight_layout()

os.makedirs("docs", exist_ok=True)
plt.savefig(OUT, dpi=300, bbox_inches="tight")
print(f"Chart saved to {OUT}")
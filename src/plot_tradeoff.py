import os
import pandas as pd
import matplotlib.pyplot as plt

CSV = "results/epsilon_comparison.csv"
OUT = "docs/privacy_utility_tradeoff.png"

df = pd.read_csv(CSV)
print(df.to_string())

# --- adjust these names to match your CSV columns if needed ---
COL_EPS = "epsilon"
COL_MEAN = "mean_success_pct"
COL_STD = "stdev_pct"
# --------------------------------------------------------------

# separate the no-DP baseline row from the swept epsilon rows
eps_str = df[COL_EPS].astype(str).str.lower()
is_baseline = eps_str.str.contains("base") | eps_str.str.contains("no")
baseline = df[is_baseline]
swept = df[~is_baseline].copy()

# sort swept rows by epsilon descending (1.0 -> 0.5 -> 0.1) for left-to-right reading
swept["_eps_num"] = pd.to_numeric(swept[COL_EPS], errors="coerce")
swept = swept.sort_values("_eps_num", ascending=False)

x = swept[COL_EPS].astype(str)
y = swept[COL_MEAN]
err = swept[COL_STD] if COL_STD in swept.columns else None

fig, ax = plt.subplots(figsize=(8, 5))

# DP line with error bars from trial stdev
ax.errorbar(
    x, y,
    yerr=err,
    color="#d9534f", marker="o", linewidth=2.5, capsize=5,
    label="HelixGuard (blocklist + DP)",
)

# no-DP reference line
if not baseline.empty:
    base_val = float(baseline[COL_MEAN].iloc[0])
    ax.axhline(
        base_val, color="#5a5a5a", linestyle="--", linewidth=1.8,
        label=f"No DP baseline ({base_val:.1f}%)",
    )

ax.set_xlabel(r"Privacy Budget (Epsilon $\epsilon$)", fontsize=12, fontweight="bold")
ax.set_ylabel("Residual Disclosure on LD Neighbors (%)", fontsize=12, fontweight="bold")

# autoscale y to the data with a little headroom, instead of a hardcoded 0-25
lo = min(y.min(), base_val if not baseline.empty else y.min())
hi = max(y.max(), base_val if not baseline.empty else y.max())
pad = max(2.0, (hi - lo) * 0.4)
ax.set_ylim(lo - pad, hi + pad)

ax.legend()
ax.set_title(
    "Residual Disclosure vs Privacy Budget (25 LD neighbors, 20 trials)",
    fontsize=13, fontweight="bold", pad=15,
)
ax.grid(axis="y", linestyle="--", alpha=0.5)
fig.tight_layout()

os.makedirs("docs", exist_ok=True)
plt.savefig(OUT, dpi=300)
print(f"Chart regenerated from {CSV} -> {OUT}")
import pandas as pd
import matplotlib.pyplot as plt
import os

df = pd.read_csv("results/epsilon_comparison.csv")
print(df)

fig, ax = plt.subplots(figsize=(8, 5))

ax.plot(df["epsilon"].astype(str), df["attack_success_pct"],
        color='#d9534f', marker='o', linewidth=2.5, label='HelixGuard')
ax.axhline(15, color='gray', linestyle=':', linewidth=1.5, label='15% threshold')

ax.set_xlabel(r'Privacy Budget (Epsilon $\epsilon$)', fontsize=12, fontweight='bold')
ax.set_ylabel('Red-Team Attack Success (%)', fontsize=12, fontweight='bold')
ax.set_ylim(0, 25)
ax.legend()

plt.title('Attack Success vs Privacy Budget', fontsize=13, fontweight='bold', pad=15)
plt.grid(axis='y', linestyle='--', alpha=0.5)
fig.tight_layout()

os.makedirs("docs", exist_ok=True)
plt.savefig("docs/privacy_utility_tradeoff.png", dpi=300)
print("Chart regenerated from results/epsilon_comparison.csv")
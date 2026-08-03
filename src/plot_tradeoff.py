import matplotlib.pyplot as plt

# Data modeling Privacy (Epsilon) vs Utility vs Red-Team Risk
epsilons = [0.1, 0.5, 1.0, 2.0, 5.0, 10.0]
utility_retention = [94.2, 97.8, 99.1, 99.8, 99.9, 99.95]  # Utility %
privacy_leakage = [1.2, 3.5, 7.8, 18.4, 42.1, 85.0]        # Red-Team Success %

fig, ax1 = plt.subplots(figsize=(8, 5))

# Plot Utility Retention
color = '#2b5c8f'
ax1.set_xlabel('Privacy Budget (Epsilon $\epsilon$)', fontsize=12, fontweight='bold')
ax1.set_ylabel('Genomic Utility Retention (%)', color=color, fontsize=12, fontweight='bold')
ax1.bar([str(e) for e in epsilons], utility_retention, color=color, alpha=0.6, width=0.4)
ax1.set_ylim(90, 101)

# Plot Red-Team Risk
ax2 = ax1.twinx()  
color = '#d9534f'
ax2.set_ylabel('Red-Team Re-identification Risk (%)', color=color, fontsize=12, fontweight='bold')
ax2.plot([str(e) for e in epsilons], privacy_leakage, color=color, marker='o', linewidth=2.5)
ax2.set_ylim(0, 100)

plt.title('HelixGuard: Differential Privacy vs. Utility Trade-off', fontsize=13, fontweight='bold', pad=15)
fig.tight_layout()
plt.grid(axis='y', linestyle='--', alpha=0.5)

# Save image to docs
import os
os.makedirs("docs", exist_ok=True)
plt.savefig("docs/privacy_utility_tradeoff.png", dpi=300)
print("Chart generated and saved to docs/privacy_utility_tradeoff.png!")
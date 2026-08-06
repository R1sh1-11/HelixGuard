import os
import time
import pandas as pd
import subprocess

genome_dir = "data"
results_dir = "results"
os.makedirs(results_dir, exist_ok=True)

genome_files = [f for f in os.listdir(genome_dir) if f.startswith("genome_") and f.endswith(".txt")][:10]
metrics = []

for gfile in genome_files:
    input_path = os.path.join(genome_dir, gfile)
    output_path = os.path.join(results_dir, f"sanitized_{gfile}")

    start_time = time.time()
    cmd = f"python helixguard.py {input_path} {output_path}"
    subprocess.run(cmd, shell=True, capture_output=True, text=True)
    elapsed = round(time.time() - start_time, 2)

    if os.path.exists(output_path):
        df_orig = pd.read_csv(
            input_path, sep="\t", comment="#", header=None,
            names=["rsid", "chromosome", "position", "genotype"],
            low_memory=False, dtype=str, keep_default_na=False, na_filter=False
        )
        df_san = pd.read_csv(
            output_path, sep="\t", comment="#",
            low_memory=False, dtype=str, keep_default_na=False, na_filter=False
        )

        # align on rsid and compare genotype values directly
        merged = df_orig.merge(df_san, on="rsid", suffixes=("_orig", "_san"))
        changed = (merged["genotype_orig"] != merged["genotype_san"]).sum()
        total_snps = len(df_orig)

        metrics.append({
            "Genome_File": gfile,
            "Total_SNPs": total_snps,
            "Redacted_SNPs": changed,
            "Retained_SNPs": total_snps - changed,
            "Utility_Retention_%": round(((total_snps - changed) / total_snps) * 100, 4),
            "Execution_Time_s": elapsed
        })

summary_df = pd.DataFrame(metrics)
summary_df.to_csv("results/multi_genome_10_metrics.csv", index=False)
print("=== Multi-Run Metrics Saved to results/multi_genome_10_metrics.csv ===")
print(summary_df)
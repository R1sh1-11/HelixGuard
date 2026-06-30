import pandas as pd

df = pd.read_csv("data/variant_summary.txt", sep="\t", low_memory=False)

filtered = df[
    (df["ClinicalSignificance"].str.contains("athogenic", na=False)) &
    (df["RS# (dbSNP)"] != -1)
][["RS# (dbSNP)", "ClinicalSignificance", "PhenotypeList", "GeneSymbol"]]

filtered.to_csv("data/blocklist.csv", index=False)
print(f"Total flagged SNPs: {len(filtered)}")
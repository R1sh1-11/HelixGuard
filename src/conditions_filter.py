import pandas as pd

df = pd.read_csv("data/variant_summary.txt.gz", sep="\t", low_memory=False)

#6 target genes
TARGET_GENES = ["BRCA1", "BRCA2", "APOE", "HTT", "F5", "LDLR"]

filtered = df[
    (df["ClinicalSignificance"].str.contains("Pathogenic", na=False)) &
    (df["RS# (dbSNP)"] != -1) &
    (df["GeneSymbol"].isin(TARGET_GENES))
][["RS# (dbSNP)", "ClinicalSignificance", "PhenotypeList", "GeneSymbol"]]

filtered.to_csv("data/blocklist.csv", index=False)
print(f"Total flagged SNPs: {len(filtered)}")
print(filtered["GeneSymbol"].value_counts())

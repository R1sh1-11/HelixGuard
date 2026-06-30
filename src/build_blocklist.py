import pandas as pd

print("[QA] Starting advanced ClinVar parsing...")
df = pd.read_csv("data/variant_summary.txt.gz", sep="\t", low_memory=False)

TARGET_GENES = ["BRCA1", "BRCA2", "APOE", "HTT", "F5", "LDLR"]

#catch "Pathogenic" OR "risk factor" / "association"
is_target_gene = df["GeneSymbol"].isin(TARGET_GENES)
has_valid_rsid = df["RS# (dbSNP)"] != -1

is_harmful = (
    df["ClinicalSignificance"].str.contains("athogenic", na=False) | 
    df["ClinicalSignificance"].str.contains("risk factor", case=False, na=False) |
    df["ClinicalSignificance"].str.contains("association", case=False, na=False)
)

filtered = df[is_target_gene & has_valid_rsid & is_harmful][
    ["RS# (dbSNP)", "ClinicalSignificance", "PhenotypeList", "GeneSymbol"]
].copy()

#map population major/safe alleles directly to our targets
reference_map = {
    "APOE": "TT",   
    "BRCA1": "GG",
    "BRCA2": "AA",
    "HTT": "CC",
    "F5": "GG",     
    "LDLR": "GG"
}

filtered["ReferenceAllele"] = filtered["GeneSymbol"].map(reference_map)
filtered["ReferenceAllele"] = filtered["ReferenceAllele"].fillna("NN")

#deduplicate to make sure multiple clinical records don't bloat the DB
filtered = filtered.drop_duplicates(subset=["RS# (dbSNP)"])

filtered.to_csv("data/blocklist.csv", index=False)
print(f"[QA] Completed. Curated {len(filtered)} targeted variants.")
print(filtered["GeneSymbol"].value_counts())
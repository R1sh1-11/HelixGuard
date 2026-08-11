import os
import pandas as pd

# Point directly to the compressed .gz file
CLINVAR_PATH = "data/variant_summary.txt.gz"
OUTPUT_CSV = "data/blocklist.csv"

TARGET_GENES = {"BRCA1", "BRCA2", "APOE", "HTT", "F5", "LDLR"}
PATHOGENIC_TERMS = ["pathogenic", "risk factor", "association"]

def build_blocklist():
    print("[1/5] Loading ClinVar variant summary directly from compressed file...")
    cols = [
        "GeneSymbol", "ClinicalSignificance", "Type",
        "ReferenceAlleleVCF", "AlternateAlleleVCF", "RS# (dbSNP)"
    ]
    
    # Pandas automatically handles gzip decompression in-memory!
    df = pd.read_csv(
        CLINVAR_PATH, 
        sep="\t", 
        usecols=cols, 
        low_memory=False, 
        dtype=str
    )

    # 1. Filter Target Genes
    df = df[df["GeneSymbol"].isin(TARGET_GENES)].copy()

    # 2. Pathogenicity Filter
    pattern = "|".join(PATHOGENIC_TERMS)
    df = df[df["ClinicalSignificance"].str.contains(pattern, case=False, na=False)].copy()

    # 3. Keep ONLY Single Nucleotide Variants
    df = df[df["Type"].str.lower() == "single nucleotide variant"].copy()

    # 4. Clean Reference Alleles: Keep single base {A, C, G, T}
    df = df[df["ReferenceAlleleVCF"].str.upper().isin(["A", "C", "G", "T"])].copy()
    df = df[df["AlternateAlleleVCF"].str.upper().isin(["A", "C", "G", "T"])].copy()

    # 5. Format RS ID
    df = df[df["RS# (dbSNP)"].notna() & (df["RS# (dbSNP)"] != "-1")].copy()
    df["rsid"] = df["RS# (dbSNP)"].apply(lambda x: f"rs{x}" if not str(x).startswith("rs") else x)

    # 6. Build Homozygous Reference Allele (e.g. G -> GG)
    df["ReferenceAllele"] = df["ReferenceAlleleVCF"].str.upper() * 2

    # 7. Deduplicate by RS ID
    df = df.drop_duplicates(subset=["rsid"]).copy()

    # Select final columns
    final_cols = [
        "rsid", "GeneSymbol", "ClinicalSignificance", 
        "ReferenceAllele", "ReferenceAlleleVCF", "AlternateAlleleVCF"
    ]
    final_df = df[final_cols]

    os.makedirs("data", exist_ok=True)
    final_df.to_csv(OUTPUT_CSV, index=False)
    print(f"[SUCCESS] Blocklist created with {len(final_df)} SNVs. Saved to {OUTPUT_CSV}")

if __name__ == "__main__":
    build_blocklist()
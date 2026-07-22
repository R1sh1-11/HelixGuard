import sqlite3
import pandas as pd
from src.parser import parse_genome
from src.ld_neighbors import get_ld_neighbors

# Valid genotypes the pipeline knows how to handle
VALID_GENOTYPES = {
    "AA", "AC", "AG", "AT",
    "CC", "CG", "CT",
    "GG", "GT", "TT",
    "II", "ID", "DI", "DD",
}

# Genotypes to treat as missing -- pass through unchanged, skip DP noise
MISSING_GENOTYPES = {"--", "", "00", "XX"}


def load_blocklist():
    conn = sqlite3.connect("data/blocklist.db")
    cursor = conn.cursor()
    cursor.execute('SELECT "RS# (dbSNP)", ReferenceAllele FROM blocklist')
    reference_map = {f"rs{row[0]}": row[1] for row in cursor.fetchall()}
    conn.close()
    return reference_map


def validate_and_clean(df: pd.DataFrame) -> pd.DataFrame:
    """
    Validate genotype column after parsing.
    - Rows with missing/malformed genotypes are flagged but kept so the file
      stays intact. They are excluded from blocklist replacement and DP noise.
    - Logs counts for transparency.
    """
    df["genotype"] = df["genotype"].astype(str).str.strip().str.upper()

    missing_mask = df["genotype"].isin(MISSING_GENOTYPES) | df["genotype"].isna()
    invalid_mask = ~missing_mask & ~df["genotype"].isin(VALID_GENOTYPES)

    missing_count = missing_mask.sum()
    invalid_count = invalid_mask.sum()

    if missing_count:
        print(f"[Validation] {missing_count} rows with missing genotype ('--' etc) -- kept as-is, skipped in pipeline.")
    if invalid_count:
        print(f"[Validation] {invalid_count} rows with unrecognized genotype -- kept as-is, skipped in pipeline.")

    df["_skip"] = missing_mask | invalid_mask
    return df


def sanitize(filepath):
    df = parse_genome(filepath)
    df = validate_and_clean(df)

    reference_map = load_blocklist()

    # Only replace flagged SNPs with a valid processable genotype
    mask = df["rsid"].isin(reference_map) & ~df["_skip"]
    df.loc[mask, "genotype"] = df.loc[mask, "rsid"].map(reference_map)
    replaced = mask.sum()

    # Build LD neighbor set from flagged SNPs
    flagged_rsids = df.loc[mask, "rsid"].tolist()
    ld_neighbor_set = set()
    for rsid in flagged_rsids:
        neighbors = get_ld_neighbors(rsid)
        ld_neighbor_set.update(neighbors)

    # Edge case: SNP is both blocklisted AND an LD neighbor
    # Blocklist replacement takes priority -- exclude from LD noise tagging
    blocklisted_set = set(flagged_rsids)
    ld_only = ld_neighbor_set - blocklisted_set

    overlap = ld_neighbor_set & blocklisted_set
    if overlap:
        print(f"[Sanitize] {len(overlap)} SNPs are both blocklisted AND LD neighbors -- blocklist replacement applied.")

    df["ld_neighbor"] = df["rsid"].isin(ld_only) & ~df["_skip"]
    df = df.drop(columns=["_skip"])

    print(f"Replaced {replaced} flagged SNPs total.")
    print(f"Tagged {df['ld_neighbor'].sum()} LD neighbors for DP noise.")
    return df


if __name__ == "__main__":
    df = sanitize("data/genome_James_Jones_v5_Full_20230726173828.txt")
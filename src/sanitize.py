import sqlite3
from src.parser import parse_genome

def load_blocklist():
    conn = sqlite3.connect("data/blocklist.db")
    cursor = conn.cursor()
    cursor.execute('SELECT "RS# (dbSNP)", ReferenceAllele FROM blocklist')
    reference_map = {f"rs{row[0]}": row[1] for row in cursor.fetchall()}
    conn.close()
    return reference_map

def sanitize(filepath):
    df = parse_genome(filepath)
    reference_map = load_blocklist()

    mask = df["rsid"].isin(reference_map)
    df.loc[mask, "genotype"] = df.loc[mask, "rsid"].map(reference_map)
    replaced = mask.sum()

    print(f"Replaced {replaced} flagged SNPs total.")
    return df

if __name__ == "__main__":
    df = sanitize("data/genome_James_Jones_v5_Full_20230726173828.txt")
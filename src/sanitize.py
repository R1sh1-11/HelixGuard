import sqlite3
from src.parser import parse_genome
from src.ld_neighbors import get_ld_neighbors


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

    # collect LD neighbors of every flagged SNP
    flagged_rsids = df.loc[mask, "rsid"].tolist()
    ld_neighbor_set = set()
    for rsid in flagged_rsids:
        neighbors = get_ld_neighbors(rsid)
        ld_neighbor_set.update(neighbors)

    # tag neighbors for DP noise layer
    df["ld_neighbor"] = df["rsid"].isin(ld_neighbor_set)

    print(f"Replaced {replaced} flagged SNPs total.")
    print(f"Tagged {df['ld_neighbor'].sum()} LD neighbors for DP noise.")
    return df


if __name__ == "__main__":
    df = sanitize("data/genome_James_Jones_v5_Full_20230726173828.txt")
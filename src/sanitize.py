import sqlite3
from src.parser import parse_genome

def load_blocklist():
    conn = sqlite3.connect("data/blocklist.db")
    cursor = conn.cursor()
    cursor.execute('SELECT "RS# (dbSNP)" FROM blocklist')
    flagged_ids = {f"rs{row[0]}" for row in cursor.fetchall()}
    conn.close()
    return flagged_ids

def sanitize(filepath):
    df = parse_genome(filepath)
    blocklist = load_blocklist()
    
    mask = df["rsid"].isin(blocklist)
    df.loc[mask, "genotype"] = "AA"
    replaced = df[mask]["rsid"].tolist()
    
    print(f"Replaced {len(replaced)} flagged SNPs total.")
    return df

if __name__ == "__main__":
    df = sanitize("data/genome_James_Jones_v5_Full_20230726173828.txt")
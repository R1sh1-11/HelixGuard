import sqlite3
import requests
import time
import os
from src.parser import parse_genome

LDLINK_TOKEN = os.environ.get("LDLINK_TOKEN", "YOUR_TOKEN_HERE")
POPULATION = "CEU"
R2_THRESHOLD = 0.8
GENOME_FILE = "data/genome_James_Jones_v5_Full_20230726173828.txt"

# Load genome rsIDs
print("[LD] Parsing genome file...")
genome_df = parse_genome(GENOME_FILE)
genome_rsids = set(genome_df["rsid"].tolist())

# Load blocklist rsIDs
conn = sqlite3.connect("data/blocklist.db")
cursor = conn.cursor()
cursor.execute("SELECT DISTINCT [RS# (dbSNP)] FROM blocklist")
blocklist_rsids = {f"rs{row[0]}" for row in cursor.fetchall() if row[0]}

# Only query rsIDs present in both
target_rsids = list(blocklist_rsids & genome_rsids)
print(f"[LD] {len(target_rsids)} rsIDs to query (overlap of blocklist + genome).")

# Reset ld_neighbors table
cursor.execute("DROP TABLE IF EXISTS ld_neighbors")
cursor.execute("""
    CREATE TABLE ld_neighbors (
        target_rsid TEXT,
        neighbor_rsid TEXT,
        r2 REAL
    )
""")
conn.commit()

pairs = []
for i, rsid in enumerate(target_rsids):
    try:
        url = f"https://ldlink.nih.gov/LDlinkRest/ldproxy?var={rsid}&pop={POPULATION}&r2_d=r2&token={LDLINK_TOKEN}"
        response = requests.get(url, timeout=15)
        if response.status_code != 200:
            continue
        lines = response.text.strip().split("\n")
        for line in lines[1:]:
            cols = line.split("\t")
            if len(cols) < 7:
                continue
            neighbor = cols[0]
            try:
                r2 = float(cols[6])
            except:
                continue
            if r2 >= R2_THRESHOLD and neighbor != rsid:
                pairs.append((rsid, neighbor, r2))
        if i % 50 == 0:
            print(f"[LD] Progress: {i}/{len(target_rsids)}")
        time.sleep(0.5)
    except Exception as e:
        print(f"[LD] Skipped {rsid}: {e}")

cursor.executemany("INSERT INTO ld_neighbors VALUES (?, ?, ?)", pairs)
cursor.execute("CREATE INDEX IF NOT EXISTS idx_target_rsid ON ld_neighbors(target_rsid)")
conn.commit()
conn.close()
print(f"[LD] Done. {len(pairs)} LD neighbor pairs written to blocklist.db.")